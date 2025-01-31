#!/usr/bin/env python
"""
cif_emdb_translator.py

This module translates mmCIF files (.cif) into .xml header files
that comply with EMDB 3.x.x.x schemas

Copyright [2014-2016] EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
"""

# pylint: disable=protected-access

__author__ = "Ardan Patwardhan, Sanja Abbott"
__email__ = "ardan@ebi.ac.uk, sanja@ebi.ac.uk"
__date__ = "2019-05-17"
__version__ = "1.0"

import os
import sys
import re
import datetime
import logging
import io
from optparse import OptionParser  # pylint: disable=deprecated-module
from lxml import etree

# Deployment paths
from wwpdb.utils.config.ConfigInfo import getSiteId
from wwpdb.utils.config.ConfigInfoApp import ConfigInfoAppEm
from mmcif.io.IoAdapterCore import IoAdapterCore
from . import emdb
from ..EmdbSchema import EmdbSchema


class Cif(object):
    """Class to represent parsed cif file conforming to needed interface"""

    def __init__(self, container):
        """
        Initialize and stash data
        :param container: DataContainer
        """
        self.__block = container[0]

        # Create a tree
        # print(dir(self.__block))
        # print self.__block.getObjNameList()

    def __getitem__(self, item):
        # Handle if abc in self.cif - key is an integer
        if type(item) == int:  # noqa: E721
            nl = self.__block.getObjNameList()
            if item >= len(nl):
                raise StopIteration
            else:
                return nl[item]

        return self.get(item)

    def get(self, catname, err=False):  # pylint: disable=unused-argument
        """Returns tuple representation of category or none like:

        a= [[('_database_2.database_id', u'PDB'), ('_database_2.database_code', u'0XXX')],
            [('_database_2.database_id', u'WWPDB'), ('_database_2.database_code', u'D_10001')],
            [('_database_2.database_id', u'EMDB'), ('_database_2.database_code', u'EMD-0000')]]

            '?' and '.' items not included
        """

        # DepUI reader use to treat a key that was not in dictionary as KeyError for [] or None for get method.
        # However, if in dictionary both methods return [] if not in model file
        dc_obj = self.__block.getObj(catname)
        if not dc_obj:
            return []

        retlist = []

        # print dir(dc_obj)
        attrlist = dc_obj.getAttributeList()
        for row in range(dc_obj.getRowCount()):
            rowlist = []
            for attr in attrlist:
                value = dc_obj.getValue(attributeName=attr, rowIndex=row)
                if value not in ["?", "."]:
                    rowlist.append(("_%s.%s" % (catname, attr), value))
            retlist.append(rowlist)

        return retlist


class CifEMDBTranslator(object):
    """Class for translating files from/to cif from EMDB XML 3.0"""

    translation_log = None  # type: TranslationLog

    class Constants(object):
        """
        There are many constants in use for the translation.
        They have been collected here for ease of use.
        """

        XML_OUT_VERSION = "3.0.10.0"
        XML_VERSION = XML_OUT_VERSION.replace('.', '_')

        # Cif categories
        CITATION = "citation"
        CITATION_AUTHOR = "citation_author"
        DATABASE_2 = "database_2"
        EM_ADMIN = "em_admin"
        EM_DEPUI = "em_depui"
        EM_EULER_ANGLE_ASSIGNMENT = "em_euler_angle_assignment"
        EM_AUTHOR_LIST = "em_author_list"
        AUDIT_AUTHOR = "audit_author"
        EM_BUFFER = "em_buffer"
        EM_BUFFER_COMPONENT = "em_buffer_component"
        EM_DB_REFERENCE = "em_db_reference"
        EM_DB_REFERENCE_AUXILIARY = "EM_DB_REFERENCE_auxiliary"
        EM_CRYSTAL_FORMATION = "em_crystal_formation"
        EM_DIFFRACTION_SHELL = "em_diffraction_shell"
        EM_DIFFRACTION_STATS = "em_diffraction_stats"
        EM_CTF_CORRECTION = "em_ctf_correction"
        EM_EMBEDDING = "em_embedding"
        EM_FIDUCIAL_MARKERS = "em_fiducial_markers"
        EM_FINAL_CLASSIFICATION = "em_final_classification"
        EM_FINAL_2D_CLASSIFICATION = "em_final_two_d_classification"
        EM_3D_RECONSTRUCTION = "em_3d_reconstruction"
        EM_FSC_CURVE = "em_fsc_curve"
        EM_SAMPLE_SUPPORT = "em_sample_support"
        EM_GRID_PRETREATMENT = "em_grid_pretreatment"
        EM_HELICAL_ENTITY = "em_helical_entity"
        EM_HIGH_PRESSURE_FREEZING = "em_high_pressure_freezing"
        EM_IMAGE_SCANS = "em_image_scans"
        EM_MAP = "em_map"
        EM_IMAGING = "em_imaging"
        EM_DIFFRACTION = "em_diffraction"
        EM_TOMOGRAPHY = "em_tomography"
        EM_3D_FITTING = "em_3d_fitting"
        EM_3D_FITTING_LIST = "em_3d_fitting_list"
        EM_ENTITY_ASSEMBLY_MOLWT = "em_entity_assembly_molwt"
        EM_ENTITY_ASSEMBLY_NATURALSOURCE = "em_entity_assembly_naturalsource"
        EM_ENTITY_ASSEMBLY_SYNTHETIC = "em_entity_assembly_synthetic"
        EM_IMAGE_PROCESSING = "em_image_processing"
        EM_IMAGE_RECORDING = "em_image_recording"
        EM_OBSOLETE = "em_obsolete"
        EM_PARTICLE_SELECTION = "em_particle_selection"
        EM_ENTITY_ASSEMBLY_RECOMBINANT = "em_entity_assembly_recombinan"
        EM_FOCUSED_ION_BEAM = "em_focused_ion_beam"
        EM_ULTRAMICROTOMY = "em_ultramicrotomy"
        EM_SHADOWING = "em_shadowing"
        EM_SOFTWARE = "em_software"
        EM_SPECIALIST_OPTICS = "em_imaging_optics"
        EM_SPECIMEN = "em_specimen"
        EM_STAINING = "em_staining"
        EM_START_MODEL = "em_start_model"
        EM_EXPERIMENT = "em_experiment"
        EM_SUPERSEDE = "em_supersede"
        EM_SUPPORT_FILM = "em_support_film"
        EM_ENTITY_ASSEMBLY = "em_entity_assembly"
        EM_SINGLE_PARTICLE_ENTITY = "em_single_particle_entity"
        EM_3D_CRYSTAL_ENTITY = "em_3d_crystal_entity"
        EM_TOMOGRAPHY_SPECIMEN = "em_tomography_specimen"
        EM_2D_CRYSTAL_ENTITY = "em_2d_crystal_entity"
        EM_VIRUS_ENTITY = "em_virus_entity"
        EM_VIRUS_NATURAL_HOST = "em_virus_natural_host"
        EM_VIRUS_SHELL = "em_virus_shell"
        EM_VITRIFICATION = "em_vitrification"
        EM_VOLUME_SELECTION = "em_volume_selection"
        ENTITY = "entity"
        ENTITY_POLY = "entity_poly"
        ENTITY_SRC_GEN = "entity_src_gen"
        ENTITY_SRC_NAT = "entity_src_nat"
        PDBX_ENTITY_SRC_SYN = "pdbx_entity_src_syn"
        EXPTL = "exptl"
        PDBX_DATABASE_STATUS = "pdbx_database_status"
        PDBX_DATABASE_RELATED = "pdbx_database_related"
        PDBX_ENTITY_NONPOLY = "pdbx_entity_nonpoly"
        PDBX_DEPOSITOR_INFO = "pdbx_struct_ref_seq_depositor_info"
        STRUCT_REF = "struct_ref"
        PDBX_OBS_SPR = "pdbx_database_PDB_obs_spr"
        PDBX_AUDIT_SUPPORT = "pdbx_audit_support"
        PDBX_CONTACT_AUTHOR = "pdbx_contact_author"
        STRUCT = "struct"
        STRUCT_KEYWORDS = "struct_keywords"
        PDBX_DICT_ITEM_MAPPING = "pdbx_dict_item_mapping"
        PDBX_AUDIT_REVISION_HISTORY = "pdbx_audit_revision_history"
        PDBX_AUDIT_REVISION_DETAILS = "pdbx_audit_revision_details"
        PDBX_AUDIT_REVISION_GROUP = "pdbx_audit_revision_group"
        PDBX_AUDIT_REVISION_CATEGORY = "pdbx_audit_revision_category"
        PDBX_AUDIT_REVISION_ITEM = "pdbx_audit_revision_item"

        # Keys
        K_EM_DIFFRACTION_STATS_ID = "em_diffraction_stats_id"
        K_SAMPLE_SUPPORT_ID = "sample_support_id"
        K_IMAGE_PROCESSING_ID = "image_processing_id"
        K_IMAGE_RECORDING_ID = "image_recording_id"
        K_IMAGING_ID = "imaging_id"
        K_3D_FITTING_ID = "3d_fitting_id"
        K_SPECIMEN_ID = "specimen_id"
        K_ENTITY_ASSEMBLY_ID = "entity_assembly_id"
        K_EM_TOMOGRAPHY_SPECIMEN_ID = "em_tomography_specimen_id"
        K_ENTITY_ID = "entity_id"
        K_ID = "id"
        K_BUFFER_ID = "buffer_id"

        # Software categories
        SOFT_CLASSIFICATION = "CLASSIFICATION"
        SOFT_CTF_CORRECTION = "CTF CORRECTION"
        SOFT_FINAL_EULER_ASSIGNMENT = "FINAL EULER ASSIGNMENT"
        SOFT_IMAGE_ACQUISITION = "IMAGE ACQUISITION"
        SOFT_INITIAL_EULER_ASSIGNMENT = "INITIAL EULER ASSIGNMENT"
        SOFT_MODEL_FITTING = "MODEL FITTING"
        SOFT_PARTICLE_SELECTION = "PARTICLE SELECTION"
        SOFT_RECONSTRUCTION = "RECONSTRUCTION"
        SOFT_MOLECULAR_REPLACEMENT = "MOLECULAR REPLACEMENT"
        SOFT_LATTICE_DISTORTION_CORRECTION = "LATTICE DISTORTION CORRECTION"
        SOFT_SYMMETRY_DETERMINATION = "SYMMETRY DETERMINATION"
        SOFT_CRYSTALLOGRAPHY_MERGING = "CRYSTALLOGRAPHY MERGING"

        # EM methods
        EMM_EC = "electronCrystallography"
        EMM_HEL = "helical"
        EMM_SP = "singleParticle"
        EMM_STOM = "subtomogramAveraging"
        EMM_TOM = "tomography"

        # Extension types
        EXT_RNA_TYPE = "rna"
        EXT_DNA_TYPE = "dna"
        EXT_PROT_TYPE = "protein_or_peptide"
        EXT_SACC_TYPE = "saccharide"
        EXT_OTHM_TYPE = "other_macromolecule"
        EXT_LIG_TYPE = "ligand"
        EXT_COMP_TYPE = "complex"
        EXT_VIR_TYPE = "virus"
        EXT_ORGCELL_TYPE = "organelle_or_cellular_component"
        EXT_TISS_TYPE = "tissue"
        EXT_CELL_TYPE = "cell"
        EXT_BASE_PREPARATION_TYPE = "base_preparation_type"
        EXT_TOMOGRAPHY_PREPARATION_TYPE = "tomography_preparation_type"
        EXT_CRYST_PREPARATION_TYPE = "crystallography_preparation_type"
        EXT_BASE_MICROSCOPY_TYPE = "base_microscopy_type"
        EXT_TOMOGRAPHY_MICROSCOPY_TYPE = "tomography_microscopy_type"
        EXT_CRYST_MICROSCOPY_TYPE = "crystallography_microscopy_type"
        EXT_SP_PROC_TYPE = "singleparticle_processing_type"
        EXT_HEL_PROC_TYPE = "helical_processing_type"
        EXT_TOM_PROC_TYPE = "tomography_processing_type"
        EXT_STOM_PROC_TYPE = "subtomogram_averaging_processing_type"
        EXT_CRYST_PROC_TYPE = "crystallography_processing_type"
        EXT_BASE_SOURCE_TYPE = "base_source_type"

        # Map types
        MAP_PRIMARY = "primary map"
        MAP_HALF = "half map"
        MAP_ADD = "additional map"
        MAP_MASK = "mask"

        # Revision history types
        PRIMARY_MAP = "Primary map"
        HALF_MAP = "Half map"
        MASK = "Mask"
        ADDITIONAL_MAP = "Additional map"
        FSC = "FSC"
        EM_METADATA = "EM metadata"
        IMAGE = "Image"
        STRUCTURE_MODEL = "Structure model"

        # Entity types
        ENT_CYCLIC_PSEUDO_PEPTIDE = "cyclic-pseudo-peptide"
        ENT_DNA = "polydeoxyribonucleotide"
        ENT_DNA_RNA = "polydeoxyribonucleotide/polyribonucleotide hybrid"
        ENT_POLYPEPTIDE_D = "polypeptide(D)"
        ENT_POLYPEPTIDE_L = "polypeptide(L)"
        ENT_RNA = "polyribonucleotide"
        ENT_PEPTIDE_NUCLEIC_ACID = "peptide nucleic acid"
        ENT_POLYSACCHARIDE_D = "polysaccharide(D)"
        ENT_POLYSACCHARIDE_L = "polysaccharide(L)"
        ENT_OTHER = "other"

        # Units
        U_ANG = u"\u212B"  # u'\u00C5'
        U_DEG = "deg"
        U_DEGF = "degrees"
        U_FIB_DOSE_RATE = "ions/(cm^2*s)"
        U_KDA_NM = "kDa/nm"
        U_KEL = "K"
        U_KVOLT = "kV"
        U_MDA = "MDa"
        U_MG_ML = "mg/mL"
        U_MMOL = "mM"
        U_NM = "nm"
        U_NMF = "nanometer"
        U_PAMP = "pA"
        U_NAMP = "nA"
        U_SEC = "s"
        U_KPA = "kPa"
        U_PERCENT = "percentage"
        U_MICROM = u"\u00B5" + "m"
        U_MM = "mm"
        U_MRAD = "mrad"
        U_EOVERANGSQR = "e/" + U_ANG + "^2"
        U_PIXEL = "pixel"
        U_EV = "eV"

        # Other constants
        CIF_EMDB_ASSOC = "associated EM volume"
        CIF_EMDB_OTHER = "other EM volume"
        CIF_EMDB_CONSENSUS = "consensus EM volume"
        CIF_EMDB_FOCUSED = "focused EM volume"
        CIF_AUTHOR_RE = re.compile(r"^([A-Za-z0-9 \'\-.]+), (([A-Z\-]+\.)*)")
        CIF_HALF_MAP_RE = re.compile(r"^D_[0-9]+\_em\-half\-volume\_P([0-9]+)\.map")
        CIF_ADD_MAP_RE = re.compile(r"^D_[0-9]+\_em\-additional\-volume\_P([0-9]+)\.map")
        CIF_EMD_ID_RE = re.compile(r"EMD\-([0-9]){4,}")
        # format: D_1000218615_em-mask-volume_P1.map.V1
        CIF_MSK_MAP_RE = re.compile(r"^D_[0-9]+\_em\-mask\-volume\_P([0-9]+)\.map")
        DA2MDA = 1.0 / 1000000.0
        PROC_SITE_CIF2XML = {"PDBE": "PDBe", "RCSB": "RCSB", "PDBJ": "PDBj", "PDBC": "PDBc"}
        AGG_STATE_CIF2XML = {
            "PARTICLE": "particle",
            "FILAMENT": "filament",
            "2D ARRAY": "twoDArray",
            "3D ARRAY": "threeDArray",
            "HELICAL ARRAY": "helicalArray",
            "TISSUE": "tissue",
            "CELL": "cell",
        }
        MAP_DATA_TYPE_CIF2XML = {
            "Image stored as signed byte": "IMAGE STORED AS SIGNED BYTE",
            "Image stored as signed integer (2 bytes)": "IMAGE STORED AS SIGNED INTEGER (2 BYTES)",
            "Image stored as floating point number (4 bytes)": "IMAGE STORED AS FLOATING POINT NUMBER (4 BYTES)",
        }

        MAP_REVISION_DETAILS_TYPE = {
            "Initial release": "INITIAL_RELEASE",
            "Coordinate replacement": "COORDINATE_REPLACEMENT",
            "Obsolete": "OBSOLETE",
            "Remediation": "REMEDIATION",
            "Data added": "DATA_ADDED",
            "Data updated": "DATA_UPDATED",
            "Data removed": "DATA_REMOVED"
        }

        MAP_REVISION_PROVIDER = {
            "author": "AUTHOR",
            "repository": "REPOSITORY"
        }

        MAP_REVISION_GROUP = {
            "Experimental data": "EXPERIMENTAL_DATA",
            "Derived data": "DERIVED_DATA",
            "Experimental summary": "EXPERIMENTAL_SUMMARY",
            "Experimental preparation": "EXPERIMENTAL_PREPARATION",
            "Database references": "DATABASE_REFERENCES",
            "Data collection": "DATA_COLLECTION",
            "Advisory": "ADVISORY",
            "Source and taxonomy": "SOURCE_AND_TAXONOMY",
            "Refinement description": "REFINEMENT_DESCRIPTION",
            "Data processing": "DATA_PROCESSING",
            "Structure summary": "STRUCTURE_SUMMARY",
            "Version format compliance": "VERSION_FORMAT_COMPLIANCE",
            "Other": "OTHER"
        }

        MAP_FILES = {
            'Primary map': 'primary_map',
            'Image': 'image',
            'Half map': 'half_map',
            'Mask': 'mask',
            'FSC': 'fsc',
            'Additional map': 'additional_map',
            'EM metadata': 'metadata',
            'Structure model': 'model'
        }

        INFO_LOG_FILE_NAME = "INFO_cifEMDBTranslation.log"
        WARN_LOG_FILE_NAME = "WARN_cifEMDBTranslation.log"
        ERR_LOG_FILE_NAME = "ERROR_cifEMDBTranslation.log"

        INFO_LOG_STRING = "info_log_string"
        WARN_LOG_STRING = "warn_log_string"
        ERROR_LOG_STRING = "error_log_string"

        INFO_LOG_LEVEL = logging.INFO
        WARN_LOG_LEVEL = logging.WARNING
        ERROR_LOG_LEVEL = logging.ERROR

        NYI = "NOT YET IMPLEMENTED"
        REQUIRED_ALERT = "PROBLEM: "
        NOT_REQUIRED_ALERT = "POTENTIAL PROBLEM: "
        INFO_ALERT = "INFO: "
        VALIDATION_ERROR = "VALIDATION ERROR "
        CHANGE_MADE = "CHANGE MADE: "
        NOT_CHANGED_FOR_NOW = "NOT CHANGED FOR NOW: "

        CENTER_NAMES_USED_FOR_AUTHORS = {
            "Accelerated Technologies Center for Gene to 3D Structure (ATCG3D)",
            "Assembly, Dynamics and Evolution of Cell-Cell and Cell-Matrix Adhesions (CELLMAT)",
            "Atoms-to-Animals: The Immune Function Network (IFN)",
            "Bacterial targets at IGS-CNRS, France (BIGS)",
            "Berkeley Structural Genomics Center (BSGC)",
            "Center for Eukaryotic Structural Genomics (CESG)",
            "Center for High-Throughput Structural Biology (CHTSB)",
            "Center for Membrane Proteins of Infectious Diseases (MPID)",
            "Center for Structural Biology of Infectious Diseases (CSBID)",
            "Center for Structural Genomics of Infectious Diseases (CSGID)",
            "Center for Structures of Membrane Proteins (CSMP)",
            "Center for the X-ray Structure Determination of Human Transporters (TransportPDB)",
            "Chaperone-Enabled Studies of Epigenetic Regulation Enzymes (CEBS)",
            "Enzyme Discovery for Natural Product Biosynthesis (NatPro)",
            "QCRG Structural Biology Consortium",
            "GPCR Network (GPCR)",
            "Integrated Center for Structure and Function Innovation (ISFI)",
            "Israel Structural Proteomics Center (ISPC)",
            "Joint Center for Structural Genomics (JCSG)",
            "Marseilles Structural Genomics Program @ AFMB (MSGP)",
            "Medical Structural Genomics of Pathogenic Protozoa (MSGPP)",
            "Membrane Protein Structural Biology Consortium (MPSBC)",
            "Membrane Protein Structures by Solution NMR (MPSbyNMR)",
            "Midwest Center for Macromolecular Research (MCMR)",
            "Midwest Center for Structural Genomics (MCSG)",
            "Mitochondrial Protein Partnership (MPP)",
            "Montreal-Kingston Bacterial Structural Genomics Initiative (BSGI)",
            "Mycobacterium Tuberculosis Structural Proteomics Project (XMTB)",
            "New York Consortium on Membrane Protein Structure (NYCOMPS)",
            "New York SGX Research Center for Structural Genomics (NYSGXRC)",
            "New York Structural GenomiX Research Consortium (NYSGXRC)",
            "New York Structural Genomics Research Consortium (NYSGRC)",
            "Northeast Structural Genomics Consortium (NESG)",
            "Nucleocytoplasmic Transport: a Target for Cellular Control (NPCXstals)",
            "Ontario Centre for Structural Proteomics (OCSP)",
            "Oxford Protein Production Facility (OPPF)",
            "Paris-Sud Yeast Structural Genomics (YSG)",
            "Partnership for Nuclear Receptor Signaling Code Biology (NHRs)",
            "Partnership for Stem Cell Biology (STEMCELL)",
            "Partnership for T-Cell Biology (TCELL)",
            "Program for the Characterization of Secreted Effector Proteins (PCSEP)",
            "Protein Structure Factory (PSF)",
            "RIKEN Structural Genomics/Proteomics Initiative (RSGI)",
            "Scottish Structural Proteomics Facility (SSPF)",
            "Seattle Structural Genomics Center for Infectious Disease (SSGCID)",
            "South Africa Structural Targets Annotation Database (SASTAD)",
            "Southeast Collaboratory for Structural Genomics (SECSG)",
            "Structural Genomics Consortium (SGC)",
            "Structural Genomics Consortium for Research on Gene Expression (SGCGES)",
            "Structural Genomics of Pathogenic Protozoa Consortium (SGPP)",
            "Structural Proteomics in Europe (SPINE)",
            "Structural Proteomics in Europe 2 (SPINE-2)",
            "Structure 2 Function Project (S2F)",
            "Structure, Dynamics and Activation Mechanisms of Chemokine Receptors (CHSAM)",
            "Structure-Function Analysis of Polymorphic CDI Toxin-Immunity Protein Complexes (UC4CDI)",
            "Structure-Function Studies of Tight Junction Membrane Proteins (TJMP)",
            "Structures of Mtb Proteins Conferring Susceptibility to Known Mtb Inhibitors (MTBI)",
            "TB Structural Genomics Consortium (TBSGC)",
            "Transcontinental EM Initiative for Membrane Protein Structure (TEMIMPS)",
            "Transmembrane Protein Center (TMPC)",
        }

        # em_admin.last_update should not be initialized twice
        MMCIF_TO_XSD = {  # pylint: disable=duplicate-key
            "_em_3d_fitting_list.accession_code": '<xs:element name="access_code"/>',
            "_em_3d_fitting_list.pdb_entry_id": '<xs:element name="access_code"/>',
            "_em_3d_fitting_list.pdb_chain_id": '<xs:element name="chain_id" type="chain_pdb_id" minOccurs="0" maxOccurs="unbounded"/>',
            "_em_3d_fitting_list.chain_id": '<xs:element name="chain_id" type="token" minOccurs="0" maxOccurs="unbounded"/>',
            "_em_3d_fitting_list.pdb_chain_residue_range": '<xs:element name="residue_range" minOccurs="0" maxOccurs="1"/>',
            "_em_3d_fitting_list.chain_residue_range": '<xs:element name="residue_range" minOccurs="0" maxOccurs="1"/>',
            "_em_3d_fitting_list.number_of_copies_in_final_model": '<xs:element name="number_of_copies_in_final_model" minOccurs="0"/>',
            "_em_3d_fitting_list.source_name": '<xs:element name="source_name" minOccurs="0" maxOccurs="1"/>',
            "_em_3d_fitting_list.type": '<xs:element name=initial_model_type" minOccurs="0" maxOccurs="1"/>',
            "_em_3d_fitting_list.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_3d_fitting_list": '<xs:element name="access_code"/>',
            # _em_3d_fitting.entry_id
            # _em_3d_fitting.id
            "_em_3d_fitting.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            # _em_3d_fitting.initial_refinement_model_id
            # _em_3d_fitting.method
            "_em_3d_fitting.overall_b_value": '<xs:element name="overall_bvalue" type="xs:float" minOccurs="0">',
            "_em_3d_fitting.ref_protocol": '<xs:element name="refinement_protocol" minOccurs="0">',
            "_em_3d_fitting.ref_space": '<xs:element name="refinement_space" type="xs:token" minOccurs="0">',
            # _em_3d_fitting.software_name
            "_em_3d_fitting.target_criteria": '<xs:element name="target_criteria" type="xs:token" minOccurs="0">',
            "_em_3d_reconstruction.num_particles": '<xs:element name="number_images_used" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_3d_reconstruction.num_class_averages": '<xs:element name="number_classes_used" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_3d_reconstruction.algorithm": '<xs:element name="algorithm" type="reconstruction_algorithm_type" minOccurs="0"/>',
            "_em_3d_reconstruction.resolution_method": '<xs:element name="resolution_method" minOccurs="0">',
            "_em_3d_reconstruction.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_3d_reconstruction.resolution": '<xs:element name="resolution" minOccurs="0">',
            "_em_buffer.pH": '<xs:element name="ph">',
            "_em_buffer.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_entity_assembly.chimera": '<xs:attribute name="chimera" type="xs:boolean" fixed="true"/>',
            "_em_entity_assembly.id": '<xs:attribute name="supramolecule_id" type="xs:positiveInteger" use="required"/>',
            "_em_entity_assembly.name": '<xs:element name="name" type="sci_name_type">',
            "_em_entity_assembly.go_id": '<xs:element name="category" minOccurs="0">',
            "_em_entity_assembly.parent_id": '<xs:element name="parent" type="xs:nonNegativeInteger">',
            "_em_entity_assembly.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_entity_assembly.entity_id_list": '<xs:element name="macromolecule_list" minOccurs="0">',
            "_em_experiment.id": '<xs:attribute name="structure_determination_id" type="xs:positiveInteger" use="required"/>',
            "_em_experiment.aggregation_state": '<xs:element name="aggregation_state">',
            "_em_experiment.reconstruction_method": '<xs:element name="method">',
            "_em_image_scans.scanner_model": '<xs:element name="scanner" minOccurs="0">',
            "_em_image_scans.used_frames_per_image": '<xs:element name="frames_per_image" type="xs:token" minOccurs="0"/>',
            "_em_image_scans.sampling_interval": '<xs:element name="sampling_interval" minOccurs="0">',
            "_em_imaging.id": '<xs:attribute name="microscopy_id" type="xs:positiveInteger" use="required"/>',
            "_em_imaging.microscope_model": '<xs:element name="microscope">',
            "_em_imaging.illumination_mode": '<xs:element name="illumination_mode">',
            "_em_imaging.mode": '<xs:element name="imaging_mode">',
            "_em_imaging.electron_source": '<xs:element name="electron_source">',
            "_em_imaging.nominal_magnification": '<xs:element name="nominal_magnification" type="allowed_magnification" minOccurs="0"/>',
            "_em_imaging.calibrated_magnification": '<xs:element name="calibrated_magnification" type="allowed_magnification" minOccurs="0"/>',
            "_em_imaging.specimen_holder_model": '<xs:element name="specimen_holder_model" minOccurs="0">',
            "_em_imaging.cryogen": '<xs:element name="cooling_holder_cryogen" minOccurs="0">',
            "_em_imaging.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_imaging.accelerating_voltage": '<xs:element name="acceleration_voltage">',
            "_em_imaging.c2_aperture_diameter": '<xs:element name="c2_aperture_diameter" minOccurs="0">',
            "_em_imaging.nominal_cs": '<xs:element name="nominal_cs" minOccurs="0">',
            "_em_imaging.nominal_defocus_min": '<xs:element name="nominal_defocus_min" minOccurs="0">',
            "_em_imaging.calibrated_defocus_min": '<xs:element name="calibrated_defocus_min" minOccurs="0">',
            "_em_imaging.nominal_defocus_max": '<xs:element name="nominal_defocus_max" minOccurs="0">',
            "_em_imaging.calibrated_defocus_max": '<xs:element name="calibrated_defocus_max" minOccurs="0">',
            "_em_imaging.recording_temperature_maximum": '<xs:element name="temperature" minOccurs="0">',
            "_em_imaging.recording_temperature_minimum": '<xs:element name="temperature" minOccurs="0">',
            "_em_sample_support.grid_type": '<xs:element name="model" type="xs:token" minOccurs="0">',
            "_em_sample_support.grid_material": '<xs:element name="material" minOccurs="0">',
            "_em_sample_support.grid_mesh_size": '<xs:element name="mesh" type="xs:positiveInteger" minOccurs="0">',
            "_em_sample_support.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_single_particle_entity.point_symmetry": '<xs:element name="point_group">',
            "_em_virus_entity.empty": '<xs:element name="virus_empty" type="xs:boolean"/>',
            "_em_virus_entity.enveloped": '<xs:element name="virus_enveloped" type="xs:boolean"/>',
            "_em_virus_entity.virus_isolate": '<xs:element name="virus_isolate">',
            "_em_virus_entity.virus_type": '<xs:element name="virus_type">',
            "_em_vitrification.instrument": '<xs:element name="instrument" minOccurs="0">',
            "_em_vitrification.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_vitrification.humidity": '<xs:element name="chamber_humidity" minOccurs="0">',
            "_em_vitrification.chamber_temperature": '<xs:element name="chamber_temperature" minOccurs="0">',
            "_em_vitrification.cryogen_name": '<xs:element name="cryogen_name">',
            "_em_2d_crystal_entity.length_a": '<xs:element name="a" type="cell_type"/>',
            "_em_2d_crystal_entity.length_b": '<xs:element name="b" type="cell_type"/>',
            "_em_2d_crystal_entity.length_c": '<xs:element name="c" type="cell_type"/>',
            "_em_2d_crystal_entity.c_sampling_length": '<xs:element name="c_sampling_length" type="cell_type" minOccurs="0"/>',
            "_em_2d_crystal_entity.angle_gamma": '<xs:element name="gamma" type="cell_angle_type"/>',
            "_em_2d_crystal_entity.alpha": '<xs:element name="alpha" type="cell_angle_type minOccurs="0"/>',
            "_em_2d_crystal_entity.beta": '<xs:element name="beta" type="cell_angle_type" minOccurs="0"/>',
            "_em_2d_crystal_entity.space_group_name_H-M": '<xs:element name="space_group" type="xs:token">',
            "_em_3d_crystal_entity.length_a": '<xs:element name="a" type="cell_type"/>',
            "_em_3d_crystal_entity.length_b": '<xs:element name="b" type="cell_type"/>',
            "_em_3d_crystal_entity.length_c": '<xs:element name="c" type="cell_type"/>',
            "_em_3d_crystal_entity.c_sampling_length": '<xs:element name="c_sampling_length" type="cell_type" minOccurs="0"/>',
            "_em_3d_crystal_entity.angle_gamma": '<xs:element name="gamma" type="cell_angle_type"/>',
            "_em_3d_crystal_entity.angle_alpha": '<xs:element name="alpha" type="cell_angle_type" minOccurs="0"/>',
            "_em_3d_crystal_entity.angle_beta": '<xs:element name="beta" type="cell_angle_type" minOccurs="0"/>',
            "_em_3d_crystal_entity.space_group_name": '<xs:element name="space_group" type="xs:token">',
            "_em_buffer_component.concentration": '<xs:element name="concentration" minOccurs="0">',
            "_em_buffer_component.concentration_units": '<xs:element name="concentration" minOccurs="0">',
            "_em_buffer_component.formula": '<xs:element name="formula" type="formula_type" minOccurs="0"/>',
            "_em_buffer_component.name": '<xs:element name="name" type="xs:token" minOccurs="0">',
            "_em_crystal_formation.time_unit": '<xs:complexType name="crystal_formation_time_type">',
            "_em_crystal_formation.lipid_protein_ratio": '<xs:element name="lipid_protein_ratio" type="xs:float" minOccurs="0"/>',
            "_em_crystal_formation.lipid_mixture": '<xs:element name="lipid_mixture" type="xs:token" minOccurs="0"/>',
            "_em_crystal_formation.instrument": '<xs:element name="instrument" type="xs:token" minOccurs="0"/>',
            "_em_crystal_formation.atmosphere": '<xs:element name="atmosphere" type="xs:token" minOccurs="0"/>',
            "_em_crystal_formation.temperature": '<xs:element name="temperature" type="crystal_formation_temperature_type" minOccurs="0"/>',
            "_em_crystal_formation.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_crystal_formation.time": '<xs:complexType name="crystal_formation_time_type">',
            "_em_ctf_correction.phase_reversal_anisotropic": '<xs:element name="anisotropic" type="xs:boolean" minOccurs="0"/>',
            "_em_ctf_correction.phase_reversal_correction_space": '<xs:element name="correction_space" type="correction_space_type" minOccurs="0"/>',
            "_em_ctf_correction.amplitude_correction_factor": '<xs:element name="factor" type="xs:float" minOccurs="0"/>',
            "_em_ctf_correction.amplitude_correction_space": '<xs:element name="correction_space" type="correction_space_type" minOccurs="0"/>',
            "_em_ctf_correction.correction_operation": '<xs:element name="correction_operation" minOccurs="0">',
            "_em_ctf_correction.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_image_recording.id": '<xs:attribute name="image_recording_id" type="xs:positiveInteger"/>',
            "_em_image_recording.detector_mode": '<xs:element name="detector_mode" minOccurs="0">',
            "_em_image_recording.num_grids_imaged": '<xs:element name="number_grids_imaged" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_image_recording.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_image_recording.num_diffraction_images": '<xs:element name="number_diffraction_images" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_image_recording.num_real_images": '<xs:element name="number_real_images" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_image_recording.film_or_detector_model": '<xs:element name="film_or_detector_model">',
            "_em_image_recording.average_exposure_time": '<xs:element name="average_exposure_time" minOccurs="0">',
            "_em_image_recording.avg_electron_dose_per_image": '<xs:element name="average_electron_dose_per_image" minOccurs="0">',
            "_em_entity_assembly_molwt.value": '<xs:complexType name="molecular_weight_type">',
            "_em_entity_assembly_naturalsource.organ": '<xs:element name="organ" type="xs:token" minOccurs="0"/>',
            "_em_entity_assembly_naturalsource.tissue": '<xs:element name="tissue" type="xs:token" minOccurs="0">',
            "_em_entity_assembly_naturalsource.cell": '<xs:element name="cell" type="xs:token" minOccurs="0">',
            "_em_entity_assembly_naturalsource.organelle": '<xs:element name="organelle" type="xs:token" minOccurs="0">',
            "_em_entity_assembly_naturalsource.cellular_location": '<xs:element name="cellular_location" type="xs:token" minOccurs="0">',
            "_em_entity_assembly_naturalsource.strain": '<xs:element name="sci_species_strain" type="xs:string" minOccurs="0" maxOccurs="1"/>',
            "_em_entity_assembly_naturalsource.ncbi_tax_id": '<xs:attribute name="database">',
            "_em_entity_assembly_naturalsource.organism": '<xs:element name="organism" type="organism_type"/>',
            "_em_entity_assembly_naturalsource.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_entity_assembly_synthetic.organ": '<xs:element name="organ" type="xs:token" minOccurs="0"/>',
            "_em_entity_assembly_synthetic.tissue": '<xs:element name="tissue" type="xs:token" minOccurs="0">',
            "_em_entity_assembly_synthetic.cell": '<xs:element name="cell" type="xs:token" minOccurs="0">',
            "_em_entity_assembly_synthetic.organelle": '<xs:element name="organelle" type="xs:token" minOccurs="0">',
            "_em_entity_assembly_synthetic.cellular_location": '<xs:element name="cellular_location" type="xs:token" minOccurs="0">',
            "_em_entity_assembly_synthetic.strain": '<xs:element name="sci_species_strain" type="xs:string" minOccurs="0" maxOccurs="1"/>',
            "_em_entity_assembly_synthetic.ncbi_tax_id": '<xs:attribute name="database">',
            "_em_entity_assembly_synthetic.organism": '<xs:element name="organism" type="organism_type"/>',
            "_em_image_processing.details": '<xs:element name="details" type="xs:token" minOccurs="0"/>',
            "_em_image_processing.image_recording_id": '<xs:element name="image_recording_id" type="xs:positiveInteger"/>',
            "_em_image_processing.id": '<xs:attribute name="image_processing_id" type="xs:positiveInteger" use="required"/>',
            "_em_imaging_optics.energyfilter_slit_width": '<xs:element name="slit_width" minOccurs="0">',
            "_em_imaging_optics.energyfilter_lower": '<xs:element name="lower_energy_threshold" minOccurs="0">',
            "_em_imaging_optics.energyfilter_upper": '<xs:element name="upper_energy_threshold" minOccurs="0">',
            "_em_imaging_optics.phase_plate": '<xs:element name="phase_plate" type="xs:token" minOccurs="0"/>',
            "_em_imaging_optics.sph_aberration_corrector": '<xs:element name="sph_aberration_corrector" type="xs:token" minOccurs="0"/>',
            "_em_imaging_optics.chr_aberration_corrector": '<xs:element name="chr_aberration_corrector" type="xs:token" minOccurs="0"/>',
            "_em_imaging_optics.energyfilter_name": '<xs:element name="name" type="xs:token" minOccurs="0">',
            "_em_imaging_optics.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_map.format": '<xs:attribute name="format" fixed="CCP4" use="required"/>',
            "_em_map.size_kb": '<xs:attribute name="size_kbytes" type="xs:positiveInteger" use="required"/>',
            "_em_map.data_type": '<xs:element name="data_type" type="map_data_type"/>',
            "_em_map.label": '<xs:element name="label" type="xs:token" minOccurs="0"/>',
            "_em_map.annotation_details": '<xs:element name="annotation_details" type="xs:string" minOccurs="0"/>',
            "_em_map.contour_level": '<xs:element name="level" type="xs:float" minOccurs="0">',
            "_em_map.contour_level_source": '<xs:element name="source" minOccurs="0">',
            "_em_map.symmetry_space_group": '<xs:complexType name="applied_symmetry_type">',
            "_em_map.dimensions_col": '<xs:element name="col" type="xs:positiveInteger"/>',
            "_em_map.dimensions_row": '<xs:element name="row" type="xs:positiveInteger"/>',
            "_em_map.dimensions_sec": '<xs:element name="sec" type="xs:positiveInteger"/>',
            "_em_map.origin_col": '<xs:element name="col" type="xs:integer"/>',
            "_em_map.origin_row": '<xs:element name="row" type="xs:integer"/>',
            "_em_map.origin_sec": '<xs:element name="sec" type="xs:integer"/>',
            "_em_map.spacing_x": '<xs:element name="x" type="xs:positiveInteger"/>',
            "_em_map.spacing_y": '<xs:element name="y" type="xs:nonNegativeInteger"/>',
            "_em_map.spacing_z": '<xs:element name="z" type="xs:nonNegativeInteger"/>',
            "_em_map.cell_a": '<xs:element name="a" type="cell_type"/>',
            "_em_map.cell_b": '<xs:element name="b" type="cell_type"/>',
            "_em_map.cell_c": '<xs:element name="c" type="cell_type"/>',
            "_em_map.cell_alpha": '<xs:element name="alpha" type="cell_angle_type" minOccurs="0"/>',
            "_em_map.cell_beta": '<xs:element name="beta" type="cell_angle_type" minOccurs="0"/>',
            "_em_map.cell_gamma": '<xs:element name="gamma" type="cell_angle_type"/>',
            "_em_map.file": '<xs:element name="file">',
            "_em_map.axis_order_fast": '<xs:element name="fast">',
            "_em_map.axis_order_medium": '<xs:element name="medium">',
            "_em_map.axis_order_slow": '<xs:element name="slow">',
            "_em_map.statistics_minimum": '<xs:element name="minimum" type="xs:float"/>',
            "_em_map.statistics_maximum": '<xs:element name="maximum" type="xs:float"/>',
            "_em_map.statistics_average": '<xs:element name="average" type="xs:float"/>',
            "_em_map.statistics_std": '<xs:element name="std" type="xs:float"/>',
            "_em_map.pixel_spacing_x": '<xs:element name="x" type="pixel_spacing_type"/>',
            "_em_map.pixel_spacing_y": '<xs:element name="y" type="pixel_spacing_type"/>',
            "_em_map.pixel_spacing_z": '<xs:element name="z" type="pixel_spacing_type"/>',
            "_em_particle_selection.reference_model": '<xs:element name="reference_model" type="xs:token" minOccurs="0"/>',
            "_em_particle_selection.method": '<xs:element name="method" type="xs:string" minOccurs="0"/>',
            "_em_particle_selection.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_particle_selection.num_particles_selected": '<xs:element name="number_selected" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_software.name": '<xs:element name="name" type="xs:token" minOccurs="0"/>',
            "_em_software.version": '<xs:element name="version" type="xs:token" minOccurs="0"/>',
            "_em_software.details": '<xs:element name="processing_details" type="xs:string" minOccurs="0"/>',
            "_em_specimen.concentration": '<xs:element name="concentration" minOccurs="0">',
            "_em_specimen.id": '<xs:element name="specimen_preparation_id" type="xs:positiveInteger"/>',
            "_em_specimen.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_staining.type": '<xs:element name="type">',
            "_em_staining.material": '<xs:element name="material" type="xs:token">',
            "_em_staining.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_support_film.thickness": '<xs:element name="film_thickness" minOccurs="0">',
            "_em_support_film.id": '<xs:attribute name="film_type_id" type="xs:positiveInteger" use="required"/>',
            "_em_support_film.material": '<xs:element name="film_material" type="xs:token" minOccurs="0">',
            "_em_support_film.topology": '<xs:element name="film_topology" minOccurs="0">',
            "_em_virus_natural_host.ncbi_tax_id": '<xs:attribute name="database">',
            "_em_virus_natural_host.organism": '<xs:element name="organism" type="organism_type">',
            "_em_virus_natural_host.strain": '<xs:element name="strain" type="xs:token" minOccurs="0"/>',
            "_em_virus_shell.triangulation": '<xs:element name="triangulation" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_virus_shell.name": '<xs:element name="name" type="xs:token" nillable="false" minOccurs="0"/>',
            "_em_virus_shell.id": '<xs:attribute name="shell_id" type="xs:positiveInteger"/>',
            "_em_virus_shell.diameter": '<xs:element name="diameter" minOccurs="0">',
            "_em_admin.current_status": '<xs:element name="code" type="code_type"/>',
            "_em_admin.deposition_site": '<xs:element name="deposition">',
            "_em_admin.deposition_date": '<xs:element name="deposition" type="xs:date">',
            "_em_admin.header_release_date": '<xs:element name="header_release" type="xs:date" minOccurs="0">',
            "_em_admin.map_release_date": '<xs:element name="map_release" type="xs:date" minOccurs="0">',
            "_em_admin.obsoleted_date": '<xs:element name="obsolete" type="xs:date" minOccurs="0">',
            "_em_admin.last_update": '<xs:element name="update" type="xs:date">',  # noqa: F601 pylint: disable=duplicate-key
            "_em_admin.title": '<xs:element name="title" type="xs:token">',
            "_em_admin.details": '<xs:element name="details" type="xs:token" minOccurs="0">',
            "_em_admin.keywords": '<xs:element name="keywords" type="xs:token" minOccurs="0">',
            "_em_admin.composite_map": '<xs:element name="composite_map" type="xs:boolean">',
            "_em_euler_angle_assignment.type": '<xs:element name="type">',
            "_em_euler_angle_assignment.proj_matching_num_projections": '<xs:element name="number_reference_projections" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_euler_angle_assignment.proj_matching_merit_function": '<xs:element name="merit_function" type="xs:token" minOccurs="0"/>',
            "_em_euler_angle_assignment.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_euler_angle_assignment.proj_matching_angular_sampling": '<xs:element name="angular_sampling" minOccurs="0">',
            "_em_author_list.identifier_ORCID": '<xs:attribute name="ORCID" type="ORCID_type"/>',
            "_em_db_reference.access_code": '<xs:element name="emdb_id" type="emdb_id_type"/>',
            "_em_db_reference.relationship": '<xs:element name="relationship" minOccurs="0">',
            "_em_db_reference.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_db_reference_auxiliary.link_type": '<xs:element name="type">',
            "_em_db_reference_auxiliary.link": '<xs:element name="link">',
            "_em_diffraction_shell.id": '<xs:attribute name="shell_id" type="xs:positiveInteger"/>',
            "_em_diffraction_shell.high_resolution": '<xs:element name="high_resolution">',
            "_em_diffraction_shell.low_resolution": '<xs:element name="low_resolution">',
            "_em_diffraction_shell.num_structure_factors": '<xs:element name="number_structure_factors" type="xs:positiveInteger"/>',
            "_em_diffraction_shell.phase_residual": '<xs:element name="phase_residual" type="xs:float"/>',
            "_em_diffraction_shell.fourier_space_coverage": '<xs:element name="fourier_space_coverage" type="xs:float">',
            "_em_diffraction_shell.multiplicity": '<xs:element name="multiplicity" type="xs:float"/>',
            "_em_diffraction_stats.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_diffraction_stats.num_intensities_measured": '<xs:element name="number_intensities_measured" type="xs:positiveInteger"/>',
            "_em_diffraction_stats.num_structure_factors": '<xs:element name="number_structure_factors" type="xs:positiveInteger"/>',
            "_em_diffraction_stats.fourier_space_coverage": '<xs:element name="fourier_space_coverage" type="xs:float"/>',
            "_em_diffraction_stats.r_sym": '<xs:element name="r_sym" type="xs:float" minOccurs="0"/>',
            "_em_diffraction_stats.r_merge": '<xs:element name="r_merge" type="xs:float"/>',
            "_em_diffraction_stats.overall_phase_error": '<xs:element name="overall_phase_error" type="xs:float" minOccurs="0"/>',
            "_em_diffraction_stats.overall_phase_residual": '<xs:element name="overall_phase_residual" type="xs:float"  minOccurs="0"/>',
            "_em_diffraction_stats.phase_error_rejection_criteria": '<xs:element name="phase_error_rejection_criteria" type="xs:token"/>',
            "_em_diffraction_stats.high_resolution": '<xs:element name="high_resolution">',
            "_em_embedding.material": '<xs:element name="material" type="xs:token">',
            "_em_embedding.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_fiducial_markers.manufacturer": '<xs:element name="manufacturer" type="xs:token" minOccurs="0">',
            "_em_fiducial_markers.diameter": '<xs:element name="diameter" type="fiducial_marker_diameter_type"/>',
            "_em_final_classification.avg_num_images_per_class": '<xs:element name="average_number_members_per_class" minOccurs="0">',
            "_em_final_classification.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_final_classification.num_classes": '<xs:element name="number_classes" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_fsc_curve.file": '<xs:element name="file">',
            "_em_fsc_curve.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_grid_pretreatment.type": '<xs:element name="type" type="xs:token" minOccurs="0"/>',
            "_em_grid_pretreatment.atmosphere": '<xs:element name="atmosphere" minOccurs="0">',
            "_em_grid_pretreatment.time": '<xs:element name="time" minOccurs="0">',
            "_em_grid_pretreatment.pressure": '<xs:element name="pressure" minOccurs="0">',
            "_em_high_pressure_freezing.instrument": '<xs:element name="instrument" minOccurs="0">',
            "_em_high_pressure_freezing.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_diffraction.camera_length": '<xs:element name="camera_length">',
            "_em_tomography.axis1_min_angle": '<xs:element name="min_angle" minOccurs="0">',
            "_em_tomography.axis2_min_angle": '<xs:element name="min_angle" minOccurs="0">',
            "_em_tomography.axis1_max_angle": '<xs:element name="max_angle" minOccurs="0">',
            "_em_tomography.axis2_max_angle": '<xs:element name="max_angle" minOccurs="0">',
            "_em_tomography.axis1_angle_increment": '<xs:element name="angle_increment" minOccurs="0">',
            "_em_tomography.axis2_angle_increment": '<xs:element name="angle_increment" minOccurs="0">',
            "_em_tomography.dual_tilt_axis_rotation": '<xs:element name="axis_rotation" fixed="90" minOccurs="0">',
            "_em_obsolete.date": '<xs:element name="date" type="xs:date"/>',
            "_em_obsolete.entry": '<xs:element name="entry" type="emdb_id_type"/>',
            "_em_obsolete.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_focused_ion_beam.instrument": '<xs:element name="instrument">',
            "_em_focused_ion_beam.ion": '<xs:element name="ion">',
            "_em_focused_ion_beam.voltage": '<xs:element name="voltage" type="fib_voltage_type"/>',
            "_em_focused_ion_beam.current": '<xs:element name="current" type="fib_current_type"/>',
            "_em_focused_ion_beam.dose_rate": '<xs:element name="dose_rate" type="fib_dose_rate_type" minOccurs="0"/>',
            "_em_focused_ion_beam.duration": '<xs:element name="duration" type="fib_duration_type"/>',
            "_em_focused_ion_beam.temperature": '<xs:element name="temperature" type="temperature_type"/>',
            "_em_focused_ion_beam.initial_thickness": '<xs:element name="initial_thickness" type="fib_initial_thickness_type">',
            "_em_focused_ion_beam.final_thickness": '<xs:element name="final_thickness" type="fib_final_thickness_type"/>',
            "_em_focused_ion_beam.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_ultramicrotomy.instrument": '<xs:element name="instrument" type="xs:token"/>',
            "_em_ultramicrotomy.temperature": '<xs:element name="temperature" type="temperature_type"/>',
            "_em_ultramicrotomy.final_thickness": '<xs:element name="final_thickness" type="ultramicrotomy_final_thickness_type"/>',
            "_em_ultramicrotomy.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_shadowing.details": '<xs:element name="details" type="xs:string" minOccurs="0">',
            "_em_shadowing.material": '<xs:element name="material" type="xs:token">',
            "_em_shadowing.angle": '<xs:element name="angle">',
            "_em_shadowing.thickness": '<xs:element name="thickness">',
            "_em_start_model.type": '<xs:attribute name="type_of_model" type="xs:token"/>',
            "_em_start_model.random_conical_tilt_num_images": '<xs:element name="number_images" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_start_model.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_start_model.orthogonal_tilt_num_images": '<xs:element name="number_images" type="xs:positiveInteger" minOccurs="0"/>',
            "_em_start_model.emdb_id": '<xs:element name="emdb_id" type="emdb_id_type" minOccurs="0"/>',
            "_em_start_model.pdb_id": '<xs:element name="pdb_id" type="pdb_code_type"/>',
            "_em_start_model.insilico_model": '<xs:element name="insilico_model" type="xs:token" minOccurs="0"/>',
            "_em_start_model.orthogonal_tilt_angle1": '<xs:element name="tilt_angle1" minOccurs="0">',
            "_em_start_model.orthogonal_tilt_angle2": '<xs:element name="tilt_angle2" minOccurs="0">',
            "_em_start_model.random_conical_tilt_angle": '<xs:element name="tilt_angle" minOccurs="0">',
            "_em_supersede.date": '<xs:element name="date" type="xs:date"/>',
            "_em_supersede.entry": '<xs:element name="entry" type="emdb_id_type"/>',
            "_em_supersede.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_tomography_specimen.cryo_protectant": '<xs:element name="cryo_protectant" type="xs:token" minOccurs="0">',
            "_em_volume_selection.num_tomograms": '<xs:element name="number_tomograms" type="xs:positiveInteger"/>',
            "_em_volume_selection.num_volumes_extracted": '<xs:element name="number_images_used" type="xs:positiveInteger"/>',
            "_em_volume_selection.reference_model": '<xs:element name="reference_model" type="xs:token" minOccurs="0">',
            "_em_volume_selection.method": '<xs:element name="method" type="xs:string" minOccurs="0">',
            "_em_volume_selection.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_pdbx_database_status.process_site": '<xs:element name="processing_site" minOccurs="0">',
            "_struct.title": '<xs:element name="title" type="xs:token">',
            "_struct_keywords.text": '<xs:element name="keywords" type="xs:string" minOccurs="0">',
            "_citation.title": '<xs:element name="title" type="xs:token"/>',
            "_citation.journal_full": '<xs:element name="journal" type="xs:token" minOccurs="0"/>',
            "_citation.journal_abbrev": '<xs:element name="journal_abbreviation" type="xs:token"/>',
            "_citation.country": '<xs:element name="country" type="xs:token" minOccurs="0"/>',
            "_citation.journal_issue": '<xs:element name="issue" type="xs:positiveInteger" minOccurs="0"/>',
            "_citation.journal_volume": '<xs:element name="volume" type="xs:string" nillable="true" minOccurs="0"/>',
            "_citation.page_first": '<xs:element name="first_page" type="page_type" nillable="false" minOccurs="0"/>',
            "_citation.page_last": '<xs:element name="last_page" type="page_type" minOccurs="0"/>',
            "_citation.year": '<xs:element name="year" minOccurs="0">',
            "_citation.language": '<xs:element name="language" type="xs:language" minOccurs="0"/>',
            "_citation.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_citation.book_title": '<xs:element name="title" type="xs:token"/>',
            "_citation.book_publisher": '<xs:element name="publisher" type="xs:token" minOccurs="0"/>',
            "_pdbx_database_related.db_name": '<xs:element name="db_name" type="token"/>',
            "_pdbx_database_related.db_id": '<xs:element name="emdb_id" type="emdb_id_type"/>',
            "_pdbx_database_related.content_type": '<xs:element name="relationship" minOccurs="0">',
            "_pdbx_database_related.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_em_entity_assembly_recombinan.strain": '<xs:element name="recombinant_strain" type="xs:token" minOccurs="0"/>',
            "_em_entity_assembly_recombinan.cell": '<xs:element name="recombinant_cell" type="xs:token" minOccurs="0"/>',
            "_em_entity_assembly_recombinan.plasmid": '<xs:element name="recombinant_plasmid" type="xs:token" minOccurs="0"/>',
            "_entity.id": '<xs:attribute name="macromolecule_id" type="xs:positiveInteger" use="required"/>',
            "_entity.pdbx_number_of_molecules": '<xs:element name="number_of_copies" type="pos_int_or_string_type" minOccurs="0"/>',
            "_entity.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_entity.pdbx_description": '<xs:element name="name" type="sci_name_type"/>',
            "_entity_src_nat.pdbx_organ": '<xs:element name="organ" type="xs:token" minOccurs="0"/>',
            "_entity_src_gen.pdbx_gene_src_organ": '<xs:element name="organ" type="xs:token" minOccurs="0"/>',
            "_entity_src_nat.tissue": '<xs:element name="tissue" type="xs:token" minOccurs="0"/>',
            "_entity_src_gen.gene_src_tissue": '<xs:element name="tissue" type="xs:token" minOccurs="0"/>',
            "_entity_src_nat.pdbx_cell": '<xs:element name="cell" type="xs:token" minOccurs="0"/>',
            "_entity_src_gen.pdbx_gene_src_cell": '<xs:element name="cell" type="xs:token" minOccurs="0"/>',
            "_entity_src_nat.pdbx_organelle": '<xs:element name="organelle" type="xs:token" minOccurs="0"/>',
            "_entity_src_gen.pdbx_gene_src_organelle": '<xs:element name="organelle" type="xs:token" minOccurs="0"/>',
            "_entity_src_nat.pdbx_cellular_location": '<xs:element name="cellular_location" type="xs:token" minOccurs="0"/>',
            "_entity_src_gen.pdbx_gene_src_cellular_location": '<xs:element name="cellular_location" type="xs:token" minOccurs="0"/>',
            "_entity.pdbx_ec": '<xs:element name="ec_number" minOccurs="0" maxOccurs="unbounded">',
            "_entity_src_gen.plasmid_name": '<xs:element name="recombinant_plasmid" type="xs:token" minOccurs="0"/>',
            "_entity_src_gen.pdbx_gene_src_cell_line": '<xs:element name="cell" type="xs:token" minOccurs="0"/>',
            "_pdbx_entity_nonpoly.comp_id": '<xs:element name="formula" type="formula_type" minOccurs="0"/>',
            "_em_helical_entity.angular_rotation_per_subunit": '<xs:element name="delta_phi">',
            "_em_helical_entity.axial_rise_per_subunit": '<xs:element name="delta_z">',
            "_em_helical_entity.axial_symmetry": '<xs:element name="axial_symmetry">',
            "_database_2.database_code": '<xs:attribute name="emdb_id" type="emdb_id_type" use="required"/>',
            "_citation.unpublished_flag": '<xs:attribute name="published" type="xs:boolean" use="required"/>',
            "_citation.pdbx_database_id_PubMed": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_citation.pdbx_database_id_DOI": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_citation.book_id_ISBN": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_citation.journal_id_ISSN": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_citation.abstract_id_CAS": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_citation.journal_id_CSD": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_citation.database_id_Medline": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_citation.journal_id_ASTM": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_em_entity_assembly_recombinan.ncbi_tax_id": '<xs:attribute name="database">',
            "_em_entity_assembly_recombinan.organism": '<xs:element name="recombinant_organism" type="organism_type">',
            "_pdbx_database_status.pdbx_annotator": '<xs:element name="annotator" minOccurs="0">',
            "_pdbx_database_PDB_obs_spr.date": '<xs:element name="date" type="xs:date"/>',
            "_pdbx_database_PDB_obs_spr.replace_pdb_id": '<xs:element name="entry" type="emdb_id_type"/>',
            "_pdbx_database_PDB_obs_spr.details": '<xs:element name="details" type="xs:string" minOccurs="0"/>',
            "_pdbx_database_PDB_obs_spr.pdb_id": '<xs:element name="entry" type="emdb_id_type"/>',
            "_pdbx_audit_support.funding_organization": '<xs:element name="funding_body" type="xs:token"/>',
            "_pdbx_audit_support.grant_number": '<xs:element name="code" type="xs:token minOccurs="0"/>',
            "_pdbx_audit_support.country": '<xs:element name="country" type="xs:token" minOccurs="0"/>',
            "_pdbx_contact_author.role": '<xs:element name="role">',
            "_pdbx_contact_author.name_salutation": '<xs:element name="title">',
            "_pdbx_contact_author.name_first": '<xs:element name="first_name" type="xs:token">',
            "_pdbx_contact_author.name_mi": '<xs:element name="middle_name">',
            "_pdbx_contact_author.name_last": '<xs:element name="last_name" type="xs:token">',
            "_pdbx_contact_author.organization_type": '<xs:element name="organization">',
            "_pdbx_contact_author.address_1": '<xs:element name="street" type="xs:string"/>',
            "_pdbx_contact_author.city": '<xs:element name="town_or_city" type="xs:token"/>',
            "_pdbx_contact_author.state_province": '<xs:element name="state_or_province" type="xs:token"/>',
            "_pdbx_contact_author.country": '<xs:element name="country" type="xs:token"/>',
            "_pdbx_contact_author.postal_code": '<xs:element name="post_or_zip_code" type="xs:token"/>',
            "_pdbx_contact_author.email": '<xs:element name="email">',
            "_pdbx_contact_author.phone": '<xs:element name="telephone" type="telephone_number_type"/>',
            "_pdbx_contact_author.fax": '<xs:element name="fax" type="telephone_number_type"/>',
            "_entity_src_nat.pdbx_ncbi_taxonomy_id": '<xs:attribute name="database">',
            "_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id": '<xs:attribute name="database">',
            "_pdbx_entity_src_syn.ncbi_taxonomy_id": '<xs:attribute name="database">',
            "_entity_poly.pdbx_seq_one_letter_code": '<xs:element name="string" type="xs:token" minOccurs="0">',
            "_struct_ref.db_name": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_struct_ref.db_code": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_struct_ref.pdbx_db_accession": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_pdbx_struct_ref_seq_depositor_info.db_name": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_pdbx_struct_ref_seq_depositor_info.db_accession": '<xs:element name="external_references" minOccurs="0" maxOccurs="unbounded">',
            "_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id": '<xs:attribute name="database">',
            "_entity_src_gen.pdbx_host_org_scientific_name": '<xs:element name="organism" type="organism_type"/>',
            "_entity.formula_weight": '<xs:element name="experimental" minOccurs="0">',
            "_entity_src_nat.common_name": '<xs:element name="organism" type="organism_type"/>',
            "_entity_src_nat.pdbx_organism_scientific": '<xs:element name="organism" type="organism_type"/>',
            "_entity_src_gen.pdbx_gene_src_scientific_name": '<xs:element name="organism" type="organism_type"/>',
            "_pdbx_entity_src_syn.organism_scientific": '<xs:element name="organism" type="organism_type"/>',
            "_entity_src_nat.strain": '<xs:element name="strain" type="xs:token" minOccurs="0"/>',
            "_entity_src_gen.gene_src_strain": '<xs:element name="strain" type="xs:token" minOccurs="0"/>',
            "_pdbx_entity_src_syn.strain": '<xs:element name="strain" type="xs:token" minOccurs="0"/>',
            "_exptl.method": '<xs:element name="method">',
            "_audit_author.identifier_ORCID": '<xs:attribute name="ORCID" type="ORCID_type"/>',
            "_pdbx_audit_revision_history.revision_date": '<xs:attribute name="date" type="xs:date" use="required"/>',
            "_pdbx_audit_revision_details and pdbx_audit_revision_history": '<xs:element name="change_list" maxOccurs="1">',
            "_pdbx_audit_revision_details.type": '<xs:element name="revision_type">',
            "_pdbx_audit_revision_details.provider": '<xs:element name="provider">',
            "_pdbx_audit_revision_details.description": '<xs:element name="description" type="xs:token" minOccurs="0"/>',
            "_pdbx_audit_revision_details.details": '<xs:element name="details" type="xs:token" minOccurs="0"/>',
            "_pdbx_audit_revision_group.group": '<xs:element name="revision_group" minOccurs="0">',
            "_pdbx_audit_revision_history.internal_part_number": '<xs:attribute name="part" type="xs:positiveInteger"/>',
            "_pdbx_audit_revision_history.data_content_type": '<xs:attribute name="revision_type" type="xs:token" use="required"/>',
            "_pdbx_audit_revision_history.part_number": '<xs:attribute name="part" type="xs:positiveInteger"/>',
            "_pdbx_audit_revision_history.type": '<xs:attribute name="revision_action" type="xs:token" use="required"/>',
            "_pdbx_audit_revision_category.category": '<xs:element name="category" type="revision_category_or_item_type" minOccurs="1" maxOccurs="unbounded"/>',
            "_pdbx_audit_revision_item.item": '<xs:element name="item" type="revision_category_or_item_type" minOccurs="1" maxOccurs="unbounded"/>'
        }

    class ALog(object):
        """Class containing one log information"""

        def __init__(self, cif_item=None, setter_func=None, xsd=None, em_for_emd=None, fmt_cif_value=None, parent_el_req=None, soft_name=None, log_text=None):
            self.cif_item = cif_item
            self.setter_func = setter_func
            self.schema_entity = xsd
            self.em_for_emd = em_for_emd
            self.fmt_cif_value = fmt_cif_value
            self.parent_el_req = parent_el_req
            self.soft_name = soft_name
            self.log_text = log_text

        def create_std_info_log_text(self, title_str):
            info1 = title_str + "The value (%s) is given to (%s)." % (self.fmt_cif_value, self.setter_func)
            info2 = "The value came from (%s) in cif." % self.cif_item
            info3 = "The value will be written for schema component (%s)" % self.schema_entity
            info4 = None
            if self.parent_el_req:
                info4 = "The parent element IS required."
            info5 = None
            if self.em_for_emd is not None:
                info5 = "The _em category for the above category is (%s)" % self.em_for_emd
            info6 = None
            if self.soft_name is not None:
                info6 = "Software name is (%s)" % self.soft_name

            self.log_text = ". ".join(filter(None, (info1, info2, info3, info4, info5, info6))) + "."

        def create_problematic_log_text(self, title_str, direct_text=None):
            if direct_text is None:
                problem1 = title_str + "Function (%s) failure" % self.setter_func
                problem2 = "The value for cif category (%s) not given" % self.cif_item
                problem3 = "The value is for schema component (%s)" % self.schema_entity
                problem4 = None
                if self.parent_el_req:
                    problem4 = "The parent element IS required."
                problem5 = None
                if self.em_for_emd is not None:
                    problem5 = "The _em category for the above category is (%s)" % self.em_for_emd
                problem6 = None
                if self.soft_name is not None:
                    problem6 = "Software name is (%s)" % self.soft_name

                self.log_text = ". ".join(filter(None, (problem1, problem2, problem3, problem4, problem5, problem6))) + "."
            else:
                self.log_text = title_str + direct_text

    class EntryLog(object):
        """Class containing translation log for an entry"""

        _info_title = "INFO: "
        _warn_title = "POTENTIAL PROBLEM: "
        _error_title = "PROBLEM: "
        _change_title = "CHANGE MADE: "
        _not_changed_for_now_title = "NOT CHANGED FOR NOW: "
        _validation_title = "VALIDATION ERROR "

        def __init__(self, entry_ID):
            self.ID = entry_ID
            # logs lists
            self.error_logs = []
            self.warn_logs = []
            self.info_logs = []

        @property
        def is_error_log_empty(self):
            if len(self.error_logs) == 0:
                return True
            return False

        @property
        def id(self):
            return self.ID

        @id.setter
        def id(self, id_value):
            self.ID = id_value

        @property
        def errors(self):
            return self.error_logs

        @property
        def warnings(self):
            return self.warn_logs

        @property
        def infos(self):
            return self.info_logs

        @property
        def info_title(self):
            return self._info_title

        @property
        def warn_title(self):
            return self._warn_title

        @property
        def error_title(self):
            return self._error_title

        @property
        def change_title(self):
            return self._change_title

        @property
        def not_changed_for_now_title(self):
            return self._not_changed_for_now_title

        @property
        def validation_title(self):
            return self._validation_title

        def add_info(self, cif_item, setter_func, xsd, em_for_emd, fmt_cif_value=None, parent_el_req=None, soft_name=None):
            info = CifEMDBTranslator.ALog(cif_item, setter_func, xsd, em_for_emd, fmt_cif_value, parent_el_req, soft_name)
            info.create_std_info_log_text("(" + self.ID + ")" + self._info_title)
            self.info_logs.append(info)

        def add_warn(self, cif_item, setter_func, xsd, em_for_emd, fmt_cif_value=None, parent_el_req=None, soft_name=None):
            warn = CifEMDBTranslator.ALog(cif_item, setter_func, xsd, em_for_emd, fmt_cif_value, parent_el_req, soft_name)
            warn.create_problematic_log_text("(" + self.ID + ")" + self._warn_title)
            self.warn_logs.append(warn)

        def add_err(self, cif_item, setter_func, xsd, em_for_emd, fmt_cif_value=None, parent_el_req=None, soft_name=None):
            err = CifEMDBTranslator.ALog(cif_item, setter_func, xsd, em_for_emd, fmt_cif_value, parent_el_req, soft_name)
            err.create_problematic_log_text("(" + self.ID + ")" + self._error_title)
            self.error_logs.append(err)

    class TranslationLog(object):
        """Container class for translation logs"""

        def __init__(self):
            self.entry_logs = []

        @property
        def logs(self):
            return self.entry_logs

    def __init__(self, info_log=None, warn_log=None, error_log=None):  # pylint: disable=unused-argument
        self.cif_file_name = None
        self.cif_file_read = False  # flag set once the cif file is read
        self.cif = None  # cif dictionary
        self.xml_out = None  # xml object representing conversion from cif
        self.map_cif_dict = {}
        self.emdb_id_u = None
        self.__show_private = False
        self.__show_log_id = False
        self.show_log_on_console = False
        self.info_log_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), self.Constants.INFO_LOG_FILE_NAME)
        self.info_log_string = None
        self.warn_log_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), self.Constants.WARN_LOG_FILE_NAME)
        self.warn_log_string = None
        self.error_log_file_name = os.path.join(os.path.dirname(os.path.abspath(__file__)), self.Constants.ERR_LOG_FILE_NAME)
        self.error_log_string = None
        self.logger = None
        self.parent_logger_level = None
        self.translation_log = self.TranslationLog()
        self.entry_in_translation_log = None
        # create_xml enables the XML out to be created - if False the creation of the output file and its validation shouldn't happen
        self.create_xml = True
        emdb.GenerateDSNamespaceDefs_ = {
            "entry_type": 'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                          'xsi:noNamespaceSchemaLocation="https://ftp.ebi.ac.uk/pub/databases/em_ebi/emdb_related/emdb-schemas/'
                          'emdb_schemas/v3/v{0}/emdb.xsd" '
                          'version="{1}"'.format(self.Constants.XML_VERSION, self.Constants.XML_OUT_VERSION)
        }
        emdb.GenerateDSNamespaceTypePrefixes_ = {}

    @property
    def is_translation_log_empty(self):
        logs = self.translation_log.logs
        if len(logs) == 0:
            return True
        else:
            for a_log in logs:
                if not a_log.is_error_log_empty:
                    return False
            return True

    @property
    def current_entry_log(self):
        return self.entry_in_translation_log

    @property
    def get_translation_log(self):
        return self.translation_log

    def get_show_log_id(self):
        return self.__show_log_id

    def set_show_log_id(self, value):
        self.__show_log_id = value

    def del_show_log_id(self):
        del self.__show_log_id

    _show_log_id = property(get_show_log_id, set_show_log_id, del_show_log_id, "_show_log_id's docstring")

    def get_show_private(self):
        return self.__show_private

    def set_show_private(self, value):
        self.__show_private = value

    def del_show_private(self):
        del self.__show_private

    _show_private = property(get_show_private, set_show_private, del_show_private, "show_private's docstring")

    def __del__(self):
        """
        Use to roll back the level of the translator parent logger
        """
        if self.parent_logger_level is not None:
            logging.getLogger().setLevel(self.parent_logger_level)
        else:
            logging.getLogger().setLevel(logging.INFO)

    def open_logs(self, open_info=True, open_warn=True, open_error=True):
        if open_info:
            self.info_log_string = self.open_log_stream(self.Constants.INFO_LOG_STRING)
        if open_warn:
            self.warn_log_string = self.open_log_stream(self.Constants.WARN_LOG_STRING)
        if open_error:
            self.error_log_string = self.open_log_stream(self.Constants.ERROR_LOG_STRING)

    def close_logs(self, close_info=True, close_warn=True, close_error=True):
        """"""
        if close_info:
            self.info_log_string.close()
        if close_warn:
            self.warn_log_string.close()
        if close_error:
            self.error_log_string.close()

    def write_to_file(self, log_file_name, log_content=None):
        if log_file_name is not None:
            log_file_hdl = open(log_file_name, "a")
            if log_content is not None:
                log_file_hdl.write(log_content)
            log_file_hdl.close()

    def write_a_logger_log(self, log_str, log_file_name):
        """

        :param log_str:
        :param log_file_name:
        """
        log_content = log_str.getvalue()
        if self.show_log_on_console:
            print("%s" % log_content)
        self.write_to_file(log_file_name, log_content)

    def write_logger_logs(self, write_error_log=False, write_warn_log=False, write_info_log=False):
        """

        :param write_error_log: A flag - if true the error log is written out
        :param write_warn_log: A flag - if true the warning log is written out
        :param write_info_log: A flag - if true the info log is written out
        """
        if write_error_log:
            if self.error_log_string is not None and not self.error_log_string.closed:
                # write the error log buffer to self.error_log_file_name
                self.write_a_logger_log(self.error_log_string, self.error_log_file_name)
                self.error_log_string.close()

        if write_warn_log:
            if self.warn_log_string is not None and not self.warn_log_string.closed:
                # write the warning log buffer to self.warn_log_file_name
                self.write_a_logger_log(self.warn_log_string, self.warn_log_file_name)
                self.warn_log_string.close()

        if write_info_log:
            if self.info_log_string is not None and not self.info_log_string.closed:
                # write the info log buffer to self.info_log_file_name
                self.write_a_logger_log(self.info_log_string, self.info_log_file_name)
                self.info_log_string.close()

    def init_logger_log(self, log_file_name_in):
        """"""
        if not os.path.exists(log_file_name_in):
            self.write_to_file(log_file_name_in)

    def initialise_logging(self, log_info=False, info_log_file_name=None, log_warn=False, warn_log_file_name=None, log_error=False, error_log_file_name=None):
        """
        Sets the console logging

        Firstly, the level of the parent logger is set. This value is used in
        the __del__ function in order to roll back the value. The parent logging level
        is then set to critical in order to produce the least information.

        The translator's logger sets up three logging files:
        1. INFO+
        2. WARN+
        3. ERROR+
        These files can be disabled on the command line.
        """
        if log_info:
            if info_log_file_name is not None:
                self.info_log_file_name = info_log_file_name
            self.init_logger_log(self.info_log_file_name)

        if log_warn:
            if warn_log_file_name is not None:
                self.warn_log_file_name = warn_log_file_name
            self.init_logger_log(self.warn_log_file_name)

        if log_error:
            if error_log_file_name is not None:
                self.error_log_file_name = error_log_file_name
            self.init_logger_log(self.error_log_file_name)

        # note the logging level of the parent logger
        self.parent_logger_level = logging.getLogger().getEffectiveLevel()

        # keep only critical information for the translator
        logging.getLogger().setLevel(60)
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        self.set_logger_logs(log_info, log_warn, log_error)

        # prevent logging from being sent to the upper logger and console
        self.logger.propagate = False

    def set_logger_logging(
        self, log_info=False, log_warn=False, log_error=False, show_log_on_console=False, info_log_file_name=None, warn_log_file_name=None, error_log_file_name=None
    ):
        """
        This method switches on any console logging
        """
        self.show_log_on_console = show_log_on_console
        if log_info or log_warn or log_error:
            self.initialise_logging(log_info, info_log_file_name, log_warn, warn_log_file_name, log_error, error_log_file_name)
        else:
            self.set_logger_logs(False, False, False)

    def get_log_level(self, log_string_name):
        """

        :param log_string:
        :return:
        """
        ret_str = ""
        if log_string_name.find(self.Constants.INFO_LOG_STRING) != -1:
            ret_str = self.Constants.INFO_LOG_LEVEL
        elif log_string_name.find(self.Constants.WARN_LOG_STRING) != -1:
            ret_str = self.Constants.WARN_LOG_LEVEL
        elif log_string_name.find(self.Constants.ERROR_LOG_STRING) != -1:
            ret_str = self.Constants.ERROR_LOG_LEVEL
        return ret_str

    def open_log_stream(self, log_string_name):
        """

        :param log_string:
        :return:
        """
        log_string = None
        if self.logger is not None:
            log_string = io.StringIO(initial_value=u"", newline=u"\n")
            log_file_hdl = logging.StreamHandler(log_string)
            log_file_hdl.name = log_string_name
            log_level = self.get_log_level(log_string_name)
            log_file_hdl.setLevel(log_level)
            log_file_hdl.setFormatter(logging.Formatter(u""))
            self.logger.addHandler(log_file_hdl)

        return log_string

    def remove_logger_hdl(self, log_string_name):
        """"""
        log_hdl = self.get_logger_handle(log_string_name)
        if log_hdl is not None:
            self.logger.removeHandler(log_hdl)

    def set_logger_logs(self, log_info, log_warn, log_err):
        """
        This function stops logging into log files depending on flags passed on

        Parameters:
        @params log_info: if True info log is recorded
        @params log_warn: if True warning log is recorded
        @params log_err: if True error log is recorded
        """
        if self.logger is not None:
            if log_info:
                self.info_log_string = self.open_log_stream(self.Constants.INFO_LOG_STRING)
            else:
                self.remove_logger_hdl(self.Constants.INFO_LOG_STRING)
                if self.info_log_string is not None:
                    self.info_log_string.close()
            if log_warn:
                self.warn_log_string = self.open_log_stream(self.Constants.WARN_LOG_STRING)
            else:
                self.remove_logger_hdl(self.Constants.WARN_LOG_STRING)
                if self.warn_log_string is not None:
                    self.warn_log_string.close()
            if log_err:
                self.error_log_string = self.open_log_stream(self.Constants.ERROR_LOG_STRING)
            else:
                self.remove_logger_hdl(self.Constants.ERROR_LOG_STRING)
                if self.error_log_string is not None:
                    self.error_log_string.close()

    def get_logger_handle(self, log_filename_or_string):
        """Returns the log handle for a log file name"""
        ret_log_hdl = None
        for log_hdl in self.logger.handlers[:]:
            if isinstance(log_hdl, logging.FileHandler):
                if log_hdl.baseFilename.find(log_filename_or_string) != -1:
                    ret_log_hdl = log_hdl
            elif isinstance(log_hdl, logging.StreamHandler):
                if log_hdl.get_name().find(log_filename_or_string) != -1:
                    ret_log_hdl = log_hdl
        return ret_log_hdl

    def log_formatted(self, log_str, txt):
        """"""
        if self.logger is not None and log_str is not None:
            if log_str.closed:
                if log_str == self.info_log_string:
                    self.info_log_string = self.open_log_stream(self.Constants.INFO_LOG_STRING)
                    log_str = self.info_log_string
                if log_str == self.warn_log_string:
                    self.warn_log_string = self.open_log_stream(self.Constants.WARN_LOG_STRING)
                    log_str = self.warn_log_string
                if log_str == self.error_log_string:
                    self.error_log_string = self.open_log_stream(self.Constants.ERROR_LOG_STRING)
                    log_str = self.error_log_string
            if log_str is not None:
                log_str.write("\n(" + self.entry_in_translation_log.id + ")" + txt + u"\n")

    def read_emd_map_v2_cif_file(self):
        """
        Read mmCIF dictionary that contains information
        about how the em categories map to the emd categories.
        Please note: the mapping is not one to one for certain categories
        """
        siteId = getSiteId()
        cIA = ConfigInfoAppEm(siteId)
        em_map_file_name = cIA.get_emd_mapping_file_path()
        io_adapt = IoAdapterCore()
        map_cat_list = io_adapt.readFile(em_map_file_name)
        cat = self.Constants.PDBX_DICT_ITEM_MAPPING
        container = map_cat_list[0]
        dc_obj = container.getObj(cat)
        if dc_obj is not None:
            for row in range(len(dc_obj.data)):
                emd_cat = dc_obj.getValueOrDefault(attributeName="item_name_dst", defaultValue="", rowIndex=row)
                if emd_cat != "":
                    em_cat = dc_obj.getValueOrDefault(attributeName="item_name_src", defaultValue="", rowIndex=row)
                    self.map_cif_dict.update({emd_cat: em_cat})

    def read_cif_in_file(self, cif_file_name):
        """
        Read cif file into the cif dictionary object

        Parameters:
        @param cif_file_name: name of cif file
        """
        const = self.Constants
        self.cif_file_name = cif_file_name
        io_util = IoAdapterCore()
        container = io_util.readFile(
            inputFilePath=cif_file_name,
            selectList=[
                const.CITATION,
                const.CITATION_AUTHOR,
                const.DATABASE_2,
                const.EM_ADMIN,
                const.EM_DEPUI,
                const.EM_EULER_ANGLE_ASSIGNMENT,
                const.AUDIT_AUTHOR,
                const.EM_AUTHOR_LIST,
                const.EM_BUFFER,
                const.EM_BUFFER_COMPONENT,
                const.EM_DB_REFERENCE,
                const.EM_DB_REFERENCE_AUXILIARY,
                const.EM_CRYSTAL_FORMATION,
                const.EM_DIFFRACTION_SHELL,
                const.EM_DIFFRACTION_STATS,
                const.EM_CTF_CORRECTION,
                const.EM_EMBEDDING,
                const.EM_FIDUCIAL_MARKERS,
                const.EM_FINAL_CLASSIFICATION,
                const.EM_FINAL_2D_CLASSIFICATION,
                const.EM_3D_RECONSTRUCTION,
                const.EM_FSC_CURVE,
                const.EM_SAMPLE_SUPPORT,
                const.EM_GRID_PRETREATMENT,
                const.EM_HELICAL_ENTITY,
                const.EM_HIGH_PRESSURE_FREEZING,
                const.EM_IMAGE_SCANS,
                const.EM_MAP,
                const.EM_IMAGING,
                const.EM_DIFFRACTION,
                const.EM_TOMOGRAPHY,
                const.EM_ENTITY_ASSEMBLY_MOLWT,
                const.EM_3D_FITTING,
                const.EM_3D_FITTING_LIST,
                const.EM_ENTITY_ASSEMBLY_NATURALSOURCE,
                const.EM_ENTITY_ASSEMBLY_SYNTHETIC,
                const.EM_IMAGE_PROCESSING,
                const.EM_IMAGE_RECORDING,
                const.EM_PARTICLE_SELECTION,
                const.EM_ENTITY_ASSEMBLY_RECOMBINANT,
                const.EM_FOCUSED_ION_BEAM,
                const.EM_ULTRAMICROTOMY,
                const.EM_SHADOWING,
                const.EM_SOFTWARE,
                const.EM_SPECIALIST_OPTICS,
                const.EM_SPECIMEN,
                const.EM_STAINING,
                const.EM_START_MODEL,
                const.EM_EXPERIMENT,
                const.EM_SUPPORT_FILM,
                const.EM_ENTITY_ASSEMBLY,
                const.EM_SINGLE_PARTICLE_ENTITY,
                const.EM_3D_CRYSTAL_ENTITY,
                const.EM_TOMOGRAPHY_SPECIMEN,
                const.EM_2D_CRYSTAL_ENTITY,
                const.EM_VIRUS_ENTITY,
                const.EM_VIRUS_NATURAL_HOST,
                const.EM_VIRUS_SHELL,
                const.EM_VITRIFICATION,
                const.EM_VOLUME_SELECTION,
                const.ENTITY,
                const.ENTITY_POLY,
                const.ENTITY_SRC_GEN,
                const.ENTITY_SRC_NAT,
                const.EXPTL,
                const.PDBX_DATABASE_RELATED,
                const.PDBX_DATABASE_STATUS,
                const.PDBX_ENTITY_NONPOLY,
                const.PDBX_DEPOSITOR_INFO,
                const.PDBX_OBS_SPR,
                const.PDBX_AUDIT_SUPPORT,
                const.PDBX_CONTACT_AUTHOR,
                const.STRUCT,
                const.STRUCT_REF,
                const.STRUCT_KEYWORDS,
                const.PDBX_ENTITY_SRC_SYN,
                const.EM_SUPERSEDE,
                const.EM_OBSOLETE,
                const.PDBX_AUDIT_REVISION_HISTORY,
                const.PDBX_AUDIT_REVISION_DETAILS,
                const.PDBX_AUDIT_REVISION_GROUP,
                const.PDBX_AUDIT_REVISION_CATEGORY,
                const.PDBX_AUDIT_REVISION_ITEM
            ],
        )
        self.cif = Cif(container)
        if container is not None or container == {}:
            self.cif_file_read = True

    def write_xml_out_file(self, xml_out_file_name):
        """
        Write out XML file.
        A translation from cif to XML needs to have taken place.

        Parameters:
        @param xml_out_file_name: name of the xml file
        """
        # self.xml_out is the xml object representing conversion from cif
        if self.xml_out is None or self.xml_out.has__content() is False:
            txt = u"There is no content to write out. No output file will be written."
            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
            self.log_formatted(self.error_log_string, "(" + self.entry_in_translation_log.id + ")" + self.Constants.REQUIRED_ALERT + txt)
            return
        # xml_out_file is the actual xml file that will be written into and saved
        xml_out_file = open(xml_out_file_name, "w") if xml_out_file_name else sys.stdout
        xml_out_file.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        self.xml_out.export(xml_out_file, 0, name_="emd")

        if xml_out_file is not sys.stdout:
            xml_out_file.close()

    def translate_cif_to_xml(self):
        """Translate cif file to EMDB xml 3.0"""

        def is_number(astr):
            try:
                float(astr)
                return True
            except ValueError:
                pass

            try:
                int(astr)
                return True
            except ValueError:
                pass

            try:
                import unicodedata

                unicodedata.numeric(astr)
                return True
            except (TypeError, ValueError):
                pass

            return False

        def cif_bool(cif_bool_value):
            """
            Convert from cif boolean values to Python

            Parameters:
            @param cif_bool_value: input boolean as specified in cif
            @return: Python boolean value
            """
            return bool(cif_bool_value in ["1", "YES"])

        def get_cif_item(cif_key, cif_category):
            """
            Helper function that returns the full key given
            the cif key and category

            @param cif_key: key for a (key,value) pair
            @param cif_category: cif category with the list of (key,value) tuples
            @return full cif key
            """
            if cif_category is not None and cif_key is not None:
                return "_" + cif_category + "." + cif_key
            else:
                return None

        def get_cif_value(cif_key, cif_category, cif_list=None):
            """
            With the assumption that the cif_category is implemented as a
            list of tuples, find value corresponding to key

            Parameters:
            @param cif_key: key for a (key,value) pair
            @param cif_category: cif category with the list of (key,value) tuples
            @param cif_list: If this is 'None' then the cif list in the wrapping function is used
            @return: value or None
            """
            if cif_category is not None and cif_key is not None:
                full_key = "_" + cif_category + "." + cif_key
                category_lists = [cif_list] if cif_list else self.cif[cif_category]
                for category_list in category_lists:
                    cif_value = [item for item in category_list if item[0] == full_key]
                    if cif_value is not None:
                        if len(cif_value) == 1 and len(cif_value[0]) == 2:
                            ret = cif_value[0][1]
                            # Handle unicode returned from mmCIF parsers, as Python 2 unidcode does not support upper()
                            if ret:
                                ret = str(ret)
                            return ret
                        else:
                            return None
            else:
                return None

        def get_xsd_for_cif_item(cif_item):
            """
            Helper function that returns the schema construct for the given cif category

            Parameters:
            @param cif_item: full name of the cif category, e.g. _em_admin.title
            """
            xsd = ""
            if cif_item in const.MMCIF_TO_XSD:  # pylint: disable=possibly-used-before-assignment
                xsd = const.MMCIF_TO_XSD[cif_item]
            else:
                txt = u"CIF item (%s) not found in the MMCIF_TO_XSD dictionary." % cif_item
                self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                self.log_formatted(self.error_log_string, "(" + self.entry_in_translation_log.id + ")" + const.REQUIRED_ALERT + txt)
            return xsd

        # def get_em_from_emd(emd_cif_item):
        #     """
        #     Every _emd has its _em category. This method returns _em category for the given _emd (cif_item)
        #
        #     Parameters:
        #     @param cif_item: full name of the cif category, e.g. _em_admin.title
        #     """
        #     if emd_cif_item.find("_emd") == -1:
        #         return None
        #     if self.map_cif_dict is not None and len(self.map_cif_dict) > 0:
        #         em_cif_item = self.map_cif_dict.get(emd_cif_item, None)
        #         if em_cif_item is None:
        #             return "Not provided"
        #         else:
        #             return em_cif_item
        #             # else:
        #             #     txt = u'Mapping from _emd to _em space does not exist for (%s).' % emd_cif_item
        #             #     self.current_entry_log.warn_logs.append(self.ALog(log_text=self.current_entry_log.warn_title + txt))
        #             #     self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)

        def is_cif_item_required(cif_item):
            """
            Helper function. Establishes if the cif item is required.
            The information is found in the value of the (cif_item, xsd) dictionary.
            The checks look for 'minOccurs="0"' for elements and 'use="required"' for attribs
            """
            if cif_item is not None:
                if cif_item.find("xs:element") != -1:
                    if cif_item.find('minOccurs="0"') == -1:
                        # not found minOccurs="0" - the elment is required
                        return True
                    else:
                        return False
                elif cif_item.find("xs:attribute") != -1:
                    if cif_item.find('use="required"') != -1:
                        # use="required" found so the attribute is required
                        return True
                    else:
                        return False
            else:
                return False

        def log_as_str(log_str, cif_item, setter_func, xsd, em_for_emd, alert_str, fmt_cif_value=None, parent_el_req=None, soft_name=None):
            """
            Formats the logging messages and uses the logging level
            """
            if self.logger is not None:
                if log_str is not None:
                    alert = "(" + self.entry_in_translation_log.id + ")" + alert_str
                    spacing = u" " * len(alert)
                    log_str.write(u"\n")
                    cif_item_txt = u""
                    xsd_txt = u"\n"
                    parent_txt = u"\n"
                    em_for_emd_txt = u"The _em category for the above category is (%s)."
                    soft_txt = u"Software name is (%s)"
                    if log_str == self.info_log_string:
                        if fmt_cif_value is not None:
                            log_str.write(alert + u"The value (%s) is given to (%s).\n" % (fmt_cif_value, setter_func.__name__))
                        cif_item_txt = u"The value came from (%s) in cif."
                        xsd_txt = u"The value will be written for schema component (%s)."
                    else:
                        log_str.write(alert + u"Function (%s) failure.\n" % setter_func.__name__)
                        cif_item_txt = u"The value for cif category (%s) not given."
                        xsd_txt = u"The value is for schema component (%s)."
                        if parent_el_req:
                            parent_txt = u"The parent element IS required."
                        else:
                            parent_txt = u"The parent element IS NOT required. This trumps the above requirement."
                    if cif_item is not None:
                        log_str.write("%s%s\n" % (spacing, (cif_item_txt % cif_item)))
                    if xsd is not None:
                        log_str.write("%s%s\n" % (spacing, (xsd_txt % xsd)))
                    if parent_el_req is not None:
                        log_str.write("%s%s\n" % (spacing, parent_txt))
                    if em_for_emd is not None:
                        log_str.write("%s%s\n" % (spacing, (em_for_emd_txt % em_for_emd)))
                    if soft_name is not None:
                        log_str.write("%s%s\n" % (spacing, (soft_txt % soft_name)))

        def log(cif_item, setter_func, fmt_cif_value=None, parent_el_req=None, soft_name=None):
            """
            Logs the outcome of setting the value for cif_item
            using setter_func

            Params:
            @param cif_item: cif_value is for this cif_item
            @param setter_func: The function used to set cif_value
            @param fmt_cif_value: The cif value is formatted if the value is set
            This is used as a flag that setting went well
            """
            xsd = None
            em_for_emd = None
            req = None
            if cif_item is not None:
                xsd = get_xsd_for_cif_item(cif_item)
                # em_for_emd = get_em_from_emd(cif_item)
            if xsd is not None:
                req = is_cif_item_required(xsd)
            elif parent_el_req is None:
                req = parent_el_req

            if fmt_cif_value is not None:
                # this is an info
                self.entry_in_translation_log.add_info(cif_item, setter_func.__name__, xsd, em_for_emd, fmt_cif_value, soft_name=soft_name)
                if self.info_log_string is None or self.info_log_string.closed:
                    self.open_logs(open_info=True, open_warn=False, open_error=False)
                log_as_str(self.info_log_string, cif_item, setter_func, xsd, em_for_emd, const.INFO_ALERT, fmt_cif_value, soft_name=soft_name)
            else:
                # this is either a warning or an error
                if req:
                    self.entry_in_translation_log.add_err(cif_item, setter_func.__name__, xsd, em_for_emd, parent_el_req=parent_el_req, soft_name=soft_name)
                    # This is an error
                    if self.error_log_string is None or self.error_log_string.closed:
                        self.open_logs(open_info=False, open_warn=False, open_error=True)
                    log_as_str(self.error_log_string, cif_item, setter_func, xsd, em_for_emd, const.REQUIRED_ALERT, parent_el_req=parent_el_req, soft_name=soft_name)
                else:
                    # This is a warning since not required
                    # self.entry_in_translation_log.add_warn(cif_item, setter_func.__name__, xsd, em_for_emd, parent_el_req=parent_el_req, soft_name=soft_name)
                    if self.warn_log_string is None or self.warn_log_string.closed:
                        self.open_logs(open_info=False, open_warn=True, open_error=False)
                        # log_as_str(self.warn_log_string, cif_item, setter_func, xsd, em_for_emd, const.NOT_REQUIRED_ALERT, parent_el_req=parent_el_req, soft_name=soft_name)

        def set_cif_value(setter_func, *cif_key_cat, **addons):
            """
            This function performs one translation of a cif value into an XML using generateDS generated code (setter_func).
            The function extends functionality of the original safe_set function

            One translation involves passing of a cif value into the setter function (setter_func).
            The cif value can be calculated before the function call and passed on. In that case the call is:

            set_cif_value(setter_func, cif_value=a_value). In this case the logger only logs INFO.

            When the cif value is not passed into the function, *cif_key_cat as two unnamed arguments have to be given so that
            can be unpacked as a cif_key, cif_cat tupple and used in getting the cif value from the get_cif_item(cif_key, cif_cat).
            In this case the call is:

            set_cif_value(setter_func, cif_key, cif_cat)

            If the list of cif categories is needed in getting the cif value in addition to the cif_key, cif_cat arguments
            the named argument cif_list must be provided:

            set_cif_value(setter_func, cif_key, cif_cat, cif_list=a_cif_list)

            If the cif value need to be converted into another type (e.g. integer) the named argument fmt is used:

            set_cif_value(setter_func, cif_key, cif_cat, fmt=int)

            Formatting may be from a string into a date in which case the call is:

            set_cif_value(setter_func, cif_key, cif_cat, fmt='date')

            If the cif value should be from a dictionary the call is:

            set_cif_value(setter_func, cif_key, cif_cat, fmt=a_dict)

            Lastly, as far as formatting is concerned, it could be a function (e.g. convert string to upper case or
            a lambda expression) so the calls might be:

            set_cif_value(setter_func, cif_key, cif_cat, fmt=str.upper)
            set_cif_value(setter_func, cif_key, cif_cat, fmt=lambda x: func(x))

            There are cases when a constructor is needed in order to set the cif value. The calls are then:

            set_cif_value(setter_func, cif_key, cif_cat, constructor=emdb.a_type)
            set_cif_value(setter_func, cif_key, cif_cat, cif_list=a_cif_list, constructor=emdb.a_type)

            If the constructor needs additional values the following named arguments can be provided:
            'units', 'type', 'ncbi', 'private', and 'res_type'. Only 'units' and 'res_type' appear simultaneously.
            The rest are never used in a conjuction with the other named arguments apart from fmt.
            The calls that add constructor's functionality are:

            set_cif_value(setter_func, cif_key, cif_cat, constructor=emdb.a_type, units=a_units)
            set_cif_value(setter_func, cif_key, cif_cat, constructor=emdb.a_type, type=a_type)
            set_cif_value(setter_func, cif_key, cif_cat, constructor=emdb.a_type, ncbi=a_ncbi_value)
            set_cif_value(setter_func, cif_key, cif_cat, constructor=emdb.a_type, private=True)
            set_cif_value(setter_func, cif_key, cif_cat, constructor=emdb.a_type, res_type='BY AUTHOR')

            The 'parent_el_req' named argument is used to override the requirement for the XML element/argument
            created in this translation by it's parent element requirement. The call can be:

            set_cif_value(setter_func, cif_key, cif_cat, parent_el_req=False)

            The 'soft_name' named argument is used to pass the name of software as it's the only way to pass it onto logger
            """
            cif_key = None
            cif_cat = None
            cif_item = None
            if len(cif_key_cat) > 1:
                cif_key, cif_cat = cif_key_cat
            cif_list = addons.get("cif_list", None)
            fmt = addons.get("fmt", None)
            constructor = addons.get("constructor", None)
            units = addons.get("units", None)
            cif_value = addons.get("cif_value", None)
            the_type = addons.get("type", None)
            ncbi = addons.get("ncbi", None)
            private = addons.get("private", None)
            parent_el_req = addons.get("parent_el_req", None)
            res_type = addons.get("res_type", None)
            soft_name = addons.get("soft_name", None)

            # if cif_value is given key, cat and list shouldn't be used
            if cif_key is not None and cif_cat is not None:
                cif_item = get_cif_item(cif_key, cif_cat)
                if cif_value is None:
                    cif_value = get_cif_value(cif_key, cif_cat, cif_list)

            if cif_value is not None and cif_value != "NULL":
                if fmt is not None:
                    if fmt == "date":
                        # pre-set the value of date in case cif_value is not in a correct format; this value will be written out
                        fmt_cif_value = datetime.datetime.strptime("1000-01-01", "%Y-%m-%d").date()
                        try:
                            fmt_cif_value = datetime.datetime.strptime(cif_value, "%Y-%m-%d").date()
                        except Exception as exp:
                            sub_txt = const.CHANGE_MADE + u"Date set to (%s) for (%s) instead of (%s)." % (fmt_cif_value, "_" + cif_cat + "." + cif_key, cif_value)
                            txt = str(exp) + u" " + sub_txt
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, "(" + self.entry_in_translation_log.id + ")" + const.REQUIRED_ALERT + txt)
                    elif isinstance(fmt, dict):
                        fmt_cif_value = fmt[cif_value]
                    else:
                        fmt_cif_value = fmt(cif_value)
                else:
                    fmt_cif_value = cif_value
                if constructor is not None:
                    if units is None and the_type is None and ncbi is None and private is None and res_type is None:
                        constructed_cif_value = constructor(valueOf_=fmt_cif_value)

                    if units is not None:
                        constructed_cif_value = constructor(valueOf_=fmt_cif_value, units=units)
                    if the_type is not None:
                        constructed_cif_value = constructor(valueOf_=fmt_cif_value, type_=the_type)
                    if ncbi is not None:
                        constructed_cif_value = constructor(valueOf_=fmt_cif_value, ncbi=ncbi)
                    if private is not None:
                        constructed_cif_value = constructor(valueOf_=fmt_cif_value, private=private)
                    if res_type is not None:
                        if units is not None:
                            constructed_cif_value = constructor(valueOf_=fmt_cif_value, units=units, res_type=res_type)
                        else:
                            constructed_cif_value = constructor(valueOf_=fmt_cif_value, res_type=res_type)
                else:
                    constructed_cif_value = fmt_cif_value
                setter_func(constructed_cif_value)  # pylint: disable=possibly-used-before-assignment
                log(cif_item, setter_func, fmt_cif_value, soft_name=soft_name)
            else:
                log(cif_item, setter_func, parent_el_req=parent_el_req, soft_name=soft_name)

        def make_list_of_dicts(cif_category, key_item, cif_list=None, min_length=3):
            """
            Create a dictionary from a cif category where
            the items are keyed by the value of key_item

            Parameters:
            @param cif_category: cif category (without leadin underscore)
            @param key_item: item label -> value becomes the key
            @param cif_list: If this is provided then this is taken as the list of cif categories
            @param min_length: by default should be 3 or greater (id and foreign key the other two)
            @return dict_in: dictionary with a list of lists
            """
            dict_in = {}
            if cif_category in self.cif or cif_list is not None:
                list_in = cif_list or self.cif[cif_category]
                for list_item in list_in:
                    if len(list_item) < min_length:
                        continue
                    key = get_cif_value(key_item, cif_category, list_item)
                    if key in dict_in:
                        dict_in[key].append(list_item)
                    else:
                        dict_in[key] = [list_item]
            return dict_in

        def make_dict(cif_category, key_item, min_length=3):
            """
            Create a dictionary from a cif category where the items are keyed by the value of key_item.
            The difference from make_list_of_dicts is that it is assumed that there is only one list.

            Parameters:
            @param cif_category: cif category (without leading underscore)
            @param key_item: item label -> value becomes the key
            @param min_length: by default should be 3 or greater
                               (id and foreign key the other two)
            @return dict_in: dictionary with a list
            """
            dict_in = {}
            if cif_category in self.cif:
                list_in = self.cif[cif_category]
                for list_item in list_in:
                    if len(list_item) < min_length:
                        continue
                    key = get_cif_value(key_item, cif_category, list_item)
                    if key is not None:
                        dict_in[key] = list_item
            return dict_in

        def set_software_list(software_category, category_dict, setter_func):
            """
            Creates a software list for the specified software_category and
            sets it using the setter_func.

            Parameters:
            @param software_category: a software category, e.g. RECONSTRUCTION
            @param category_dict: dictionary of list of lists keyed by software_category
            @param setter_func: setter_func for adding software list
            XSD: <xs:complexType name="software_list_type"> is
            .. a sequence of 1 element
            """

            def set_software_type(software_category, soft, soft_in):
                """
                XSD: <xs:element name="software" type="software_type" minOccurs="0" maxOccurs="unbounded"/> has
                .. 3 elements
                """

                def set_el_name(soft, soft_in):
                    """
                    XSD: <xs:element name="name" type="xs:token" minOccurs="0"/>
                    CIF: _em_software.name
                    """
                    set_cif_value(soft.set_name, "name", const.EM_SOFTWARE, cif_list=soft_in, soft_name=software_category)

                def set_el_version(soft, soft_in):
                    """
                    XSD: <xs:element name="version" type="xs:token" minOccurs="0"/>
                    CIF: _em_software.version
                    """
                    set_cif_value(soft.set_version, "version", const.EM_SOFTWARE, cif_list=soft_in)

                def set_el_processing_details(soft, soft_in):
                    """
                    XSD: <xs:element name="processing_details" type="xs:string" minOccurs="0"/>
                    CIF: _em_software.details
                    """
                    set_cif_value(soft.set_processing_details, "details", const.EM_SOFTWARE, cif_list=soft_in)

                # element 1
                set_el_name(soft, soft_in)
                # element 2
                set_el_version(soft, soft_in)
                # element 3
                set_el_processing_details(soft, soft_in)

            if software_category in category_dict:
                soft_list = emdb.software_list_type()
                for soft_in in category_dict[software_category]:
                    soft = emdb.software_type()
                    set_software_type(software_category, soft, soft_in)
                    if soft.has__content():
                        soft_list.add_software(soft)
                if soft_list.has__content():
                    setter_func(soft_list)

        def assert_get_value(key, dic):
            """
            Return dict[value] but throw an exception if key is not found

            @param key: key of (key,value) pair
            @param dic: dictionary
            """
            if key not in dic:
                raise KeyError("Key %s not found in dictionary %s" % (key, dic))
            return dic[key]

        def format_author(auth_in):
            """
            Convert author from CIF format Smith, J.S. to pubmed/emdb format Smith JS
            allowed_authors = ["Ashish",
                               "Garcia-Moreno E., B.",
                               "van 't Hag, L.",
                               "Nur 'Izzah, N.",
                               "Ihsanawati",
                               "Preeti",
                               "Morigen",
                               "Nolte-'T Hoen, E.N.M."

            Parameters:
            @param auth_in: string author name in CIF format
            @return auth_out: author in EMDB format
            """
            auth_out = auth_in
            if auth_in not in const.CENTER_NAMES_USED_FOR_AUTHORS:
                auth_match = re.match(const.CIF_AUTHOR_RE, auth_in)
                if auth_match is not None and not auth_in.isspace():
                    match_groups = auth_match.groups()
                    auth_out = "%s %s" % (match_groups[0].replace(".", ""), match_groups[1].replace(".", ""))
                else:
                    auth_out = ""
                    txt = u"Author name: (%s) is not in a required CIF format." % auth_in
                    self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                    self.log_formatted(self.error_log_string, "(" + self.entry_in_translation_log.id + ")" + const.REQUIRED_ALERT + txt)
            return auth_out

        def set_admin_type(admin):
            """
            XSD: <xs:complexType name="admin_type"> is
            ..a sequence of 15 elements
            """

            def set_version_type(admin_status):
                """
                XSD: <xs:complexType name="version_type"> is
                ..a sequence of 5 elements
                """

                def set_el_date(admin_status):
                    """
                    XSD: <xs:element name="date" minOccurs="0">
                    CIF: _em_admin.last_update
                    """
                    set_cif_value(admin_status.set_date, "last_update", const.EM_ADMIN, fmt="date")

                def set_el_code(admin_status):
                    """
                    XSD: <xs:element name="code" type="code_type">
                    CIF: _em_admin.current_status HPUB
                    """
                    set_cif_value(admin_status.set_code, "current_status", const.EM_ADMIN, constructor=emdb.code_type)

                def set_el_processing_site(admin_status):
                    """
                    XSD: <xs:element name="processing_site" minOccurs="0">
                    CIF: _pdbx_database_status.process_site
                    NOT YET IMPLEMENTED
                    Values to expect: PDBe, PDBj, RCSB, PDBc
                    """
                    set_cif_value(admin_status.set_processing_site, "process_site", const.PDBX_DATABASE_STATUS, fmt=const.PROC_SITE_CIF2XML)

                def set_el_annotator(admin_status):
                    """
                    XSD: <xs:element name="annotator" minOccurs="0">
                    CIF: _pdbx_database_status.pdbx_annotator
                    """
                    if self.__show_private:
                        set_cif_value(admin_status.set_annotator, "pdbx_annotator", const.PDBX_DATABASE_STATUS, constructor=emdb.annotatorType, private="true")

                def set_el_details(admin_status):
                    """
                    XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                    CIF: _em_admin.details
                    """
                    set_cif_value(admin_status.set_details, "details", const.EM_ADMIN)

                # element 1
                set_el_date(admin_status)
                # element 2
                set_el_code(admin_status)
                # element 3q
                set_el_processing_site(admin_status)
                # element 4
                set_el_annotator(admin_status)
                # element 5
                set_el_details(admin_status)

            def set_attr_composite_map(admin):
                """
                XSD: <xs:attribute name="composite_map" type="xs:boolean">
                CIF: _em_admin.composite_map
                """
                composite = get_cif_value("composite_map", const.EM_ADMIN)
                if composite is not None and composite == "YES":
                    set_cif_value(admin.set_composite_map, cif_value=True)

            def set_el_status_history_list():
                """
                XSD: <xs:element name="status_history_list" type="version_list_type" minOccurs="0">
                CIF: doesn't exist yet
                TO_IMPLEMENT when the OneDep is changed: admin.set_status_history_list(...)
                """

            def set_el_current_status(admin):
                """
                XSD: <xs:element name="current_status" type="version_type">
                """
                status = emdb.version_type()
                set_version_type(status)
                admin.set_current_status(status)

            def set_el_revision_history(admin):
                """
                List of entry's changes. The fisrt revision is from the initial release
                XSD: <xs:element minOccurs="0" name="revision_history">
                CIF: pdbx_audit_revision_history
                """

                def set_revisions_type(revisions, revision_history_in):
                    """
                    XSD: <xs:element name="revision" type="revision_history_type" maxOccurs="unbounded"/>
                    CIF: pdbx_audit_revision_history
                    """

                    def set_revision_type(revision, revision_list, revision_categories, revision_items):
                        """
                        XSD: <xs:complexType name="revision_history_type"> is
                        ... a sequence of 1 element and two attributes
                        """

                        def set_attr_version(revision, revision_list):
                            """
                            XSD: <xs:attribute name="version" type="revision_history_version_type" use="required"/>
                            CIF: (pdbx_audit_revision_history.major_revision).(pdbx_audit_revision_history.minor_revision)
                            """
                            major_revision = get_cif_value("major_revision", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_list[0])
                            minor_revision = get_cif_value("minor_revision", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_list[0])
                            revision_full = ".".join([major_revision, minor_revision])
                            set_cif_value(revision.set_version, cif_value=revision_full)

                        def set_attr_date(revision, revision_list):
                            """
                            XSD: <xs:attribute name="date" type="xs:date" use="required"/>
                            CIF: pdbx_audit_revision_history.revision_date
                            """
                            set_cif_value(revision.set_date, "revision_date", const.PDBX_AUDIT_REVISION_HISTORY,
                                          cif_list=revision_list[0], fmt="date")

                        def set_el_change_list(revision, revision_list):
                            """
                            XSD: <xs:element name="change_list" maxOccurs="1">
                            CIF: pdbx_audit_revision_details and pdbx_audit_revision_history
                            """

                            def set_change_list_type(change_list, revision_list):
                                """
                                <xs:element name="change_list" maxOccurs="1"> is
                                ... a sequence of revision_change_sub_group substitutions
                                """

                                def set_base_revision_change_type(base_content_type, revision_detail_in, revision_gr_in):
                                    """
                                    <xs:complexType name="base_revision_change_type"> is
                                    ... a sequence of 5 elements
                                    """

                                    def set_el_revision_type(base_content_type, revision_detail_in):
                                        """
                                        XSD: <xs:element name="revision_type">
                                        CIF: pdbx_audit_revision_details.type
                                        """
                                        set_cif_value(base_content_type.set_revision_type, "type", const.PDBX_AUDIT_REVISION_DETAILS, cif_list=revision_detail_in, fmt=const.MAP_REVISION_DETAILS_TYPE)

                                    def set_el_provider(base_content_type, revision_detail_in):
                                        """
                                        XSD: <xs:element name="provider">
                                        CIF: pdbx_audit_revision_details.provider
                                        """
                                        set_cif_value(base_content_type.set_provider, "provider", const.PDBX_AUDIT_REVISION_DETAILS, cif_list=revision_detail_in, fmt=const.MAP_REVISION_PROVIDER)

                                    def set_el_description(base_content_type, revision_detail_in):
                                        """
                                        XSD: <xs:element name="description" type="xs:token" minOccurs="0"/>
                                        CIF: pdbx_audit_revision_details.description
                                        """
                                        set_cif_value(base_content_type.set_description, "description", const.PDBX_AUDIT_REVISION_DETAILS, cif_list=revision_detail_in)

                                    def set_el_details(base_content_type, revision_detail_in):
                                        """
                                        XSD: <xs:element name="details" type="xs:token" minOccurs="0"/>
                                        CIF: pdbx_audit_revision_details.details
                                        """
                                        set_cif_value(base_content_type.set_details, "details", const.PDBX_AUDIT_REVISION_DETAILS, cif_list=revision_detail_in)

                                    def set_el_revision_group(base_content_type, revision_gr_in):
                                        """
                                        XSD: <xs:element name="revision_group" minOccurs="0">
                                        CIF: pdbx_audit_revision_group.group
                                        """
                                        set_cif_value(base_content_type.set_revision_group, "group", const.PDBX_AUDIT_REVISION_GROUP, cif_list=revision_gr_in, fmt=const.MAP_REVISION_GROUP)

                                    if revision_detail_in is not None:
                                        # element 1
                                        set_el_revision_type(base_content_type, revision_detail_in)
                                        # element 2
                                        set_el_provider(base_content_type, revision_detail_in)
                                        # element 3
                                        set_el_description(base_content_type, revision_detail_in)
                                        # element 4
                                        set_el_details(base_content_type, revision_detail_in)

                                    if revision_gr_in is not None:
                                        # element 5
                                        set_el_revision_group(base_content_type, revision_gr_in)

                                def set_part_revision_change_type(part_content_type, revision_in, revision_detail_in):
                                    """
                                    <xs:complexType name="part_revision_change_type"> is
                                    .. an extenstion of <xs:extension base="base_revision_change_type"> and
                                    .. has 1 extra attribute
                                    """
                                    def set_attr_part(part_content_type, revision_in):
                                        """
                                        XSD: <xs:attribute name="part" type="xs:positiveInteger"/>
                                        CIF: pdbx_audit_revision_history.internal_part_number = 2
                                        """
                                        set_cif_value(part_content_type.set_part, "part_number", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_in)

                                    set_base_revision_change_type(part_content_type, revision_detail_in, revision_gr_in)
                                    # attribute 1
                                    set_attr_part(part_content_type, revision_in)

                                def set_categories_and_items(categories, items, revision_history_in,
                                                             revision_categories, revision_items):
                                    """
                                    <xs:complexType name="metadata_revision_type"> extends <xs:extension base="base_revision_change_type"> and
                                    ... has 2 more elements
                                    """

                                    def set_revision_category_or_item_type(category_or_item, revision_in, details_in):
                                        """
                                        <xs:complexType name="revision_category_or_item_type">
                                        ... has 3 attributes
                                        """
                                        def set_attrib_revision_type(category_or_item, revision_in):
                                            """
                                            XSD: <xs:attribute name="revision_type" type="xs:token" use="required"/>
                                            CIF: pdbx_audit_revision_history.data_content_type
                                            """
                                            set_cif_value(category_or_item.set_revision_type, "data_content_type", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_in, fmt=const.MAP_FILES)

                                        def set_attrib_part(category_or_item, revision_in):
                                            """
                                            XSD: <xs:attribute name="part" type="xs:positiveInteger"/>
                                            CIF: pdbx_audit_revision_history.part_number
                                            """
                                            set_cif_value(category_or_item.set_part, "part_number", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_in)

                                        def set_attrib_revision_action(category_or_item, details_in):
                                            """
                                            XSD: <xs:attribute name="revision_action" type="xs:token" use="required"/>
                                            CIF: pdbx_audit_revision_details.type
                                            """
                                            set_cif_value(category_or_item.set_revision_action, "type", const.PDBX_AUDIT_REVISION_DETAILS, cif_list=details_in)

                                        # attribute 1
                                        set_attrib_revision_type(category_or_item, revision_in)
                                        # attribute 2
                                        set_attrib_part(category_or_item, revision_in)
                                        # attribute 3
                                        set_attrib_revision_action(category_or_item, details_in)

                                    def set_el_category(cat, revision_cat_in, revision_in, details_in):
                                        """
                                        <xs:element name="category" type="revision_category_or_item_type" minOccurs="1" maxOccurs="unbounded"/>
                                        XSD: <xs:element name="category" type="revision_category_or_item_type" minOccurs="1" maxOccurs="unbounded"/>
                                        CIF: pdbx_audit_revision_category.category
                                        """
                                        set_revision_category_or_item_type(cat, revision_in, details_in)
                                        set_cif_value(cat.set_valueOf_, "category", const.PDBX_AUDIT_REVISION_CATEGORY, cif_list=revision_cat_in)

                                    def set_el_item(item, revision_it_in, revision_in, details_in):
                                        """
                                        <xs:element name="item" type="revision_category_or_item_type" minOccurs="1" maxOccurs="unbounded"/>
                                        XSD: <xs:element name="item" type="revision_category_or_item_type" minOccurs="1" maxOccurs="unbounded"/>
                                        CIF: pdbx_audit_revision_item.item
                                        """
                                        set_revision_category_or_item_type(item, revision_in, details_in)
                                        set_cif_value(item.set_valueOf_, "item", const.PDBX_AUDIT_REVISION_ITEM, cif_list=revision_it_in)

                                    # element 1
                                    for revision_cat_in in revision_categories:
                                        cat = emdb.revision_category_or_item_type()
                                        cat_ordinal = get_cif_value("revision_ordinal", const.PDBX_AUDIT_REVISION_CATEGORY, cif_list=revision_cat_in)
                                        revision_in = revision_history_in.get(cat_ordinal)
                                        details_in = revision_details_in.get(cat_ordinal)
                                        set_el_category(cat, revision_cat_in, revision_in, details_in)
                                        categories.add_category(cat)

                                    # element 2
                                    for revision_it_in in revision_items:
                                        item = emdb.revision_category_or_item_type()
                                        item_ordinal = get_cif_value("revision_ordinal", const.PDBX_AUDIT_REVISION_ITEM, cif_list=revision_it_in)
                                        revision_in = revision_history_in.get(item_ordinal)
                                        details_in = revision_details_in.get(item_ordinal)
                                        set_el_item(item, revision_it_in, revision_in, details_in)
                                        items.add_item(item)

                                revision_details_in = make_dict(const.PDBX_AUDIT_REVISION_DETAILS, "revision_ordinal")
                                revision_group_in = make_dict(const.PDBX_AUDIT_REVISION_GROUP, "revision_ordinal")
                                revision_history_in = make_dict(const.PDBX_AUDIT_REVISION_HISTORY, "ordinal")

                                for revision_in in revision_list:
                                    history_ordinal = get_cif_value("ordinal", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_in)
                                    revision_detail_in = revision_details_in.get(history_ordinal)
                                    revision_gr_in = revision_group_in.get(history_ordinal)
                                    data_content_type = get_cif_value("data_content_type", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_in)
                                    if data_content_type == const.PRIMARY_MAP:
                                        primary_map = emdb.base_revision_change_type()
                                        primary_map.original_tagname_ = "primary_map"
                                        set_base_revision_change_type(primary_map, revision_detail_in, revision_gr_in)
                                        if primary_map.has__content():
                                            change_list.add_revision_change_sub_group(primary_map)
                                    elif data_content_type == const.IMAGE:
                                        image = emdb.part_revision_change_type()
                                        image.original_tagname_ = "image"
                                        set_part_revision_change_type(image, revision_in, revision_detail_in)
                                        if image.has__content():
                                            change_list.add_revision_change_sub_group(image)
                                    elif data_content_type == const.HALF_MAP:
                                        half_map = emdb.part_revision_change_type()
                                        half_map.original_tagname_ = "half_map"
                                        set_part_revision_change_type(half_map, revision_in, revision_detail_in)
                                        if half_map.has__content():
                                            change_list.add_revision_change_sub_group(half_map)
                                    elif data_content_type == const.MASK:
                                        mask = emdb.part_revision_change_type()
                                        mask.original_tagname_ = "mask"
                                        set_part_revision_change_type(mask, revision_in, revision_detail_in)
                                        if mask.has__content():
                                            change_list.add_revision_change_sub_group(mask)
                                    elif data_content_type == const.ADDITIONAL_MAP:
                                        add_map = emdb.part_revision_change_type()
                                        add_map.original_tagname_ = "additional_map"
                                        set_part_revision_change_type(add_map, revision_in, revision_detail_in)
                                        if add_map.has__content():
                                            change_list.add_revision_change_sub_group(add_map)
                                    elif data_content_type == const.FSC:
                                        fsc = emdb.part_revision_change_type()
                                        fsc.original_tagname_ = "fsc"
                                        set_part_revision_change_type(fsc, revision_in, revision_detail_in)
                                        if fsc.has__content():
                                            change_list.add_revision_change_sub_group(fsc)
                                    elif data_content_type == const.STRUCTURE_MODEL:
                                        model = emdb.part_revision_change_type()
                                        model.original_tagname_ = "model"
                                        set_part_revision_change_type(model, revision_in, revision_detail_in)
                                        if model.has__content():
                                            change_list.add_revision_change_sub_group(model)
                                    elif data_content_type == const.EM_METADATA:
                                        # metadata can only have once instance and one of categories and items
                                        metadata = emdb.metadata_revision_type()
                                        categories = emdb.categoriesType()
                                        items = emdb.itemsType()
                                        metadata.original_tagname_ = "metadata"
                                        set_base_revision_change_type(metadata, revision_detail_in, revision_gr_in)
                                        set_categories_and_items(categories, items, revision_history_in,
                                                                 revision_categories, revision_items)
                                        if categories.has__content():
                                            metadata.set_categories(categories)
                                        if items.has__content():
                                            metadata.set_items(items)
                                        if metadata.has__content():
                                            change_list.add_revision_change_sub_group(metadata)

                            change_list = emdb.change_listType()
                            set_change_list_type(change_list, revision_list)
                            revision.set_change_list(change_list)

                        # attribute 1
                        set_attr_version(revision, revision_list)
                        # attribute 2
                        set_attr_date(revision, revision_list)
                        # element 1
                        set_el_change_list(revision, revision_list)

                    def set_revisons(revision_dict, rev_list, revisions_in, cif_cat):
                        """
                        Helper function
                        """
                        for revision, revision_list in rev_list.items(): # history
                            revisions = []
                            for a_revision in revision_list: # history list
                                a_ordinal = get_cif_value("ordinal", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=a_revision)
                                for rev in revisions_in.values(): # cat or items
                                    revision_ordinal = get_cif_value("revision_ordinal", cif_cat, cif_list=rev)
                                    if a_ordinal == revision_ordinal:
                                        major_revision = get_cif_value("major_revision", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=a_revision)
                                        minor_revision = get_cif_value("minor_revision", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=a_revision)
                                        revision_full = ".".join([major_revision, minor_revision])
                                        if revision == revision_full:
                                            revisions.append(rev)
                                            # break
                            revision_dict[revision] = revisions

                    revision_num = ""
                    revision_lists = {}
                    revision_list = []

                    for revision_in in revision_history_in.values():
                        major_revision = get_cif_value("major_revision", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_in)
                        minor_revision = get_cif_value("minor_revision", const.PDBX_AUDIT_REVISION_HISTORY, cif_list=revision_in)
                        revision_full = ".".join([major_revision, minor_revision])
                        if revision_num == revision_full:
                            revision_list.append(revision_in)
                        else:
                            if revision_num != "":
                                revision_lists.update({revision_num: revision_list})
                            revision_num = revision_full
                            revision_list = []
                            revision_list.append(revision_in)
                        revision_lists.update({revision_full: revision_list})

                    # get categories and items into revisions
                    revisions_categories = {key: [] for key in revision_lists}
                    revisions_items = {key: [] for key in revision_lists}

                    revision_category_in = make_dict(const.PDBX_AUDIT_REVISION_CATEGORY, "ordinal")
                    revision_item_in = make_dict(const.PDBX_AUDIT_REVISION_ITEM, "ordinal")

                    set_revisons(revisions_categories, revision_lists, revision_category_in, const.PDBX_AUDIT_REVISION_CATEGORY)
                    set_revisons(revisions_items, revision_lists, revision_item_in, const.PDBX_AUDIT_REVISION_ITEM)

                    for revision_num, revision_list in revision_lists.items():
                        revision = emdb.revision_history_type()
                        set_revision_type(revision, revision_list, revisions_categories[revision_num], revisions_items[revision_num])
                        revisions.add_revision(revision)

                revision_history_in = make_dict(const.PDBX_AUDIT_REVISION_HISTORY, "ordinal")

                revisions = emdb.revision_historyType()
                set_revisions_type(revisions, revision_history_in)
                if revisions.has__content():
                    admin.set_revision_history(revisions)

            def set_el_sites(admin):
                """
                Deposition and processing sites
                XSD: <xs:element name="sites">
                """

                def set_sites_type(dep_proc_sites):
                    """
                    XSD: <xs:element name="sites"> has
                    .. 2 elements
                    """

                    def set_el_deposition(dep_proc_sites):
                        """
                        XSD: <xs:element name="deposition">
                        CIF: _em_admin.deposition_site PDBE
                        """
                        set_cif_value(dep_proc_sites.set_deposition, "deposition_site", const.EM_ADMIN, fmt=const.PROC_SITE_CIF2XML)

                    def set_el_last_processing(dep_proc_sites):
                        """
                        XSD: <xs:element name="last_processing">
                        CIF: _pdbx_database_status.process_site PDBE
                        """
                        set_cif_value(dep_proc_sites.set_last_processing, "process_site", const.PDBX_DATABASE_STATUS, fmt=const.PROC_SITE_CIF2XML)

                    # element 1
                    set_el_deposition(dep_proc_sites)
                    # element 2
                    set_el_last_processing(dep_proc_sites)

                dep_proc_sites = emdb.sitesType()
                set_sites_type(dep_proc_sites)
                admin.set_sites(dep_proc_sites)

            def set_el_key_dates(admin):
                """
                XSD: <xs:element name="key_dates">
                """

                def set_key_dates_type(key_dates):
                    """
                    XSD: <xs:element name="key_dates"> is a
                    ..sequence of 5 elements
                    """

                    def set_el_deposition(key_dates):
                        """
                        XSD: <xs:element name="deposition" type="xs:date">
                        CIF: _em_admin.deposition_date 2016-10-02
                        """
                        set_cif_value(key_dates.set_deposition, "deposition_date", const.EM_ADMIN, fmt="date")

                    def set_el_header_release(key_dates):
                        """
                        XSD: <xs:element name="header_release" type="xs:date" minOccurs="0">
                        CIF: _em_admin.header_release_date 2016-11-09
                        """
                        set_cif_value(key_dates.set_header_release, "header_release_date", const.EM_ADMIN, fmt="date")

                    def set_el_map_release(key_dates):
                        """
                        XSD: <xs:element name="map_release" type="xs:date" minOccurs="0">
                        CIF: _em_admin.map_release_date ?
                        """
                        set_cif_value(key_dates.set_map_release, "map_release_date", const.EM_ADMIN, fmt="date")

                    def set_el_obsolete(key_dates):
                        """
                        XSD: <xs:element name="obsolete" type="xs:date" minOccurs="0">
                        CIF: _em_admin.obsoleted_date ?
                        """
                        set_cif_value(key_dates.set_obsolete, "obsoleted_date", const.EM_ADMIN, fmt="date")

                    def set_el_update(key_dates):
                        """
                        XSD: <xs:element name="update" type="xs:date">
                        CIF: _em_admin.last_update 2016-11-09
                        """
                        set_cif_value(key_dates.set_update, "last_update", const.EM_ADMIN, fmt="date")

                    # element 1
                    set_el_deposition(key_dates)
                    # element 2
                    set_el_header_release(key_dates)
                    # element 3
                    set_el_map_release(key_dates)
                    # element 4
                    set_el_obsolete(key_dates)
                    # element 5
                    set_el_update(key_dates)

                key_dates = emdb.key_datesType()
                set_key_dates_type(key_dates)
                admin.set_key_dates(key_dates)

            def set_el_obsolete_list(admin):
                """
                This list contains old entries that have been replaced
                because of this newer entry.
                XSD: <xs:element name="obsolete_list" minOccurs="0"> has
                .. 1 element <xs:element name="entry" type="supersedes_type" maxOccurs="unbounded">
                CIF: _pdbx_database_PDB_obs_spr
                """

                def set_obsolete_list_type(obs_list, obsolete_in):
                    """
                    XSD: <xs:element name="obsolete_list" minOccurs="0"> has
                    .. 1 element of supersedes_type
                    """

                    def set_supersedes_type(obs_entry, obs_in):
                        """
                        XSD: <xs:complexType name="supersedes_type"> has
                        .. 3 elements
                        """

                        def set_el_date(obs_entry, obs_in):
                            """
                            XSD: <xs:element name="date" type="xs:date"/>
                            CIF: _em_obsolete.date
                            """
                            set_cif_value(obs_entry.set_date, "date", const.EM_OBSOLETE, cif_list=obs_in, fmt="date")

                        def set_el_entry(obs_entry, obs_in):
                            """
                            XSD: <xs:element name="entry" type="emdb_id_type"/>
                            CIF: _em_obsolete.entry
                            """
                            set_cif_value(obs_entry.set_entry, "entry", const.EM_OBSOLETE, cif_list=obs_in)

                        def set_el_details(obs_entry, obs_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _em_obsolete.details
                            """
                            set_cif_value(obs_entry.set_details, "details", const.EM_OBSOLETE, cif_list=obs_in)

                        # element 1
                        set_el_date(obs_entry, obs_in)
                        # element 2
                        set_el_entry(obs_entry, obs_in)
                        # element 3
                        set_el_details(obs_entry, obs_in)

                    for _obs_id, obs_in in obsolete_in.items():
                        obs_entry = emdb.supersedes_type()
                        set_supersedes_type(obs_entry, obs_in)
                        obs_list.add_entry(obs_entry)

                obsolete_in = make_dict(const.EM_OBSOLETE, "id")
                obs_list = emdb.obsolete_listType()
                set_obsolete_list_type(obs_list, obsolete_in)
                if obs_list.has__content():
                    admin.set_obsolete_list(obs_list)

            def set_el_superseded_by_list(admin):
                """
                If this element appears means that the current entry is obsoleted.
                The newer entry which replaces this entry is listed here.
                XSD: <xs:element name="superseded_by_list" minOccurs="0"> has
                .. 1 element <xs:element name="entry" type="supersedes_type" maxOccurs="unbounded">
                CIF: _em_supersede
                """

                def set_supersede_by_list(supersed_list, supr_in):
                    """
                    XSD: <xs:element name="superseded_by_list" minOccurs="0"> has
                    .. 1 element of supersedes_type
                    """

                    def set_supersede_entry(spr_entry, spr_in):
                        """
                        XSD: <xs:complexType name="supersedes_type"> has
                        .. 3 elements
                        """

                        def set_el_date(spr_entry, spr_in):
                            """
                            XSD: <xs:element name="date" type="xs:date"/>
                            CIF: _em_supersede.date
                            """
                            set_cif_value(spr_entry.set_date, "date", const.EM_SUPERSEDE, cif_list=spr_in, fmt="date")

                        def set_el_entry(spr_entry, spr_in):
                            """
                            XSD: <xs:element name="entry" type="emdb_id_type"/>
                            CIF: _em_supersede.entry
                            """
                            set_cif_value(spr_entry.set_entry, "entry", const.EM_SUPERSEDE, cif_list=spr_in)

                        def set_el_details(spr_entry, spr_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _em_supersede.details
                            """
                            set_cif_value(spr_entry.set_details, "details", const.EM_SUPERSEDE, cif_list=spr_in)

                        # element 1
                        set_el_date(spr_entry, spr_in)
                        # element 2
                        set_el_entry(spr_entry, spr_in)
                        # element 3
                        set_el_details(spr_entry, spr_in)

                    for _spr_id, spr_in in supr_in.items():
                        spr_entry = emdb.supersedes_type()
                        set_supersede_entry(spr_entry, spr_in)
                        supersed_list.add_entry(spr_entry)

                supr_in = make_dict(const.EM_SUPERSEDE, "id")
                supersede_list = emdb.superseded_by_listType()
                set_supersede_by_list(supersede_list, supr_in)
                if supersede_list.has__content():
                    admin.set_superseded_by_list(supersede_list)

            def set_el_grant_support(admin, aud_sup_in):
                """
                XSD: <xs:element name="grant_support" minOccurs="0">
                CIF: _pdbx_audit_support
                """

                def set_grant_support_type(grant_support, aud_sup_in):
                    """
                    XSD: <xs:element name="grant_support" minOccurs="0"> has
                    .. 1 element of grant_reference_type
                    """

                    def set_grant_reference_type(grant_ref, aud_sup):
                        """
                        XSD: <xs:element name="grant_reference" type="grant_reference_type" maxOccurs="unbounded"/> has
                        .. 3 elements
                        """

                        def set_el_funding_body(grant_ref, aud_sup, parent_req=False):
                            """
                            XSD: <xs:element name="funding_body" type="xs:token"/>
                            CIF: _pdbx_audit_support.funding_organization
                            """
                            set_cif_value(grant_ref.set_funding_body, "funding_organization", const.PDBX_AUDIT_SUPPORT, cif_list=aud_sup, parent_el_req=parent_req)

                        def set_el_code(grant_ref, aud_sup, parent_req=False):
                            """
                            XSD: <xs:element name="code" type="xs:token minOccures="0"/>
                            CIF: _pdbx_audit_support.grant_number
                            """
                            set_cif_value(grant_ref.set_code, "grant_number", const.PDBX_AUDIT_SUPPORT, cif_list=aud_sup, parent_el_req=parent_req)

                        def set_el_country(grant_ref, aud_sup, parent_req=False):
                            """
                            XSD: <xs:element name="country" type="xs:token" minOccurs="0"/>
                            CIF: _pdbx_audit_support.country
                            """
                            set_cif_value(grant_ref.set_country, "country", const.PDBX_AUDIT_SUPPORT, cif_list=aud_sup, parent_el_req=parent_req)

                        # element 1
                        set_el_funding_body(grant_ref, aud_sup, parent_req=False)
                        # element 2
                        set_el_code(grant_ref, aud_sup, parent_req=False)
                        # element 3
                        set_el_country(grant_ref, aud_sup, parent_req=False)

                    for _aud_sup_key, aud_sup in aud_sup_in.items():
                        el_funding_body = get_cif_value("funding_organization", const.PDBX_AUDIT_SUPPORT, cif_list=aud_sup)
                        # el_code = get_cif_value('grant_number', const.PDBX_AUDIT_SUPPORT, cif_list=aud_sup)
                        # el_country = get_cif_value('country', const.PDBX_AUDIT_SUPPORT, cif_list=aud_sup)
                        if el_funding_body is not None:  # or el_code is not None or el_country is not None:
                            grant_ref = emdb.grant_reference_type()
                            set_grant_reference_type(grant_ref, aud_sup)
                            if grant_ref.has__content():
                                grant_support.add_grant_reference(grant_ref)

                grant_support = emdb.grant_supportType()
                set_grant_support_type(grant_support, aud_sup_in)
                if grant_support.has__content():
                    admin.set_grant_support(grant_support)

            def set_el_contact_author(admin, contact_auth_in):
                """
                XSD: <xs:element name="contact_author" maxOccurs="unbounded" minOccurs="0">
                CIF: _pdbx_contact_author
                """

                def set_contact_author_type(cont_author, contact_auth_in):
                    """
                    <xs:element name="contact_author" maxOccurs="unbounded" minOccurs="0"> has
                    .. a base of contact_details_type and
                    .. 1 attribute
                    """

                    def set_contact_details_type(cont_author, contact_auth_in, parent_req):
                        """
                        XSD: <xs:complexType name="contact_details_type"> has
                        .. 14 elements
                        @param parent_req: the requirement for the elements depend on
                                the parent element requirement
                        """

                        def set_el_role(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="role">
                            CIF: _pdbx_contact_author.role
                            """
                            set_cif_value(cont_author.set_role, "role", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, fmt=str.upper, parent_el_req=parent_req)

                        def set_el_title(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="title">
                            CIF: _pdbx_contact_author.name_salutation
                            """
                            set_cif_value(cont_author.set_title, "name_salutation", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, fmt=str.upper, parent_el_req=parent_req)

                        def set_el_first_name(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="first_name" type="xs:token">
                            CIF: _pdbx_contact_author.name_first
                            """
                            set_cif_value(cont_author.set_first_name, "name_first", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        def set_el_middle_name(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="middle_name">
                            CIF: _pdbx_contact_author.name_mi
                            """
                            set_cif_value(cont_author.set_middle_name, "name_mi", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, fmt=str.upper, parent_el_req=parent_req)

                        def set_el_last_name(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="last_name" type="xs:token">
                            CIF: _pdbx_contact_author.name_last
                            """
                            set_cif_value(cont_author.set_last_name, "name_last", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        def set_el_organization(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="organization">
                            CIF: _pdbx_contact_author.organization_type
                            """
                            set_cif_value(
                                cont_author.set_organization,
                                "organization_type",
                                const.PDBX_CONTACT_AUTHOR,
                                cif_list=contact_auth_in,
                                constructor=emdb.organizationType,
                                type="",
                                fmt=str.upper,
                                parent_el_req=parent_req,
                            )

                        def set_el_street(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="street" type="xs:string">
                            CIF: _pdbx_contact_author.address_1
                            CIF: _pdbx_contact_author.address_2
                            CIF: _pdbx_contact_author.address_3
                            """
                            address_1 = get_cif_value("address_1", const.PDBX_CONTACT_AUTHOR, contact_auth_in)
                            address_2 = get_cif_value("address_2", const.PDBX_CONTACT_AUTHOR, contact_auth_in)
                            address_3 = get_cif_value("address_3", const.PDBX_CONTACT_AUTHOR, contact_auth_in)
                            full_address = None
                            if address_1 is not None and address_2 is not None and address_3 is not None:
                                full_address = ", ".join(filter(None, (address_1, address_2, address_3)))
                            set_cif_value(
                                cont_author.set_street, "address_1", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, cif_value=full_address, parent_el_req=parent_req
                            )

                        def set_el_town_or_city(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="town_or_city" type="xs:token">
                            CIF: _pdbx_contact_author.city
                            """
                            set_cif_value(cont_author.set_town_or_city, "city", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        def set_el_state_or_province(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="state_or_province" type="xs:token">
                            CIF: _pdbx_contact_author.state_province
                            """
                            set_cif_value(cont_author.set_state_or_province, "state_province", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        def set_el_country(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="country" type="xs:token">
                            CIF: _pdbx_contact_author.country
                            """
                            set_cif_value(cont_author.set_country, "country", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        def set_el_post_or_zip_code(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="post_or_zip_code" type="xs:token">
                            CIF: _pdbx_contact_author.postal_code
                            """
                            set_cif_value(cont_author.set_country, "postal_code", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        def set_el_email(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="email">
                            CIF: _pdbx_contact_author.email
                            """
                            set_cif_value(cont_author.set_email, "email", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        def set_el_telephone(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="telephone" type="telephone_number_type">
                            CIF: _pdbx_contact_author.phone
                            """
                            set_cif_value(cont_author.set_telephone, "phone", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        def set_el_fax(cont_author, contact_auth_in, parent_req):
                            """
                            XSD: <xs:element name="fax" type="telephone_number_type">
                            CIF: _pdbx_contact_author.fax
                            """
                            set_cif_value(cont_author.set_fax, "fax", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in, parent_el_req=parent_req)

                        # element 1
                        set_el_role(cont_author, contact_auth_in, parent_req)
                        # element 2
                        set_el_title(cont_author, contact_auth_in, parent_req)
                        # element 3
                        set_el_first_name(cont_author, contact_auth_in, parent_req)
                        # element 4
                        set_el_middle_name(cont_author, contact_auth_in, parent_req)
                        # element 5
                        set_el_last_name(cont_author, contact_auth_in, parent_req)
                        # element 6
                        set_el_organization(cont_author, contact_auth_in, parent_req)
                        # element 7
                        set_el_street(cont_author, contact_auth_in, parent_req)
                        # element 8
                        set_el_town_or_city(cont_author, contact_auth_in, parent_req)
                        # element 9
                        set_el_state_or_province(cont_author, contact_auth_in, parent_req)
                        # element 10
                        set_el_country(cont_author, contact_auth_in, parent_req)
                        # element 11
                        set_el_post_or_zip_code(cont_author, contact_auth_in, parent_req)
                        # element 12
                        set_el_email(cont_author, contact_auth_in, parent_req)
                        # element 13
                        set_el_telephone(cont_author, contact_auth_in, parent_req)
                        # element 14
                        set_el_fax(cont_author, contact_auth_in, parent_req)

                    def set_attr_private(cont_author):
                        """
                        XSD: <xs:attribute fixed="true" name="private" use="required"/>
                        """
                        cont_author.set_private("true")

                    # base
                    set_contact_details_type(cont_author, contact_auth_in, False)
                    # attribute 1
                    set_attr_private(cont_author)

                el_role = get_cif_value("role", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_title = get_cif_value("name_salutation", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_first_name = get_cif_value("name_first", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_middle_name = get_cif_value("name_mi", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_last_name = get_cif_value("name_last", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_organization = get_cif_value("organization_type", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_address_1 = get_cif_value("address_1", const.PDBX_CONTACT_AUTHOR, contact_auth_in)
                el_address_2 = get_cif_value("address_2", const.PDBX_CONTACT_AUTHOR, contact_auth_in)
                el_address_3 = get_cif_value("address_3", const.PDBX_CONTACT_AUTHOR, contact_auth_in)
                el_town_or_city = get_cif_value("city", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_state_or_province = get_cif_value("state_province", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_country = get_cif_value("country", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_post_or_zip_code = get_cif_value("postal_code", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_email = get_cif_value("email", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_telephone = get_cif_value("phone", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                el_fax = get_cif_value("fax", const.PDBX_CONTACT_AUTHOR, cif_list=contact_auth_in)
                contact_author_type_list = [
                    el_role,
                    el_title,
                    el_first_name,
                    el_middle_name,
                    el_last_name,
                    el_organization,
                    el_address_1,
                    el_address_2,
                    el_address_3,
                    el_town_or_city,
                    el_state_or_province,
                    el_country,
                    el_post_or_zip_code,
                    el_email,
                    el_telephone,
                    el_fax,
                ]
                if any(x is not None for x in contact_author_type_list):
                    cont_author = emdb.contact_authorType()
                    set_contact_author_type(cont_author, contact_auth_in)
                    if cont_author.has__content():
                        admin.add_contact_author(cont_author)

            def set_el_title(admin):
                """
                XSD: <xs:element name="title" type="xs:token">
                The value is _struct.title if _em_depui.same_title_as_pdb.
                Otherwise it is _em_admin.title
                CIF: _em_depui.same_title_as_pdb  YES/NO
                    YES: CIF: _struct.title
                    NO: CIF: _em_admin.title
                """
                set_cif_value(admin.set_title, "title", const.EM_ADMIN)

            def set_el_authors_list(admin):
                """
                XSD: <xs:element name="authors_list_type"> is a sequence of elements <author>
                CIF: _em_depui.same_authors_as_pdb YES/NO
                    YES: CIF: _em_author_list
                        _em_author_list.ordinal
                        _em_author_list.author
                        _em_author_list.identifier_ORCID
                        1 'Test, T.'     0000-0002-5251-4674
                        2 'Benton, D.J.' 0000-0001-6748-9339
                    NO or ./?: _audit_author
                        _audit_author.address
                        _audit_author.name
                        _audit_author.pdbx_ordinal
                        _audit_author.identifier_ORCID
                        ? 'First, A.'           1  0000-0002-5251-4674
                        ? 'Second-Second, B.' 2  0000-0001-6748-9339
                """

                def set_authors_list_type(authors_list, authors_in):
                    """
                    XSD: <xs:element name="authors_list"> has
                     ... 1 element of author_ORCID_type
                    """

                    def set_author_orcid_type(author_with_ORCID, auth_in):
                        """
                        XSD: <xs:complexType name="author_ORCID_type"> extends author_type and hashttps://rcsbpdb.atlassian.net/browse/DAOTHER-2725has
                        ... 1 attribute
                        XSD: <xs:attribute name="ORCID" type="ORCID_type"/>
                            CIF: _audit_author.identifier_ORCID  ? 'First, A.'           1  0000-0002-5251-4674
                            CIF: _em_author_list.ordinal            1
                                 _em_author_list.author             'Turner, J.'
                                 _em_author_list.identifier_ORCID   0000-0002-5251-4674
                        """

                        set_cif_value(author_with_ORCID.set_ORCID, "identifier_ORCID", const.EM_AUTHOR_LIST, cif_list=auth_in)
                        author = get_cif_value("author", const.EM_AUTHOR_LIST, cif_list=auth_in)

                        fmt_auth = format_author(author)
                        if fmt_auth != "":
                            author_with_ORCID.set_valueOf_(fmt_auth)
                        else:
                            txt = u"Author (%s) is not added to the list of authors as the format is wrong." % author
                            self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                            self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)

                    for _auth_id, auth_in in authors_in.items():
                        author_with_orcid = emdb.author_ORCID_type()
                        set_author_orcid_type(author_with_orcid, auth_in)
                        authors_list.add_author(author_with_orcid)

                authors_in = make_dict(const.EM_AUTHOR_LIST, "ordinal", 2)

                authors_list = emdb.authors_listType()
                set_authors_list_type(authors_list, authors_in)
                admin.set_authors_list(authors_list)

            def set_el_details():
                """
                EMDB administration details
                XSD: <xs:element name="details" type="xs:token" minOccurs="0">
                CIF: _em_admin.details
                Deprecated (2014-10-21)
                """

            def set_el_keywords(admin, key_words):
                """
                XSD: <name="z" type="xs:string" xs:element minOccurs="0">
                DEPRECATED 2014-10-21
                """
                set_cif_value(admin.set_keywords, "text", const.STRUCT_KEYWORDS, cif_list=key_words)

            def set_el_replace_existing_entry():
                """
                XSD: <xs:element name="replace_existing_entry" type="xs:boolean">
                CIF: _em_admin.replace_existing_entry_flag NO
                DEPRECATED (2014-10-21)
                """

            # attribute 1
            set_attr_composite_map(admin)
            # element 1
            set_el_status_history_list()
            # element 2
            set_el_current_status(admin)
            # element 3
            set_el_revision_history(admin)
            # element 4
            set_el_sites(admin)
            # element 5
            set_el_key_dates(admin)
            # element 6
            set_el_obsolete_list(admin)
            # element 7
            set_el_superseded_by_list(admin)
            # element 8
            aud_sup_in = make_dict(const.PDBX_AUDIT_SUPPORT, "ordinal", 2)
            set_el_grant_support(admin, aud_sup_in)
            # element 9
            # set_el_microscopy_center(admin)
            # element 10
            if self.__show_private:
                contact_auth_in = make_dict(const.PDBX_CONTACT_AUTHOR, "id")
                set_el_contact_author(admin, contact_auth_in)
            # element 11
            set_el_title(admin)
            # element 12
            set_el_authors_list(admin)
            # element 13
            set_el_details()
            # element 14
            keywords_in = make_dict(const.STRUCT_KEYWORDS, "entry_id")
            for key_words in keywords_in.values():
                pdbx_keywords = get_cif_value("text", const.STRUCT_KEYWORDS, cif_list=key_words)
                if pdbx_keywords is not None:
                    set_el_keywords(admin, key_words)
            # element 15
            set_el_replace_existing_entry()

        def set_crossreferences_type(cross_references):
            """
            Sets <xs:complexType name="crossreferences_type"> as
            ...a sequence of 4 elements
            elements 2 and 3 require CIF value for _em_db_reference.db_name
            """

            def set_el_citation_list(cross_references):
                """
                XSD: <xs:element name="citation_list">
                """

                def set_citation_list_type(citation_list):
                    """
                    XSD: <xs:element name="citation_list"> has 2 elements:
                        1. <xs:element name="primary_citation">
                        2. <xs:element name="secondary_citation" maxOccurs="unbounded" minOccurs="0">
                        CIF: _citation
                        CIF: _citation_author
                        _citation_author.citation_id:
                            - for the primary citation: 'primary'
                            - for the secondary citations can be anything.
                              However, they are usually numbered
                    """

                    def set_el_external_references(pub, cite_in, cite_ref_type_list):
                        """
                        XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                        CIF: _citation.pdbx_database_id_PubMed   ?
                        CIF: _citation.pdbx_database_id_DOI      ?
                        CIF: _citation.book_id_ISBN              ?
                        CIF: _citation.journal_id_ISSN           ?
                        CIF: _citation.abstract_id_CAS           ?
                        CIF: _citation.journal_id_CSD            0353
                        CIF: _citation.database_id_Medline       ?
                        CIF: _citation.journal_id_ASTM           ?
                        Citation references
                        """
                        init_txt = u"Citation references:"
                        refs = []
                        for cite_ref_type in cite_ref_type_list:
                            cite_ref_value = get_cif_value(cite_ref_type[0], const.CITATION, cite_in)
                            # cite_ref_item = get_cif_item(cite_ref_type[0], const.CITATION)
                            if cite_ref_value is not None:
                                ext_ref = emdb.external_referencesType()
                                ext_ref.set_type(cite_ref_type[1])
                                if cite_ref_type[1] == "DOI":  # 'doi':
                                    cite_ref_value = "doi:%s" % cite_ref_value
                                ext_ref.set_valueOf_(cite_ref_value)
                                refs.append(cite_ref_type[1])
                                pub.add_external_references(ext_ref)
                        if refs is not []:
                            ref_txt = ", ".join(refs)
                            txt = init_txt + ref_txt
                            self.current_entry_log.info_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.info_title + txt))
                            self.log_formatted(self.info_log_string, const.INFO_ALERT + txt)

                    def create_auth_dict(auth_dict):
                        """
                        Creates a dictionary of authors where key is citation id
                        CIF: _citation_author.citation_id
                        CIF: _citation_author.name
                        CIF: _citation_author.ordinal
                        primary  1  'Fitzgerald, P.M.D.'
                        primary  2  'McKeever, B.M.'
                        2        1  'Navia, M.A.'
                        """
                        cite_authors_in = self.cif.get(const.CITATION_AUTHOR)
                        if cite_authors_in != []:
                            for cite_author in cite_authors_in:
                                auth_id = get_cif_value("citation_id", const.CITATION_AUTHOR, cite_author)
                                if auth_id not in auth_dict:
                                    auth_dict[auth_id] = []
                                # add the author and their ordinal id
                                auth_name = get_cif_value("name", const.CITATION_AUTHOR, cite_author)
                                auth_ordinal_id = get_cif_value("ordinal", const.CITATION_AUTHOR, cite_author)
                                auth_orcid = get_cif_value("identifier_ORCID", const.CITATION_AUTHOR, cite_author)
                                auth_dict[auth_id].append((auth_name, auth_ordinal_id, auth_orcid))
                            # Sort the author lists according to ordinal
                            for auth_id in auth_dict:
                                auth_dict[auth_id].sort(key=lambda item: int(item[1]))
                        else:
                            txt = u"CIF category (%s) missing." % const.CITATION_AUTHOR
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, "(" + self.entry_in_translation_log.id + ")" + const.REQUIRED_ALERT + txt)

                    def set_el_author(pub, cite_id_in, auth_dict):
                        """
                        XSD: <xs:element name="author" type="author_order_type" maxOccurs="unbounded"/>
                        Citation authors
                        """
                        if cite_id_in in auth_dict:
                            if len(auth_dict[cite_id_in]) > 0:
                                for auth_in in auth_dict[cite_id_in]:
                                    author = emdb.author_order_type(valueOf_=format_author(auth_in[0]), ORCID=auth_in[2], order=int(auth_in[1]))
                                    if author.has__content():
                                        pub.add_author(author)
                            else:
                                txt = u"No authors for citation id (%s) found. At least one is required." % cite_id_in
                                self.current_entry_log.error_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                )
                                self.log_formatted(self.error_log_string, "(" + self.entry_in_translation_log.id + ")" + const.REQUIRED_ALERT + txt)

                    def set_journal_citation(jrnl, cite_in, cite_id_in, auth_dict, cite_ref_type_list, jrnl_abbrev_in):
                        """
                        <xs:element name="journal_citation" substitutionGroup="citation_type"> has
                        .. 1 attribute and
                        .. a sequence of 13 elements.
                        """

                        def set_att_published(jrnl, jrnl_abbrev_in):
                            """
                            <xs:attribute name="published" type="xs:boolean" use="required"/>
                            No proper flag to distinguish between journal and book depositions.
                            Use the journal abbreviation as an implicit flag
                            """
                            # ja_item = get_cif_item('journal_abbrev', const.CITATION)
                            value_given = None
                            if jrnl_abbrev_in.lower() in ["to be published", "suppressed"]:
                                value_given = False
                            else:
                                value_given = True
                            if value_given is not None:
                                jrnl.set_published(value_given)
                                txt = u"The value (%s) is given to (journal_citation.set_published)." % value_given
                                self.current_entry_log.info_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.info_title + txt))
                                self.log_formatted(self.info_log_string, const.INFO_ALERT + txt)

                        def set_el_title(jrnl, cite_in):
                            """
                            XSD: <xs:element name="title" type="xs:token"/>
                            CIF: _citation.title
                            """
                            set_cif_value(jrnl.set_title, "title", const.CITATION, cif_list=cite_in)

                        def set_el_journal(jrnl, cite_in):
                            """
                            XSD: <xs:element name="journal" type="xs:token" minOccurs="0"/>
                            CIF: _citation.journal_full
                            """
                            set_cif_value(jrnl.set_journal, "journal_full", const.CITATION, cif_list=cite_in)

                        def set_el_journal_abbreviation(jrnl, cite_in):
                            """
                            XSD: <xs:element name="journal_abbreviation" type="xs:token">
                            CIF: _citation.journal_abbrev 'To Be Published'
                            _citation.journal_abbrev is not required in mmCIF dict!!! It is in XSD
                            """
                            set_cif_value(jrnl.set_journal_abbreviation, "journal_abbrev", const.CITATION, cif_list=cite_in)

                        def set_el_country(jrnl, cite_in):
                            """
                            XSD: <xs:element name="country" type="xs:token" minOccurs="0"/>
                            CIF: _citation.country ?
                            """
                            set_cif_value(jrnl.set_country, "country", const.CITATION, cif_list=cite_in)

                        def set_el_issue(jrnl, cite_in):
                            """
                            XSD: <xs:element name="issue" type="xs:positiveInteger"
                                  minOccurs="0"/>
                            CIF: _citation.journal_issue ?
                            """
                            set_cif_value(jrnl.set_issue, "journal_issue", const.CITATION, cif_list=cite_in)

                        def set_el_volume(jrnl, cite_in):
                            """
                            XSD: <xs:element name="volume" type="xs:string" nillable="true" minOccurs="0"/>
                            CIF: _citation.journal_volume ?
                            """
                            set_cif_value(jrnl.set_volume, "journal_volume", const.CITATION, cif_list=cite_in)

                        def set_el_first_page(jrnl, cite_in):
                            """
                            XSD: <xs:element name="first_page" type="page_type" nillable="false" minOccurs="0"/>
                            CIF: _citation.page_first ?
                            """
                            set_cif_value(jrnl.set_first_page, "page_first", const.CITATION, cif_list=cite_in)

                        def set_el_last_page(jrnl, cite_in):
                            """
                            XSD: <xs:element name="last_page" type="page_type" minOccurs="0"/>
                            CIF: _citation.page_last ?
                            """
                            set_cif_value(jrnl.set_last_page, "page_last", const.CITATION, cif_list=cite_in)

                        def set_el_year(jrnl, cite_in):
                            """
                            XSD: <xs:element name="year" minOccurs="0">
                            CIF: _citation.year ?
                            """
                            set_cif_value(jrnl.set_year, "year", const.CITATION, cif_list=cite_in)

                        def set_el_language(jrnl, cite_in):
                            """
                            XSD: <xs:element name="language" type="xs:language" minOccurs="0"/>
                            CIF: _citation.language ?
                            """
                            set_cif_value(jrnl.set_language, "language", const.CITATION, cif_list=cite_in)

                        def set_el_details(jrnl, cite_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _citation.details ?
                            """
                            set_cif_value(jrnl.set_details, "details", const.CITATION, cif_list=cite_in)

                        jrnl.original_tagname_ = "journal_citation"
                        # attribute 1
                        set_att_published(jrnl, jrnl_abbrev_in)
                        # element 1
                        set_el_author(jrnl, cite_id_in, auth_dict)
                        # element 2
                        set_el_title(jrnl, cite_in)
                        # element 3
                        set_el_journal(jrnl, cite_in)
                        # element 4
                        set_el_journal_abbreviation(jrnl, cite_in)
                        # element 5
                        set_el_country(jrnl, cite_in)
                        # element 6
                        set_el_issue(jrnl, cite_in)
                        # element 7
                        set_el_volume(jrnl, cite_in)
                        # element 8
                        set_el_first_page(jrnl, cite_in)
                        # element 9
                        set_el_last_page(jrnl, cite_in)
                        # element 10
                        set_el_year(jrnl, cite_in)
                        # element 11
                        set_el_language(jrnl, cite_in)
                        # element 12
                        set_el_external_references(jrnl, cite_in, cite_ref_type_list)
                        # element 13
                        set_el_details(jrnl, cite_in)

                    def set_non_journal_citation(non_jrnl, cite_in, cite_id_in, auth_dict, cite_ref_type_list):
                        """
                        <xs:element name="non_journal_citation" substitutionGroup="citation_type"> has
                        ..1 attribute and
                        ..a sequence of 14 elements
                        """

                        def set_attr_published(non_jrnl, cite_in):
                            """
                            XSD: <xs:attribute name="published" type="xs:boolean" use="required"/>
                            CIF: _citation.unpublished_flag ?
                            """
                            unpublished_flag = get_cif_value("unpublished_flag", const.CITATION, cite_in)
                            # uf_item = get_cif_item('unpublished_flag', const.CITATION)
                            value_given = None
                            if unpublished_flag is None or unpublished_flag == "Y":
                                value_given = False
                            else:
                                value_given = True
                            if value_given is not None:
                                non_jrnl.set_published(value_given)
                                txt = u"The value (%s) is given to (non_journal_citation.set_published)." % value_given
                                self.current_entry_log.info_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.info_title + txt))
                                self.log_formatted(self.info_log_string, const.INFO_ALERT + txt)

                        def set_el_editor():
                            """
                            XSD: <xs:element name="editor" type="author_order_type" maxOccurs="unbounded" minOccurs="0"/>
                            CIF: TO BE IMPLEMENTED: add editor???
                            """

                        def set_el_non_journal_citation(non_jrnl, cite_in):
                            """
                            XSD: <xs:element name="title" type="xs:token">
                            CIF: _citation.book_title ?
                            """
                            set_cif_value(non_jrnl.set_title, "book_title", const.CITATION, cif_list=cite_in)

                        def set_el_thesis_title():
                            """
                            XSD: <xs:element name="thesis_title" type="xs:token" minOccurs="0">
                            Deprecated (2014-10-21)
                            """

                        def set_el_chapter_title():
                            """
                            XSD: <xs:element name="chapter_title" type="xs:token" minOccurs="0"/>
                            CIF: TO BE IMPLEMENTED
                            """

                        def set_el_volume():
                            """
                            XSD: <xs:element name="volume" type="xs:string" minOccurs="0"/>
                            CIF: TO BE IMPLEMENTED:
                            """

                        def set_el_publisher(non_jrnl, cite_in):
                            """
                            XSD: <xs:element name="publisher" type="xs:token" minOccurs="0"/>
                            CIF: _citation.book_publisher ?
                            """
                            set_cif_value(non_jrnl.set_publisher, "book_publisher", const.CITATION, cif_list=cite_in)

                        def set_el_publisher_location(non_jrnl, cite_in):
                            """
                            XSD: <xs:element name="publisher_location" type="xs:token" minOccurs="0">
                            CIF: _citation.book_publisher_city ?
                            CIF: _citation.country ?
                            Add city and country together
                            """
                            city = get_cif_value("book_publisher_city", const.CITATION, cite_in)
                            country = get_cif_value("country", const.CITATION, cite_in)
                            location = ""
                            if city is not None:
                                if country is not None:
                                    location = " ".join([city, country])
                                else:
                                    location = city
                            else:
                                if country is not None:
                                    location = country

                            set_cif_value(non_jrnl.set_publisher_location, "country", const.CITATION, cif_value=location)

                        def set_el_first_page(non_jrnl, cite_in):
                            """
                            XSD: <xs:element name="first_page" type="page_type" minOccurs="0"/>
                            CIF: _citation.page_first ?
                            """
                            set_cif_value(non_jrnl.set_first_page, "page_first", const.CITATION, cif_list=cite_in, fmt=const.CITATION)

                        def set_el_last_page(non_jrnl, cite_in):
                            """
                            XSD: <xs:element name="last_page" type="page_type" minOccurs="0"/>
                            CIF: _citation.page_last ?
                            """
                            set_cif_value(non_jrnl.set_last_page, "page_last", const.CITATION, cif_list=cite_in)

                        def set_el_year(non_jrnl, cite_in):
                            """
                            XSD: <xs:element name="year">
                            CIF: _citation.year ?
                            """
                            set_cif_value(non_jrnl.set_year, "year", const.CITATION, cif_list=cite_in)

                        def set_el_language(non_jrnl, cite_in):
                            """
                            XSD: <xs:element name="language" type="xs:language" minOccurs="0"/>
                            CIF: _citation.language ?
                            """
                            set_cif_value(non_jrnl.set_language, "language", const.CITATION, cif_list=cite_in)

                        def set_el_details(non_jrnl, cite_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _citation.details ?
                            """
                            set_cif_value(non_jrnl.set_details, "details", const.CITATION, cif_list=cite_in)

                        non_jrnl.original_tagname_ = "non_journal_citation"
                        # attribute 1
                        set_attr_published(non_jrnl, cite_in)
                        # element 1
                        set_el_author(non_jrnl, cite_id_in, auth_dict)
                        # element 2
                        set_el_editor()
                        # element 3
                        set_el_non_journal_citation(non_jrnl, cite_in)
                        # element 4
                        set_el_thesis_title()
                        # element 5
                        set_el_chapter_title()
                        # element 6
                        set_el_volume()
                        # element 7
                        set_el_publisher(non_jrnl, cite_in)
                        # element 8
                        set_el_publisher_location(non_jrnl, cite_in)
                        # element 9
                        set_el_first_page(non_jrnl, cite_in)
                        # element 10
                        set_el_last_page(non_jrnl, cite_in)
                        # element 11
                        set_el_year(non_jrnl, cite_in)
                        # element 12
                        set_el_language(non_jrnl, cite_in)
                        # element 13
                        set_el_external_references(non_jrnl, cite_in, cite_ref_type_list)
                        # element 14
                        set_el_details(non_jrnl, cite_in)

                    auth_dict = {}
                    create_auth_dict(auth_dict)
                    # Citations
                    cite_ref_type_list = [
                        ["pdbx_database_id_PubMed", "PUBMED"],
                        ["pdbx_database_id_DOI", "DOI"],
                        ["book_id_ISBN", "ISBN"],
                        ["journal_id_ISSN", "ISSN"],
                        ["abstract_id_CAS", "CAS"],
                        ["journal_id_CSD", "CSD"],
                        ["database_id_Medline", "MEDLINE"],
                        ["journal_id_ASTM", "ASTM"],
                    ]
                    # get all citation values from cif
                    cite_list_in = self.cif.get(const.CITATION)
                    if cite_list_in != []:
                        # set a flag for primary citations
                        any_primary_citations = False
                        # get all citation author values from cif
                        for cite_in in cite_list_in:
                            # CIF:_citation.id
                            cite_id_in = get_cif_value(const.K_ID, const.CITATION, cite_in)
                            citation = None
                            if cite_id_in == "primary":
                                if not any_primary_citations:
                                    any_primary_citations = True
                                citation = emdb.primary_citationType()
                                citation_list.set_primary_citation(citation)
                            else:
                                citation = emdb.secondary_citationType()
                                citation_list.add_secondary_citation(citation)
                            if citation is not None:
                                # Is this a book (non-journal) or journal citation?
                                # These are the two citation_types
                                jrnl_abbrev_in = get_cif_value("journal_abbrev", const.CITATION, cite_in)
                                if jrnl_abbrev_in is not None:
                                    jrnl = emdb.journal_citation()
                                    set_journal_citation(jrnl, cite_in, cite_id_in, auth_dict, cite_ref_type_list, jrnl_abbrev_in)
                                    citation.set_citation_type(jrnl)
                                else:
                                    non_jrnl = emdb.non_journal_citation()
                                    set_non_journal_citation(non_jrnl, cite_in, cite_id_in, auth_dict, cite_ref_type_list)
                                    citation.set_citation_type(non_jrnl)
                            else:
                                txt = u"Citations cannot be set. The value for (_citation_author.citation_id) should be either (primary) or a number."
                                self.current_entry_log.error_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                )
                                self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                        if not any_primary_citations:
                            txt = u"No primary citations given."
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                    else:
                        txt = u"CIF category (%s) missing." % const.CITATION
                        self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                citation_list = emdb.citation_listType()
                set_citation_list_type(citation_list)
                cross_references.set_citation_list(citation_list)

            def set_el_emdb_list(cross_references, x_ref_dict_in):
                """
                <xs:element name="emdb_list"
                    type="emdb_cross_reference_list_type"
                    minOccurs="0">
                """

                def set_emdb_cross_ref_list_type(emdb_ref_list, x_ref_dict_in):
                    """
                    XSD: <xs:element name="emdb_list" type="emdb_cross_reference_list_type" minOccurs="0">
                    ..has only one element of
                    <xs:element name="emdb_reference" type="emdb_cross_reference_type" maxOccurs="unbounded"/>
                    ..is a sequence of 3 elements

                    EMDB entries related to this entry
                    Related EMDB and PDB entries - for some reason the entries
                    up to now have the EMDB related in PDBX_DATABASE_RELATED
                    CIF:     EMDB . EMD-4134 'other EM volume'
                    CIF: _pdbx_database_related
                    """

                    def set_ref_el_emdb_id(emdb_ref, emdb_ref_in):
                        """
                        XSD: <xs:element name="emdb_id" type="emdb_id_type"/>
                        CIF: _em_db_reference.access_code (MANDATORY)
                        """
                        acc_code = get_cif_value("access_code", const.EM_DB_REFERENCE, emdb_ref_in)
                        if acc_code is not None:
                            acc_code = acc_code.strip()
                            acc_code_match = re.match(const.CIF_EMD_ID_RE, acc_code)
                            if acc_code_match is not None:  # a match
                                set_cif_value(emdb_ref.set_emdb_id, "access_code", const.EM_DB_REFERENCE, cif_list=emdb_ref_in)
                            else:
                                # acc_code is not in a format of EMD-xxxx...x
                                corrected = False
                                txt = (
                                    u"The value for (_em_db_reference.access_code) is (%s) and it is in a wrong format. If a new value is given the message follows." % acc_code
                                )
                                self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                                self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)
                                if (len(acc_code) == 4 or len(acc_code) == 5) and acc_code.isdigit():
                                    # acc_code is a 4 or 5-digit number - add EMD- to it
                                    correct_acc_code = "EMD-" + acc_code
                                    corrected = True
                                if acc_code.find("EMDB") != -1:
                                    # acc_code is in a EMDB-xxxx format - remove B
                                    correct_acc_code = acc_code.replace("B", "")
                                    corrected = True
                                if corrected:
                                    set_cif_value(emdb_ref.set_emdb_id, "access_code", const.EM_DB_REFERENCE, cif_list=emdb_ref_in, cif_value=correct_acc_code)  # pylint: disable=possibly-used-before-assignment
                                    txt = u"emdb_id is set to (%s) as (_em_db_reference.access_code) is (%s)." % (correct_acc_code, acc_code)
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                        else:
                            txt = u"Cannot set crossreference emdb_id as the required value for (_em_db_reference.access_code) is not given."
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                    def set_ref_el_relationship(emdb_ref, emdb_ref_in):
                        """
                        XSD: <xs:element name="relationship" minOccurs="0">
                        CIF: _em_db_reference.relationship ?
                        """
                        rel_in = get_cif_value("relationship", const.EM_DB_REFERENCE, emdb_ref_in)
                        # rel_item = get_cif_item('relationship', const.EM_DB_REFERENCE)
                        txt = None
                        if rel_in == "IN FRAME":
                            emdb_ref.set_relationship(emdb.relationshipType(in_frame="FULLOVERLAP"))
                            txt = u"The value (FULLOVERLAP) is given to (emdb_ref.set_relationship)."
                        else:
                            emdb_ref.set_relationship(emdb.relationshipType(other="unknown"))
                            txt = u"The value (unknown) is given to (emdb_ref.set_relationship)."
                        if txt is not None:
                            self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                            self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)

                    def set_ref_el_details(emdb_ref, emdb_ref_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _em_db_reference.details ?
                        """
                        set_cif_value(emdb_ref.set_details, "details", const.EM_DB_REFERENCE, cif_list=emdb_ref_in)

                    def set_rel_el_emdb_id(cross_ref, rel_in):
                        """
                        XSD: <xs:element name="emdb_id" type="emdb_id_type"/>
                        CIF: _pdbx_database_related.db_id
                        """
                        db_id = get_cif_value("db_id", const.PDBX_DATABASE_RELATED, rel_in)
                        if db_id is not None:
                            db_id_match = re.match(const.CIF_EMD_ID_RE, db_id)
                            if db_id_match is not None:  # a match
                                set_cif_value(cross_ref.set_emdb_id, "db_id", const.PDBX_DATABASE_RELATED, cif_list=rel_in)
                            else:
                                # db_id is not in a format of EMD-xxxx
                                corrected = False
                                txt = u"The value for (_pdbx_database_related.db_id) (%s) is in a wrong format. If a new value is given the message follows." % db_id
                                self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                                self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)
                                if (len(db_id) == 4 or len(db_id) == 5) and db_id.isdigit():
                                    # db_id is a 4 or 5-digit number - add EMD- to it
                                    correct_db_id = "EMD-" + db_id
                                    corrected = True
                                if db_id.find("EMDB") != -1:
                                    # db_id is in a EMDB-xxxx...x format - remove B
                                    correct_db_id = db_id.replace("B", "")
                                    corrected = True
                                if db_id.find("D_") != -1:
                                    # db_id is given as e.g. D_1000232117; should be EMD-xxxx
                                    txt = u"emdb_id cannot be set as the value for (_pdbx_database_related.db_id) is given as (%s)." % db_id
                                    self.current_entry_log.error_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                    )
                                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                                if corrected:
                                    set_cif_value(cross_ref.set_emdb_id, "db_id", const.PDBX_DATABASE_RELATED, cif_list=rel_in, cif_value=correct_db_id)  # pylint: disable=possibly-used-before-assignment
                                    txt = u"emdb_id is set to (%s) as the value for (_pdbx_database_related.db_id) is given as (%s)." % (correct_db_id, db_id)
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                        else:
                            set_cif_value(cross_ref.set_emdb_id, "db_id", const.PDBX_DATABASE_RELATED, cif_list=rel_in)

                    def set_rel_el_relationship(cross_ref, rel_in):
                        """
                        XSD: <xs:element name="relationship" minOccurs="0">
                        CIF: _pdbx_database_related.content_type
                        """
                        content_type = get_cif_value("content_type", const.PDBX_DATABASE_RELATED, rel_in)
                        txt = None
                        if content_type == const.CIF_EMDB_ASSOC:
                            cross_ref.set_relationship(emdb.relationshipType(other=const.CIF_EMDB_ASSOC))
                            txt = u"The value (%s) is given to (cross_ref.set_relationship)." % const.CIF_EMDB_ASSOC
                        elif content_type == const.CIF_EMDB_OTHER:
                            cross_ref.set_relationship(emdb.relationshipType(other=const.CIF_EMDB_OTHER))
                            txt = u"The value (%s) is given to (cross_ref.set_relationship)." % const.CIF_EMDB_OTHER
                        elif content_type == const.CIF_EMDB_CONSENSUS:
                            cross_ref.set_relationship(emdb.relationshipType(other=const.CIF_EMDB_CONSENSUS))
                            txt = u"The value (%s) is given to (cross_ref.set_relationship)." % const.CIF_EMDB_CONSENSUS
                        elif content_type == const.CIF_EMDB_FOCUSED:
                            cross_ref.set_relationship(emdb.relationshipType(other=const.CIF_EMDB_FOCUSED))
                            txt = u"The value (%s) is given to (cross_ref.set_relationship)." % const.CIF_EMDB_FOCUSED
                        else:
                            txt = u"No value is given to (cross_ref.set_relationship)."
                        if txt is not None:
                            self.current_entry_log.info_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.info_title + txt))
                            self.log_formatted(self.info_log_string, const.INFO_ALERT + txt)

                    def set_rel_el_details(cross_ref, rel_in):
                        """
                        XSD: <xs:element name="details" minOccurs="0" type="xs:string"/>
                        CIF: _pdbx_database_related.details
                        """
                        set_cif_value(cross_ref.set_details, "details", const.PDBX_DATABASE_RELATED, cif_list=rel_in)

                    rel_entries_dict_in = make_list_of_dicts(const.PDBX_DATABASE_RELATED, "db_name")
                    if "EMDB" in x_ref_dict_in:
                        emdb_ref_list_in = x_ref_dict_in["EMDB"]
                        for emdb_ref_in in emdb_ref_list_in:
                            emdb_ref = emdb.emdb_cross_reference_type()
                            set_ref_el_emdb_id(emdb_ref, emdb_ref_in)
                            set_ref_el_relationship(emdb_ref, emdb_ref_in)
                            set_ref_el_details(emdb_ref, emdb_ref_in)

                            if emdb_ref.has__content():
                                emdb_ref_list.add_emdb_reference(emdb_ref)

                    if "EMDB" in rel_entries_dict_in:
                        emdb_rel_list_in = rel_entries_dict_in["EMDB"]
                        for rel_in in emdb_rel_list_in:
                            em_id = get_cif_value("db_id", const.PDBX_DATABASE_RELATED, rel_in)
                            db2_in = assert_get_value(const.DATABASE_2, self.cif)
                            dict_db2_in = {t[0]: t[1] for t in db2_in}
                            emdb_id = dict_db2_in[('_database_2.database_id', 'EMDB')][1]
                            if em_id != emdb_id:
                                cross_ref = emdb.emdb_cross_reference_type()
                                set_rel_el_emdb_id(cross_ref, rel_in)
                                set_rel_el_relationship(cross_ref, rel_in)
                                set_rel_el_details(cross_ref, rel_in)

                                if cross_ref.has__content():
                                    emdb_ref_list.add_emdb_reference(cross_ref)

                emdb_ref_list = emdb.emdb_cross_reference_list_type()
                set_emdb_cross_ref_list_type(emdb_ref_list, x_ref_dict_in)
                if emdb_ref_list.has__content():
                    cross_references.set_emdb_list(emdb_ref_list)

            def set_el_pdb_list(cross_references, x_ref_dict_in):
                """
                XSD: <xs:element name="pdb_list" type="pdb_cross_reference_list_type" minOccurs="0">
                """

                def set_pdb_cross_ref_list_type(pdb_ref_list, x_ref_dict_in):
                    """
                    XSD: <xs:element name="pdb_list" type="pdb_cross_reference_list_type" minOccurs="0">
                    ..has only one element of
                    <xs:element name="pdb_reference" type="pdb_cross_reference_type" maxOccurs="unbounded">
                    ..is a sequence of 3 elements
                    """

                    def set_el_pdb_id(pdb_ref, pdb_ref_in):
                        """
                        XSD: <xs:element name="emdb_id" type="emdb_id_type"/>
                        CIF: _em_db_reference.access_code
                        """
                        set_cif_value(pdb_ref.set_pdb_id, "access_code", const.EM_DB_REFERENCE, cif_list=pdb_ref_in, fmt=str.lower, parent_el_req=False)

                    def set_el_relationship(pdb_ref, pdb_ref_in):
                        """
                        XSD: <xs:element name="relationship" minOccurs="0">
                        CIF: _em_db_reference.relationship ?
                        """
                        rel_in = get_cif_value("relationship", const.EM_DB_REFERENCE, pdb_ref_in)
                        # rel_item = get_cif_item('relationship', const.EM_DB_REFERENCE)
                        txt = None
                        if rel_in == "IN FRAME":
                            pdb_ref.set_relationship(emdb.relationshipType(in_frame="FULLOVERLAP"))
                            txt = u"The value (FULLOVERLAP) is given to (pdb_ref.set_relationship)."
                        else:
                            pdb_ref.set_relationship(emdb.relationshipType(other="unknown"))
                            txt = u"The value (unknown) is given to (pdb_ref.set_relationship)."
                        if pdb_ref.has__content():
                            pdb_ref_list.add_pdb_reference(pdb_ref)
                            if txt is not None:
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)

                    def set_el_details(pdb_ref, pdb_ref_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _em_db_reference.details ?
                        """
                        set_cif_value(pdb_ref.set_details, "details", const.EM_DB_REFERENCE, cif_list=pdb_ref_in)

                    pdb_ref_list_in = x_ref_dict_in["PDB"]
                    for pdb_ref_in in pdb_ref_list_in:
                        el_emdb_id = get_cif_value("access_code", const.EM_DB_REFERENCE, cif_list=pdb_ref_in)
                        el_details = get_cif_value("details", const.EM_DB_REFERENCE)
                        if any(x is not None for x in [el_emdb_id, el_details]):
                            pdb_ref = emdb.pdb_cross_reference_type()
                            # element 1
                            set_el_pdb_id(pdb_ref, pdb_ref_in)
                            # element 2
                            set_el_relationship(pdb_ref, pdb_ref_in)
                            # element 3
                            set_el_details(pdb_ref, pdb_ref_in)

                if "PDB" in x_ref_dict_in:
                    pdb_ref_list = emdb.pdb_cross_reference_list_type()
                    set_pdb_cross_ref_list_type(pdb_ref_list, x_ref_dict_in)
                    if pdb_ref_list.has__content():
                        cross_references.set_pdb_list(pdb_ref_list)

            def set_el_other_db_list(cross_references, other_db_ref_dict_in):
                """
                <xs:element name="other_db_list"
                    type="other_cross_reference_list_type"
                    minOccurs="0">
                """

                def set_other_db_cross_ref_list_type(other_db_ref_list, other_db_ref_dict_in):
                    """
                    XSD: <xs:element name="other_db_list" type="other_db_cross_reference_list_type" minOccurs="0">
                    ..has only one element of
                    <xs:element name="db_reference" type="db_cross_reference_type" maxOccurs="unbounded"/>
                    ..is a sequence of 4 elements
                    """

                    def set_ref_el_db_name(other_db_ref, other_db_ref_in):
                        """
                        XSD: <xs:element name="db_name" type="token"/>
                        CIF: _pdbx_database_related.db_name (MANDATORY)
                        """
                        set_cif_value(other_db_ref.set_db_name, "db_name", const.PDBX_DATABASE_RELATED, cif_list=other_db_ref_in)

                    def set_ref_el_accession_id(other_db_ref, other_db_ref_in):
                        """
                        XSD: <xs:element name="accession_id" type="token"/>
                        CIF: _pdbx_database_related.db_id (MANDATORY)
                        """
                        set_cif_value(other_db_ref.set_accession_id, "db_id", const.PDBX_DATABASE_RELATED, cif_list=other_db_ref_in)

                    def set_ref_el_content_type(other_db_ref, other_db_ref_in):
                        """
                        XSD: <xs:element name="content_type" minOccurs="0" type="xs:string"/>
                        CIF: _pdbx_database_related.details
                        """
                        set_cif_value(other_db_ref.set_content_type, "content_type", const.PDBX_DATABASE_RELATED, cif_list=other_db_ref_in)

                    def set_ref_el_details(other_db_ref, other_db_ref_in):
                        """
                        XSD: <xs:element name="details" minOccurs="0" type="xs:string"/>
                        CIF: _pdbx_database_related.details
                        """
                        set_cif_value(other_db_ref.set_details, "details", const.PDBX_DATABASE_RELATED, cif_list=other_db_ref_in)

                    if "SASBDB" in other_db_ref_dict_in:
                        other_db_ref_list_in = other_db_ref_dict_in["SASBDB"]
                        for other_db_ref_in in other_db_ref_list_in:
                            other_db_ref = emdb.other_db_cross_reference_type()
                            set_ref_el_db_name(other_db_ref, other_db_ref_in)
                            set_ref_el_accession_id(other_db_ref, other_db_ref_in)
                            set_ref_el_content_type(other_db_ref, other_db_ref_in)
                            set_ref_el_details(other_db_ref, other_db_ref_in)

                            if other_db_ref.has__content():
                                other_db_ref_list.add_db_reference(other_db_ref)

                other_db_ref_list = emdb.other_db_cross_reference_list_type()
                set_other_db_cross_ref_list_type(other_db_ref_list, other_db_ref_dict_in)
                if other_db_ref_list.has__content():
                    cross_references.set_other_db_list(other_db_ref_list)

            def set_el_auxiliary_link_list(cross_references):
                """
                <xs:element name="auxiliary_link_list" minOccurs="0">
                ..a list of one and only one element
                """

                def set_aux_link_type(aux, aux_in):
                    """
                    <xs:complexType name="auxiliary_link_type">
                    ..a sequence of 3 elements
                    """

                    def set_el_type(aux, aux_in):
                        """
                        XSD: <xs:element name="type">
                        CIF: _em_db_reference_auxiliary.link_type
                        """
                        set_cif_value(aux.set_type, "type", const.EM_DB_REFERENCE_AUXILIARY, cif_list=aux_in)

                    def set_el_link(aux, aux_in):
                        """
                        XSD: <xs:element name="link">
                        CIF: _em_db_reference_auxiliary.link
                        """
                        set_cif_value(aux.set_link, "link", const.EM_DB_REFERENCE_AUXILIARY, cif_list=aux_in)

                    def set_el_details():
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: NOT YET IMPLEMENTED
                        """

                    # element 1
                    set_el_type(aux, aux_in)
                    # element 2
                    set_el_link(aux, aux_in)
                    # element 3
                    set_el_details()

                aux_link_list = emdb.auxiliary_link_listType()
                aux_link_list_in = self.cif.get(const.EM_DB_REFERENCE_AUXILIARY)
                for aux_in in aux_link_list_in:
                    aux = emdb.auxiliary_link_type()
                    set_aux_link_type(aux, aux_in)
                    if aux.has__content():
                        aux_link_list.add_auxiliary_link(aux)
                if aux_link_list.has__content():
                    cross_references.set_auxiliary_link_list(aux_link_list)

            # element 1
            set_el_citation_list(cross_references)
            x_ref_dict_in = make_list_of_dicts(const.EM_DB_REFERENCE, "db_name")
            # element 2
            set_el_emdb_list(cross_references, x_ref_dict_in)
            # element 3
            set_el_pdb_list(cross_references, x_ref_dict_in)
            # element 4
            other_db_ref_dict_in = make_list_of_dicts(const.PDBX_DATABASE_RELATED, "db_name")
            set_el_other_db_list(cross_references, other_db_ref_dict_in)
            # element 5
            set_el_auxiliary_link_list(cross_references)

        def set_sample_type(sample):
            """
            Sets <xs:element name="sample" type="sample_type"> as
            ..a sequence of 3 elements
            """

            def set_el_name():
                """
                XSD: <xs:element name="name" type="xs:token">
                CIF: _em_entity_assembly.name 'Israeli acute paralysis virus'
                This is set at the end as it is the name given to
                the supramolecule that has parent id = 0
                """

            def set_base_source_type(src, cif_category, src_in):
                """
                Parameters:
                @param src: either an object for supramolecule or macromolecule
                @param cif_category: _entity_src_nat, _entity_src_gen, _pdbx_entity_src_syn or _em_entity_assembly_naturalsource
                @param src_in:
                XSD: <xs:complexType name="base_source_type"> has
                    .. 2 attribute and
                    .. 3 elements
                """

                def set_attr_database(src, cif_category, src_in):
                    """
                    XSD: <xs:attribute name="database">; <xs:enumeration value="NCBI"/>
                    CIF: _entity_src_nat.pdbx_ncbi_taxonomy_id
                    CIF: _entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id
                    CIF: _pdbx_entity_src_syn.ncbi_taxonomy_id
                    CIF: _em_entity_assembly_naturalsource.ncbi_tax_id
                    """
                    tax_id_in = None
                    a_dict = {
                        const.ENTITY_SRC_NAT: "pdbx_ncbi_taxonomy_id",
                        const.ENTITY_SRC_GEN: "pdbx_gene_src_ncbi_taxonomy_id",
                        const.PDBX_ENTITY_SRC_SYN: "ncbi_taxonomy_id",
                        const.EM_ENTITY_ASSEMBLY_NATURALSOURCE: "ncbi_tax_id"
                    }
                    if cif_category is not None:
                        cif_key = a_dict.get(cif_category, None)
                        tax_id_in = get_cif_value(cif_key, cif_category, src_in)
                        if cif_key is None or cif_key.isspace():
                            txt = u"Cannot set the database attribute as the cif category (%s) is not one of: (%s)." % (cif_category, a_dict)
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                        else:
                            set_cif_value(src.set_database, cif_key, cif_category, cif_list=src_in, cif_value="NCBI")

                    return tax_id_in

                def set_attr_synthetically_produced(src_in):
                    """
                    XSD: <xs:attribute name="synthetically_produced">/>
                    CIF: _em_entity_assembly_source.source YES
                    """
                    source_dict = make_dict(const.EM_ENTITY_ASSEMBLY, "id")
                    for i in source_dict.keys():
                        assembly_type = get_cif_value("source", const.EM_ENTITY_ASSEMBLY, source_dict[i])
                        assembly_id = get_cif_value("id", const.EM_ENTITY_ASSEMBLY, source_dict[i])
                        if assembly_type == "SYNTHETIC":
                            nat_src_id = get_cif_value("entity_assembly_id", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, src_in)
                            if assembly_id == nat_src_id:
                                set_cif_value(src.set_synthetically_produced, cif_value=True)

                def set_el_organism(src, tax_id_in, cif_category, src_in):
                    """
                    XSD: <xs:element name="organism" type="organism_type">
                    CIF: _entity_src_nat.pdbx_organism_scientific
                    CIF: _entity_src_gen.pdbx_gene_src_scientific_name
                    CIF: _pdbx_entity_src_syn.organism_scientific
                    For some entries _entity_src_nat.common_name is not given but _entity_src_nat.pdbx_organism_scientific is
                    CIF: _entity_src_nat.pdbx_organism_scientific
                    CIF: _em_entity_assembly_naturalsource.organism
                    """
                    a_dict = {
                        const.ENTITY_SRC_NAT: "pdbx_organism_scientific",
                        const.ENTITY_SRC_GEN: "pdbx_gene_src_scientific_name",
                        const.PDBX_ENTITY_SRC_SYN: "organism_scientific",
                        const.EM_ENTITY_ASSEMBLY_NATURALSOURCE: "organism",
                    }
                    if cif_category is not None:
                        cif_key = a_dict.get(cif_category, None)
                        if cif_key is not None:
                            common_name = get_cif_value(cif_key, cif_category, cif_list=src_in)
                            if cif_category == const.ENTITY_SRC_NAT:
                                if common_name is None:
                                    cif_key = "pdbx_organism_scientific"
                                    org_sci_name = get_cif_value(cif_key, cif_category, cif_list=src_in)
                                    if org_sci_name is not None and not org_sci_name.isspace():
                                        txt = (
                                            u"The value for (_entity_src_nat.common_name) is not given so the value for (_entity_src_nat.pdbx_organism_scientific): (%s) is used."
                                            % org_sci_name
                                        )
                                        self.current_entry_log.warn_logs.append(
                                            self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                        )
                                        self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                                    else:
                                        org_sci_name = "Unspecified"
                                    if tax_id_in is not None and not tax_id_in.isspace():
                                        set_cif_value(
                                            src.set_organism,
                                            cif_key,
                                            cif_category,
                                            cif_list=src_in,
                                            constructor=emdb.organism_type,
                                            ncbi=tax_id_in,
                                            cif_value=org_sci_name,
                                            parent_el_req=False,
                                        )
                                    else:
                                        set_cif_value(
                                            src.set_organism, cif_key, cif_category, cif_list=src_in, constructor=emdb.organism_type, cif_value=org_sci_name, parent_el_req=False
                                        )
                                else:
                                    # common_name is not None
                                    if tax_id_in is not None:
                                        set_cif_value(
                                            src.set_organism,
                                            cif_key,
                                            cif_category,
                                            cif_list=src_in,
                                            constructor=emdb.organism_type,
                                            ncbi=tax_id_in,
                                            cif_value=common_name,
                                            parent_el_req=False,
                                        )
                                    else:
                                        set_cif_value(
                                            src.set_organism, cif_key, cif_category, cif_list=src_in, constructor=emdb.organism_type, cif_value=common_name, parent_el_req=False
                                        )
                            else:
                                if tax_id_in is not None:
                                    set_cif_value(src.set_organism, cif_key, cif_category, cif_list=src_in, constructor=emdb.organism_type, ncbi=tax_id_in, parent_el_req=False)
                                else:
                                    set_cif_value(src.set_organism, cif_key, cif_category, cif_list=src_in, constructor=emdb.organism_type, parent_el_req=False)
                        else:
                            txt = u"Cannot set the organism element as cif item (%s) is not of of (%s)." % (cif_key, a_dict)
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                def set_el_strain(src, cif_category, src_in):
                    """
                    XSD: <xs:element name="strain" type="xs:token" minOccurs="0"/>
                    XPath: /element(*,base_source_type)/strain
                    CIF: _entity_src_nat.strain
                    CIF: _entity_src_gen.gene_src_strain
                    CIF: _pdbx_entity_src_syn.strain
                    CIF: _em_entity_assembly_naturalsource.strain
                    """
                    a_dict = {const.ENTITY_SRC_NAT: "strain", const.ENTITY_SRC_GEN: "gene_src_strain", const.PDBX_ENTITY_SRC_SYN: "strain", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE: "strain"}
                    if cif_category is not None:
                        cif_key = a_dict.get(cif_category, None)
                        if cif_key is not None:
                            set_cif_value(src.set_strain, cif_key, cif_category, cif_list=src_in)
                        else:
                            txt = u"Cannot set the strain element as cif item (%s) is not one in (%s)." % (cif_key, a_dict)
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                def set_el_synonym_organism():
                    """
                    XSD: <xs:element name="synonym_organism" type="xs:token" minOccurs="0">
                    Deprecated (2014-10-21)
                    """

                def set_el_details(src, cif_category, src_in):
                    """
                    XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    XPath: /element(*,base_source_type)/details
                    CIF: _em_entity_assembly_naturalsource.details
                    """
                    det = get_cif_value("details", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, src_in)
                    if det is not None:
                        set_cif_value(src.set_details, "details", cif_category, cif_list=src_in)

                # attribute 1
                tax_id_in = set_attr_database(src, cif_category, src_in)
                # attribute 2
                set_attr_synthetically_produced(src_in)
                # element 1
                set_el_organism(src, tax_id_in, cif_category, src_in)
                # element 2
                set_el_strain(src, cif_category, src_in)
                # element 3
                set_el_synonym_organism()
                # element 4
                set_el_details(src, cif_category, src_in)

            def get_rec_exp_dict(ent_id_in, src_dicts, cif_cat_in, is_supramolecule):
                """
                Replace above
                Helper function: It returns a dictionary needed for recombinant expression containing a cif category depending on several requirements.

                SUPRAMOLECULE are all map entries. Cif category _em_entity_assembly:
                source molecule: em_entity_assembly_naturalsource - should always be present.
                if _em_entity_assembly.source = RECOMBINANT;  then get expression system from _em_entity_assembly_recombinan

                MACROMOLECULE

                from entity and its children.
                 If _entity.type is 'polymer' than _entity.src_method gives the return dictionary as:
                                1. nat; The return category is _entity_src_nat
                                2. man; The return category is _entity_src_gen
                                3. syn; The return category is _entity_src_syn

                if _entity_type != polymer then no source information will be present.

                Details for both SUPRAMOLECULE and MACROMOLECULE  are required by the depUI on separate pages so you will always get data for both (map+model) or just supramolecule (map only).

                Parameters:
                @param ent_id_in:
                @param src_dicts
                @param is_supramolecule: True for supramolecules; False for macromolecules
                """
                if is_supramolecule:
                    src = get_cif_value("source", const.EM_ENTITY_ASSEMBLY, cif_cat_in)
                    if src == "RECOMBINANT":
                        if src_dicts.get("rec_exp_dict_in", None) is None:
                            txt = u"(_em_entity_assembly_recombinan) category missing for creating recombinant expression for the map and model supramolecule entry (%s)." % ent_id_in
                            self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                            self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)
                        else:
                            return src_dicts["rec_exp_dict_in"]
                else:  # macromolecules
                    entity_dict = make_dict(const.ENTITY, "id")
                    if ent_id_in in entity_dict:
                        ent_type = get_cif_value("type", const.ENTITY, entity_dict[ent_id_in])
                        if ent_type == "polymer":
                            ent_src_method = get_cif_value("src_method", const.ENTITY, entity_dict[ent_id_in])
                            if ent_src_method == "nat":
                                pass  # CHECK THIS
                            elif ent_src_method == "man":
                                if src_dicts.get("ent_src_gen_dict", None) is None:
                                    txt = (
                                        u"(_entity_src_gen) category missing for creating recombinant expression for the map and model entry (%s) where (_em_entity_assembly.source) is (RECOMBINANT) and polymer is (man)."
                                        % ent_id_in
                                    )
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)
                                else:
                                    return src_dicts["ent_src_gen_dict"]
                            if ent_src_method == "syn":
                                if src_dicts.get("ent_src_syn_dict", None) is None:
                                    txt = (
                                        u"(_entity_src_syn) category missing for creating recombinant expression for the map and model entry (%s) where (_em_entity_assembly.source) is (RECOMBINANT) and polymer is (man)."
                                        % ent_id_in
                                    )
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)
                                else:
                                    return src_dicts["ent_src_syn_dict"]

            def set_recombinant_source_type(r_exp, rec_exp_in):
                """
                XSD: <xs:complexType name="recombinant_source_type"> has
                .. 1 attribute and
                .. a sequence of 5 elements
                """

                def set_attr_database(r_exp):
                    """
                    XSD: <xs:attribute name="database"> has
                    the value of "NCBI"
                    CIF: _em_entity_assembly_recombinan.ncbi_tax_id
                    """
                    set_cif_value(r_exp.set_database, "ncbi_tax_id", const.EM_ENTITY_ASSEMBLY_RECOMBINANT, cif_list=rec_exp_in, cif_value="NCBI")

                def set_el_recombinant_organism(r_exp, rec_exp_in):
                    """
                    XSD: <xs:element name="recombinant_organism" type="organism_type">
                    CIF: _em_entity_assembly_recombinan.organism 'Escherichia coli'
                    """
                    rec_exp_dict = dict(rec_exp_in)
                    if "_em_entity_assembly_recombinan.ncbi_tax_id" in rec_exp_dict:
                        tax_id = get_cif_value("ncbi_tax_id", const.EM_ENTITY_ASSEMBLY_RECOMBINANT, rec_exp_in)
                        if tax_id is not None or tax_id.isspace():
                            set_cif_value(
                                r_exp.set_recombinant_organism, "organism", const.EM_ENTITY_ASSEMBLY_RECOMBINANT, cif_list=rec_exp_in, constructor=emdb.organism_type, ncbi=tax_id
                            )
                        else:
                            txt = u"The value for (_em_entity_assembly_recombinan.ncbi_tax_id) is not given. It is required for setting the recombinant expression organism."
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                    elif "_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id" in rec_exp_dict:
                        tax_id = get_cif_value("pdbx_host_org_ncbi_taxonomy_id", const.ENTITY_SRC_GEN, rec_exp_in)
                        if tax_id is not None or tax_id.isspace():
                            set_cif_value(
                                r_exp.set_recombinant_organism,
                                "pdbx_host_org_scientific_name",
                                const.ENTITY_SRC_GEN,
                                cif_list=rec_exp_in,
                                constructor=emdb.organism_type,
                                ncbi=tax_id,
                            )
                        else:
                            txt = u"The value for (_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id) is not given. It is required for setting the recombinant expression organism."
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                    elif "_pdbx_entity_src_syn.ncbi_taxonomy_id" in rec_exp_dict:
                        pass
                    else:
                        # shouldn't happen
                        txt = u"Cannot set recombinant organism. Information is missing in (%s)." % rec_exp_dict
                        self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                def set_el_recombinant_strain(r_exp, rec_exp_in):
                    """
                    XSD: <xs:element name="recombinant_strain" type="xs:token" minOccurs="0"/>
                    CIF: _em_entity_assembly_recombinan.strain ?
                    """
                    set_cif_value(r_exp.set_recombinant_strain, "strain", const.EM_ENTITY_ASSEMBLY_RECOMBINANT, cif_list=rec_exp_in)

                def set_el_recombinant_cell(r_exp, rec_exp_in):
                    """
                    XSD: <xs:element name="recombinant_cell" type="xs:token" minOccurs="0">
                    CIF: _em_entity_assembly_recombinan.cell ?
                    """
                    set_cif_value(r_exp.set_recombinant_cell, "cell", const.EM_ENTITY_ASSEMBLY_RECOMBINANT, cif_list=rec_exp_in)

                def set_el_recombinant_plasmid(r_exp, rec_exp_in):
                    """
                    XSD: <xs:element name="recombinant_plasmid" type="xs:token" minOccurs="0"/>
                    CIF: _em_entity_assembly_recombinan.plasmid 'pHis-Parallel 1'
                    """
                    set_cif_value(r_exp.set_recombinant_plasmid, "plasmid", const.EM_ENTITY_ASSEMBLY_RECOMBINANT, cif_list=rec_exp_in)

                def set_el_recombinant_synonym_organism():
                    """
                    <xs:element name="synonym_organism" type="xs:token" minOccurs="0">
                    Deprecated (2014-10-21)
                    """

                # attribute 1
                set_attr_database(r_exp)
                # element 1
                set_el_recombinant_organism(r_exp, rec_exp_in)
                # element 2
                set_el_recombinant_strain(r_exp, rec_exp_in)
                # element 3
                set_el_recombinant_cell(r_exp, rec_exp_in)
                # element 4
                set_el_recombinant_plasmid(r_exp, rec_exp_in)
                # element 5
                set_el_recombinant_synonym_organism()

            def set_el_supramolecule_list(sample):
                """
                XSD: <xs:element name="supramolecule_list" maxOccurs="1">
                """

                def set_sup_mol_weight(supra_mol, mol_wt_in):
                    """
                    Set molecular weight object of supra_mol from values in sup_in

                    Parameters:
                    @param supra_mol: supramolecule object with a set_molecular_weight method
                    @param mol_wt_in: cif em_entity_assembly_molwt
                    """

                    def set_molecular_weight_type(mol_weight):
                        """
                        XSD: <xs:complexType name="molecular_weight_type"> is
                        .. a sequence of 3 elements
                        CIF: _em_entity_assembly_molwt.value
                        """
                        wt_units_in = get_cif_value("units", const.EM_ENTITY_ASSEMBLY_MOLWT, mol_wt_in)
                        if wt_units_in == "MEGADALTONS":
                            wt_units = const.U_MDA
                        elif wt_units_in == "KILODALTONS/NANOMETER":
                            wt_units = const.U_KDA_NM
                        else:
                            wt_units = const.U_MDA

                        exprtl = get_cif_value("experimental", const.EM_ENTITY_ASSEMBLY_MOLWT, mol_wt_in)
                        if exprtl == "YES":
                            # element 1
                            # XSD: <xs:element name="experimental" minOccurs="0">
                            set_cif_value(mol_weight.set_experimental, "value", const.EM_ENTITY_ASSEMBLY_MOLWT, cif_list=mol_wt_in, constructor=emdb.experimentalType, units=wt_units)
                        else:
                            # element 2
                            # XSD: <xs:element name="theoretical" minOccurs="0">
                            set_cif_value(mol_weight.set_theoretical, "value", const.EM_ENTITY_ASSEMBLY_MOLWT, cif_list=mol_wt_in, constructor=emdb.experimentalType, units=wt_units)
                            # element 3
                            # XSD: <xs:element name="method" type="xs:token" minOccurs="0"/>
                            # CIF: method doesn't exist in the dictionary!!!!
                            # mol_weight.set_method()

                    mol_weight = emdb.molecular_weight_type()
                    set_molecular_weight_type(mol_weight)
                    if mol_weight.has__content():
                        supra_mol.set_molecular_weight(mol_weight)

                def set_sup_mol_base_source(sup_mol, cif_category, src_in):
                    """
                    Creates the base source element for supramolecules

                    Parameters:
                    @param sup_mol: supramolecule source object of supramolecule source of
                    emdb.{virus,cell,tissue,organelle,complex,sample}_source_type
                    @param cif_category: contains natural source info (em_entity_assembly_naturalsource)
                    @param src_in: source in cif: natural from em_entity_assembly_naturalsource
                    XSD: <xs:complexType name="base_source_type">
                    """
                    set_base_source_type(sup_mol, cif_category, src_in)

                def set_sup_mol_nat_src(nat_src, sup_mol, cif_category, src_dict_in, flags_dict):
                    """
                    Set natural source for a supramolecule

                    Parameters:
                    @param nat_src: natural_source object of different type,
                           depending on the supramolecule
                    @param sup_mol: supramolecule object
                    @param cif_category: contains natural source info
                           (em_entity_assembly_naturalsource for supramolecules)
                    @param src_dict_in: source dictionary keyed by supramolecule id,
                                        each value is a list of natural sources
                    @param flags_dict: a dictionary containing boolean values for the following:
                           add_nat_src: flag for adding or setting natural source
                           (only true used at the moment)
                           the following true if required:
                           add_organ, add_tissue, add_cell, add_organelle, add_cellular_location
                    Supramolecule natural source is
                    .. an extension of base="base_source_type" and
                    .. a sequence of 5 possible elements
                    """

                    def set_el_organ(nat_src, cif_category, sup_mol_nat_src_in):
                        """
                        XSD: <xs:element name="organ" type="xs:token" minOccurs="0"/>
                        CIF: _em_entity_assembly_naturalsource.organ .
                        """
                        set_cif_value(nat_src.set_organ, "organ", cif_category, cif_list=sup_mol_nat_src_in)

                    def set_el_tissue(nat_src, cif_category, sup_mol_nat_src_in):
                        """
                        XSD: <xs:element name="tissue" type="xs:token" minOccurs="0">
                        CIF: _em_entity_assembly_naturalsource.tissue ?
                        """
                        set_cif_value(nat_src.set_tissue, "tissue", cif_category, cif_list=sup_mol_nat_src_in)

                    def set_el_cell(nat_src, cif_category, sup_mol_nat_src_in):
                        """
                        XSD: <xs:element name="cell" type="xs:token" minOccurs="0">
                        CIF: _em_entity_assembly_naturalsource.cell    ?
                        """
                        set_cif_value(nat_src.set_cell, "cell", cif_category, cif_list=sup_mol_nat_src_in)

                    def set_el_organelle(nat_src, cif_category, sup_mol_nat_src_in):
                        """
                        XSD: <xs:element name="organelle" type="xs:token" minOccurs="0">
                        CIF: _em_entity_assembly_naturalsource.organelle ?
                        """
                        set_cif_value(nat_src.set_organelle, "organelle", cif_category, cif_list=sup_mol_nat_src_in)

                    def set_el_cellular_location(nat_src, cif_category, sup_mol_nat_src_in):
                        """
                        XSD: <xs:element name="cellular_location" type="xs:token" minOccurs="0">
                        CIF: _em_entity_assembly_naturalsource.cellular_location ?
                        """
                        set_cif_value(nat_src.set_cellular_location, "cellular_location", cif_category, cif_list=sup_mol_nat_src_in)

                    for sup_mol_nat_src_in in src_dict_in:
                        set_sup_mol_base_source(nat_src, cif_category, sup_mol_nat_src_in)

                        if flags_dict["add_organ"]:
                            set_el_organ(nat_src, cif_category, sup_mol_nat_src_in)
                        if flags_dict["add_tissue"]:
                            set_el_tissue(nat_src, cif_category, sup_mol_nat_src_in)
                        if flags_dict["add_cell"]:
                            set_el_cell(nat_src, cif_category, sup_mol_nat_src_in)
                        if flags_dict["add_organelle"]:
                            set_el_organelle(nat_src, cif_category, sup_mol_nat_src_in)
                        if flags_dict["add_cellular_location"]:
                            set_el_cellular_location(nat_src, cif_category, sup_mol_nat_src_in)

                        if nat_src.has__content():
                            if flags_dict["add_nat_src"]:
                                sup_mol.add_natural_source(nat_src)
                            else:
                                sup_mol.set_natural_source(nat_src)
                            break

                def set_sup_mol_rec_exp(sup_mol, sup_mol_id_in, rec_exp_dict_in, add_rec_exp=False, virus=False):
                    """
                    Set recombinant expression for a supramolecule

                    Parameters:
                    @param sup_mol: supramolecule object
                    @param sup_mol_id_in: cif id for supramolecule
                    @param rec_exp_dict_in: recombinant expression dictionary keyed by supramolecule id
                           each value is a list of recombinant expression systems for the supramolecule
                    @param add_rec_exp: boolean -if true add multiple recombinant expression systems instead of just setting one
                    @param virus: boolean - if true use host_system instead of recombinant_expression
                    """
                    if rec_exp_dict_in is not None:
                        if sup_mol_id_in in rec_exp_dict_in:
                            for rec_exp_in in rec_exp_dict_in[sup_mol_id_in]:
                                r_exp = emdb.recombinant_source_type()
                                set_recombinant_source_type(r_exp, rec_exp_in)
                                if r_exp.has__content():
                                    if add_rec_exp is False:
                                        if virus is False:
                                            sup_mol.set_recombinant_expression(r_exp)
                                        else:
                                            sup_mol.set_host_system(r_exp)
                                        break
                                    else:
                                        sup_mol.add_recombinant_expression(r_exp)

                def set_base_sup_mol(sup_mol, sup_in, sup_mol_id_in, sample=None):
                    """
                    Set base parameters of base_supramolecule_type

                    Parameters:
                    @param sup_mol: object wrapping supramolecule element
                                    - this object will be updated
                    @param sup_in: cif em_entity_assembly category dictionary
                    @param sup_mol_id_in: cif id of supramolecule
                    XSD: <xs:complexType name="base_supramolecule_type"> has
                    .. 1 attribute and
                    .. a sequence of 9 elements
                    """

                    def set_attr_id(sup_mol, sup_mol_id_in):
                        """
                        XSD: <xs:attribute name="id" type="xs:positiveInteger" use="required"/>
                        CIF: _em_entity_assembly.id 1
                        """
                        set_cif_value(sup_mol.set_supramolecule_id, "id", const.EM_ENTITY_ASSEMBLY, cif_value=sup_mol_id_in, fmt=int)

                    def set_el_name(sup_mol, sup_in):
                        """
                        XSD: <xs:element name="name" type="sci_name_type">
                        CIF: _em_entity_assembly.name 'Israeli acute paralysis virus'
                        """
                        if isinstance(sup_mol, emdb.virus_supramolecule_type):
                            virus_name = get_cif_value("name", const.EM_ENTITY_ASSEMBLY, cif_list=sup_in)
                            if virus_name is None:
                                virus_name_from_nat_source = get_cif_value("organism", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=sup_in)
                                if virus_name_from_nat_source is not None:
                                    set_cif_value(
                                        sup_mol.set_name,
                                        "organism",
                                        const.EM_ENTITY_ASSEMBLY_NATURALSOURCE,
                                        cif_list=sup_in,
                                        constructor=emdb.sci_name_type,
                                        cif_value=virus_name_from_nat_source,
                                    )
                                    txt = u"(_em_entity_assembly.name) is not given so the value for (_em_entity_assembly_naturalsource.organism) is used."
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                                else:
                                    set_cif_value(sup_mol.set_name, "organism", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=sup_in, constructor=emdb.sci_name_type, cif_value="Unspecified")
                                    txt = u"The values for (_em_entity_assembly.name) and (_em_entity_assembly_naturalsource.organism) are not given. The value of Unspecified is used."
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                            else:
                                set_cif_value(sup_mol.set_name, "name", const.EM_ENTITY_ASSEMBLY, cif_list=sup_in, constructor=emdb.sci_name_type)
                        else:
                            set_cif_value(sup_mol.set_name, "name", const.EM_ENTITY_ASSEMBLY, cif_list=sup_in, constructor=emdb.sci_name_type)

                    def set_el_category(sup_mol, sup_in):
                        """
                        XSD: <xs:element name="category" minOccurs="0">
                        CIF: _em_entity_assembly.go_id ?

                        At the moment only the GO category is being captured; others:
                        ARBITRARY DEFINITION and PROTEIN ONTOLOGY are not supported
                        """
                        set_cif_value(sup_mol.set_category, "go_id", const.EM_ENTITY_ASSEMBLY, cif_list=sup_in, constructor=emdb.categoryType, type="GO")

                    def set_el_parent(sup_mol, sup_in):
                        """
                        XSD: <xs:element name="parent" type="xs:nonNegativeInteger">
                        CIF: _em_entity_assembly.parent_id 0

                        This is where the sample name is set if the parent is 0
                        XSD: name="sample" type="sample_type" element 1:
                        XSD: <xs:element name="name" type="sci_name_type">
                        CIF: _em_entity_assembly.name 'Israeli acute paralysis virus'
                        """
                        set_cif_value(sup_mol.set_parent, "parent_id", const.EM_ENTITY_ASSEMBLY, cif_list=sup_in, fmt=int)

                        # The name given to the supramolecule that has parent id = 0
                        parent_id = get_cif_value("parent_id", const.EM_ENTITY_ASSEMBLY, sup_in)
                        if parent_id == "0":
                            set_cif_value(sample.set_name, "name", const.EM_ENTITY_ASSEMBLY, cif_list=sup_in, constructor=emdb.sci_name_type)

                    def set_el_macromolecule_list(sup_mol, sup_in):
                        """
                        XSD: <xs:element name="macromolecule_list" minOccurs="0">
                        CIF: _em_entity_assembly.entity_id_list '1, 2, 3, 4'
                        """
                        macro_list_str_in = get_cif_value("entity_id_list", const.EM_ENTITY_ASSEMBLY, sup_in)
                        id_list_item = get_cif_item("entity_id_list", const.EM_ENTITY_ASSEMBLY)
                        if macro_list_str_in is not None:
                            macro_list_in = macro_list_str_in.rstrip().rstrip(",").split(",")
                            if macro_list_in is not None and len(macro_list_in) > 0:
                                macro_list = emdb.macromolecule_list_type()
                                for m_in in macro_list_in:
                                    a_macromol = emdb.macromoleculeType(macromolecule_id=int(m_in))
                                    a_macromol.original_tagname_ = "macromolecule"
                                    macro_list.add_macromolecule(a_macromol)
                                txt = u"Macromolecule (%s) added to the list of macromolecules." % int(m_in)  # pylint: disable=undefined-loop-variable
                                self.current_entry_log.info_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.info_title + txt))
                                self.log_formatted(self.info_log_string, const.INFO_ALERT + txt)
                                if macro_list.has__content():
                                    sup_mol.set_macromolecule_list(macro_list)
                            else:
                                txt = u"No macromolecule found for (%s)." % id_list_item
                                self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                                self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)

                    def set_el_details(sup_mol, sup_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _em_entity_assembly.details ?
                        """
                        set_cif_value(sup_mol.set_details, "details", const.EM_ENTITY_ASSEMBLY, cif_list=sup_in)

                    def set_el_number_of_copies():
                        """
                        XSD: <xs:element name="number_of_copies" type="pos_int_or_string_type" minOccurs="0">
                        CIF: This information doesn't exist in the cif files
                             nor is captured in OneDep
                             (unless the entry has a full overlap PDB model)
                        """

                    def set_el_oligomeric_state():
                        """
                        XSD: <xs:element name="oligomeric_state" type="pos_int_or_string_type" minOccurs="0">
                            Deprecated. Oligomeric state in parent,
                            or if sample, the oligomeric state of sample.
                        """

                    def set_el_external_references():
                        """
                        XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0" >
                        CIF: DOESN'T EXIST - TALK TO ARDAN?
                        """

                    def set_el_recombinant_exp_flag():
                        """
                        XSD: <xs:element name="recombinant_exp_flag" type="xs:boolean" maxOccurs="1" minOccurs="0" >
                        Deprecated 2014/12/3
                        """

                    # attribute 1
                    set_attr_id(sup_mol, sup_mol_id_in)
                    # element 1
                    set_el_name(sup_mol, sup_in)
                    # element 2
                    set_el_category(sup_mol, sup_in)
                    # element 3
                    set_el_parent(sup_mol, sup_in)
                    # element 4
                    set_el_macromolecule_list(sup_mol, sup_in)
                    # element 5
                    set_el_details(sup_mol, sup_in)
                    # element 6
                    set_el_number_of_copies()
                    # element 7
                    set_el_oligomeric_state()
                    # element 8
                    set_el_external_references()
                    # element 9
                    set_el_recombinant_exp_flag()

                def set_supramolecule_list(sup_list, sup_mol_dicts):
                    """
                    XSD: <xs:element name="supramolecule_list" maxOccurs="1"> is a
                    ..a sequence of <xs:element ref="supramolecule" maxOccurs="unbounded"/>

                    """

                    def set_complex_supramolecule_type(complex_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample):
                        """
                        XSD: <xs:element name="complex_supramolecule" substitutionGroup="supramolecule" type="complex_supramolecule_type"/> is
                        .. an extension of base="base_supramolecule_type" with
                        .. 1 attribute and
                        .. a sequence of 4 elements
                        """

                        def set_attr_chimera(complex_sup_mol, sup_in):
                            """
                            XSD: <xs:attribute fixed="true" name="chimera" type="xs:boolean"/>
                            CIF: _em_entity_assembly.chimera ?/YES (NO cannot be given)
                            """
                            set_cif_value(complex_sup_mol.set_chimera, "chimera", const.EM_ENTITY_ASSEMBLY, cif_list=sup_in, fmt=bool)

                        def set_el_natural_source(complex_sup_mol, sup_mol_id_in, sup_mol_dicts):
                            """
                            XSD: <xs:element name="natural_source" type="complex_source_type" minOccurs="0"/> is
                            .. an extension of base="base_source_type" and
                            .. a sequence of organ, tissue, cell, organelle and cellular location
                            """
                            nat_src_dict_in = sup_mol_dicts["nat_src_dict_in"]
                            if sup_mol_id_in in nat_src_dict_in:
                                sup_mol_dict_in = nat_src_dict_in[sup_mol_id_in]
                                cmpx_dict = {"add_nat_src": True, "add_organ": True, "add_tissue": True, "add_cell": True, "add_organelle": True, "add_cellular_location": True}
                                complex_natural_source_type_list = []
                                attr_ncbis = []
                                el_organisms = []
                                el_strains = []
                                for src_in in sup_mol_dict_in:
                                    attr_ncbi = get_cif_value("ncbi_tax_id", const.EM_VIRUS_NATURAL_HOST, cif_list=src_in)
                                    attr_ncbis.extend([attr_ncbi])
                                    el_organism = get_cif_value("organism", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_organisms.extend([el_organism])
                                    el_strain = get_cif_value("strain", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_strains.extend([el_strain])
                                complex_natural_source_type_list.extend(attr_ncbis)
                                complex_natural_source_type_list.extend(el_organisms)
                                complex_natural_source_type_list.extend(el_strains)
                                if any(x is not None for x in complex_natural_source_type_list):
                                    cns = emdb.complex_source_type()
                                    set_sup_mol_nat_src(cns, complex_sup_mol, const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, sup_mol_dict_in, cmpx_dict)

                        def set_el_recombinant_expression(complex_sup_mol, sup_mol_id_in, rec_exp_dict_in):
                            """
                            XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" maxOccurs="unbounded" minOccurs="0">
                            """
                            set_sup_mol_rec_exp(complex_sup_mol, sup_mol_id_in, rec_exp_dict_in, add_rec_exp=True)

                        def set_el_molecular_weight(complex_sup_mol, sup_mol_id_in, s_mol_wt_dict_in):
                            """
                            XSD: <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0" />
                            """
                            if sup_mol_id_in in s_mol_wt_dict_in:
                                set_sup_mol_weight(complex_sup_mol, s_mol_wt_dict_in[sup_mol_id_in])

                        def set_el_ribosome_details():
                            """
                            XSD:<xs:element name="ribosome-details" type="xs:string" minOccurs="0">
                            Deprecated (2014/11/17)
                            """
                            # if sup_type == 'RIBOSOME':
                            #    complex_sup_mol.set_ribosome_details('RIBOSOME')

                        # set up the supramolecule specific tagname explicitly
                        # as DSgenerate doesn't provide it
                        complex_sup_mol.original_tagname_ = "complex_supramolecule"
                        # the extension
                        set_base_sup_mol(complex_sup_mol, sup_in, sup_mol_id_in, sample)
                        # attribute 1
                        set_attr_chimera(complex_sup_mol, sup_in)
                        # element 1
                        set_el_natural_source(complex_sup_mol, sup_mol_id_in, sup_mol_dicts)
                        # element 2
                        rec_exp_dict_in = get_rec_exp_dict(sup_mol_id_in, sup_mol_dicts, sup_in, is_supramolecule=True)
                        set_el_recombinant_expression(complex_sup_mol, sup_mol_id_in, rec_exp_dict_in)
                        # element 3
                        set_el_molecular_weight(complex_sup_mol, sup_mol_id_in, sup_mol_dicts["s_mol_wt_dict_in"])
                        # element 4
                        set_el_ribosome_details()

                    def set_virus_supramolecule_type(virus_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample):
                        """
                        XSD: <xs:element name="virus_supramolecule" substitutionGroup="supramolecule" type="virus_supramolecule_type"/>
                        ..an extension of base="base_supramolecule_type" with
                        ..a sequence of 14 elements
                        """

                        def set_el_sci_species_name(virus_sup_mol, nat_src_in):
                            """
                            XSD: <xs:element name="sci_species_name" type="virus_species_name_type" minOccurs="0">
                            """

                            def set_virus_species_name_type(virus_name):
                                """
                                XSD: <xs:element name="sci_species_name" type="virus_species_name_type" minOccurs="0"> has
                                CIF: _em_entity_assembly_naturalsource.organism 'Oryctolagus cuniculus'
                                .. a xs:token with 1 attribute
                                """

                                def set_attr_ncbi(virus_name, nat_src_in):
                                    """
                                    XSD: <xs:attribute name="ncbi" type="xs:positiveInteger"/>
                                    CIF: _em_entity_assembly_naturalsource.ncbi_tax_id
                                    """
                                    tax_id = get_cif_value("ncbi_tax_id", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, nat_src_in)
                                    if tax_id is None:
                                        # this shouldn't happen - it's mandatory
                                        tax_id = 32644  # ID for unknown
                                        set_cif_value(virus_name.set_ncbi, "ncbi_tax_id", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=nat_src_in, cif_value=tax_id)
                                        txt = u"The value for (_em_entity_assembly_naturalsource.ncbi_tax_id) is not given. The value set is (%s) for unknown." % tax_id
                                        self.current_entry_log.warn_logs.append(
                                            self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                        )
                                        self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                                    else:
                                        set_cif_value(virus_name.set_ncbi, "ncbi_tax_id", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=nat_src_in, fmt=int)

                                # attribute 1
                                set_attr_ncbi(virus_name, nat_src_in)
                                set_cif_value(virus_name.set_valueOf_, "organism", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=nat_src_in)
                                # additionally, set the supramolecule name
                                set_cif_value(virus_sup_mol.set_name, "organism", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=nat_src_in, constructor=emdb.sci_name_type)

                            virus_name = emdb.virus_species_name_type()
                            set_virus_species_name_type(virus_name)
                            if virus_name.has__content():
                                virus_sup_mol.set_sci_species_name(virus_name)

                        def set_el_sci_species_strain(virus_sup_mol, nat_src_in):
                            """
                            XSD: <xs:element name="sci_species_strain" type="xs:string" maxOccurs="1" minOccurs="0">
                            XPath: /element(*,virus_supramolecule_type)/sci_species_strain
                            CIF: _em_entity_assembly_naturalsource.strain ?
                            """
                            set_cif_value(virus_sup_mol.set_sci_species_strain, "strain", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=nat_src_in)

                        def set_el_natural_host(virus_sup_mol, sup_mol_id_in, virus_nat_host_dict_in):
                            """
                            XSD: <xs:element name="natural_host" type="virus_host_type" minOccurs="0" maxOccurs="unbounded"/>
                            """

                            def set_virus_natural_host_type(nat_host, nat_host_in):
                                """
                                XSD: <xs:element name="natural_host" type="virus_host_type" minOccurs="0" maxOccurs="unbounded"/>
                                CIF: _em_virus_natural_host.id                     1
                                CIF: _em_virus_natural_host.entity_assembly_id   1
                                .. an extension of base_source_type that has
                                .. 1 attribute and
                                .. a sequence of 3 elements
                                """

                                def set_attr_database(nat_host, nat_host_in):
                                    """
                                    XSD: <xs:attribute name="database">
                                    CIF: _em_virus_natural_host.ncbi_tax_id 7460
                                    """
                                    set_cif_value(nat_host.set_database, "ncbi_tax_id", const.EM_VIRUS_NATURAL_HOST, cif_list=nat_host_in, cif_value="NCBI")

                                def set_el_organism(org, nat_host_in):
                                    """
                                    XSD: <xs:element name="organism" type="organism_type"/> has
                                    .. 1 attribute
                                    """

                                    def set_attr_ncbi(org, nat_host_in):
                                        """
                                        XSD: <xs:attribute name="ncbi" type="xs:positiveInteger"/>
                                        CIF: _em_virus_natural_host.ncbi_tax_id
                                        """
                                        set_cif_value(org.set_ncbi, "ncbi_tax_id", const.EM_VIRUS_NATURAL_HOST, cif_list=nat_host_in, fmt=int, parent_el_req=False)

                                    def set_organism(org, nat_host_in):
                                        """
                                        XSD: <xs:element name="organism" type="organism_type">
                                        CIF: _em_virus_natural_host.organism 'Apis mellifera'
                                        """
                                        set_cif_value(org.set_valueOf_, "organism", const.EM_VIRUS_NATURAL_HOST, cif_list=nat_host_in, parent_el_req=False)

                                    # attribute 1
                                    set_attr_ncbi(org, nat_host_in)
                                    # element
                                    set_organism(org, nat_host_in)

                                def set_el_strain(nat_host, nat_host_in):
                                    """
                                    XSD: <xs:element name="strain" type="xs:token" minOccurs="0"/>
                                    CIF: _em_virus_natural_host.strain ?
                                    """
                                    set_cif_value(nat_host.set_strain, "strain", const.EM_VIRUS_NATURAL_HOST, cif_list=nat_host_in, parent_el_req=False)

                                def set_el_synonym_organism():
                                    """
                                    XSD: <xs:element name="synonym_organism" type="xs:token" minOccurs="0">
                                    Deprecated (2014-10-21)
                                    """

                                # attribute 1
                                set_attr_database(nat_host, nat_host_in)
                                # element 1
                                org = emdb.organism_type()
                                set_el_organism(org, nat_host_in)
                                nat_host.set_organism(org)
                                # element 2
                                set_el_strain(nat_host, nat_host_in)
                                # element 3
                                set_el_synonym_organism()

                            if sup_mol_id_in in virus_nat_host_dict_in:
                                nat_host_list_in = virus_nat_host_dict_in[sup_mol_id_in]
                                for nat_host_in in nat_host_list_in:
                                    attr_ncbi = get_cif_value("ncbi_tax_id", const.EM_VIRUS_NATURAL_HOST, cif_list=nat_host_in)
                                    organism = get_cif_value("organism", const.EM_VIRUS_NATURAL_HOST, cif_list=nat_host_in)
                                    if any(x is not None for x in [attr_ncbi, organism]):
                                        nat_host = emdb.virus_host_type()
                                        set_virus_natural_host_type(nat_host, nat_host_in)
                                        if nat_host.has__content():
                                            virus_sup_mol.add_natural_host(nat_host)

                        def set_el_host_system(virus_sup_mol, sup_mol_id_in, rec_exp_dict_in):
                            """
                            XSD: <xs:element name="host_system" type="recombinant_source_type">
                            """
                            set_sup_mol_rec_exp(virus_sup_mol, sup_mol_id_in, rec_exp_dict_in, virus=True)

                        def set_el_molecular_weight(virus_sup_mol, sup_mol_id_in, s_mol_wt_dict_in):
                            """
                            XSD:  <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0"/>
                            """
                            if sup_mol_id_in in s_mol_wt_dict_in:
                                set_sup_mol_weight(virus_sup_mol, s_mol_wt_dict_in[sup_mol_id_in])

                        def set_el_virus_shell(virus_sup_mol, sup_mol_id_in, virus_shell_dict_in):
                            """
                            XSD:  <xs:element name="virus_shell" maxOccurs="unbounded" minOccurs="0"> has
                            .. 1 attribute and
                            .. a sequence of 3 elements
                            """

                            def set_attr_id(virus_shell, vs_in):
                                """
                                XSD: <xs:attribute name="shell_id" type="xs:positiveInteger"/>
                                CIF: _em_virus_shell.id
                                """
                                set_cif_value(virus_shell.set_shell_id, const.K_ID, const.EM_VIRUS_SHELL, cif_list=vs_in, fmt=int)

                            def set_el_name(virus_shell, vs_in):
                                """
                                XSD: <xs:element name="name" type="xs:token" nillable="false" minOccurs="0"/>
                                CIF: _em_virus_shell.name
                                """
                                set_cif_value(virus_shell.set_name, "name", const.EM_VIRUS_SHELL, cif_list=vs_in)

                            def set_el_diameter(virus_shell, vs_in):
                                """
                                XSD: <xs:element name="diameter" minOccurs="0">
                                CIF: _em_virus_shell.diameter
                                """
                                set_cif_value(
                                    virus_shell.set_diameter, "diameter", const.EM_VIRUS_SHELL, cif_list=vs_in, constructor=emdb.diameterType, fmt=float, units=const.U_ANG
                                )

                            def set_el_triangulation(virus_shell, vs_in):
                                """
                                XSD: <xs:element name="triangulation" type="xs:positiveInteger" minOccurs="0">
                                CIF: _em_virus_shell.triangulation
                                """
                                set_cif_value(virus_shell.set_triangulation, "triangulation", const.EM_VIRUS_SHELL, cif_list=vs_in, fmt=int)

                            if sup_mol_id_in in virus_shell_dict_in:
                                virus_shell_in = virus_shell_dict_in[sup_mol_id_in]
                                for vs_in in virus_shell_in:
                                    virus_shell = emdb.virus_shellType()
                                    # attribute 1
                                    set_attr_id(virus_shell, vs_in)
                                    # element 1
                                    set_el_name(virus_shell, vs_in)
                                    # element 2
                                    set_el_diameter(virus_shell, vs_in)
                                    # element 3
                                    set_el_triangulation(virus_shell, vs_in)

                                    if virus_shell.has__content():
                                        virus_sup_mol.add_virus_shell(virus_shell)

                        def set_el_virus_type(virus_sup_mol, virus_in):
                            """
                            XSD:  <xs:element name="virus_type">; base xs:token with restriction
                            CIF: _em_virus_entity.virus_type PRION
                            """
                            set_cif_value(virus_sup_mol.set_virus_type, "virus_type", const.EM_VIRUS_ENTITY, cif_list=virus_in)

                        def set_el_virus_isolate(virus_sup_mol, virus_in):
                            """
                            XSD: <xs:element name="virus_isolate">; base xs:token with restriction
                            CIF: _em_virus_entity.virus_isolate
                            """
                            virus_iso = get_cif_value("virus_isolate", const.EM_VIRUS_ENTITY, virus_in)
                            if virus_iso is None:
                                virus_iso = "OTHER"  # default value
                                set_cif_value(virus_sup_mol.set_virus_isolate, "virus_isolate", const.EM_VIRUS_ENTITY, cif_list=virus_in, cif_value=virus_iso)
                                txt = u"The value for (_em_virus_entity.virus_isolate) is not given. Set to (OTHER)."
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                            else:
                                set_cif_value(virus_sup_mol.set_virus_isolate, "virus_isolate", const.EM_VIRUS_ENTITY, cif_list=virus_in)

                        def set_el_virus_enveloped(virus_sup_mol, virus_in):
                            """
                            XSD:  <xs:element name="virus_enveloped" type="xs:boolean"/>
                            CIF: _em_virus_entity.enveloped YES
                            """
                            virus_env = get_cif_value("enveloped", const.EM_VIRUS_ENTITY, virus_in)
                            if virus_env is None:
                                virus_env = "False"  # default value
                                set_cif_value(virus_sup_mol.set_virus_enveloped, "enveloped", const.EM_VIRUS_ENTITY, cif_list=virus_in, cif_value=virus_env)
                                txt = u"The value for (_em_virus_entity.enveloped) is not given. Set to (False)."
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                            else:
                                set_cif_value(virus_sup_mol.set_virus_enveloped, "enveloped", const.EM_VIRUS_ENTITY, cif_list=virus_in, fmt=cif_bool)

                        def set_el_virus_empty(virus_sup_mol, virus_in):
                            """
                            XSD:  <xs:element name="virus_empty" type="xs:boolean"/>
                            CIF: _em_virus_entity.empty YES
                            """
                            virus_empty = get_cif_value("empty", const.EM_VIRUS_ENTITY, virus_in)
                            if virus_empty is None:
                                virus_empty = "True"  # default value
                                set_cif_value(virus_sup_mol.set_virus_empty, "empty", const.EM_VIRUS_ENTITY, cif_list=virus_in, cif_value=virus_empty)
                                txt = u"The value for (_em_virus_entity.enveloped) is not given. Set to (True)."
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                            else:
                                set_cif_value(virus_sup_mol.set_virus_empty, "empty", const.EM_VIRUS_ENTITY, cif_list=virus_in, fmt=cif_bool)

                        def set_el_syn_species_name():
                            """
                            XSD: <xs:element name="syn_species_name" type="xs:string" maxOccurs="1" minOccurs="0">
                            Deprecated (2014/11/17)
                            """

                        def set_el_sci_species_serotype():
                            """
                            XSD: <xs:element name="sci_species_serotype" type="xs:string" maxOccurs="1" minOccurs="0">
                            Deprecated (2014/11/17)
                            """

                        def set_el_sci_species_serocomplex():
                            """
                            XSD: <xs:element name="sci_species_serocomplex" type="xs:string" maxOccurs="1" minOccurs="0">
                            Deprecated (2014/11/17)
                            """

                        def set_el_sci_species_subspecies():
                            """
                            XSD: <xs:element name="sci_species_subspecies" type="xs:string" maxOccurs="1" minOccurs="0">
                            Deprecated (2014/11/17)
                            """

                        # set up the supramolecule specific tagname explicitly
                        # as DSgenerate doesn't provide it
                        virus_sup_mol.original_tagname_ = "virus_supramolecule"
                        set_base_sup_mol(virus_sup_mol, sup_in, sup_mol_id_in, sample)
                        if sup_mol_id_in in nat_src_dict_in:
                            nat_src_list_in = nat_src_dict_in[sup_mol_id_in]
                            len_nat_src_list_in = len(nat_src_list_in)
                            if len_nat_src_list_in > 0:
                                if len_nat_src_list_in > 1:
                                    txt = u"Only the first row of the (%s) (_em_entity_assembly_naturalsource) rows for supramolecule (%s)." % (len_nat_src_list_in, sup_mol_id_in)
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)
                                nat_src_in = nat_src_list_in[0]
                                # element 1
                                set_el_sci_species_name(virus_sup_mol, nat_src_in)
                                # element 2
                                set_el_sci_species_strain(virus_sup_mol, nat_src_in)
                            else:
                                txt = u"Empty natural source category for virus!"
                                self.current_entry_log.error_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                )
                                self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                        # element 3
                        virus_nat_host_dict_in = sup_mol_dicts["virus_nat_host_dict_in"]
                        set_el_natural_host(virus_sup_mol, sup_mol_id_in, virus_nat_host_dict_in)
                        # element 4
                        rec_exp_dict_in = get_rec_exp_dict(sup_mol_id_in, sup_mol_dicts, sup_in, is_supramolecule=True)
                        set_el_host_system(virus_sup_mol, sup_mol_id_in, rec_exp_dict_in)
                        # element 5
                        s_mol_wt_dict_in = sup_mol_dicts["s_mol_wt_dict_in"]
                        set_el_molecular_weight(virus_sup_mol, sup_mol_id_in, s_mol_wt_dict_in)
                        # element 6
                        virus_shell_dict_in = sup_mol_dicts["virus_shell_dict_in"]
                        set_el_virus_shell(virus_sup_mol, sup_mol_id_in, virus_shell_dict_in)

                        virus_in = {}
                        virus_dict_in = sup_mol_dicts["virus_dict_in"]
                        if virus_dict_in is not None:
                            if sup_mol_id_in in virus_dict_in:
                                virus_in = virus_dict_in[sup_mol_id_in]
                                if virus_in is not None:
                                    # element 7
                                    set_el_virus_type(virus_sup_mol, virus_in)
                                    # element 8
                                    set_el_virus_isolate(virus_sup_mol, virus_in)
                                    # element 9
                                    set_el_virus_enveloped(virus_sup_mol, virus_in)
                                    # element 10
                                    set_el_virus_empty(virus_sup_mol, virus_in)
                            else:
                                txt = u"Cannot set virus type. This supramolecule with id=(%s) is not in the (_em_virus_entity) category: (%s)." % (sup_mol_id_in, virus_dict_in)
                                self.current_entry_log.error_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                )
                                self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                        else:
                            txt = u"Cannot set virus type. The (_em_virus_entity) category does not exist in cif."
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                        # element 11
                        set_el_syn_species_name()
                        # element 12
                        set_el_sci_species_serotype()
                        # element 13
                        set_el_sci_species_serocomplex()
                        # element 14
                        set_el_sci_species_subspecies()

                    def set_orgorcell_supmol_type(org_or_cell_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample):
                        """
                        XSD: <xs:element name="organelle_or_cellular_component_supramolecule" substitutionGroup="supramolecule" type="organelle_or_cellular_component_supramolecule_type"/>
                        XSD: type="organelle_or_cellular_component_supramolecule_type" has
                        .. 3 additional elements to the base ones
                        """

                        def set_el_natural_source(org_or_cell_sup_mol, nat_src_dict_in):
                            """
                            XSD: <xs:element name="natural_source">
                            """
                            if sup_mol_id_in in nat_src_dict_in:
                                sup_mol_dict_in = nat_src_dict_in[sup_mol_id_in]
                                org_dict = {"add_nat_src": True, "add_organ": True, "add_tissue": True, "add_cell": True, "add_organelle": True, "add_cellular_location": True}
                                organelle_natural_source_type_list = []
                                attr_ncbis = []
                                el_organisms = []
                                el_strains = []
                                el_organs = []
                                el_tissues = []
                                el_cells = []
                                el_organelles = []
                                el_cellular_locations = []
                                for src_in in sup_mol_dict_in:
                                    attr_ncbi = get_cif_value("ncbi_tax_id", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    attr_ncbis.extend([attr_ncbi])
                                    el_organism = get_cif_value("organism", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_organisms.extend([el_organism])
                                    el_strain = get_cif_value("strain", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_strains.extend([el_strain])
                                    el_organ = get_cif_value("organ", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_organs.extend([el_organ])
                                    el_tissue = get_cif_value("tissue", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_tissues.extend([el_tissue])
                                    el_cell = get_cif_value("cell", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_cells.extend([el_cell])
                                    el_organelle = get_cif_value("organelle", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_organelles.extend([el_organelle])
                                    el_cellular_location = get_cif_value("cellular_location", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_cellular_locations.extend([el_cellular_location])
                                organelle_natural_source_type_list.extend(attr_ncbis)
                                organelle_natural_source_type_list.extend(el_organisms)
                                organelle_natural_source_type_list.extend(el_strains)
                                organelle_natural_source_type_list.extend(el_organs)
                                organelle_natural_source_type_list.extend(el_tissues)
                                organelle_natural_source_type_list.extend(el_cells)
                                organelle_natural_source_type_list.extend(el_organelles)
                                organelle_natural_source_type_list.extend(el_cellular_locations)
                                if any(x is not None for x in organelle_natural_source_type_list):
                                    cns = emdb.organelle_source_type()
                                    set_sup_mol_nat_src(cns, org_or_cell_sup_mol, const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, sup_mol_dict_in, org_dict)

                        def set_el_molecular_weight(org_or_cell_sup_mol, sup_mol_id_in, s_mol_wt_dict_in):
                            """
                            XSD: <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0"/>
                            """
                            if sup_mol_id_in in s_mol_wt_dict_in:
                                set_sup_mol_weight(org_or_cell_sup_mol, s_mol_wt_dict_in[sup_mol_id_in])

                        def set_el_recombinant_expression(org_or_cell_sup_mol, sup_mol_id_in, rec_exp_dict_in):
                            """
                            XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" maxOccurs="1" minOccurs="0" >
                            """
                            set_sup_mol_rec_exp(org_or_cell_sup_mol, sup_mol_id_in, rec_exp_dict_in)

                        # set up the supramolecule specific tagname explicitly
                        # as DSgenerate doesn't provide it
                        org_or_cell_sup_mol.original_tagname_ = "organelle_or_cellular_component_supramolecule"
                        set_base_sup_mol(org_or_cell_sup_mol, sup_in, sup_mol_id_in, sample)
                        # element 1
                        set_el_natural_source(org_or_cell_sup_mol, nat_src_dict_in)
                        # element 2
                        set_el_molecular_weight(org_or_cell_sup_mol, sup_mol_id_in, s_mol_wt_dict_in)
                        # element 3
                        rec_exp_dict_in = get_rec_exp_dict(sup_mol_id_in, sup_mol_dicts, sup_in, is_supramolecule=True)
                        set_el_recombinant_expression(org_or_cell_sup_mol, sup_mol_id_in, rec_exp_dict_in)

                    def set_tissue_supramolecule_type(tissue_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample):
                        """
                        XSD: <xs:element name="tissue_supramolecule" substitutionGroup="supramolecule" type="tissue_supramolecule_type"/>
                        XSD: type="tissue_supramolecule_type" has
                        .. 1 additional element to the base ones
                        """

                        def set_el_natural_source(tissue_sup_mol, sup_mol_id_in, sup_mol_dicts):
                            """
                            XSD: <xs:element name="natural_source">
                            """
                            nat_src_dict_in = sup_mol_dicts["nat_src_dict_in"]
                            if sup_mol_id_in in nat_src_dict_in:
                                sup_mol_dict_in = nat_src_dict_in[sup_mol_id_in]
                                tiss_dict = {"add_nat_src": True, "add_organ": True, "add_tissue": True, "add_cell": False, "add_organelle": False, "add_cellular_location": False}
                                tissue_natural_source_type_list = []
                                attr_ncbis = []
                                el_organisms = []
                                el_strains = []
                                el_organs = []
                                el_tissues = []
                                for src_in in sup_mol_dict_in:
                                    attr_ncbi = get_cif_value("ncbi_tax_id", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    attr_ncbis.extend([attr_ncbi])
                                    el_organism = get_cif_value("organism", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_organisms.extend([el_organism])
                                    el_strain = get_cif_value("strain", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_strains.extend([el_strain])
                                    el_organ = get_cif_value("organ", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_organs.extend([el_organ])
                                    el_tissue = get_cif_value("tissue", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_tissues.extend([el_tissue])
                                tissue_natural_source_type_list.extend(attr_ncbis)
                                tissue_natural_source_type_list.extend(el_organisms)
                                tissue_natural_source_type_list.extend(el_strains)
                                tissue_natural_source_type_list.extend(el_organs)
                                tissue_natural_source_type_list.extend(el_tissues)
                                if any(x is not None for x in tissue_natural_source_type_list):
                                    cns = emdb.tissue_source_type()
                                    set_sup_mol_nat_src(cns, tissue_sup_mol, const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, sup_mol_dict_in, tiss_dict)

                        # set up the supramolecule specific tagname explicitly
                        # as DSgenerate doesn't provide it
                        tissue_sup_mol.original_tagname_ = "tissue_supramolecule"
                        set_base_sup_mol(tissue_sup_mol, sup_in, sup_mol_id_in, sample)
                        # element 1
                        set_el_natural_source(tissue_sup_mol, sup_mol_id_in, sup_mol_dicts)

                    def set_cell_supramolecule_type(cell_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample):
                        """
                        XSD: <xs:element name="cell_supramolecule" substitutionGroup="supramolecule" type="cell_supramolecule_type"/>
                        XSD: type="cell_supramolecule_type" has
                        .. 1 additional element to the base ones
                        """

                        def set_el_natural_source(cell_sup_mol, sup_mol_id_in, sup_mol_dicts):
                            """
                            XSD: <xs:element name="natural_source">
                            """
                            nat_src_dict_in = sup_mol_dicts["nat_src_dict_in"]
                            if sup_mol_id_in in nat_src_dict_in:
                                sup_mol_dict_in = nat_src_dict_in[sup_mol_id_in]
                                cell_dict = {"add_nat_src": True, "add_organ": True, "add_tissue": True, "add_cell": True, "add_organelle": False, "add_cellular_location": False}
                                cell_natural_source_type_list = []
                                attr_ncbis = []
                                el_organisms = []
                                el_strains = []
                                el_organs = []
                                el_tissues = []
                                el_cells = []
                                for src_in in sup_mol_dict_in:
                                    attr_ncbi = get_cif_value("ncbi_tax_id", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    attr_ncbis.extend([attr_ncbi])
                                    el_organism = get_cif_value("organism", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_organisms.extend([el_organism])
                                    el_strain = get_cif_value("strain", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_strains.extend([el_strain])
                                    el_organ = get_cif_value("organ", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_organs.extend([el_organ])
                                    el_tissue = get_cif_value("tissue", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_tissues.extend([el_tissue])
                                    el_cell = get_cif_value("cell", const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, cif_list=src_in)
                                    el_cells.extend([el_cell])
                                cell_natural_source_type_list.extend(attr_ncbis)
                                cell_natural_source_type_list.extend(el_organisms)
                                cell_natural_source_type_list.extend(el_strains)
                                cell_natural_source_type_list.extend(el_organs)
                                cell_natural_source_type_list.extend(el_tissues)
                                cell_natural_source_type_list.extend(el_cells)
                                if any(x is not None for x in cell_natural_source_type_list):
                                    cns = emdb.cell_source_type()
                                    set_sup_mol_nat_src(cns, cell_sup_mol, const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, sup_mol_dict_in, cell_dict)

                        # set up the supramolecule specific tagname explicitly
                        # as DSgenerate doesn't provide it
                        cell_sup_mol.original_tagname_ = "cell_supramolecule"
                        set_base_sup_mol(cell_sup_mol, sup_in, sup_mol_id_in, sample)
                        # element 1
                        set_el_natural_source(cell_sup_mol, sup_mol_id_in, sup_mol_dicts)

                    sup_list_in = self.cif.get(const.EM_ENTITY_ASSEMBLY, None)
                    if sup_list_in is not None:
                        for sup_in in sup_list_in:
                            # get id from CIF: _em_entity_assembly.id 1
                            sup_mol_id_in = get_cif_value(const.K_ID, const.EM_ENTITY_ASSEMBLY, sup_in)
                            if sup_mol_id_in is not None:
                                # CIF:  _em_entity_assembly.type - can be:
                                # (RIBOSOME or COMPLEX), VIRUS, ORGANELLE OR CELLULAR COMPONENT,
                                # TISSUE, CELL
                                sup_type = get_cif_value("type", const.EM_ENTITY_ASSEMBLY, sup_in)
                                if sup_type is not None:
                                    if sup_type in ["RIBOSOME", "COMPLEX"]:
                                        complex_sup_mol = emdb.complex_supramolecule_type()
                                        set_complex_supramolecule_type(complex_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample)
                                        if complex_sup_mol.has__content():
                                            sup_list.add_supramolecule(complex_sup_mol)
                                    elif sup_type == "VIRUS":
                                        virus_sup_mol = emdb.virus_supramolecule_type()
                                        set_virus_supramolecule_type(virus_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample)
                                        if virus_sup_mol.has__content():
                                            sup_list.add_supramolecule(virus_sup_mol)
                                    elif sup_type == "ORGANELLE OR CELLULAR COMPONENT":
                                        org_or_cell_sup_mol = emdb.organelle_or_cellular_component_supramolecule_type()
                                        set_orgorcell_supmol_type(org_or_cell_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample)
                                        if org_or_cell_sup_mol.has__content():
                                            sup_list.add_supramolecule(org_or_cell_sup_mol)
                                    elif sup_type == "TISSUE":
                                        tissue_sup_mol = emdb.tissue_supramolecule_type()
                                        set_tissue_supramolecule_type(tissue_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample)
                                        if tissue_sup_mol.has__content():
                                            sup_list.add_supramolecule(tissue_sup_mol)
                                    elif sup_type == "CELL":
                                        cell_sup_mol = emdb.cell_supramolecule_type()
                                        set_cell_supramolecule_type(cell_sup_mol, sup_in, sup_mol_id_in, sup_mol_dicts, sample)
                                        if cell_sup_mol.has__content():
                                            sup_list.add_supramolecule(cell_sup_mol)
                                    else:
                                        txt = u"Supramolecule type not implemented. (_em_entity_assembly.type) is (%s)" % sup_type
                                        self.current_entry_log.error_logs.append(
                                            self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                        )
                                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                                else:
                                    # supramolecule type is None
                                    txt = (
                                        u"Supramolecule type for supramolecule with id=(%s) is not given in the cif file. This supramolecule cannot be written into the output XML file."
                                        % sup_mol_id_in
                                    )
                                    self.current_entry_log.error_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                    )
                                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                            else:
                                # supramolecule id is None!
                                txt = u"Supramolecule id (_em_entity_assembly.id) missing in (%s)." % sup_in
                                self.current_entry_log.error_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                )
                                self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                    else:
                        txt = u"CIF category (_em_entity_assembly) is missing."
                        self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                # Create a natural source dictionary with the supramolecule as the key
                # order by _em_entity_assembly_naturalsource.entity_assembly_id
                nat_src_dict_in = make_list_of_dicts(const.EM_ENTITY_ASSEMBLY_NATURALSOURCE, const.K_ENTITY_ASSEMBLY_ID)
                # Create a recombinant expression dictionary with the supramolecule as the key
                rec_exp_dict_in = make_list_of_dicts(const.EM_ENTITY_ASSEMBLY_RECOMBINANT, const.K_ENTITY_ASSEMBLY_ID)
                # Molecular weight dictionary with supramolecule as the key
                # ordered by em_entity_assembly_molwt.entity_assembly_id
                s_mol_wt_dict_in = make_dict(const.EM_ENTITY_ASSEMBLY_MOLWT, const.K_ENTITY_ASSEMBLY_ID)
                # Virus dictionary with supramolecule as key
                # ordred by _em_virus_entity.entity_assembly_id
                virus_dict_in = make_dict(const.EM_VIRUS_ENTITY, const.K_ENTITY_ASSEMBLY_ID)
                # Virus natural host dictionary with supramolecule as key -> dictionary of lists
                # ordered by _em_virus_natural_host.entity_assembly_id
                virus_nat_host_dict_in = make_list_of_dicts(const.EM_VIRUS_NATURAL_HOST, const.K_ENTITY_ASSEMBLY_ID)
                # Virus shell dictionary with supramolecule ID as key -> dictionary of lists
                virus_shell_dict_in = make_list_of_dicts(const.EM_VIRUS_SHELL, const.K_ENTITY_ASSEMBLY_ID)
                sup_mol_dicts = {
                    "nat_src_dict_in": nat_src_dict_in,
                    "rec_exp_dict_in": rec_exp_dict_in,
                    "s_mol_wt_dict_in": s_mol_wt_dict_in,
                    "virus_dict_in": virus_dict_in,
                    "virus_nat_host_dict_in": virus_nat_host_dict_in,
                    "virus_shell_dict_in": virus_shell_dict_in,
                }

                sup_list = emdb.supramolecule_listType()
                set_supramolecule_list(sup_list, sup_mol_dicts)
                sample.set_supramolecule_list(sup_list)

            def set_el_macromolecule_list(sample):
                """
                XSD: <xs:element name="macromolecule_list" type="macromolecule_list_type" minOccurs="0">
                CIF: _entity_poly.type (shows the types of micromolecules):
                    cyclic-pseudo-peptide, other, peptide nucleic acid,
                    polydeoxyribonucleotide
                    polydeoxyribonucleotide/polyribonucleotide hybrid
                    polypeptide(D), polypeptide(L), polyribonucleotide
                    polysaccharide(D), polysaccharide(L)
                """

                def set_mol_weight(mol, ent_in):
                    """
                    Set molecular weight object of mol from values in ent_in

                    Parameters:
                    @param mol: molecule object with a set_molecular_weight method
                    @param ent_in: cif entity
                    XSD: <xs:complexType name="molecular_weight_type"> is
                    .. a sequence of 3 of elements
                    """

                    def set_el_experimental(mol_weight, ent_in):
                        """
                        XSD: <xs:element name="experimental" minOccurs="0">
                        CIF: _entity.formula_weight 28088.863
                        """
                        theo_wt_in = get_cif_value("formula_weight", const.ENTITY, ent_in)
                        if theo_wt_in is not None:
                            set_cif_value(
                                mol_weight.set_theoretical,
                                "formula_weight",
                                const.ENTITY,
                                cif_list=ent_in,
                                constructor=emdb.theoreticalType,
                                units=const.U_MDA,
                                fmt=lambda x: float(x) * const.DA2MDA,
                            )
                            return True
                        else:
                            return False

                    def set_el_theoretical(mol_weight, ent_in):
                        """
                        XSD: <xs:element name="theoretical" minOccurs="0">
                        CIF: _entity.pdbx_formula_weight_exptl
                        """
                        exp_wt_in = get_cif_value("pdbx_formula_weight_exptl", const.ENTITY, ent_in)
                        if exp_wt_in is not None:
                            set_cif_value(
                                mol_weight.set_experimental,
                                "pdbx_formula_weight_exptl",
                                const.ENTITY,
                                cif_list=ent_in,
                                constructor=emdb.experimentalType,
                                units=const.U_MDA,
                                fmt=lambda x: float(x) * const.DA2MDA,
                            )
                            return True
                        else:
                            return False

                    def set_el_method():
                        """
                        XSD:  <xs:element name="method" type="xs:token" minOccurs="0"/>
                        CIF: ????
                        code ???
                        """

                    mol_weight = emdb.molecular_weight_type()
                    # element 1
                    ex_value_set = set_el_experimental(mol_weight, ent_in)
                    # element 2
                    th_value_set = set_el_theoretical(mol_weight, ent_in)
                    # element 3
                    set_el_method()
                    if ex_value_set or th_value_set:
                        mol.set_molecular_weight(mol_weight)

                def set_mol_base_source(macromol_src, cif_category, src_in):
                    """
                    Creates the base source element for macromolecules

                    Parameters:
                    @param macromol_src: macromolecule source object of emdb.molecule_source_type()
                    @param cif_category: contains natural source info (entity_src_nat, entity_src_gen or pdbx_entity_src_syn (not yet used in mmCIF))
                    @param src_in: source in cif: natural, genetically modified or synthetic
                    XSD: <xs:complexType name="base_source_type">
                    """
                    set_base_source_type(macromol_src, cif_category, src_in)

                def make_mol_src(cif_category, src_in, flag_dict):
                    """
                    Make and return molecule source element

                    Parameters:
                    @param cif_category: contains the source info
                                         (entity_src_nat, entity_src_gen or pdbx_entity_src_syn)
                    @param src_in: source in cif: natural, genetically modified or synthetic
                    @param flag_dict:add_organ, add_tissue, add_cell,
                                     add_organelle, add_cellular_location - true if required
                    @return: molecule source object of emdb.molecule_source_type
                    XSD: <xs:complexType name="macromolecule_source_type"> has
                    .. a base of base_source_type
                    .. a sequence of 5 elements
                    """

                    def set_el_organ(macromol_src, cif_category, src_in):
                        """
                        XSD: <xs:element name="organ" type="xs:token" minOccurs="0"/>
                        CIF: _entity_src_nat.pdbx_organ
                        CIF: _entity_src_gen.pdbx_gene_src_organ
                        """
                        a_dict = {const.ENTITY_SRC_NAT: "pdbx_organ", const.ENTITY_SRC_GEN: "pdbx_gene_src_organ", const.PDBX_ENTITY_SRC_SYN: None}
                        cif_key = a_dict.get(cif_category, None)
                        if cif_key is not None:
                            set_cif_value(macromol_src.set_organ, cif_key, cif_category, cif_list=src_in)

                    def set_el_tissue(macromol_src, cif_category, src_in):
                        """
                        XSD: <xs:element name="tissue" type="xs:token" minOccurs="0">
                        CIF: _entity_src_nat.tissue
                        CIF: _entity_src_gen.gene_src_tissue
                        """
                        a_dict = {const.ENTITY_SRC_NAT: "tissue", const.ENTITY_SRC_GEN: "gene_src_tissue", const.PDBX_ENTITY_SRC_SYN: None}
                        cif_key = a_dict.get(cif_category, None)
                        if cif_key is not None:
                            set_cif_value(macromol_src.set_tissue, cif_key, cif_category, cif_list=src_in)

                    def set_el_cell(macromol_src, cif_category, src_in):
                        """
                        XSD: <xs:element name="cell" type="xs:token" minOccurs="0">
                        CIF: _entity_src_nat.pdbx_cell
                        CIF: _entity_src_gen.pdbx_gene_src_cell
                        """
                        a_dict = {const.ENTITY_SRC_NAT: "pdbx_cell", const.ENTITY_SRC_GEN: "pdbx_gene_src_cell", const.PDBX_ENTITY_SRC_SYN: None}
                        cif_key = a_dict.get(cif_category, None)
                        if cif_key is not None:
                            set_cif_value(macromol_src.set_cell, cif_key, cif_category, cif_list=src_in)

                    def set_el_organelle(macromol_src, cif_category, src_in):
                        """
                        XSD: <xs:element name="organelle" type="xs:token" minOccurs="0">
                        CIF: _entity_src_nat.pdbx_organelle
                        CIF: _entity_src_gen.pdbx_gene_src_organelle
                        """
                        a_dict = {const.ENTITY_SRC_NAT: "pdbx_organelle", const.ENTITY_SRC_GEN: "pdbx_gene_src_organelle", const.PDBX_ENTITY_SRC_SYN: None}
                        cif_key = a_dict.get(cif_category, None)
                        if cif_key is not None:
                            set_cif_value(macromol_src.set_organelle, cif_key, cif_category, cif_list=src_in)

                    def set_el_cellular_location(macromol_src, cif_category, src_in):
                        """
                        XSD: <xs:element name="cellular_location" type="xs:token" minOccurs="0">
                        CIF: _entity_src_nat.pdbx_cellular_location
                        CIF: _entity_src_gen.pdbx_gene_src_cellular_location
                        """
                        a_dict = {const.ENTITY_SRC_NAT: "pdbx_cellular_location", const.ENTITY_SRC_GEN: "pdbx_gene_src_cellular_location", const.PDBX_ENTITY_SRC_SYN: None}
                        cif_key = a_dict.get(cif_category, None)
                        if cif_key is not None:
                            set_cif_value(macromol_src.set_cellular_location, cif_key, cif_category, cif_list=src_in)

                    macromol_src = emdb.macromolecule_source_type()
                    # base
                    set_mol_base_source(macromol_src, cif_category, src_in)

                    if flag_dict["add_organ"]:
                        # element 1
                        set_el_organ(macromol_src, cif_category, src_in)
                    if flag_dict["add_tissue"]:
                        # element 2
                        set_el_tissue(macromol_src, cif_category, src_in)
                    if flag_dict["add_cell"]:
                        # element 3
                        set_el_cell(macromol_src, cif_category, src_in)
                    if flag_dict["add_organelle"]:
                        # element 4
                        set_el_organelle(macromol_src, cif_category, src_in)
                    if flag_dict["add_cellular_location"]:
                        # element 5
                        set_el_cellular_location(macromol_src, cif_category, src_in)

                    return macromol_src

                def set_source_in_base_mol(mol, ent_src_dict=None, cif_category=None):
                    """
                    Set base parameters of base_macromolecule_type

                    Parameters:
                    @param mol: object wrapping macromolecule element - this will be updated
                    @param ent_src_dict: dictionary of cif natural,
                        genetically modified or synthetic sources keyed by entity_id
                    @param cif_category: contains the source info
                    XSD: <xs:element name="natural_source" type="molecule_source_type" minOccurs="0"/>
                    """
                    if ent_src_dict is not None and cif_category is not None:
                        mol_src_dict = {"add_organ": True, "add_tissue": True, "add_cell": True, "add_organelle": True, "add_cellular_location": True}
                        nat_src = make_mol_src(cif_category, ent_src_dict, mol_src_dict)
                        if nat_src.has__content():
                            mol.set_natural_source(nat_src)

                def set_base_mol_with_dict(mol, ent_in, ent_id_in, mol_dicts):
                    """
                    Set the dictionary and cif category for base_macromolecule_type

                    Parameters:
                    @param mol: object wrapping macromolecule element - this will be updated
                    @param ent_in: cif entity category dictionary
                    @param ent_id_in: cif entity_id
                    @param mol_dicts: a dictionary of the following dictionaries:
                    entSrcDict - natural source dictionary
                    ent_src_gen_dict - genetically modified source dictionary
                    ent_src_syn_dict - synthetic dictionary source dictionary

                    Molecule source is set depending in which dictionary the entry ID is.
                    XSD: <xs:complexType name="base_macromolecule_type"> has
                    .. 3 attributes and
                    .. 7 elements
                    """

                    def set_attr_id(mol, ent_in):
                        """
                        XSD: <xs:attribute name="macromolecule_id" type="xs:positiveInteger" use="required"/>
                        CIF: _entity.id 1
                        """
                        set_cif_value(mol.set_macromolecule_id, "id", const.ENTITY, cif_list=ent_in, fmt=int)

                    def set_attr_mutant():
                        """
                        XSD: <xs:attribute name="mutant" type="xs:boolean"/>
                        CIF: _entity.pdbx_mutation ?
                        N-M code ????
                        """

                    def set_attr_chimera():
                        """
                        XSD: <xs:attribute name="chimera" type="xs:boolean"/>
                        CIF: ?????
                        N-M code ????
                        """

                    def set_el_name(mol, ent_in):
                        """
                        XSD: <xs:element name="name" type="sci_name_type">
                        CIF: _entity.pdbx_description uL2
                        """
                        set_cif_value(mol.set_name, "pdbx_description", const.ENTITY, cif_list=ent_in, constructor=emdb.sci_name_type)

                    def set_el_natural_source(mol, ent_id_in, mol_dicts):
                        """
                        XSD: <xs:element name="natural_source" type="macromolecule_source_type" minOccurs="0"/>
                        """
                        if ent_id_in in mol_dicts["ent_src_nat_dict"]:
                            ent_src_nat_dict = mol_dicts["ent_src_nat_dict"]
                            the_dict = ent_src_nat_dict[ent_id_in]
                            set_source_in_base_mol(mol, ent_src_dict=the_dict, cif_category=const.ENTITY_SRC_NAT)
                        elif ent_id_in in mol_dicts["ent_src_gen_dict"]:
                            ent_src_gen_dict = mol_dicts["ent_src_gen_dict"]
                            the_dict = ent_src_gen_dict[ent_id_in]
                            set_source_in_base_mol(mol, ent_src_dict=the_dict, cif_category=const.ENTITY_SRC_GEN)
                        elif ent_id_in in mol_dicts["ent_src_syn_dict"]:
                            ent_src_syn_dict = mol_dicts["ent_src_syn_dict"]
                            the_dict = ent_src_syn_dict[ent_id_in]
                            # not implemented in mmCif yet
                            set_source_in_base_mol(mol, ent_src_dict=the_dict, cif_category=const.PDBX_ENTITY_SRC_SYN)
                        else:
                            # shouldn't happen but it is at the moment
                            pass

                    def set_el_molecular_weight(mol, ent_in):
                        """
                        XSD: <xs:element name="molecular_weight" type="molecular_weight_type" minOccurs="0">
                        CIF: _entity.formula_weight
                        """
                        set_mol_weight(mol, ent_in)

                    def set_el_details(mol, ent_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _entity.details
                        """
                        set_cif_value(mol.set_details, "details", const.ENTITY, cif_list=ent_in)

                    def set_el_number_of_copies(mol, ent_in):
                        """
                        XSD: <xs:element name="number_of_copies" type="pos_int_or_string_type" minOccurs="0"/>
                        CIF: _entity.pdbx_number_of_molecules
                        """
                        set_cif_value(mol.set_number_of_copies, "pdbx_number_of_molecules", const.ENTITY, cif_list=ent_in, fmt=int)

                    def set_el_oligomeric_state():
                        """
                        XSD: <xs:element name="oligomeric_state" type="pos_int_or_string_type" minOccurs="0"/>
                        CIF: ?????? Deprecated?
                        """

                    def set_el_recombinant_exp_flag():
                        """
                        XSD: <xs:element name="recombinant_exp_flag" type="xs:boolean" maxOccurs="1" minOccurs="0">
                        Deprecated 2014/12/3
                        """

                    # attribute 1
                    set_attr_id(mol, ent_in)
                    # attribute 2
                    set_attr_mutant()
                    # attribute 3
                    set_attr_chimera()
                    # element 1
                    set_el_name(mol, ent_in)
                    # element 2
                    set_el_natural_source(mol, ent_id_in, mol_dicts)
                    # element 3
                    set_el_molecular_weight(mol, ent_in)
                    # element 4
                    set_el_details(mol, ent_in)
                    # element 1
                    set_el_number_of_copies(mol, ent_in)
                    # element 6
                    set_el_oligomeric_state()
                    # element 7
                    set_el_recombinant_exp_flag()

                def set_mol_seq(mol, ent_poly_in, ent_id_in, ent_ref_dict):
                    """
                    Set the sequence object of mol using entity_poly as input,
                    and also the external references using pdbx_struct_ref_seq_depositor_info

                    Parameters:
                    @param mol: molecule object with a set_sequence method
                    @param ent_poly_in: cif entity_poly
                    @param ent_id_in: cif entity id
                    @param ent_ref_dict: dictionary of pdbx_struct_ref_seq_depositor_info keyed by entity_id (dictionary of lists)
                    XSD: <xs:element name="sequence"> is
                    .. a sequence of 3 elements
                    """

                    def set_el_string(seq, ent_poly_in):
                        """
                        XSD: <xs:element name="string" type="xs:token" minOccurs="0">
                        CIF: _entity_poly.pdbx_seq_one_letter_code
                        MGRVIRGQRKGAGSVFRAHVKHRKGAARLRAVDFAERHGYIKGIVKDIIHDPGRGAPLAKVVFRDPYRFKKRTE
                        """
                        set_cif_value(seq.set_string, "pdbx_seq_one_letter_code", const.ENTITY_POLY, cif_list=ent_poly_in)

                    def set_el_discrepancy_list():
                        """
                        XSD: <xs:element name="discrepancy_list" minOccurs="0">
                        CIF: ??
                        """

                    def set_el_external_references(seq, ent_id_in, ent_ref_dict):
                        """
                        XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                        """

                        def set_external_references_type(cross_ref, rel_in):
                            """
                            XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0"> has
                            .. 1 attribute
                            .. an extension of xs:token
                            FOr Map-model entries:
                            CIF: _struct_ref.db_name UNP
                            CIF: _struct_ref.db_code ?
                            For map only entries:
                            CIF: _pdbx_struct_ref_seq_depositor_info.db_name UNP
                            CIF: _pdbx_struct_ref_seq_depositor_info.db_accession ?
                            """
                            db2_in = assert_get_value(const.DATABASE_2, self.cif)
                            db_id_dict = make_list_of_dicts(const.DATABASE_2, "database_id", db2_in, 2)
                            if "PDB" in db_id_dict:
                                dict_db2_in = {t[0]: t[1] for t in db2_in}
                                pdb_id = dict_db2_in[('_database_2.database_id', 'PDB')][1]
                                if pdb_id is not None:
                                    db_code = get_cif_value("db_code", const.STRUCT_REF, rel_in)
                                    if db_code is not None:
                                        db_in = get_cif_value("db_name", const.STRUCT_REF, rel_in)
                                        if db_in == "UNP":
                                            set_cif_value(cross_ref.set_type, "db_name", const.STRUCT_REF, cif_list=rel_in, cif_value="UNIPROTKB")
                                            set_cif_value(cross_ref.set_valueOf_, "pdbx_db_accession", const.STRUCT_REF, cif_list=rel_in)
                                        elif db_in == "GB":
                                            set_cif_value(cross_ref.set_type, "db_name", const.STRUCT_REF, cif_list=rel_in, cif_value="GENBANK")
                                            set_cif_value(cross_ref.set_valueOf_, "db_code", const.STRUCT_REF, cif_list=rel_in)
                                        elif db_in is not None and db_in != "PDB":
                                            set_cif_value(cross_ref.set_type, "db_name", const.STRUCT_REF, cif_list=rel_in)
                                            set_cif_value(cross_ref.set_valueOf_, "db_code", const.STRUCT_REF, cif_list=rel_in)
                            if "PDB" not in db_id_dict:
                                db_code = get_cif_value("db_accession", const.PDBX_DEPOSITOR_INFO, rel_in)
                                if db_code is not None:
                                    db_in = get_cif_value("db_name", const.PDBX_DEPOSITOR_INFO, rel_in)
                                    if db_in == "UNP":
                                        set_cif_value(cross_ref.set_type, "db_name", const.PDBX_DEPOSITOR_INFO, cif_list=rel_in, cif_value="UNIPROTKB")
                                        set_cif_value(cross_ref.set_valueOf_, "db_accession", const.PDBX_DEPOSITOR_INFO, cif_list=rel_in)
                                    elif db_in == "GB":
                                        set_cif_value(cross_ref.set_type, "db_name", const.PDBX_DEPOSITOR_INFO, cif_list=rel_in, cif_value="GENBANK")
                                        set_cif_value(cross_ref.set_valueOf_, "db_accession", const.PDBX_DEPOSITOR_INFO, cif_list=rel_in)
                                    elif db_in is not None and db_in != "PDB":
                                        set_cif_value(cross_ref.set_type, "db_name", const.PDBX_DEPOSITOR_INFO, cif_list=rel_in)
                                        set_cif_value(cross_ref.set_valueOf_, "db_accession", const.PDBX_DEPOSITOR_INFO, cif_list=rel_in)

                        if ent_id_in in ent_ref_dict:
                            ent_ref_list_in = ent_ref_dict[ent_id_in]
                            for rel_in in ent_ref_list_in:
                                cross_ref = emdb.external_referencesType()
                                set_external_references_type(cross_ref, rel_in)
                                if cross_ref.has__content():
                                    seq.add_external_references(cross_ref)

                    seq = emdb.sequenceType()
                    # element 1
                    set_el_string(seq, ent_poly_in)
                    # element 2
                    set_el_discrepancy_list()
                    # element 3
                    set_el_external_references(seq, ent_id_in, ent_ref_dict)

                    if seq.has__content():
                        mol.set_sequence(seq)

                def set_mol_rec_exp(mol, ent_id_in, mol_rec_exp_dict):
                    """
                    Set recombinant expression for a macromolecule

                    Parameters:
                    @param mol: macromolecule object that will be updated with recombinant source info
                    @param ent_id_in: cif id for entity corresponding to macromolecule
                    @param mol_rec_exp_dict: dictionary for this macromolecule
                    XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" minOccurs="0"/>
                    XSD: <xs:complexType name="recombinant_source_type"> has
                    .. 1  attribute and
                    .. 5 elements
                    """
                    if mol_rec_exp_dict is not None:
                        if ent_id_in in mol_rec_exp_dict:
                            r_exp_in = mol_rec_exp_dict[ent_id_in]
                            r_exp = emdb.recombinant_source_type()
                            set_recombinant_source_type(r_exp, r_exp_in)
                            if r_exp.has__content():
                                mol.set_recombinant_expression(r_exp)
                    else:
                        txt = u"The dictionary for recombinant expression for macromolecule id (%s) not found." % ent_id_in
                        self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                        self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)

                def set_mol_syn_src(syn_src, mol, cif_category, src_dict_in, flags_dict):
                    """
                    Set synthetic source for a macromolecule

                    Parameters:
                    @param syn_src: synthetic_source object of different type,
                           depending on the macromolecule
                    @param mol: macromolecule object
                    @param cif_category: contains synthetic source info
                           (pdbx_entity_src_syn for macromolecule)
                    @param src_dict_in: source dictionary keyed by macromolecule id,
                                        each value is a list of synthetic sources
                    @param flags_dict: a dictionary containing boolean values for the following:
                           add_syn_src: flag for adding or setting synthetic source
                           (only true used at the moment)
                           the following true if required:
                           add_organ, add_tissue, add_cell, add_organelle, add_cellular_location
                    Macromolecule synthetic source is
                    .. an extension of base="base_source_type" and
                    .. a sequence of 5 possible elements
                    """

                    def set_el_organ(syn_src, cif_category, src_dict_in):
                        """
                        XSD: <xs:element name="organ" type="xs:token" minOccurs="0"/>
                        CIF: _em_entity_assembly_naturalsource.organ .
                        """
                        set_cif_value(syn_src.set_organ, "organ", cif_category, cif_list=src_dict_in)

                    def set_el_tissue(syn_src, cif_category, src_dict_in):
                        """
                        XSD: <xs:element name="tissue" type="xs:token" minOccurs="0">
                        CIF: _em_entity_assembly_naturalsource.tissue ?
                        """
                        set_cif_value(syn_src.set_tissue, "tissue", cif_category, cif_list=src_dict_in)

                    def set_el_cell(syn_src, cif_category, src_dict_in):
                        """
                        XSD: <xs:element name="cell" type="xs:token" minOccurs="0">
                        CIF: _em_entity_assembly_naturalsource.cell    ?
                        """
                        set_cif_value(syn_src.set_cell, "cell", cif_category, cif_list=src_dict_in)

                    def set_el_organelle(syn_src, cif_category, src_dict_in):
                        """
                        XSD: <xs:element name="organelle" type="xs:token" minOccurs="0">
                        CIF: _em_entity_assembly_naturalsource.organelle ?
                        """
                        set_cif_value(syn_src.set_organelle, "organelle", cif_category, cif_list=src_dict_in)

                    def set_el_cellular_location(syn_src, cif_category, src_dict_in):
                        """
                        XSD: <xs:element name="cellular_location" type="xs:token" minOccurs="0">
                        CIF: _em_entity_assembly_naturalsource.cellular_location ?
                        """
                        set_cif_value(syn_src.set_cellular_location, "cellular_location", cif_category, cif_list=src_dict_in)

                    for mol_syn_src_in in src_dict_in:
                        set_base_source_type(mol, cif_category, mol_syn_src_in)

                        if flags_dict["add_organ"]:
                            set_el_organ(syn_src, cif_category, mol_syn_src_in)
                        if flags_dict["add_tissue"]:
                            set_el_tissue(syn_src, cif_category, mol_syn_src_in)
                        if flags_dict["add_cell"]:
                            set_el_cell(syn_src, cif_category, mol_syn_src_in)
                        if flags_dict["add_organelle"]:
                            set_el_organelle(syn_src, cif_category, mol_syn_src_in)
                        if flags_dict["add_cellular_location"]:
                            set_el_cellular_location(syn_src, cif_category, mol_syn_src_in)

                        if syn_src.has__content():
                            if flags_dict["add_syn_src"]:
                                mol.add_synthetic_source(syn_src)
                            else:
                                mol.set_synthetic_source(syn_src)
                            break

                def set_rna_macromolecule_type(rna_mol, ent_in, ent_id_in, ent_poly_in, src_dicts):
                    """
                    XSD: <xs:complexType name="rna_macromolecule_type"> has
                    .. a base of base_macromolecule_type and
                    .. a sequence of 5 elements
                    """

                    def set_el_sequence(rna_mol, ent_poly_in, ent_id_in, src_dict):
                        """
                        XSD:  <xs:element name="sequence">
                        """
                        set_mol_seq(rna_mol, ent_poly_in, ent_id_in, src_dict)

                    def set_el_classification():
                        """
                        XSD: <xs:element name="classification" minOccurs="0">
                        CIF: ??
                        """

                    def set_el_structure():
                        """
                        XSD: <xs:element name="structure" type="xs:token" minOccurs="0">
                        CIF: ??
                        """

                    def set_el_synthetic_flag():
                        """
                        XSD: <xs:element name="synthetic_flag" type="xs:boolean" minOccurs="0">
                        CIF: ??
                        """

                    def set_el_ec_number(rna_mol, ent_in):
                        """
                        XSD:  <xs:element name="ec_number" maxOccurs="unbounded" minOccurs="0">
                        CIF: _entity.pdbx_ec ? (Enzyme classification. Format: "EC 3.4.11.4")
                        """
                        ec_num_in = get_cif_value("pdbx_ec", const.ENTITY, cif_list=ent_in)
                        if ec_num_in is not None:
                            if ec_num_in.find(",") != -1:
                                # there is more than one EC number: split it and write the first bit
                                ec_num = ec_num_in.split(",")[0]
                                set_cif_value(rna_mol.add_ec_number, "pdbx_ec", const.ENTITY, cif_list=ent_in, cif_value=ec_num)
                                txt = u"(_entity.pdbx_ec) is given (%s) but the value (%s) is set." % (ec_num_in, ec_num)
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                            else:
                                set_cif_value(rna_mol.add_ec_number, "pdbx_ec", const.ENTITY, cif_list=ent_in)

                    ent_ref_dict = src_dicts["ent_ref_dict"]
                    # set up the macromolecule specific tagname explicitly
                    # as DSgenerate doesn't provide it
                    rna_mol.original_tagname_ = "rna"
                    # base
                    set_base_mol_with_dict(rna_mol, ent_in, ent_id_in, src_dicts)
                    # element 1
                    set_el_sequence(rna_mol, ent_poly_in, ent_id_in, ent_ref_dict)
                    # element 2
                    set_el_classification()
                    # element 3
                    set_el_structure()
                    # element 4
                    set_el_synthetic_flag()
                    # element 5
                    set_el_ec_number(rna_mol, ent_in)

                def set_dna_macromolecule_type(dna_mol, ent_in, ent_id_in, ent_poly_in, src_dicts):
                    """
                    XSD: <xs:complexType name="dna_macromolecule_type"> has
                    .. a base of base_macromolecule_type and
                    .. a sequence of 5 elements
                    """

                    def set_el_sequence(dna_mol, ent_poly_in, ent_id_in, src_dict):
                        """
                        XSD: <xs:element name="sequence">
                        CIF: in ent_ref_dict
                        """
                        set_mol_seq(dna_mol, ent_poly_in, ent_id_in, src_dict)

                    def set_el_classification(dna_mol):
                        """
                        XSD: <xs:element name="classification" minOccurs="0">
                        CIF:
                        """
                        set_cif_value(dna_mol.set_classification, cif_value="DNA")

                    def set_el_structure():
                        """
                        XSD: <xs:element name="structure" type="xs:token" minOccurs="0">
                        Deprecated 2014/11/29
                        """

                    def set_el_synthetic_flag():
                        """
                        XSD: <xs:element name="synthetic_flag" type="xs:boolean" minOccurs="0">
                        Deprecated 2014/11/29
                        """

                    def set_el_synthetic_source(dna_mol, ent_id_in, syn_src_dict_in):
                        """
                        XSD: <xs:element name="synthetic_source" type="macromolecule_source_type" minOccurs="0"/>
                        CIF:
                        """
                        if ent_id_in in syn_src_dict_in:
                            syn_mol_dict_in = syn_src_dict_in[ent_id_in]
                            dna_dict = {"add_syn_src": True, "add_organ": True, "add_tissue": True, "add_cell": True, "add_organelle": True, "add_cellular_location": True}
                            dna_synthetic_source_type_list = []
                            attr_ncbis = []
                            el_organisms = []
                            el_strains = []
                            el_organs = []
                            el_tissues = []
                            el_cells = []
                            el_organelles = []
                            el_cellular_locations = []
                            for syn_src_in in syn_mol_dict_in:
                                attr_ncbi = get_cif_value("ncbi_tax_id", const.PDBX_ENTITY_SRC_SYN, cif_list=syn_src_in)
                                attr_ncbis.extend([attr_ncbi])
                                el_organism = get_cif_value("organism", const.PDBX_ENTITY_SRC_SYN, cif_list=syn_src_in)
                                el_organisms.extend([el_organism])
                                el_strain = get_cif_value("strain", const.PDBX_ENTITY_SRC_SYN, cif_list=syn_src_in)
                                el_strains.extend([el_strain])
                                el_organ = get_cif_value("organ", const.PDBX_ENTITY_SRC_SYN, cif_list=syn_src_in)
                                el_organs.extend([el_organ])
                                el_tissue = get_cif_value("tissue", const.PDBX_ENTITY_SRC_SYN, cif_list=syn_src_in)
                                el_tissues.extend([el_tissue])
                                el_cell = get_cif_value("cell", const.PDBX_ENTITY_SRC_SYN, cif_list=syn_src_in)
                                el_cells.extend([el_cell])
                                el_organelle = get_cif_value("organelle", const.PDBX_ENTITY_SRC_SYN, cif_list=syn_src_in)
                                el_organelles.extend([el_organelle])
                                el_cellular_location = get_cif_value("cellular_location", const.PDBX_ENTITY_SRC_SYN, cif_list=syn_src_in)
                                el_cellular_locations.extend([el_cellular_location])
                            dna_synthetic_source_type_list.extend(attr_ncbis)
                            dna_synthetic_source_type_list.extend(el_organisms)
                            dna_synthetic_source_type_list.extend(el_strains)
                            dna_synthetic_source_type_list.extend(el_organs)
                            dna_synthetic_source_type_list.extend(el_tissues)
                            dna_synthetic_source_type_list.extend(el_cells)
                            dna_synthetic_source_type_list.extend(el_organelles)
                            dna_synthetic_source_type_list.extend(el_cellular_locations)
                            if any(x is not None for x in dna_synthetic_source_type_list):
                                dna_st = emdb.macromolecule_source_type()
                                set_mol_syn_src(dna_st, dna_mol, const.PDBX_ENTITY_SRC_SYN, syn_src_dict_in, dna_dict)

                    entity_ref_dict = src_dicts["ent_ref_dict"]
                    syn_src_dict_in = src_dicts["ent_src_syn_dict"]
                    dna_mol.original_tagname_ = "dna"
                    # base
                    set_base_mol_with_dict(dna_mol, ent_in, ent_id_in, src_dicts)
                    # element 1
                    set_el_sequence(dna_mol, ent_poly_in, ent_id_in, entity_ref_dict)
                    # element 2
                    set_el_classification(dna_mol)
                    # element 3
                    set_el_structure()
                    # element 4
                    set_el_synthetic_flag()
                    # element 5
                    set_el_synthetic_source(dna_mol, ent_id_in, syn_src_dict_in)

                def set_protein_or_peptide_mol_type(p_mol, ent_in, ent_id_in, ent_poly_in, src_dicts):
                    """
                    XSD: <xs:complexType name="protein_or_peptide_macromolecule_type"> has
                    .. a base of base_macromolecule_type and
                    .. a sequence of 4 elements
                    """

                    def set_el_recombinant_expression(p_mol, ent_id_in, rec_exp_dict_in):
                        """
                        XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" minOccurs="0"/>
                        """
                        set_mol_rec_exp(p_mol, ent_id_in, rec_exp_dict_in)

                    def set_el_enantiomer(p_mol, ent_type_in):
                        """
                        XSD: <xs:element name="enantiomer">
                        """
                        if ent_type_in == const.ENT_POLYPEPTIDE_D:
                            set_cif_value(p_mol.set_enantiomer, cif_value="DEXTRO")
                        else:
                            set_cif_value(p_mol.set_enantiomer, cif_value="LEVO")

                    def set_el_sequence(p_mol, ent_poly_in, ent_id_in, src_dict):
                        """
                        XSD: <xs:element name="sequence">
                        """
                        set_mol_seq(p_mol, ent_poly_in, ent_id_in, src_dict)

                    def set_el_ec_number(p_mol, ent_in):
                        """
                        XSD: <xs:element name="ec_number" maxOccurs="unbounded" minOccurs="0">
                        CIF: _entity.pdbx_ec ? (Enzyme classification. Format: "EC 3.4.11.4")
                        """
                        ec_num_in = get_cif_value("pdbx_ec", const.ENTITY, cif_list=ent_in)
                        if ec_num_in is not None:
                            if ec_num_in.find(",") != -1:
                                # there is more than one EC number: split it and write the first bit
                                ec_num = ec_num_in.split(",")[0]
                                set_cif_value(p_mol.add_ec_number, "pdbx_ec", const.ENTITY, cif_list=ent_in, cif_value=ec_num)
                                txt = u"(_entity.pdbx_ec) is given: (%s) but the value: (%s) is set instead." % (ec_num_in, ec_num)
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                            else:
                                set_cif_value(p_mol.add_ec_number, "pdbx_ec", const.ENTITY, cif_list=ent_in)

                    ent_ref_dict = src_dicts["ent_ref_dict"]
                    # set up the macromolecule specific tagname explicitly
                    # as DSgenerate doesn't provide it
                    p_mol.original_tagname_ = "protein_or_peptide"
                    # base
                    set_base_mol_with_dict(p_mol, ent_in, ent_id_in, src_dicts)
                    # element 1
                    rec_exp_dict_in = get_rec_exp_dict(ent_id_in, src_dicts, ent_poly_in, is_supramolecule=False)
                    set_el_recombinant_expression(p_mol, ent_id_in, rec_exp_dict_in)
                    # element 2
                    set_el_enantiomer(p_mol, ent_type_in)  # pylint: disable=possibly-used-before-assignment
                    # element 3
                    set_el_sequence(p_mol, ent_poly_in, ent_id_in, ent_ref_dict)
                    # element 4
                    set_el_ec_number(p_mol, ent_in)

                def set_saccharide_mol_type(poly_mol, ent_in, ent_id_in, src_dicts):
                    """
                    <xs:complexType name="saccharide_macromolecule_type"> has
                    .. a base of base_macromolecule_type and
                    .. a sequence of 3 elements
                    """

                    def set_el_enantiomer(poly_mol):
                        """
                        XSD: <xs:element name="enantiomer">
                        """
                        if ent_type_in == const.ENT_POLYSACCHARIDE_D:
                            set_cif_value(poly_mol.set_enantiomer, cif_value="DEXTRO")
                        else:
                            set_cif_value(poly_mol.set_enantiomer, cif_value="LEVO")

                    def set_el_formula():
                        """
                        XSD: <xs:element name="formula" type="formula_type" minOccurs="0"/>
                        CIF: ???
                        """

                    def set_el_external_references():
                        """
                        XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                        CIF: ???
                        """

                    # set up the macromolecule specific tagname explicitly
                    # as DSgenerate doesn't provide it
                    poly_mol.original_tagname_ = "saccharide"
                    # base
                    set_base_mol_with_dict(poly_mol, ent_in, ent_id_in, src_dicts)
                    # element 1
                    set_el_enantiomer(poly_mol)
                    # element 2
                    set_el_formula()
                    # element 3
                    set_el_external_references()

                def set_other_macromolecule_type(other_mol, ent_in, ent_id_in, ent_poly_in, src_dicts):
                    """
                    XSD: <xs:complexType name="other_macromolecule_type"> has
                    .. a base of base_macromolecule_type and
                    .. a sequence of 5 elements
                    """

                    def set_el_sequence(other_mol, ent_poly_in, ent_id_in, src_dict):
                        """
                        XSD: <xs:element name="sequence" minOccurs="0">
                        """
                        set_mol_seq(other_mol, ent_poly_in, ent_id_in, src_dict)

                    def set_el_classification():
                        """
                        XSD: <xs:element name="classification" type="xs:token">
                        """
                        other_mol.set_classification(ent_type_in)

                    def set_el_recombinant_expression(other_mol, ent_id_in, rec_exp_dict_in):
                        """
                        XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" minOccurs="0"/>
                        """
                        set_mol_rec_exp(other_mol, ent_id_in, rec_exp_dict_in)

                    def set_el_structure():
                        """
                        XSD: <xs:element name="structure" type="xs:token">
                        Deprecated 2015/11/16; Added only for transfer of 1.9 data
                        """

                    def set_el_synthetic_flag():
                        """
                        XSD: <xs:element name="synthetic_flag" type="xs:boolean">
                        Deprecated 2015/11/16; Added only for transfer of 1.9 data
                        """

                    ent_ref_dict = src_dicts["ent_ref_dict"]
                    # set up the macromolecule specific tagname explicitly
                    # as DSgenerate doesn't provide it
                    other_mol.original_tagname_ = "other_macromolecule"
                    # base
                    set_base_mol_with_dict(other_mol, ent_in, ent_id_in, src_dicts)
                    # element 1
                    set_el_sequence(other_mol, ent_poly_in, ent_id_in, ent_ref_dict)
                    # element 2
                    set_el_classification()
                    # element 3
                    rec_exp_dict_in = get_rec_exp_dict(ent_id_in, src_dicts, ent_poly_in, is_supramolecule=False)
                    set_el_recombinant_expression(other_mol, ent_id_in, rec_exp_dict_in)
                    # element 4
                    set_el_structure()
                    # element 5
                    set_el_synthetic_flag()

                def set_ligand_macromolecule_type(lig_mol, ent_in, ent_id_in, src_dicts):
                    """
                    Ligands are non-polymers and do not feature in any source categories

                    XSD: <xs:complexType name="ligand_macromolecule_type"> has
                    .. a base of base_macromolecule_type and
                    .. a sequence of 3 elements
                    """

                    def set_el_formula(lig_mol, ent_non_poly_dict, ent_id_in):
                        """
                        XSD: <xs:element name="formula" minOccurs="0">
                        CIF: _pdbx_entity_nonpoly.comp_id MG
                        """
                        ent_non_poly_in = ent_non_poly_dict[ent_id_in]
                        if ent_non_poly_in is not None:
                            set_cif_value(lig_mol.set_formula, "comp_id", const.PDBX_ENTITY_NONPOLY, cif_list=ent_non_poly_in)

                    def set_el_external_references():
                        """
                        XSD: <xs:element name="external_references" maxOccurs="unbounded" minOccurs="0">
                        CIF:
                        """

                    def set_el_recombinant_expression():
                        """
                        XSD: <xs:element name="recombinant_expression" type="recombinant_source_type" minOccurs="0"/>
                        CIF: ??
                        There is no set_number_of_copies so why call this? Commented out for now
                        safe_set(get_cif_value('pdbx_number_of_molecules',
                                const.ENTITY, ent_in), mol.set_number_of_copies, int)
                        """

                    # set up the macromolecule specific tagname explicitly
                    # as DSgenerate doesn't provide it
                    lig_mol.original_tagname_ = "ligand"
                    # base
                    set_base_mol_with_dict(lig_mol, ent_in, ent_id_in, src_dicts)
                    # element 1
                    set_el_formula(lig_mol, ent_non_poly_dict, ent_id_in)
                    # element 2
                    set_el_external_references()
                    # element 3
                    set_el_recombinant_expression()

                mol_list = emdb.macromolecule_list_type()

                # Create dictionaries keyed by entity_id
                ent_poly_dict = make_dict(const.ENTITY_POLY, const.K_ENTITY_ID)
                ent_non_poly_dict = make_dict(const.PDBX_ENTITY_NONPOLY, const.K_ENTITY_ID)
                ent_src_nat_dict = make_dict(const.ENTITY_SRC_NAT, const.K_ENTITY_ID)
                ent_src_gen_dict = make_dict(const.ENTITY_SRC_GEN, const.K_ENTITY_ID)
                ent_src_syn_dict = make_dict(const.PDBX_ENTITY_SRC_SYN, const.K_ENTITY_ID)
                ent_ref_dict = make_list_of_dicts(const.PDBX_DEPOSITOR_INFO, const.K_ENTITY_ID)
                struct_ref_dict = make_list_of_dicts(const.STRUCT_REF, const.K_ENTITY_ID)
                ent_ref_dict.update(struct_ref_dict)
                src_dicts = {"ent_src_nat_dict": ent_src_nat_dict, "ent_src_gen_dict": ent_src_gen_dict, "ent_src_syn_dict": ent_src_syn_dict, "ent_ref_dict": ent_ref_dict}
                entity_list_in = self.cif.get(const.ENTITY, None)
                for ent_in in entity_list_in:
                    ent_id_in = get_cif_value(const.K_ID, const.ENTITY, ent_in)
                    if ent_id_in in ent_poly_dict:
                        ent_poly_in = ent_poly_dict[ent_id_in]
                        ent_type_in = get_cif_value("type", const.ENTITY_POLY, ent_poly_in)
                        if ent_type_in == const.ENT_RNA:
                            rna_mol = emdb.rna_macromolecule_type()
                            set_rna_macromolecule_type(rna_mol, ent_in, ent_id_in, ent_poly_in, src_dicts)
                            mol_list.add_macromolecule(rna_mol)
                        elif ent_type_in == const.ENT_DNA:
                            dna_mol = emdb.dna_macromolecule_type()
                            set_dna_macromolecule_type(dna_mol, ent_in, ent_id_in, ent_poly_in, src_dicts)
                            mol_list.add_macromolecule(dna_mol)
                        elif ent_type_in in [const.ENT_POLYPEPTIDE_L, const.ENT_POLYPEPTIDE_D]:
                            p_mol = emdb.protein_or_peptide_macromolecule_type()
                            set_protein_or_peptide_mol_type(p_mol, ent_in, ent_id_in, ent_poly_in, src_dicts)
                            mol_list.add_macromolecule(p_mol)
                        elif ent_type_in in [const.ENT_POLYSACCHARIDE_D, const.ENT_POLYSACCHARIDE_L]:
                            poly_mol = emdb.saccharide_macromolecule_type()
                            set_saccharide_mol_type(poly_mol, ent_in, ent_id_in, src_dicts)
                            mol_list.add_macromolecule(poly_mol)
                        elif ent_type_in in [const.ENT_CYCLIC_PSEUDO_PEPTIDE, const.ENT_DNA_RNA, const.ENT_OTHER, const.ENT_PEPTIDE_NUCLEIC_ACID]:
                            other_mol = emdb.other_macromolecule_type()
                            set_other_macromolecule_type(other_mol, ent_in, ent_id_in, ent_poly_in, src_dicts)
                            mol_list.add_macromolecule(other_mol)
                        else:
                            txt = u"Entity poly type (%s) not recognized. It needs changing to an allowed type." % ent_type_in
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                    elif ent_id_in in ent_non_poly_dict:
                        # Entity non-poly
                        lig_mol = emdb.ligand_macromolecule_type()
                        set_ligand_macromolecule_type(lig_mol, ent_in, ent_id_in, src_dicts)
                        mol_list.add_macromolecule(lig_mol)

                    if mol_list.has__content():
                        sample.set_macromolecule_list(mol_list)

            # element 1
            set_el_name()
            # element 2
            set_el_supramolecule_list(sample)
            # element 3
            set_el_macromolecule_list(sample)

        def set_base_specimen_preparation(sp_id_in, spec_prep_in, specimen):
            """
            Set base_preparation_type elements

            Parameters:
            @param sp_id_in: specimen id
            @param spec_prep_in: cif values for specific specimen preparation
            @param specimen: object for specimen preparation:
                            {tomography_preparation, single_particle_preparation,
                            subtomogram_averaging_preparation, helixal_preparation,
                            crystalography_preparation}
            XSD: <xs:element name="specimen_preparation" type="base_preparation_type"> has
            .. 1 attribute and
            .. 8 elements
            """

            def set_attr_id(specimen, spec_prep_in):
                """
                XSD: <xs:element name="specimen_preparation_id" type="xs:positiveInteger"/>
                CIF: _em_specimen.id
                """
                set_cif_value(specimen.set_preparation_id, const.K_ID, const.EM_SPECIMEN, cif_list=spec_prep_in, fmt=int)

            def set_el_concentration(specimen, spec_prep_in):
                """
                XSD: <xs:element name="concentration" minOccurs="0">
                CIF: _em_specimen.concentration
                """
                set_cif_value(specimen.set_concentration, "concentration", const.EM_SPECIMEN, cif_list=spec_prep_in, constructor=emdb.concentrationType, units=const.U_MG_ML)

            def set_el_buffer(specimen, sp_id_in, buff_dict_in, buff_comp_dict_in):
                """
                XSD:  <xs:element name="buffer" type="buffer_type" minOccurs="0">
                CIF: _
                """

                def set_buffer_type(buff, buff_in, buff_comp_dict_in):
                    """
                    XSD <xs:complexType name="buffer_type"> has
                    .. 3 elements
                    """

                    def set_el_ph(buff, buff_in):
                        """
                        XSD: <xs:element name="ph">
                        CIF: _em_buffer.pH 7.4
                        """
                        set_cif_value(buff.set_ph, "pH", const.EM_BUFFER, cif_list=buff_in, fmt=float, parent_el_req=False)

                    def set_el_component(buff, buff_comp_dict_in):
                        """
                        XSD: <xs:element name="component" maxOccurs="unbounded">
                        CIF: _em_buffer_component
                        """

                        def set_buffer_component_type(buff_comp, buff_comp_in):
                            """
                            XSD: <xs:complexType name="buffer_component_type"> has
                            .. 3 elements
                            """

                            def set_el_concentration(buff_comp, buff_comp_in):
                                """
                                XSD: <xs:element name="concentration" minOccurs="0">
                                CIF: _em_buffer_component.concentration
                                CIF: _em_buffer_component.concentration_units
                                """
                                conc_units = get_cif_value("concentration_units", const.EM_BUFFER_COMPONENT, buff_comp_in)
                                if conc_units is not None:
                                    set_cif_value(
                                        buff_comp.set_concentration,
                                        "concentration",
                                        const.EM_BUFFER_COMPONENT,
                                        cif_list=buff_comp_in,
                                        constructor=emdb.concentrationType,
                                        fmt=float,
                                        units=conc_units,
                                    )
                                else:
                                    txt = u"The value for (_em_buffer_component.concentration_units) is missing. Buffer concentration will not be set."
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)

                            def set_el_formula(buff_comp, buff_comp_in):
                                """
                                XSD: <xs:element name="formula" minOccurs="0">
                                CIF: _em_buffer_component.formula KOAc
                                """
                                set_cif_value(buff_comp.set_formula, "formula", const.EM_BUFFER_COMPONENT, cif_list=buff_comp_in)

                            def set_el_name(buff_comp, buff_comp_in):
                                """
                                XSD: <xs:element name="name" type="xs:token" minOccurs="0">
                                CIF: _em_buffer_component.name 'Potassium acetate'
                                """
                                set_cif_value(buff_comp.set_name, "name", const.EM_BUFFER_COMPONENT, cif_list=buff_comp_in)

                            # element 1
                            set_el_concentration(buff_comp, buff_comp_in)
                            # element 2
                            set_el_formula(buff_comp, buff_comp_in)
                            # element 3
                            set_el_name(buff_comp, buff_comp_in)

                        buff_id_in = get_cif_value(const.K_ID, const.EM_BUFFER, buff_in)
                        if buff_id_in in buff_comp_dict_in:
                            for buff_comp_in in buff_comp_dict_in[buff_id_in]:
                                buff_comp = emdb.buffer_component_type()
                                set_buffer_component_type(buff_comp, buff_comp_in)
                                buff.add_component(buff_comp)

                    def set_el_details(buff, buff_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                        CIF: _em_buffer.details ?
                        """
                        set_cif_value(buff.set_details, "details", const.EM_BUFFER, cif_list=buff_in)

                    # element 1
                    set_el_ph(buff, buff_in)
                    # element 2
                    set_el_component(buff, buff_comp_dict_in)
                    # element 3
                    set_el_details(buff, buff_in)

                if sp_id_in in buff_dict_in:
                    buff_in = buff_dict_in[sp_id_in]
                    buff_list = []
                    el_ph = get_cif_value("pH", const.EM_BUFFER, cif_list=buff_in)
                    buff_id_in = get_cif_value(const.K_ID, const.EM_BUFFER, buff_in)
                    if buff_id_in in buff_comp_dict_in:
                        for buff_comp_in in buff_comp_dict_in[buff_id_in]:
                            conc_units = get_cif_value("concentration_units", const.EM_BUFFER_COMPONENT, buff_comp_in)
                            el_formula = get_cif_value("formula", const.EM_BUFFER_COMPONENT, cif_list=buff_comp_in)
                            el_name = get_cif_value("name", const.EM_BUFFER_COMPONENT, cif_list=buff_comp_in)
                            if any(x is not None for x in [conc_units, el_formula, el_name]):
                                buff_list.extend([conc_units, el_formula, el_name])
                    el_details = get_cif_value("details", const.EM_BUFFER, cif_list=buff_in)
                    buff_list.extend([el_ph, el_details])
                    if any(x is not None for x in buff_list):
                        buff = emdb.buffer_type()
                        set_buffer_type(buff, buff_in, buff_comp_dict_in)
                        specimen.set_buffer(buff)

            def set_el_staining(specimen, sp_id_in, stain_dict_in):
                """
                XSD: <xs:element name="staining" minOccurs="0"> has
                .. 3 elements
                """

                def set_el_type(stain, stain_in):
                    """
                    XSD: <xs:element name="type">
                    CIF: _em_staining.type
                    This is a required cif item. If not given the value 'NEGATIVE' is set
                    """
                    staining_value = get_cif_value("type", const.EM_STAINING, cif_list=stain_in)
                    if staining_value is not None:
                        set_cif_value(stain.set_type, "type", const.EM_STAINING, cif_list=stain_in)
                    else:
                        set_cif_value(stain.set_type, "type", const.EM_STAINING, cif_list=stain_in, cif_value="NEGATIVE")
                        txt = u"No value is found for (_em_staining.type). The value (NEGATIVE) is given."
                        self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                        self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)

                def set_el_material(stain, stain_in):
                    """
                    XSD: <xs:element name="material" type="xs:token">
                    CIF: _em_staining.material
                    """
                    set_cif_value(stain.set_material, "material", const.EM_STAINING, cif_list=stain_in)

                def set_el_details(stain, stain_in):
                    """
                    XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                    CIF: _em_staining.details
                    """
                    set_cif_value(stain.set_details, "details", const.EM_STAINING, cif_list=stain_in)

                if sp_id_in in stain_dict_in:
                    stain_in = stain_dict_in[sp_id_in]
                    stain = emdb.stainingType()
                    # element 1
                    set_el_type(stain, stain_in)
                    # element 2
                    set_el_material(stain, stain_in)
                    # element 3
                    set_el_details(stain, stain_in)
                    if stain.has__content():
                        specimen.set_staining(stain)

            def set_el_sugar_embedding(specimen, sp_id_in, embed_dict_in):
                """
                XSD: <xs:element name="sugar_embedding" minOccurs="0">
                """

                def set_sugar_embedding_type(embed, embed_in):
                    """
                    XSD: <xs:element name="sugar_embedding" minOccurs="0"> has
                    .. 2 elements
                    """

                    def set_el_material(embed, embed_in):
                        """
                        XSD: <xs:element name="material" type="xs:token">
                        CIF: _em_embedding.material 'tannin and glucose'
                        """
                        set_cif_value(embed.set_material, "material", const.EM_EMBEDDING, cif_list=embed_in, parent_el_req=False)

                    def set_el_details(embed, embed_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0" >
                        CIF: _em_embedding.details ?
                        """
                        set_cif_value(embed.set_details, "details", const.EM_EMBEDDING, cif_list=embed_in)

                    # element 1
                    set_el_material(embed, embed_in)
                    # element 2
                    set_el_details(embed, embed_in)

                if sp_id_in in embed_dict_in:
                    embed_in = embed_dict_in[sp_id_in]
                    el_material = get_cif_value("material", const.EM_EMBEDDING, cif_list=embed_in)
                    el_details = get_cif_value("details", const.EM_EMBEDDING, cif_list=embed_in)
                    if any(x is not None for x in [el_material, el_details]):
                        embed = emdb.sugar_embeddingType()
                        set_sugar_embedding_type(embed, embed_in)
                        if embed.has__content():
                            specimen.set_sugar_embedding(embed)

            def set_el_shadowing(specimen, sp_id_in, shadow_dict_in):
                """
                XSD: <xs:element name="shadowing" minOccurs="0">
                """

                def set_shadowing_type(shadow):
                    """
                    XSD: <xs:element name="shadowing" minOccurs="0"> is
                    .. a sequence of 4 elements
                    """

                    def set_el_material(shadow, shadow_in):
                        """
                        XSD: <xs:element name="material" type="xs:token">
                        CIF: _em_shadowing.material 'Platinum'
                        """
                        set_cif_value(shadow.set_material, "material", const.EM_SHADOWING, cif_list=shadow_in, parent_el_req=False)

                    def set_el_angle(shadow, shadow_in):
                        """
                        XSD: <xs:element name="angle"> has
                        .. a base allowed_angle_shadowing and
                        .. 1 attribute
                        CIF: _em_shadowing.angle 20
                        IS THIS CORRECT??????
                        """
                        set_cif_value(
                            shadow.set_angle, "angle", const.EM_SHADOWING, cif_list=shadow_in, constructor=emdb.angleType, fmt=float, units=const.U_DEG, parent_el_req=False
                        )

                    def set_el_thickness(shadow, shadow_in):
                        """
                        XSD: <xs:element name="thickness">
                        CIF: _em_shadowing.thickness ?
                        """
                        set_cif_value(
                            shadow.set_thickness,
                            "thickness",
                            const.EM_SHADOWING,
                            cif_list=shadow_in,
                            constructor=emdb.thicknessType,
                            fmt=float,
                            units=const.U_NM,
                            parent_el_req=False,
                        )

                    def set_el_details(shadow, shadow_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0" >
                        CIF: _em_shadowing.details 'rotary shadowing'
                        """
                        set_cif_value(shadow.set_details, "details", const.EM_SHADOWING, cif_list=shadow_in)

                    # element 1
                    set_el_material(shadow, shadow_in)  # pylint: disable=possibly-used-before-assignment
                    # element 2
                    set_el_angle(shadow, shadow_in)
                    # element 3
                    set_el_thickness(shadow, shadow_in)
                    # element 4
                    set_el_details(shadow, shadow_in)

                if sp_id_in in shadow_dict_in:
                    shadow_in = shadow_dict_in[sp_id_in]
                    el_material = get_cif_value("material", const.EM_SHADOWING, cif_list=shadow_in)
                    el_angle = get_cif_value("angle", const.EM_SHADOWING, cif_list=shadow_in)
                    el_thickness = get_cif_value("thickness", const.EM_SHADOWING, cif_list=shadow_in)
                    el_details = get_cif_value("details", const.EM_SHADOWING, cif_list=shadow_in)
                    if any(x is not None for x in [el_material, el_angle, el_thickness, el_details]):
                        shadow = emdb.shadowingType()
                        set_shadowing_type(shadow)
                        if shadow.has__content():
                            specimen.set_shadowing(shadow)

            def set_el_grid(specimen, sp_id_in, grid_dict_in, film_dict_in, pretreat_dict_in):
                """
                XSD: <xs:element name="grid" type="grid_type">
                """

                def set_grid_type(grid, film_dict_in, pretreat_dict_in):
                    """
                    XSD: <xs:element name="grid" type="grid_type"> is
                    .. a sequence of 6 elements
                    CIF: _em_sample_support.id 1
                    """

                    def set_el_model(grid, grid_in):
                        """
                        XSD: <xs:element name="model" type="xs:token" minOccurs="0">
                        CIF: _em_sample_support.grid_type 'Quantifoil R2/2'
                        """
                        set_cif_value(grid.set_model, "grid_type", const.EM_SAMPLE_SUPPORT, cif_list=grid_in)

                    def set_el_material(grid, grid_in):
                        """
                        XSD: <xs:element name="material" minOccurs="0">
                        CIF: _em_sample_support.grid_material COPPER
                        """
                        set_cif_value(grid.set_material, "grid_material", const.EM_SAMPLE_SUPPORT, cif_list=grid_in)

                    def set_el_mesh(grid, grid_in):
                        """
                        XSD: <xs:element name="mesh" type="xs:positiveInteger" minOccurs="0">
                        CIF: _em_sample_support.grid_mesh_size 400
                        """
                        set_cif_value(grid.set_mesh, "grid_mesh_size", const.EM_SAMPLE_SUPPORT, cif_list=grid_in, fmt=int)

                    def set_el_support_film(grid, grid_id_in, film_dict_in):
                        """
                        XSD: <xs:element name="support_film" type="film_type" maxOccurs="unbounded">
                        """

                        def set_film_type(film, film_in):
                            """
                            XSD: <xs:complexType name="film_type"> has
                            .. 1 attribute and
                            .. 3 elements
                            """

                            def set_attr_id(film, film_in):
                                """
                                XSD: <xs:attribute name="film_type_id" type="xs:positiveInteger" use="required"/>
                                CIF: _em_support_film.id 1
                                """
                                set_cif_value(film.set_film_type_id, "id", const.EM_SUPPORT_FILM, cif_list=film_in, fmt=int)

                            def set_el_film_material(film, film_in):
                                """
                                XSD: <xs:element name="film_material" type="xs:token" minOccurs="0">
                                CIF: _em_support_film.material CARBON
                                """
                                set_cif_value(film.set_film_material, "material", const.EM_SUPPORT_FILM, cif_list=film_in)

                            def set_el_film_topology(film, film_in):
                                """
                                XSD: <xs:element name="film_topology">
                                CIF: _em_support_film.topology CONTINUOUS
                                """
                                set_cif_value(film.set_film_topology, "topology", const.EM_SUPPORT_FILM, cif_list=film_in)

                            def set_el_film_thickness(film, film_in):
                                """
                                XSD: <xs:element name="film_thickness" minOccurs="0">
                                CIF: _em_support_film.thickness 50
                                mmCIF dict: Thickness of the support film, in Angstrom
                                """
                                # film_thickness = get_cif_value('thickness', const.EM_SUPPORT_FILM, cif_list=film_in)
                                # if film_thickness is not None:
                                #     fl_film_thickness = float(film_thickness) * 0.1
                                #     if fl_film_thickness < 5:
                                #         txt = u'The value for (_em_support_film.thickness) is (%s) angstroms. The lowest value should be 5.0 nm.' % film_thickness
                                #         self.current_entry_log.warn_logs.append(self.ALog(log_text=self.current_entry_log.not_changed_for_now_title + txt))
                                #         self.log_formatted(self.warn_log_string, const.NOT_CHANGED_FOR_NOW + txt)
                                #     elif fl_film_thickness > 50:
                                #         txt = u'The value for (_em_support_film.thickness) is (%s) angstroms. The highest value is 50.0 nm.'% film_thickness
                                #         self.current_entry_log.warn_logs.append(self.ALog(log_text=self.current_entry_log.not_changed_for_now_title + txt))
                                #         self.log_formatted(self.warn_log_string, const.NOT_CHANGED_FOR_NOW + txt)
                                set_cif_value(
                                    film.set_film_thickness,
                                    "thickness",
                                    const.EM_SUPPORT_FILM,
                                    cif_list=film_in,
                                    constructor=emdb.film_thicknessType,
                                    fmt=lambda x: float(x) * 0.1,
                                    units=const.U_NM,
                                )

                            # attribute 1
                            set_attr_id(film, film_in)
                            # element 1
                            set_el_film_material(film, film_in)
                            # element 2
                            set_el_film_topology(film, film_in)
                            # element 3
                            set_el_film_thickness(film, film_in)

                        if grid_id_in in film_dict_in:
                            for film_in in film_dict_in[grid_id_in]:
                                film = emdb.film_type()
                                set_film_type(film, film_in)
                                grid.add_support_film(film)

                    def set_el_pretreatment(grid, grid_id_in, pretreat_dict_in):
                        """
                        XSD: <xs:element name="pretreatment" type="grid_pretreatment_type" minOccurs="0"/>
                        """

                        def set_grid_pretreatment_type(pretreat):
                            """
                            XSD: <xs:complexType name="grid_pretreatment_type"> is
                            .. a sequence of 4 elements
                            """

                            def set_el_type(pretreat, pretreat_in):
                                """
                                XSD: <xs:element name="type" type="xs:token" minOccurs="0"/>
                                CIF: _em_grid_pretreatment.type 'GLOW DISCHARGE'
                                """
                                set_cif_value(pretreat.set_type, "type", const.EM_GRID_PRETREATMENT, cif_list=pretreat_in, parent_el_req=False)

                            def set_el_time(pretreat, pretreat_in):
                                """
                                XSD: <xs:element name="time" minOccurs="0">
                                CIF: _em_grid_pretreatment.time ? or 60
                                """
                                set_cif_value(
                                    pretreat.set_time, "time", const.EM_GRID_PRETREATMENT, cif_list=pretreat_in, constructor=emdb.timeType, fmt=int, units=const.U_SEC
                                )

                            def set_el_atmosphere(pretreat, pretreat_in):
                                """
                                XSD: <xs:element name="atmosphere" minOccurs="0">
                                CIF: _em_grid_pretreatment.atmosphere ? or AIR
                                """
                                currently_allawed_values = ["AIR", "AMYLAMINE", "NITROGEN", "OTHER"]
                                atm_value = get_cif_value("atmosphere", const.EM_GRID_PRETREATMENT, pretreat_in)
                                if atm_value is not None:
                                    atm_value_up = atm_value.upper()
                                    if atm_value_up not in currently_allawed_values:
                                        # the value is not in the currently allowed values - set it to OTHER
                                        atm_value_up = "OTHER"
                                        set_cif_value(
                                            pretreat.set_atmosphere, "atmosphere", const.EM_GRID_PRETREATMENT, cif_list=pretreat_in, fmt=str.upper, cif_value=atm_value_up
                                        )
                                        txt = u"The value (%s) for (_em_grid_pretreatment.atmosphere) is not allowed so it is changed to (%s)." % (atm_value, atm_value_up)
                                        self.current_entry_log.warn_logs.append(
                                            self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                        )
                                        self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                                    else:
                                        set_cif_value(pretreat.set_atmosphere, "atmosphere", const.EM_GRID_PRETREATMENT, cif_list=pretreat_in, fmt=str.upper)

                            def set_el_pressure(pretreat, pretreat_in):
                                """
                                XSD: <xs:element name="pressure" minOccurs="0">
                                CIF: _em_grid_pretreatment.pressure ? or 101325
                                """
                                set_cif_value(
                                    pretreat.set_pressure,
                                    "pressure",
                                    const.EM_GRID_PRETREATMENT,
                                    cif_list=pretreat_in,
                                    constructor=emdb.pressureType,
                                    fmt=lambda x: float(x) * 0.001,
                                    units=const.U_KPA,
                                )

                            # element 1
                            set_el_type(pretreat, pretreat_in)  # pylint: disable=possibly-used-before-assignment
                            # element 2
                            set_el_time(pretreat, pretreat_in)
                            # element 3
                            set_el_atmosphere(pretreat, pretreat_in)
                            # element 4
                            set_el_pressure(pretreat, pretreat_in)

                        if grid_id_in in pretreat_dict_in:
                            pretreat_in = pretreat_dict_in[grid_id_in]
                            el_type = get_cif_value("type", const.EM_GRID_PRETREATMENT, cif_list=pretreat_in)
                            el_time = get_cif_value("pretreat_time", const.EM_GRID_PRETREATMENT, cif_list=pretreat_in)
                            el_atmosphere = get_cif_value("atmosphere", const.EM_GRID_PRETREATMENT, pretreat_in)
                            el_pressure = get_cif_value("pressure", const.EM_GRID_PRETREATMENT, cif_list=pretreat_in)
                            if any(x is not None for x in [el_type, el_time, el_atmosphere, el_pressure]):
                                pretreat = emdb.grid_pretreatment_type()
                                set_grid_pretreatment_type(pretreat)
                                grid.set_pretreatment(pretreat)

                    def set_el_details(grid, grid_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                        CIF: _em_sample_support.details ?
                        """
                        set_cif_value(grid.set_details, "details", const.EM_SAMPLE_SUPPORT, cif_list=grid_in)

                    grid_id_in = get_cif_value(const.K_ID, const.EM_SAMPLE_SUPPORT, grid_in)  # pylint: disable=possibly-used-before-assignment
                    # element 1
                    set_el_model(grid, grid_in)
                    # element 2
                    set_el_material(grid, grid_in)
                    # element 3
                    set_el_mesh(grid, grid_in)
                    # element 4
                    set_el_support_film(grid, grid_id_in, film_dict_in)
                    # element 5
                    set_el_pretreatment(grid, grid_id_in, pretreat_dict_in)
                    # element 6
                    set_el_details(grid, grid_in)

                if sp_id_in in grid_dict_in:
                    grid_in = grid_dict_in[sp_id_in]
                    grid = emdb.grid_type()
                    set_grid_type(grid, film_dict_in, pretreat_dict_in)
                    if grid.has__content():
                        specimen.set_grid(grid)

            def set_el_vitrification(specimen, sp_id_in, vitr_dict_in):
                """
                XSD: <xs:element name="vitrification" type="vitrification_type" minOccurs="0">
                """

                def set_vitrification_type(vitr, vitr_in):
                    """
                    XSD: <xs:complexType name="vitrification_type"> is
                    .. a sequence of 7 elements
                    _em_vitrification or _emd_specimen_vitrification?
                    They are the same; XSD references the former,
                    here the latter is used.
                    """

                    def set_el_cryogen_name(vitr, vitr_in):
                        """
                        XSD: <xs:element name="cryogen_name">
                        CIF: _em_vitrification.cryogen_name ETHANE
                        """
                        set_cif_value(vitr.set_cryogen_name, "cryogen_name", const.EM_VITRIFICATION, cif_list=vitr_in)

                    def set_el_chamber_humidity(vitr, vitr_in):
                        """
                        XSD: <xs:element name="chamber_humidity" minOccurs="0">
                        CIF: _em_vitrification.humidity 100
                        """
                        set_cif_value(
                            vitr.set_chamber_humidity, "humidity", const.EM_VITRIFICATION, cif_list=vitr_in, constructor=emdb.chamber_humidityType, units=const.U_PERCENT
                        )

                    def set_el_chamber_temperature(vitr, vitr_in):
                        """
                        XSD: <xs:element name="chamber_temperature" minOccurs="0">
                        CIF: _em_vitrification.chamber_temperature 277
                        """
                        chamber_temp = get_cif_value("chamber_temperature", const.EM_VITRIFICATION, vitr_in)
                        if chamber_temp is not None:
                            fl_chamber_temp = float(chamber_temp)
                            if fl_chamber_temp < 85:
                                txt = u"The value given for (_em_vitrification.chamber_temperature) is (%s) K. The lowest value should be 85.0 K." % fl_chamber_temp
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.not_changed_for_now_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.NOT_CHANGED_FOR_NOW + txt)
                            elif fl_chamber_temp > 300:
                                txt = u"The value given for (_em_vitrification.chamber_temperature) is (%s) K. The highest value should be 300.0 K." % fl_chamber_temp
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.not_changed_for_now_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.NOT_CHANGED_FOR_NOW + txt)

                        set_cif_value(
                            vitr.set_chamber_temperature,
                            "chamber_temperature",
                            const.EM_VITRIFICATION,
                            cif_list=vitr_in,
                            constructor=emdb.chamber_temperatureType,
                            units=const.U_KEL,
                        )

                    def set_el_instrument(vitr, vitr_in):
                        """
                        XSD: <xs:element name="instrument" minOccurs="0">
                        CIF: _em_vitrification.instrument 'FEI VITROBOT MARK III'
                        """
                        details_txt = u""
                        instrument = get_cif_value("instrument", const.EM_VITRIFICATION, cif_list=vitr_in)
                        allowed_instruments = {
                            "EMS-002 RAPID IMMERSION FREEZER",
                            "FEI VITROBOT MARK I",
                            "FEI VITROBOT MARK II",
                            "FEI VITROBOT MARK III",
                            "FEI VITROBOT MARK IV",
                            "GATAN CRYOPLUNGE 3",
                            "HOMEMADE PLUNGER",
                            "LEICA EM CPC",
                            "LEICA EM GP",
                            "LEICA KF80",
                            "LEICA PLUNGER",
                            "REICHERT-JUNG PLUNGER",
                            "SPOTITON",
                            "OTHER",
                        }
                        if instrument is not None:
                            if instrument in allowed_instruments:
                                set_cif_value(vitr.set_instrument, "instrument", const.EM_VITRIFICATION, cif_list=vitr_in)
                            else:
                                # write OTHER and add note in details
                                set_cif_value(vitr.set_instrument, "instrument", const.EM_VITRIFICATION, cif_list=vitr_in, cif_value="OTHER")
                                details_txt = (
                                    u"The value given for _em_vitrification.instrument is %s. This is not in a list of allowed values %s so OTHER is written into the XML file."
                                    % (instrument, allowed_instruments)
                                )
                        return details_txt

                    def set_el_details(vitr, vitr_in, details_txt):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                        CIF: _em_vitrification.details
                        """
                        all_details = ""
                        current_details = get_cif_value("details", const.EM_VITRIFICATION, cif_list=vitr_in)
                        if current_details is not None:
                            all_details = ". ".join((current_details, details_txt))
                        else:
                            all_details = details_txt
                        if all_details != "":
                            set_cif_value(vitr.set_details, "details", const.EM_VITRIFICATION, cif_list=vitr_in, cif_value=all_details)

                    def set_el_timed_resolved_state():
                        """
                        XSD: <xs:element name="timed_resolved_state" type="xs:token" minOccurs="0">
                        CIF: Deprecated (2014-10-21)
                        """

                    def set_el_method():
                        """
                        XSD: <xs:element name="method" type="xs:string" minOccurs="0">
                        CIF:Deprecated (2014/11/13)
                        """

                    # element 1
                    set_el_cryogen_name(vitr, vitr_in)
                    # element 2
                    set_el_chamber_humidity(vitr, vitr_in)
                    # element 3
                    set_el_chamber_temperature(vitr, vitr_in)
                    # element 4
                    details_txt = set_el_instrument(vitr, vitr_in)
                    # element 5
                    set_el_details(vitr, vitr_in, details_txt)
                    # element 6
                    set_el_timed_resolved_state()
                    # element 7
                    set_el_method()

                if sp_id_in in vitr_dict_in:
                    vitr_in = vitr_dict_in[sp_id_in]
                    vitr = emdb.vitrification_type()
                    set_vitrification_type(vitr, vitr_in)
                    if vitr.has__content():
                        specimen.set_vitrification(vitr)

            def set_el_details(specimen, spec_prep_in):
                """
                XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                CIF: _em_specimen.details
                """
                set_cif_value(specimen.set_details, "details", const.EM_SPECIMEN, cif_list=spec_prep_in)

            vitr_dict_in = make_dict(const.EM_VITRIFICATION, const.K_SPECIMEN_ID)
            stain_dict_in = make_dict(const.EM_STAINING, const.K_SPECIMEN_ID)
            embed_dict_in = make_dict(const.EM_EMBEDDING, const.K_SPECIMEN_ID)
            shadow_dict_in = make_dict(const.EM_SHADOWING, const.K_SPECIMEN_ID)
            grid_dict_in = make_dict(const.EM_SAMPLE_SUPPORT, const.K_SPECIMEN_ID)
            buff_dict_in = make_dict(const.EM_BUFFER, const.K_SPECIMEN_ID)
            film_dict_in = make_list_of_dicts(const.EM_SUPPORT_FILM, const.K_SAMPLE_SUPPORT_ID)
            pretreat_dict_in = make_dict(const.EM_GRID_PRETREATMENT, const.K_SAMPLE_SUPPORT_ID)
            buff_comp_dict_in = make_list_of_dicts(const.EM_BUFFER_COMPONENT, const.K_BUFFER_ID)

            sp_id_in = get_cif_value(const.K_ID, const.EM_SPECIMEN)

            # attribute 1
            set_attr_id(specimen, spec_prep_in)
            # element 1
            set_el_concentration(specimen, spec_prep_in)
            # element 2
            set_el_buffer(specimen, sp_id_in, buff_dict_in, buff_comp_dict_in)
            # element 3
            set_el_staining(specimen, sp_id_in, stain_dict_in)
            # element 4
            set_el_sugar_embedding(specimen, sp_id_in, embed_dict_in)
            # element 5
            set_el_shadowing(specimen, sp_id_in, shadow_dict_in)
            # element 6
            set_el_grid(specimen, sp_id_in, grid_dict_in, film_dict_in, pretreat_dict_in)
            # element 7
            set_el_vitrification(specimen, sp_id_in, vitr_dict_in)
            # element 8
            set_el_details(specimen, spec_prep_in)

        def set_tom_prep_specifics(sp_id_in, tom_prep):
            """
            Set additional elements for tomography specimen preparation

            Parameters:
            @param sp_id_in: specimen id
            @param tom_prep: object for tomography_preparation_type
            XSD: <xs:complexType name="tomography_preparation_type"> has
            .. a base (see: set_base_specimen_preparation) and
            .. a sequence of 5 elements
            CIF:  _em_tomography_specimen.id
            """

            def set_el_fiducial_markers_list(tom_prep, tom_id_in, fid_dict_in):
                """
                XSD: <xs:element name="fiducial_markers_list" minOccurs="0">
                CIF: _em_tomography_specimen.fiducial_markers YES
                """

                def set_fiducial_marker_type(fid, fid_in):
                    """
                    XSD: <xs:complexType name="fiducial_marker_type"> has
                    .. 3 elements
                    """

                    def set_el_fiducial_type():
                        """
                        XSD: <xs:element name="fiducial_type" type="xs:token" minOccurs="0"/>
                        CIF:?????
                        """

                    def set_el_manufacturer(fid, fid_in):
                        """
                        XSD: <xs:element name="manufacturer" type="xs:token" minOccurs="0">
                        CIF: _em_fiducial_markers.manufacturer 'nanoprobes'
                        """
                        set_cif_value(fid.set_manufacturer, "manufacturer", const.EM_FIDUCIAL_MARKERS, cif_list=fid_in)

                    def set_el_diameter(fid, fid_in):
                        """
                        XSD: <xs:element name="diameter" type="fiducial_marker_diameter_type"/>
                        CIF:  _em_fiducial_markers.diameter 14
                        """
                        set_cif_value(fid.set_diameter, "diameter", const.EM_FIDUCIAL_MARKERS, cif_list=fid_in, constructor=emdb.fiducial_marker_diameter_type, units=const.U_NMF)

                    # element 1
                    set_el_fiducial_type()
                    # element 2
                    set_el_manufacturer(fid, fid_in)
                    # element 3
                    set_el_diameter(fid, fid_in)

                fs_present = get_cif_value("fiducial_markers", const.EM_TOMOGRAPHY_SPECIMEN, tom_prep_in)  # pylint: disable=possibly-used-before-assignment
                if tom_id_in in fid_dict_in and fs_present == "YES":
                    fid_list = emdb.fiducial_markers_listType()
                    fid_list_in = fid_dict_in[tom_id_in]
                    for fid_in in fid_list_in:
                        fid = emdb.fiducial_marker_type()
                        set_fiducial_marker_type(fid, fid_in)
                        if fid.has__content():
                            fid_list.add_fiducial_marker(fid)
                    if fid_list.has__content():
                        tom_prep.set_fiducial_markers_list(fid_list)

            def set_el_high_pressure_freezing(tom_prep, tom_id_in, h_pfdict_in):
                """
                XSD: <xs:element name="high_pressure_freezing" minOccurs="0">
                CIF: _em_tomography_specimen.high_pressure_freezing YES
                """

                def set_high_pressure_freezing_type(h_pf, h_pf_in):
                    """
                    XSD: <xs:element name="high_pressure_freezing" minOccurs="0"> has
                    .. 2 elements
                    """

                    def set_el_instrument(h_pf, h_pf_in):
                        """
                        XSD: <xs:element name="instrument">
                        CIF: _em_high_pressure_freezing.instrument 'Leica EM HP100'
                        """
                        details_txt = u""
                        instrument = get_cif_value("instrument", const.EM_HIGH_PRESSURE_FREEZING, cif_list=h_pf_in)
                        allowed_instruments = {"BAL-TEC HPM 010", "EMS-002 RAPID IMMERSION FREEZER", "LEICA EM HPM100", "LEICA EM PACT", "LEICA EM PACT2", "OTHER"}
                        if instrument is not None:
                            if instrument in allowed_instruments:
                                set_cif_value(h_pf.set_instrument, "instrument", const.EM_HIGH_PRESSURE_FREEZING, cif_list=h_pf_in, parent_el_req=False)
                            else:
                                # write OTHER and add note in details
                                set_cif_value(h_pf.set_instrument, "instrument", const.EM_HIGH_PRESSURE_FREEZING, cif_list=h_pf_in, cif_value="OTHER", parent_el_req=False)
                                details_txt = (
                                    u"The value given for _em_high_pressure_freezing.instrument is %s. This is not in a list of allowed values %s so OTHER is written into the XML file."
                                    % (instrument, allowed_instruments)
                                )
                        return details_txt

                    def set_el_details(h_pf, h_pf_in, details_txt):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _em_high_pressure_freezing.details
                        """
                        all_details = ""
                        current_details = get_cif_value("details", const.EM_HIGH_PRESSURE_FREEZING, cif_list=h_pf_in)
                        if current_details is not None:
                            all_details = ". ".join((current_details, details_txt))
                        else:
                            all_details = details_txt
                        set_cif_value(h_pf.set_details, "details", const.EM_HIGH_PRESSURE_FREEZING, cif_list=h_pf_in, cif_value=all_details, parent_el_req=False)

                    # element 1
                    details_txt = set_el_instrument(h_pf, h_pf_in)
                    # element 2
                    set_el_details(h_pf, h_pf_in, details_txt)

                hpf_present = get_cif_value("high_pressure_freezing", const.EM_TOMOGRAPHY_SPECIMEN, tom_prep_in)
                if tom_id_in in h_pfdict_in and hpf_present == "YES":
                    h_pf_in = h_pfdict_in[tom_id_in]
                    el_instrument = get_cif_value("instrument", const.EM_HIGH_PRESSURE_FREEZING, cif_list=h_pf_in)
                    el_details = get_cif_value("details", const.EM_HIGH_PRESSURE_FREEZING, cif_list=h_pf_in)
                    if any(x is not None for x in [el_instrument, el_details]):
                        h_pf = emdb.high_pressure_freezingType()
                        set_high_pressure_freezing_type(h_pf, h_pf_in)
                        if h_pf.has__content():
                            tom_prep.set_high_pressure_freezing(h_pf)

            def set_el_embedding_material():
                """
                XSD: <xs:element name="embedding_material" minOccurs="0" type="xs:token">
                CIF: ??
                """

            def set_el_cryo_protectant(tom_prep, tom_prep_in):
                """
                XSD: <xs:element name="cryo_protectant" minOccurs="0" type="xs:token">
                CIF: _em_tomography_specimen.cryo_protectant '2% glycerol'
                """
                set_cif_value(tom_prep.set_cryo_protectant, "cryo_protectant", const.EM_TOMOGRAPHY_SPECIMEN, cif_list=tom_prep_in)

            def set_el_sectioning(tom_prep, tom_prep_in, u_tome_dict_in):
                """
                XSD: <xs:element name="sectioning" minOccurs="0"> is a choice of
                .. 3 elements
                CIF: _em_tomography_specimen.sectioning
                {ULTRAMICROTOMY, FOCUSED ION BEAM, NO SECTIONING }
                """

                def set_ultramicrotomy_type(u_tome):
                    """
                    XSD: <xs:element name="ultramicrotomy"> has
                    .. 4 elements
                    """

                    def set_el_instrument(u_tome, u_tome_in):
                        """
                        XSD: <xs:element name="instrument" type="xs:token">
                        CIF: _em_ultramicrotomy.instrument 'Leica EM UC7'
                        """
                        set_cif_value(u_tome.set_instrument, "instrument", const.EM_ULTRAMICROTOMY, cif_list=u_tome_in)

                    def set_el_temperature(u_tome, u_tome_in):
                        """
                        XSD: <xs:element name="temperature" type="temperature_type"/>
                        CIF: _em_ultramicrotomy.temperature 100
                        """
                        set_cif_value(
                            u_tome.set_temperature, "temperature", const.EM_ULTRAMICROTOMY, cif_list=u_tome_in, constructor=emdb.temperature_type, units=const.U_KEL
                        )

                    def set_el_final_thickness(u_tome, u_tome_in):
                        """
                        XSD: <xs:element name="final_thickness" type="ultramicrotomy_final_thickness_type"/>
                        CIF: _em_ultramicrotomy.final_thickness 60
                        """
                        set_cif_value(
                            u_tome.set_final_thickness,
                            "final_thickness",
                            const.EM_ULTRAMICROTOMY,
                            cif_list=u_tome_in,
                            constructor=emdb.ultramicrotomy_final_thickness_type,
                            units=const.U_NM,
                        )

                    def set_el_details(u_tome, u_tome_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _em_ultramicrotomy.details ?
                        """
                        set_cif_value(u_tome.set_details, "details", const.EM_ULTRAMICROTOMY, cif_list=u_tome_in)

                    # element 1
                    set_el_instrument(u_tome, u_tome_in)  # Will be inherited from caller scope.  #  pylint: disable=used-before-assignment
                    # element 2
                    set_el_temperature(u_tome, u_tome_in)
                    # element 3
                    set_el_final_thickness(u_tome, u_tome_in)
                    # element 4
                    set_el_details(u_tome, u_tome_in)

                def set_focused_ion_beam_type(fib):
                    """
                    XSD: <xs:element name="focused_ion_beam"> has
                    .. 10 elements
                    """

                    def set_el_instrument(fib, fib_in):
                        """
                        XSD: <xs:element name="instrument">
                        CIF: _em_focused_ion_beam.instrument 'FEI Quanta FIB'
                        """
                        details_txt = u""
                        instrument = get_cif_value("instrument", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in)
                        allowed_instruments = {"DB235", "OTHER"}
                        if instrument is not None:
                            if instrument in allowed_instruments:
                                set_cif_value(fib.set_instrument, "instrument", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in)
                            else:
                                # write OTHER and add note in details
                                set_cif_value(fib.set_instrument, "instrument", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in, cif_value="OTHER")
                                details_txt = (
                                    u"The value given for _em_focused_ion_beam.instrument is %s. This is not in a list of allowed values %s so OTHER is written into the XML file."
                                    % (instrument, allowed_instruments)
                                )
                        return details_txt

                    def set_el_ion(fib, fib_in, details_txt):
                        """
                        XSD: <xs:element name="ion">
                        CIF: _em_focused_ion_beam.ion 'gallium ion'
                        """
                        ion = get_cif_value("ion", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in)
                        allowed_ions = {"GALLIUM+", "OTHER"}
                        if ion is not None:
                            if ion in allowed_ions:
                                set_cif_value(fib.set_ion, "ion", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in)
                            else:
                                set_cif_value(fib.set_ion, "ion", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in, cif_value="OTHER")
                                add_text = (
                                    " The value given for _em_focused_ion_beam.ions is %s. This is not in a list of allowed values %s so OTHER is written into the XML file."
                                    % (ion, allowed_ions)
                                )
                                details_txt = u"".join(add_text)
                        return details_txt

                    def set_el_voltage(fib, fib_in):
                        """
                        XSD: <xs:element name="voltage" type="fib_voltage_type"/>
                        CIF: _em_focused_ion_beam.voltage  30
                        """
                        set_cif_value(fib.set_voltage, "voltage", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in, constructor=emdb.fib_voltage_type, units=const.U_KVOLT)

                    def set_el_current(fib, fib_in):
                        """
                        XSD: <xs:element name="current" type="fib_current_type"/>
                        CIF: _em_focused_ion_beam.current
                        """
                        set_cif_value(fib.set_current, "current", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in, constructor=emdb.fib_current_type, units=const.U_NAMP)

                    def set_el_dose_rate(fib, fib_in):
                        """
                        XSD: <xs:element name="dose_rate" type="fib_dose_rate_type" minOccurs="0/>
                        CIF: _em_focused_ion_beam.dose_rate
                        """
                        set_cif_value(
                            fib.set_dose_rate, "dose_rate", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in, constructor=emdb.fib_dose_rate_type, units=const.U_FIB_DOSE_RATE
                        )

                    def set_el_duration(fib, fib_in):
                        """
                        XSD: <xs:element name="duration" type="fib_duration_type"/>
                        CIF: _em_focused_ion_beam.duration
                        """
                        set_cif_value(fib.set_duration, "duration", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in, constructor=emdb.fib_duration_type, units=const.U_SEC)

                    def set_el_temperature(fib, fib_in):
                        """
                        XSD: <xs:element name="temperature" type="temperature_type"/>
                        CIF: _em_focused_ion_beam.temperature 100
                        """
                        set_cif_value(
                            fib.set_temperature, "temperature", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in, constructor=emdb.temperature_type, units=const.U_KEL
                        )

                    def set_el_initial_thickness(fib, fib_in):
                        """
                        XSD: <xs:element name="initial_thickness" type="fib_initial_thickness_type">
                        CIF: _em_focused_ion_beam.initial_thickness
                        """
                        set_cif_value(
                            fib.set_initial_thickness,
                            "initial_thickness",
                            const.EM_FOCUSED_ION_BEAM,
                            cif_list=fib_in,
                            constructor=emdb.fib_initial_thickness_type,
                            units=const.U_NM,
                        )

                    def set_el_final_thickness(fib, fib_in):
                        """
                        XSD: <xs:element name="final_thickness" type="fib_final_thickness_type"/>
                        CIF: _em_focused_ion_beam.final_thickness
                        """
                        set_cif_value(
                            fib.set_final_thickness,
                            "final_thickness",
                            const.EM_FOCUSED_ION_BEAM,
                            cif_list=fib_in,
                            constructor=emdb.fib_final_thickness_type,
                            units=const.U_NM,
                        )

                    def set_el_details(fib, fib_in, details_txt):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _em_focused_ion_beam.details
                        """
                        all_details = ""
                        current_details = get_cif_value("details", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in)
                        if current_details is not None:
                            all_details = ". ".join((current_details, details_txt))
                        else:
                            all_details = details_txt
                        set_cif_value(fib.set_details, "details", const.EM_FOCUSED_ION_BEAM, cif_list=fib_in, cif_value=all_details)

                    # element 1
                    details_txt = set_el_instrument(fib, fib_in)  # Will be inherited from caller scope.  #  pylint: disable=used-before-assignment
                    # element 2
                    set_el_ion(fib, fib_in, details_txt)
                    # element 3
                    set_el_voltage(fib, fib_in)
                    # element 4
                    set_el_current(fib, fib_in)
                    # element 5
                    set_el_dose_rate(fib, fib_in)
                    # element 6
                    set_el_duration(fib, fib_in)
                    # element 7
                    set_el_temperature(fib, fib_in)
                    # element 8
                    set_el_initial_thickness(fib, fib_in)
                    # element 9
                    set_el_final_thickness(fib, fib_in)
                    # element 10
                    set_el_details(fib, fib_in, details_txt)

                sec_in = get_cif_value("sectioning", const.EM_TOMOGRAPHY_SPECIMEN, tom_prep_in)
                if sec_in is not None:
                    if sec_in == "ULTRAMICROTOMY" and tom_id_in in u_tome_dict_in:  # Ordering sets this  #  pylint: disable=possibly-used-before-assignment
                        u_tome_in = u_tome_dict_in[tom_id_in]
                        u_tome = emdb.ultramicrotomyType()
                        set_ultramicrotomy_type(u_tome)
                        if u_tome.has__content():
                            tom_prep.set_sectioning(emdb.sectioningType(ultramicrotomy=u_tome))
                    elif sec_in == "FOCUSED ION BEAM" and tom_id_in in fib_dict_in:
                        fib_in = fib_dict_in[tom_id_in]  # fib_in must be set before calling set_focused_ion_beam_type
                        fib = emdb.focused_ion_beamType()
                        set_focused_ion_beam_type(fib)
                        if fib.has__content():
                            tom_prep.set_sectioning(emdb.sectioningType(focused_ion_beam=fib))
                    elif sec_in == "NO SECTIONING":
                        # XSD: choice 3: <xs:element name="other_sectioning" type="xs:string"/>
                        tom_prep.set_sectioning(emdb.sectioningType(other_sectioning=sec_in))

            # Create dictionaries where the tomography preparation id is key
            tom_prep_dict_in = make_dict(const.EM_TOMOGRAPHY_SPECIMEN, const.K_SPECIMEN_ID)
            fid_dict_in = make_list_of_dicts(const.EM_FIDUCIAL_MARKERS, const.K_EM_TOMOGRAPHY_SPECIMEN_ID)
            h_pfdict_in = make_dict(const.EM_HIGH_PRESSURE_FREEZING, const.K_EM_TOMOGRAPHY_SPECIMEN_ID)
            u_tome_dict_in = make_dict(const.EM_ULTRAMICROTOMY, const.K_EM_TOMOGRAPHY_SPECIMEN_ID)
            fib_dict_in = make_dict(const.EM_FOCUSED_ION_BEAM, const.K_EM_TOMOGRAPHY_SPECIMEN_ID)
            if sp_id_in in tom_prep_dict_in:
                tom_prep_in = tom_prep_dict_in[sp_id_in]
                tom_id_in = get_cif_value(const.K_ID, const.EM_TOMOGRAPHY_SPECIMEN, tom_prep_in)
                # element 1
                set_el_fiducial_markers_list(tom_prep, tom_id_in, fid_dict_in)
                # element 2
                set_el_high_pressure_freezing(tom_prep, tom_id_in, h_pfdict_in)
                # element 3
                set_el_embedding_material()
                # element 4
                set_el_cryo_protectant(tom_prep, tom_prep_in)
                # element 5
                set_el_sectioning(tom_prep, tom_prep_in, u_tome_dict_in)

        def set_el_cryst_prep_specifics(sp_id_in, cryst_prep):
            """
            Set additional elements for electron crystalography specimen preparation

            Parameters:
            @param sp_id_in:
            """

            def set_crystal_formation_type(cryst):
                """
                XSD: <xs:element name="crystal_formation"> has
                .. 7 elements
                """

                def set_el_lipid_protein_ratio(cryst, cryst_in):
                    """
                    XSD: <xs:element name="lipid_protein_ratio" type="xs:float" minOccurs="0"/>
                    CIF: _em_crystal_formation.lipid_protein_ratio 5.0
                    """
                    set_cif_value(cryst.set_lipid_protein_ratio, "lipid_protein_ratio", const.EM_CRYSTAL_FORMATION, cif_list=cryst_in, fmt=float)

                def set_el_lipid_mixture(cryst, cryst_in):
                    """
                    XSD: <xs:element name="lipid_mixture" type="xs:token" minOccurs="0"/>
                    CIF: _em_crystal_formation.lipid_mixture
                    """
                    set_cif_value(cryst.set_lipid_mixture, "lipid_mixture", const.EM_CRYSTAL_FORMATION, cif_list=cryst_in)

                def set_el_instrument(cryst, cryst_in):
                    """
                    XSD: <xs:element name="instrument" minOccurs="0">
                    CIF: _em_crystal_formation.instrument 'Langmuir trough'
                    """
                    set_cif_value(cryst.set_instrument, "instrument", const.EM_CRYSTAL_FORMATION, cif_list=cryst_in)

                def set_el_atmosphere(cryst, cryst_in):
                    """
                    XSD: <xs:element name="atmosphere" type="xs:token" minOccurs="0"/>
                    CIF: _em_crystal_formation.atmosphere
                    """
                    set_cif_value(cryst.set_atmosphere, "atmosphere", const.EM_CRYSTAL_FORMATION, cif_list=cryst_in)

                def set_el_temperature(cryst, cryst_in):
                    """
                    XSD: <xs:element name="temperature" type="crystal_formation_temperature_type" minOccurs="0"/>
                    CIF: _em_crystal_formation.temperature 298
                    """
                    set_cif_value(
                        cryst.set_temperature,
                        "temperature",
                        const.EM_CRYSTAL_FORMATION,
                        cif_list=cryst_in,
                        constructor=emdb.crystal_formation_temperature_type,
                        fmt=float,
                        units=const.U_KEL,
                    )

                def set_el_time(cryst, cryst_in):
                    """
                    XSD: <xs:element name="time" type="crystal_formation_time_type" minOccurs="0"/>
                    """

                    def set_crystal_formation_time_type(cryst_tm):
                        """
                        XSD: <xs:complexType name="crystal_formation_time_type">
                        CIF: _em_crystal_formation.time  50
                        CIF: _em_crystal_formation.time_unit DAY
                        """
                        tm_in = get_cif_value("time", const.EM_CRYSTAL_FORMATION, cryst_in)
                        tm_units_in = get_cif_value("time_unit", const.EM_CRYSTAL_FORMATION, cryst_in)
                        if tm_in is not None and tm_units_in is not None:
                            # write the value for time only if both values are given
                            set_cif_value(cryst_tm.set_valueOf_, "time", const.EM_CRYSTAL_FORMATION, cif_list=cryst_in, fmt=float)
                            set_cif_value(cryst_tm.set_units, "time_unit", const.EM_CRYSTAL_FORMATION, cif_list=cryst_in)
                            cryst.set_time(cryst_tm)

                    cryst_tm = emdb.crystal_formation_time_type()
                    set_crystal_formation_time_type(cryst_tm)

                def set_el_details(cryst, cryst_in):
                    """
                    XSD:<xs:element name="details" type="xs:string" minOccurs="0">
                    CIF: _em_crystal_formation.details
                    """
                    set_cif_value(cryst.set_details, "details", const.EM_CRYSTAL_FORMATION, cif_list=cryst_in)

                # element 1
                set_el_lipid_protein_ratio(cryst, cryst_in)  # Will never happen but should be passed to function.  #  pylint: disable=possibly-used-before-assignment
                # element 2
                set_el_lipid_mixture(cryst, cryst_in)
                # element 3
                set_el_instrument(cryst, cryst_in)
                # element 4
                set_el_atmosphere(cryst, cryst_in)
                # element 5
                set_el_temperature(cryst, cryst_in)
                # element 6
                set_el_time(cryst, cryst_in)
                # element 7
                set_el_details(cryst, cryst_in)

            cryst_dict_in = make_dict(const.EM_CRYSTAL_FORMATION, const.K_SPECIMEN_ID)

            if sp_id_in in cryst_dict_in:
                cryst_in = cryst_dict_in[sp_id_in]  # must be set before set_crystal_formation_type
                cryst = emdb.crystal_formationType()
                set_crystal_formation_type(cryst)
                if cryst.has__content():
                    cryst_prep.set_crystal_formation(cryst)

        def set_classification(ip_id_in, final_class_dict_in, cat_soft_dict_in):
            """
            Classification for image processing

            @param ip_id_in: image processing id
            @param final_class_dict_in: dictionary for final classification categories
            @return clas: classification type object
            XSD: <xs:complexType name="classification_type">
            """

            def set_classification_type(clas, f_c_in, cat_soft_dict_in):
                """
                XSD: <xs:complexType name="classification_type"> has
                .. 4 elements
                """

                def set_el_number_classes(clas, f_c_in):
                    """
                    XSD: <xs:element name="number_classes" type="xs:positiveInteger" minOccurs="0">
                    CIF: _em_final_classification.num_classes 200
                    """
                    set_cif_value(clas.set_number_classes, "num_classes", const.EM_FINAL_CLASSIFICATION, cif_list=f_c_in, fmt=int)

                def set_el_av_num_members_per_class(clas, f_c_in):
                    """
                    XSD: <xs:element name="average_number_members_per_class" minOccurs="0">
                    CIF: _em_final_classification.avg_num_images_per_class 75
                    """
                    set_cif_value(clas.set_average_number_members_per_class, "avg_num_images_per_class", const.EM_FINAL_CLASSIFICATION, cif_list=f_c_in, fmt=float)

                def set_el_software_list(clas, cat_soft_dict_in):
                    """
                    XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                    """
                    set_software_list(const.SOFT_CLASSIFICATION, cat_soft_dict_in, clas.set_software_list)

                def set_el_details(clas, f_c_in):
                    """
                    XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                    CIF: _em_final_classification.details
                    """
                    set_cif_value(clas.set_details, "details", const.EM_FINAL_CLASSIFICATION, cif_list=f_c_in)

                # element 1
                set_el_number_classes(clas, f_c_in)
                # element 2
                set_el_av_num_members_per_class(clas, f_c_in)
                # element 3
                set_el_software_list(clas, cat_soft_dict_in)
                # element 4
                set_el_details(clas, f_c_in)

            f_c_in = final_class_dict_in[ip_id_in]
            clas = emdb.classification_type()
            set_classification_type(clas, f_c_in, cat_soft_dict_in)
            return clas

        def set_struct_determination_type(struct_det):
            """
            XSD: <xs:complexType name="structure_determination_type"> has
            .. 1 attribute and
            .. a sequence of 6 elements
            """

            def set_attr_id(struct_det):
                """
                XSD: <xs:attribute name="structure_determination_id" type="xs:positiveInteger" use="required"/>
                CIF: _em_experiment.id 1
                """
                set_cif_value(struct_det.set_structure_determination_id, "id", const.EM_EXPERIMENT, fmt=int)

            def set_el_method(struct_det):
                """
                XSD: <xs:element name="method">
                CIF: _exptl.method 'X-RAY DIFFRACTION'
                CIF: _em_experiment.reconstruction_method 'SINGLE PARTICLE'
                # Structure determination - assumes only one!
                XSD-DOC: We are assuming one method per map - is that OK?
                """
                em_method_dict = {
                    "SINGLE PARTICLE": const.EMM_SP,
                    "SUBTOMOGRAM AVERAGING": const.EMM_STOM,
                    "TOMOGRAPHY": const.EMM_TOM,
                    "HELICAL": const.EMM_HEL,
                    "CRYSTALLOGRAPHY": const.EMM_EC,
                }
                # em_method is so fundamental - if it fails, let it. WHY?????
                em_method = None
                method = get_cif_value("method", const.EXPTL)
                if method == "ELECTRON CRYSTALLOGRAPHY":
                    em_method = const.EMM_EC
                    set_cif_value(struct_det.set_method, "method", const.EXPTL, cif_value=em_method)
                else:
                    metd = get_cif_value("reconstruction_method", const.EM_EXPERIMENT)
                    metd_item = get_cif_item("reconstruction_method", const.EM_EXPERIMENT)
                    if metd is not None:
                        em_method = em_method_dict.get(metd)
                    else:
                        txt = u"The value for (%s) is not given." % metd_item
                        self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                    if em_method is not None:
                        set_cif_value(struct_det.set_method, "reconstruction_method", const.EM_EXPERIMENT, cif_value=em_method)
                    else:
                        txt = (
                            u"(%s) is not a recognised structure determination method. Recognised methods are: SINGLE PARTICLE, SUBTOMOGRAM AVERAGING, TOMOGRAPHY, HELICAL, CRYSTALLOGRAPHY."
                            % metd_item
                        )
                        self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                return em_method

            def set_el_aggregation_state(struct_det):
                """
                XSD: <xs:element name="aggregation_state">
                CIF: _em_experiment.aggregation_state PARTICLE
                """
                set_cif_value(struct_det.set_aggregation_state, "aggregation_state", const.EM_EXPERIMENT, fmt=const.AGG_STATE_CIF2XML)

            def set_el_mol_and_complexes():
                """
                XSD: <xs:element name="macromolecules_and_complexes" type="macromolecules_and_complexes_type" minOccurs="0">
                CIF: _em_specimen.macromolecules_and_complexes DOESN'T EXIST YET
                """

            def set_el_specimen_prep_list(struct_det, em_method):
                """
                XSD: <xs:element name="specimen_preparation_list"> has
                .. is a sequence of 1 element - "specimen_preparation"
                """

                def set_spec_prep_list_type(spec_prep_list, em_method):
                    """
                    XSD: <xs:element name="specimen_preparation" type="base_preparation_type"> has
                    .. 1 attribute and
                    .. 8 elements
                    CIF: _em_specimen.id  1
                    """
                    spec_prep_list_in = assert_get_value(const.EM_SPECIMEN, self.cif)

                    for spec_prep_in in spec_prep_list_in:
                        sp_id_in = get_cif_value(const.K_ID, const.EM_SPECIMEN, spec_prep_in)
                        if em_method == const.EMM_SP:
                            single_part_prep = emdb.single_particle_preparation_type()
                            single_part_prep.original_tagname_ = "single_particle_preparation"
                            set_base_specimen_preparation(sp_id_in, spec_prep_in, single_part_prep)
                            spec_prep_list.add_specimen_preparation(single_part_prep)
                        elif em_method == const.EMM_STOM:
                            subtom_prep = emdb.subtomogram_averaging_preparation_type()
                            subtom_prep.original_tagname_ = "subtomogram_averaging_preparation"
                            set_base_specimen_preparation(sp_id_in, spec_prep_in, subtom_prep)
                            spec_prep_list.add_specimen_preparation(subtom_prep)
                        elif em_method == const.EMM_HEL:
                            hel_prep = emdb.helical_preparation_type()
                            hel_prep.original_tagname_ = "helical_preparation"
                            set_base_specimen_preparation(sp_id_in, spec_prep_in, hel_prep)
                            spec_prep_list.add_specimen_preparation(hel_prep)
                        elif em_method == const.EMM_TOM:
                            tom_prep = emdb.tomography_preparation_type()
                            tom_prep.original_tagname_ = "tomography_preparation"
                            set_base_specimen_preparation(sp_id_in, spec_prep_in, tom_prep)
                            set_tom_prep_specifics(sp_id_in, tom_prep)
                            spec_prep_list.add_specimen_preparation(tom_prep)
                        elif em_method == const.EMM_EC:
                            cryst_prep = emdb.crystallography_preparation_type()
                            cryst_prep.original_tagname_ = "crystallography_preparation"
                            set_base_specimen_preparation(sp_id_in, spec_prep_in, cryst_prep)
                            set_el_cryst_prep_specifics(sp_id_in, cryst_prep)
                            spec_prep_list.add_specimen_preparation(cryst_prep)
                        else:
                            txt = u"Unknown EM method: (%s). The specimen preparation section, therefore, cannot be set." % em_method
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                            break

                spec_prep_list = emdb.specimen_preparation_listType()
                set_spec_prep_list_type(spec_prep_list, em_method)
                struct_det.set_specimen_preparation_list(spec_prep_list)

            def set_el_microscopy_list(struct_det, microscopy_list):
                """
                XSD: <xs:element name="microscopy_list">
                """

                def get_tilt_axis(ts_in, axis1=True):
                    """
                    Get axis min, max and inc from a EM_TOMOGRAPHY
                    element and return an axis element

                    Parameters:
                    @param ts_in: cif dictionary item for EM_TOMOGRAPHY
                                  containing axis parameters
                    @param axis1: boolean whether this is axis1 or axis2
                    @return: axis XML element that can be added to tilt series
                    XSD: <xs:complexType name="axis_type">
                    """

                    def set_axis_type(axis, min_angle, max_angle, angle_inc, ts_in):
                        """
                        <xs:complexType name="axis_type"> has
                        .. 3 elements
                        """

                        def set_el_min_angle(axis, min_angle, ts_in):
                            """
                            XSD: <xs:element name="min_angle" minOccurs="0">
                            CIF: _em_tomography.axis1_min_angle
                            CIF: _em_tomography.axis2_min_angle
                            """
                            set_cif_value(axis.set_min_angle, min_angle, const.EM_TOMOGRAPHY, cif_list=ts_in, constructor=emdb.min_angleType, units=const.U_DEG)

                        def set_el_max_angle(axis, max_angle, ts_in):
                            """
                            XSD: <xs:element name="max_angle" minOccurs="0">
                            CIF: _em_tomography.axis1_max_angle
                            CIF: _em_tomography.axis2_max_angle
                            """
                            set_cif_value(axis.set_max_angle, max_angle, const.EM_TOMOGRAPHY, cif_list=ts_in, constructor=emdb.max_angleType, units=const.U_DEG)

                        def set_el_angle_increment(axis, angle_inc, ts_in):
                            """
                            XSD: <xs:element name="angle_increment" minOccurs="0">
                            CIF: _em_tomography.axis1_angle_increment
                            CIF: _em_tomography.axis2_angle_increment
                            """
                            set_cif_value(axis.set_angle_increment, angle_inc, const.EM_TOMOGRAPHY, cif_list=ts_in, constructor=emdb.max_angleType, units=const.U_DEG)

                        # element 1
                        set_el_min_angle(axis, min_angle, ts_in)
                        # element 2
                        set_el_max_angle(axis, max_angle, ts_in)
                        # element 3
                        set_el_angle_increment(axis, angle_inc, ts_in)

                    if axis1 is True:
                        axis_prefix = "axis1"
                    else:
                        axis_prefix = "axis2"

                    min_angle = axis_prefix + "_min_angle"
                    max_angle = axis_prefix + "_max_angle"
                    angle_inc = axis_prefix + "_angle_increment"

                    axis = emdb.axis_type()
                    set_axis_type(axis, min_angle, max_angle, angle_inc, ts_in)

                    return axis

                def set_special_optics(mic_id, mic):
                    """
                    Set special optics for microscopy

                    Parameters:
                    @param mic_id:
                    @param mic: microscopy object
                    XSD: <xs:complexType name="specialist_optics_type">
                    """

                    def set_specialist_optics_type(sp_op, sp_op_in):
                        """
                        XSD: <xs:complexType name="specialist_optics_type"> has
                        .. a sequence of 5 elements
                        """

                        def set_el_phase_plate(sp_op, sp_op_in):
                            """
                            XSD: <xs:element name="phase_plate" type="xs:token" minOccurs="0"/>
                            CIF: _em_imaging_optics.phase_plate 'Zernike phase plate'
                            """
                            set_cif_value(sp_op.set_phase_plate, "phase_plate", const.EM_SPECIALIST_OPTICS, cif_list=sp_op_in)

                        def set_el_sph_aberration_corrector(sp_op, sp_op_in):
                            """
                            XSD: <xs:element name="sph_aberration_corrector" type="xs:token" minOccurs="0"/>
                            CIF: _em_imaging_optics.sph_aberration_corrector
                            """
                            set_cif_value(sp_op.set_sph_aberration_corrector, "sph_aberration_corrector", const.EM_SPECIALIST_OPTICS, cif_list=sp_op_in)

                        def set_el_chr_aberration_corrector(sp_op, sp_op_in):
                            """
                            XSD: <xs:element name="chr_aberration_corrector" type="xs:token" minOccurs="0"/>
                            CIF: _em_imaging_optics.chr_aberration_corrector
                            """
                            set_cif_value(sp_op.set_chr_aberration_corrector, "chr_aberration_corrector", const.EM_SPECIALIST_OPTICS, cif_list=sp_op_in)

                        def set_el_energy_filter(sp_op, sp_op_in):
                            """
                            XSD: <xs:element name="energy_filter" minOccurs="0">
                            """

                            def set_energy_filter_type(eng_flt, sp_op_in):
                                """
                                XSD: <xs:element name="energy_filter" minOccurs="0"> has
                                .. 3 elements
                                """

                                def set_el_name(eng_flt, sp_op_in):
                                    """
                                    XSD: <xs:element name="name" type="xs:token" minOccurs="0">
                                    CIF: _em_imaging_optics.energyfilter_name 'FEI'
                                    """
                                    set_cif_value(eng_flt.set_name, "energyfilter_name", const.EM_SPECIALIST_OPTICS, cif_list=sp_op_in)

                                def set_el_slit_width(eng_flt, sp_op_in):
                                    """
                                    XSD: <xs:element name="slit_width" minOccurs="0">
                                    CIF: _em_imaging_optics.energyfilter_slit_width
                                    """
                                    set_cif_value(
                                        eng_flt.set_slit_width,
                                        "energyfilter_slit_width",
                                        const.EM_SPECIALIST_OPTICS,
                                        cif_list=sp_op_in,
                                        constructor=emdb.slit_widthType,
                                        units=const.U_EV,
                                    )

                                def set_el_lower_energy_threshold(eng_flt, sp_op_in):
                                    """
                                    XSD: <xs:element name="lower_energy_threshold" minOccurs="0">
                                    CIF: _em_imaging_optics.energyfilter_lower 0
                                    """
                                    eng_flt_low = get_cif_value("energyfilter_lower", const.EM_SPECIALIST_OPTICS, sp_op_in)
                                    if eng_flt_low is not None:
                                        if eng_flt_low.lstrip("-").lstrip("+").isdigit():
                                            set_cif_value(
                                                eng_flt.set_lower_energy_threshold,
                                                "energyfilter_lower",
                                                const.EM_SPECIALIST_OPTICS,
                                                cif_list=sp_op_in,
                                                constructor=emdb.lower_energy_thresholdType,
                                                units=const.U_EV,
                                            )
                                        else:
                                            # should be a float
                                            txt = u"The value for (_em_imaging_optics.energyfilter_lower) should not be: (%s)." % eng_flt_low
                                            self.current_entry_log.error_logs.append(
                                                self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                            )
                                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                                def set_el_upper_energy_threshold(eng_flt, sp_op_in):
                                    """
                                    XSD: <xs:element name="upper_energy_threshold" minOccurs="0">
                                    CIF: _em_imaging_optics.energyfilter_upper  15
                                    """
                                    eng_flt_uppr = get_cif_value("energyfilter_upper", const.EM_SPECIALIST_OPTICS, sp_op_in)
                                    if eng_flt_uppr is not None:
                                        if eng_flt_uppr.lstrip("-").lstrip("+").isdigit():
                                            set_cif_value(
                                                eng_flt.set_upper_energy_threshold,
                                                "energyfilter_upper",
                                                const.EM_SPECIALIST_OPTICS,
                                                cif_list=sp_op_in,
                                                constructor=emdb.upper_energy_thresholdType,
                                                units=const.U_EV,
                                            )
                                        else:
                                            # should be a float
                                            txt = u"The value for (_em_imaging_optics.energyfilter_upper) should not be: (%s)." % eng_flt_uppr
                                            self.current_entry_log.error_logs.append(
                                                self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                            )
                                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                                # element 1
                                set_el_name(eng_flt, sp_op_in)
                                # element 2
                                set_el_slit_width(eng_flt, sp_op_in)
                                # element 3
                                set_el_lower_energy_threshold(eng_flt, sp_op_in)
                                # element 3
                                set_el_upper_energy_threshold(eng_flt, sp_op_in)

                            eng_flt = emdb.energy_filterType()
                            set_energy_filter_type(eng_flt, sp_op_in)
                            if eng_flt.has__content():
                                sp_op.set_energy_filter(eng_flt)

                        def set_el_details(sp_op, sp_op_in):
                            """
                            XSD: <xs:element name="detals" type="xs:string" minOccurs="0"/>
                            CIF: _em_imaging_optics.details - NOT YET IMPLEMENTED IN mmcif dictionary
                            """
                            set_cif_value(sp_op.set_details, "details", const.EM_SPECIALIST_OPTICS, cif_list=sp_op_in)

                        # element 1
                        set_el_phase_plate(sp_op, sp_op_in)
                        # element 2
                        set_el_sph_aberration_corrector(sp_op, sp_op_in)
                        # element 3
                        set_el_chr_aberration_corrector(sp_op, sp_op_in)
                        # element 4
                        set_el_energy_filter(sp_op, sp_op_in)
                        # element 5
                        set_el_details(sp_op, sp_op_in)

                    # Create a specialist optics dictionary that has microscopy id as the key
                    sp_op_dict_in = make_dict(const.EM_SPECIALIST_OPTICS, const.K_IMAGING_ID)
                    if mic_id in sp_op_dict_in:
                        sp_op_in = sp_op_dict_in[mic_id]
                        sp_op = emdb.specialist_optics_type()
                        set_specialist_optics_type(sp_op, sp_op_in)
                        mic.set_specialist_optics(sp_op)

                def set_image_recording(mic_id, mic):
                    """
                    Set image recording for microscopy

                    Parameters:
                    @param mic: microscopy object
                    @param: mic_id: microscopy id
                    XSD: <xs:element name="image_recording">
                    """

                    def set_image_recording_type(im_rec, im_rec_in, im_dig_dict_in):
                        """
                        XSD: <xs:element name="image_recording"> has
                        .. 1 attribute and
                        .. 12 elements
                        """

                        def set_attr_id(im_rec, im_rec_in):
                            """
                            XSD: <xs:attribute name="image_recording_id" type="xs:positiveInteger"/>
                            CIF: _em_image_recording.id 1
                            """
                            set_cif_value(im_rec.set_image_recording_id, const.K_ID, const.EM_IMAGE_RECORDING, cif_list=im_rec_in, fmt=int)

                        def set_el_film_or_detector_model(im_rec, im_rec_in):
                            """
                            XSD: <xs:element name="film_or_detector_model">
                            CIF: _em_image_recording.film_or_detector_model 'AGFA SCIENTA FILM'
                            """
                            # film_or_detector_model has an attribute describing
                            # the type of detector - this does not seem to have been
                            # implemented in D&A
                            set_cif_value(
                                im_rec.set_film_or_detector_model,
                                "film_or_detector_model",
                                const.EM_IMAGE_RECORDING,
                                cif_list=im_rec_in,
                                constructor=emdb.film_or_detector_modelType,
                            )

                        def set_el_detector_mode(im_rec, im_rec_in):
                            """
                            XSD: <xs:element name="detector_mode" minOccurs="0">
                            CIF: _em_image_recording.detector_mode
                                {COUNTING,INTEGRATING,OTHER,SUPER-RESOLUTION}
                            """
                            set_cif_value(im_rec.set_detector_mode, "detector_mode", const.EM_IMAGE_RECORDING, cif_list=im_rec_in)

                        def set_el_digitization_details(im_rec, im_rec_in, im_dig_dict_in):
                            """
                            XSD: <xs:element name="digitization_details">
                            CIF: em_image_scans
                            """

                            def set_digitization_details_type(im_dig, im_dig_in):
                                """
                                XSD: <xs:element name="digitization_details"> has
                                .. 4 elements
                                """

                                def set_el_scanner(im_dig, im_dig_in):
                                    """
                                    XSD: <xs:element name="scanner" minOccurs="0">
                                    CIF: _em_image_scans.scanner_model 'EIKONIX IEEE 488'
                                    """
                                    set_cif_value(im_dig.set_scanner, "scanner_model", const.EM_IMAGE_SCANS, cif_list=im_dig_in)

                                def set_el_dimensions(im_dig, im_dig_in):
                                    """
                                    XSD: <xs:element name="dimensions">
                                    CIF: _em_image_scans.dimension_width 1034 (in pixels)
                                    CIF: _em_image_scans.dimension_height 1034 (in pixels)
                                    """
                                    width = get_cif_value("dimension_width", const.EM_IMAGE_SCANS, im_dig_in)
                                    height = get_cif_value("dimension_height", const.EM_IMAGE_SCANS, im_dig_in)
                                    if width is not None and height is not None:
                                        im_dig.set_dimensions(
                                            emdb.dimensionsType(
                                                width=emdb.widthType(valueOf_=int(width), units=const.U_PIXEL), height=emdb.heightType(valueOf_=int(height), units=const.U_PIXEL)
                                            )
                                        )
                                        txt = u"(im_dig.set_dimensions) set with width (%s) and height (%s)." % (width, height)
                                        self.current_entry_log.info_logs.append(
                                            self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.info_title + txt)
                                        )
                                        self.log_formatted(self.info_log_string, const.INFO_ALERT + txt)
                                    elif width is not None or height is not None:
                                        txt = u"(im_dig.set_dimensions) cannot be set since the width and height are not given. Both (_em_image_scans.dimension_height) and (_em_image_scans.dimension_height) values need to be given."
                                        self.current_entry_log.error_logs.append(
                                            self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                        )
                                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                                def set_el_sampling_interval(im_dig, im_dig_in):
                                    """
                                    XSD: <xs:element name="sampling_interval" minOccurs="0">
                                    CIF: _em_image_scans.sampling_interval 1.7(in microns)
                                    """
                                    set_cif_value(
                                        im_dig.set_sampling_interval,
                                        "sampling_interval",
                                        const.EM_IMAGE_SCANS,
                                        cif_list=im_dig_in,
                                        constructor=emdb.sampling_intervalType,
                                        fmt=float,
                                        units=const.U_MICROM,
                                    )

                                def set_el_frames_per_image(im_dig, im_dig_in):
                                    """
                                    XSD: <xs:element name="frames_per_image" type="xs:positiveInteger" minOccurs="0">
                                    CIF: _em_image_scans.used_frames_per_image 2-10
                                    """
                                    set_cif_value(im_dig.set_frames_per_image, "used_frames_per_image", const.EM_IMAGE_SCANS, cif_list=im_dig_in)

                                # element 1
                                set_el_scanner(im_dig, im_dig_in)
                                # element 2
                                set_el_dimensions(im_dig, im_dig_in)
                                # element 3
                                set_el_sampling_interval(im_dig, im_dig_in)
                                # element 4
                                set_el_frames_per_image(im_dig, im_dig_in)

                            im_rec_id_in = get_cif_value(const.K_ID, const.EM_IMAGE_RECORDING, im_rec_in)
                            if im_rec_id_in in im_dig_dict_in:
                                im_dig_in = im_dig_dict_in[im_rec_id_in]
                                im_dig = emdb.digitization_detailsType()
                                set_digitization_details_type(im_dig, im_dig_in)
                                im_rec.set_digitization_details(im_dig)

                        def set_el_number_grids_imaged(im_rec, im_rec_in):
                            """
                            XSD: <xs:element name="number_grids_imaged" type="xs:positiveInteger" minOccurs="0"/>
                            CIF: _em_image_recording.num_grids_imaged 10
                            """
                            set_cif_value(im_rec.set_number_grids_imaged, "num_grids_imaged", const.EM_IMAGE_RECORDING, cif_list=im_rec_in, fmt=int)

                        def set_el_number_real_images(im_rec, im_rec_in):
                            """
                            XSD: <xs:element name="number_real_images" type="xs:positiveInteger" minOccurs="0"/>
                            CIF: _em_image_recording.number_real_images 10
                            """
                            set_cif_value(im_rec.set_number_real_images, "num_real_images", const.EM_IMAGE_RECORDING, cif_list=im_rec_in, fmt=int)

                        def set_el_number_diffr_images(im_rec, im_rec_in):
                            """
                            XSD: <xs:element name="number_diffraction_images" type="xs:positiveInteger" minOccurs="0"/>
                            CIF: _em_image_recording.num_diffraction_images 10
                            """
                            set_cif_value(im_rec.set_number_diffraction_images, "num_diffraction_images", const.EM_IMAGE_RECORDING, cif_list=im_rec_in, fmt=int)

                        def set_el_average_exposure_time(im_rec, im_rec_in):
                            """
                            XSD: <xs:element name="average_exposure_time" minOccurs="0">
                            CIF: _em_image_recording.average_exposure_time 2.0
                            """
                            set_cif_value(
                                im_rec.set_average_exposure_time,
                                "average_exposure_time",
                                const.EM_IMAGE_RECORDING,
                                cif_list=im_rec_in,
                                constructor=emdb.average_exposure_timeType,
                                fmt=float,
                                units=const.U_SEC,
                            )

                        def set_el_av_el_dose_per_image(im_rec, im_rec_in):
                            """
                            XSD: <xs:element name="average_electron_dose_per_image" minOccurs="0">
                            CIF: _em_image_recording.avg_electron_dose_per_image
                            """
                            set_cif_value(
                                im_rec.set_average_electron_dose_per_image,
                                "avg_electron_dose_per_image",
                                const.EM_IMAGE_RECORDING,
                                cif_list=im_rec_in,
                                constructor=emdb.average_electron_dose_per_imageType,
                                fmt=float,
                                units=const.U_EOVERANGSQR,
                            )

                        def set_el_detector_distance():
                            """
                            XSD: <xs:element name="detector_distance" type="xs:string" minOccurs="0"/>
                            CIF: detectorDistance in emdb.xsd
                            Detector distance in not there in the cif!
                            """

                        def set_el_details(im_rec, im_rec_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _em_image_recording.details
                            """
                            set_cif_value(im_rec.set_details, "details", const.EM_IMAGE_RECORDING, cif_list=im_rec_in)

                        def set_el_od_range():
                            """
                            XSD: <xs:element name="od_range" type="xs:float" minOccurs="0">
                            CIF:?? _em_image_scans.od_range??
                            """

                        def set_el_bits_per_pixel():
                            """
                            XSD: <xs:element name="bits_per_pixel" type="xs:float" minOccurs="0">
                            CIF:??
                            """

                        # attribute 1
                        set_attr_id(im_rec, im_rec_in)
                        # element 1
                        set_el_film_or_detector_model(im_rec, im_rec_in)
                        # element 2
                        set_el_detector_mode(im_rec, im_rec_in)
                        # element 3
                        set_el_digitization_details(im_rec, im_rec_in, im_dig_dict_in)
                        # element 4
                        set_el_number_grids_imaged(im_rec, im_rec_in)
                        # element 5
                        set_el_number_real_images(im_rec, im_rec_in)
                        # element 6
                        set_el_number_diffr_images(im_rec, im_rec_in)
                        # element 7
                        set_el_average_exposure_time(im_rec, im_rec_in)
                        # element 8
                        set_el_av_el_dose_per_image(im_rec, im_rec_in)
                        # element 9
                        set_el_detector_distance()
                        # element 10
                        set_el_details(im_rec, im_rec_in)
                        # element 11
                        set_el_od_range()
                        # element 12
                        set_el_bits_per_pixel()

                    im_rec_dict_in = make_list_of_dicts(const.EM_IMAGE_RECORDING, const.K_IMAGING_ID)
                    im_dig_dict_in = make_dict(const.EM_IMAGE_SCANS, const.K_IMAGE_RECORDING_ID)

                    if mic_id not in im_rec_dict_in:
                        txt = u"No value for (_em_image_recording) found for microscope: (%s)." % mic_id
                        self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                    else:
                        im_rec_list_in = im_rec_dict_in[mic_id]
                        im_rec_list = emdb.image_recording_listType()
                        for im_rec_in in im_rec_list_in:
                            im_rec = emdb.image_recordingType()
                            set_image_recording_type(im_rec, im_rec_in, im_dig_dict_in)
                            if im_rec.has__content():
                                im_rec_list.add_image_recording(im_rec)

                        if im_rec_list.has__content():
                            mic.set_image_recording_list(im_rec_list)

                def set_base_microscopy(mic_in, mic, mic_id):
                    """
                    Set all elements for base microscopy

                    Parameters:
                    @param mic_in: cif dictionary item for EM_IMAGING
                    @param mic: microscopy object
                    @param: mic_id: microscopy id
                    XSD: <xs:complexType name="base_microscopy_type"> has
                    .. 1 attribute
                    .. 26 elements
                    """

                    def set_attr_id(mic):
                        """
                        XSD: <xs:attribute name="microscopy_id" type="xs:positiveInteger" use="required"/>
                        CIF: _em_imaging.id 1
                        """
                        set_cif_value(mic.set_microscopy_id, "id", const.EM_IMAGING, fmt=int)

                    def set_el_specimen_preparations():
                        """
                        XSD: <xs:element name="specimen_preparations" minOccurs="0">
                        CIF: not in cif. needs to be harvested from different cif sources
                        """

                    def set_el_microscope(mic, mic_in):
                        """
                        XSD: <xs:element name="microscope">
                        CIF: _em_imaging.microscope_model 'FEI MORGAGNI'
                        """
                        set_cif_value(mic.set_microscope, "microscope_model", const.EM_IMAGING, cif_list=mic_in)

                    def set_el_illumination_mode(mic, mic_in):
                        """
                        XSD: <xs:element name="illumination_mode">
                        CIF: _em_imaging.illumination_mode 'FLOOD BEAM'
                        """
                        set_cif_value(mic.set_illumination_mode, "illumination_mode", const.EM_IMAGING, cif_list=mic_in)

                    def set_el_imaging_mode(mic, mic_in):
                        """
                        XSD: <xs:element name="imaging_mode">
                        CIF: _em_imaging.mode 'BRIGHT FIELD'
                        """
                        set_cif_value(mic.set_imaging_mode, "mode", const.EM_IMAGING, cif_list=mic_in)

                    def set_el_electron_source(mic, mic_in):
                        """
                        XSD: <xs:element name="electron_source">
                        CIF: _em_imaging.electron_source 'FIELD EMISSION GUN'
                        """
                        set_cif_value(mic.set_electron_source, "electron_source", const.EM_IMAGING, cif_list=mic_in)

                    def set_el_acceleration_voltage(mic, mic_in):
                        """
                        XSD: <xs:element name="acceleration_voltage">
                        CIF: _em_imaging.accelerating_voltage 300
                        """
                        set_cif_value(
                            mic.set_acceleration_voltage,
                            "accelerating_voltage",
                            const.EM_IMAGING,
                            cif_list=mic_in,
                            constructor=emdb.acceleration_voltageType,
                            fmt=int,
                            units=const.U_KVOLT,
                        )

                    def set_el_c2_aperture_diameter(mic, mic_in):
                        """
                        XSD: <xs:element name="c2_aperture_diameter" minOccurs="0">
                        CIF: _em_imaging.c2_aperture_diameter 100
                        """
                        set_cif_value(
                            mic.set_c2_aperture_diameter,
                            "c2_aperture_diameter",
                            const.EM_IMAGING,
                            cif_list=mic_in,
                            constructor=emdb.c2_aperture_diameterType,
                            fmt=float,
                            units=const.U_MICROM,
                        )

                    def set_el_nominal_cs(mic, mic_in):
                        """
                        XSD: <xs:element name="nominal_cs" minOccurs="0">
                        CIF: _em_imaging.nominal_cs 2.7
                        """
                        set_cif_value(mic.set_nominal_cs, "nominal_cs", const.EM_IMAGING, cif_list=mic_in, constructor=emdb.nominal_csType, fmt=float, units=const.U_MM)

                    def set_el_nominal_defocus_min(mic, mic_in):
                        """
                        XSD: <xs:element name="nominal_defocus_min" minOccurs="0">
                        CIF: _em_imaging.nominal_defocus_min  1200 (in nm)
                        """
                        set_cif_value(
                            mic.set_nominal_defocus_min,
                            "nominal_defocus_min",
                            const.EM_IMAGING,
                            cif_list=mic_in,
                            constructor=emdb.nominal_defocus_minType,
                            units=const.U_MICROM,
                            fmt=lambda x: float(x) * 0.001,
                        )

                    def set_el_calibrated_defocus_min(mic, mic_in):
                        """
                        XSD: <xs:element name="calibrated_defocus_min" minOccurs="0">
                        CIF: _em_imaging.calibrated_defocus_min 1200 (in nm)
                        """
                        set_cif_value(
                            mic.set_calibrated_defocus_min,
                            "calibrated_defocus_min",
                            const.EM_IMAGING,
                            cif_list=mic_in,
                            constructor=emdb.calibrated_defocus_minType,
                            units=const.U_MICROM,
                            fmt=lambda x: float(x) * 0.001,
                        )

                    def set_el_nominal_defocus_max(mic, mic_in):
                        """
                        XSD: <xs:element name="nominal_defocus_max" minOccurs="0">
                        CIF: _em_imaging.nominal_defocus_max 5000 (in nm)
                        """
                        nom_def = get_cif_value("nominal_defocus_max", const.EM_IMAGING, cif_list=mic_in)
                        if nom_def is not None:
                            fl_nom_fel = float(nom_def) * 0.001
                            if fl_nom_fel < -20:
                                txt = u"_em_imaging.nominal_defocus_max (%s) is less than -20." % fl_nom_fel
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.not_changed_for_now_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.NOT_CHANGED_FOR_NOW + txt)
                            elif fl_nom_fel > 20:
                                txt = u"_em_imaging.nominal_defocus_max (%s) is larger than 20." % fl_nom_fel
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.not_changed_for_now_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.NOT_CHANGED_FOR_NOW + txt)
                            set_cif_value(
                                mic.set_nominal_defocus_max,
                                "nominal_defocus_max",
                                const.EM_IMAGING,
                                cif_list=mic_in,
                                constructor=emdb.nominal_defocus_maxType,
                                units=const.U_MICROM,
                                fmt=lambda x: float(x) * 0.001,
                            )

                    def set_el_calibrated_defocus_max(mic, mic_in):
                        """
                        XSD: <xs:element name="calibrated_defocus_max" minOccurs="0">
                        CIF: _  5000 (in nm)
                        """
                        cal_def = get_cif_value("calibrated_defocus_max", const.EM_IMAGING, cif_list=mic_in)
                        if cal_def is not None:
                            fl_cal_fel = float(cal_def) * 0.001
                            if fl_cal_fel < -20:
                                txt = u"The value given to (_em_imaging.calibrated_defocus_max): (%s) is less than -20." % fl_cal_fel
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.not_changed_for_now_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.NOT_CHANGED_FOR_NOW + txt)
                            elif fl_cal_fel > 20:
                                txt = u"The value given to (_em_imaging.calibrated_defocus_max): (%s) is larger than 20." % fl_cal_fel
                                self.current_entry_log.warn_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.not_changed_for_now_title + txt)
                                )
                                self.log_formatted(self.warn_log_string, const.NOT_CHANGED_FOR_NOW + txt)
                            set_cif_value(
                                mic.set_calibrated_defocus_max,
                                "calibrated_defocus_max",
                                const.EM_IMAGING,
                                cif_list=mic_in,
                                constructor=emdb.calibrated_defocus_maxType,
                                units=const.U_MICROM,
                                fmt=lambda x: float(x) * 0.001,
                            )

                    def set_el_nominal_magnification(mic, mic_in):
                        """
                        XSD: <xs:element name="nominal_magnification" type="allowed_magnification" minOccurs="0">
                        CIF: _em_imaging.nominal_magnification 50000
                        """
                        set_cif_value(mic.set_nominal_magnification, "nominal_magnification", const.EM_IMAGING, cif_list=mic_in, fmt=float)

                    def set_el_calibrated_magnification(mic, mic_in):
                        """
                        XSD: <xs:element name="calibrated_magnification" type="allowed_magnification" minOccurs="0">
                        CIF: _em_imaging.calibrated_magnification 50230
                        """
                        set_cif_value(mic.set_calibrated_magnification, "calibrated_magnification", const.EM_IMAGING, cif_list=mic_in, fmt=float)

                    def set_el_specimen_holder_model(mic, mic_in):
                        """
                        XSD: <xs:element name="specimen_holder_model" minOccurs="0">
                        CIF: _em_imaging.specimen_holder_model 'FEI TITAN KRIOS AUTOGRID HOLDER'
                        """
                        set_cif_value(mic.set_specimen_holder_model, "specimen_holder_model", const.EM_IMAGING, cif_list=mic_in)

                    def set_el_cooling_holder_cryogen(mic, mic_in):
                        """
                        XSD: <xs:element name="cooling_holder_cryogen" minOccurs="0">
                        CIF: _em_imaging.cryogen HELIUM
                        """
                        set_cif_value(mic.set_cooling_holder_cryogen, "cryogen", const.EM_IMAGING, cif_list=mic_in)

                    def set_el_temperature(mic, mic_in):
                        """
                        XSD: <xs:element name="temperature" minOccurs="0">
                        CIF: _em_imaging.recording_temperature_maximum 70 (deg K)
                        CIF: _em_imaging.recording_temperature_minimum 70 (deg K)
                        """
                        temp_max_in = get_cif_value("recording_temperature_maximum", const.EM_IMAGING, mic_in)
                        temp_min_in = get_cif_value("recording_temperature_minimum", const.EM_IMAGING, mic_in)
                        if temp_max_in is not None or temp_min_in is not None:
                            temp = emdb.temperatureType()
                            set_cif_value(
                                temp.set_temperature_max, "recording_temperature_maximum", const.EM_IMAGING, cif_list=mic_in, constructor=emdb.temperature_type, fmt=float, units=const.U_KEL
                            )
                            set_cif_value(
                                temp.set_temperature_min, "recording_temperature_minimum", const.EM_IMAGING, cif_list=mic_in, constructor=emdb.temperature_type, fmt=float, units=const.U_KEL
                            )
                            mic.set_temperature(temp)

                    def set_el_alignment_procedure(mic, mic_in):
                        """
                        XSD: <xs:element name="alignment_procedure" minOccurs="0">
                        CIF: _em_imaging.alignment_procedure {BASIC,COMA FREE,NONE,OTHER,ZEMLIN TABLEAU}
                        CIF: _em_imaging.residual_tilt (Residual tilt, in milliradians)
                        """
                        align_proc = get_cif_value("alignment_procedure", const.EM_IMAGING, mic_in)
                        if align_proc is not None:
                            ali = emdb.alignment_procedureType()
                            tilt = get_cif_value("residual_tilt", const.EM_IMAGING, mic_in)
                            if align_proc == "NONE":
                                ali_none = emdb.noneType()
                                ali.set_none(ali_none)
                            elif align_proc == "BASIC":
                                ali_basic = emdb.basicType()
                                if tilt is not None:
                                    ali_basic.set_residual_tilt(emdb.residual_tilt_type(valueOf_=float(tilt), units=const.U_MRAD))
                                ali.set_basic(ali_basic)
                            elif align_proc == "ZEMLIN TABLEAU":
                                ali_zem = emdb.zemlin_tableauType()
                                ali.set_zemlin_tableau(ali_zem)
                            elif align_proc == "COMA FREE":
                                ali_cf = emdb.coma_freeType()
                                if tilt is not None:
                                    ali_cf.set_residual_tilt(emdb.residual_tilt_type(valueOf_=float(tilt), units=const.U_MRAD))
                                ali.set_coma_free(ali_cf)
                            elif align_proc == "OTHER":
                                ali_other = emdb.otherType()
                                ali.set_other(ali_other)
                            if ali.has__content():
                                mic.set_alignment_procedure(ali)

                    def set_el_specialist_optics(mic, mic_id):
                        """
                        XSD: <xs:element name="specialist_optics" type="specialist_optics_type" minOccurs="0"/>
                        CIF: em_imaging_optics
                        """
                        set_special_optics(mic_id, mic)

                    def set_el_software_list():
                        """
                        XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                        CIF:
                        """

                    def set_el_details(mic, mic_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                        CIF: _em_imaging.details
                        """
                        set_cif_value(mic.set_details, "details", const.EM_IMAGING, cif_list=mic_in)

                    def set_el_date():
                        """
                        XSD: <xs:element name="date" type="xs:date" minOccurs="0">
                        Deprecated (2014-10-21)
                        """

                    def set_el_image_recording_list(mic_id, mic):
                        """
                        XSD: <xs:element name="image_recording_list">
                        """
                        set_image_recording(mic_id, mic)

                    def set_el_specimen_holder():
                        """
                        XSD: <xs:element name="specimen_holder" type="xs:string" minOccurs="0">
                        Deprecated (2014/11/12)
                        """

                    def set_el_tilt_angle_min():
                        """
                        XSD: <xs:element name="tilt_angle_min" minOccurs="0">
                        Deprecated (2014/11/14)
                        """

                    def set_el_tilt_angle_max():
                        """
                        XSD: <xs:element name="tilt_angle_max" minOccurs="0">
                        Deprecated (2014/11/14)
                        """

                    # attribute 1
                    set_attr_id(mic)
                    # element 1
                    set_el_specimen_preparations()
                    # element 2
                    set_el_microscope(mic, mic_in)
                    # element 3
                    set_el_illumination_mode(mic, mic_in)
                    # element 4
                    set_el_imaging_mode(mic, mic_in)
                    # element 5
                    set_el_electron_source(mic, mic_in)
                    # element 6
                    set_el_acceleration_voltage(mic, mic_in)
                    # element 7
                    set_el_c2_aperture_diameter(mic, mic_in)
                    # element 8
                    set_el_nominal_cs(mic, mic_in)
                    # element 9
                    set_el_nominal_defocus_min(mic, mic_in)
                    # element 10
                    set_el_calibrated_defocus_min(mic, mic_in)
                    # element 11
                    set_el_nominal_defocus_max(mic, mic_in)
                    # element 12
                    set_el_calibrated_defocus_max(mic, mic_in)
                    # element 13
                    set_el_nominal_magnification(mic, mic_in)
                    # element 14
                    set_el_calibrated_magnification(mic, mic_in)
                    # element 15
                    set_el_specimen_holder_model(mic, mic_in)
                    # element 16
                    set_el_cooling_holder_cryogen(mic, mic_in)
                    # element 17
                    set_el_temperature(mic, mic_in)
                    # element 18
                    set_el_alignment_procedure(mic, mic_in)
                    # element 19
                    set_el_specialist_optics(mic, mic_id)
                    # element 20
                    set_el_software_list()
                    # element 21
                    set_el_details(mic, mic_in)
                    # element 22
                    set_el_date()
                    # element 23
                    set_el_image_recording_list(mic_id, mic)
                    # element 24
                    set_el_specimen_holder()
                    # element 25
                    set_el_tilt_angle_min()
                    # element 26
                    set_el_tilt_angle_max()

                def set_tilt_parameters(mic_id, mic):
                    """
                    Set parameters for tilt series in tomography and EC microscopy.
                    There can be more than one tilt series per microscopy.

                    Parameters:
                    @param mic: microscopy object
                    @param: mic_id: microscopy id
                    XSD: <xs:element name="tilt_series" type="tilt_series_type" maxOccurs="unbounded">
                    """

                    def set_tilt_series_type(tilt_ser, ts_in):
                        """
                        XSD: <xs:complexType name="tilt_series_type"> has
                        .. 3 elements
                        """

                        def set_el_axis1(tilt_ser, ts_in):
                            """
                            XSD: <xs:element name="axis1" type="axis_type"/>
                            CIF:
                            """
                            axis1 = get_tilt_axis(ts_in, axis1=True)
                            if axis1.has__content():
                                set_cif_value(tilt_ser.set_axis1, axis1)

                        def set_el_axis2(tilt_ser, ts_in):
                            """
                            XSD: <xs:element name="axis2" minOccurs="0"> has
                            .. a base of axis_type and
                            .. 1 element <xs:element name="axis_rotation" minOccurs="0">???
                            CIF:
                            """
                            axis2 = get_tilt_axis(ts_in, axis1=False)
                            if axis2.has__content():
                                set_cif_value(tilt_ser.set_axis2, axis2)

                        def set_el_axis_rotation(tilt_ser, ts_in):
                            """
                            XSD: <xs:element name="axis_rotation" fixed="90" minOccurs="0">
                            CIF: _em_tomography.dual_tilt_axis_rotation
                            """
                            set_cif_value(
                                tilt_ser.set_axis_rotation,
                                "dual_tilt_axis_rotation",
                                const.EM_TOMOGRAPHY,
                                cif_list=ts_in,
                                constructor=emdb.axis_rotationType,
                                units=const.U_DEG,
                            )

                        # element 1
                        set_el_axis1(tilt_ser, ts_in)
                        # element 2
                        set_el_axis2(tilt_ser, ts_in)
                        # element 3
                        set_el_axis_rotation(tilt_ser, ts_in)

                    # Create list of tomography tilt series that has microscopy id as a key.
                    tilt_dict_in = make_list_of_dicts(const.EM_TOMOGRAPHY, const.K_IMAGING_ID)
                    if mic_id in tilt_dict_in:
                        tilt_list_in = tilt_dict_in[mic_id]
                        for ts_in in tilt_list_in:
                            tilt_ser = emdb.tilt_series_type()
                            set_tilt_series_type(tilt_ser, ts_in)
                            if tilt_ser.has__content():
                                mic.add_tilt_series(tilt_ser)

                def set_cryst_mic_specifics(mic_id, cryst_mic):
                    """
                    Method that sets the crystalograhy specific microscopy elements

                    Parameters:
                    @param mic_id: microscopy id
                    @param cryst_mic: crystallography microscopy object
                    XSD: <xs:complexType name="crystallography_microscopy_type"> has
                    .. a base (see: set_base_microscopy) and
                    .. 1 element and
                    .. a choice between 2 elements
                    """

                    def set_el_camera_length(cryst_mic, cry_mic_in):
                        """
                        XSD: <xs:element name="camera_length">
                        CIF: _em_diffraction.camera_length
                        """
                        set_cif_value(
                            cryst_mic.set_camera_length,
                            "camera_length",
                            const.EM_DIFFRACTION,
                            cif_list=cry_mic_in,
                            constructor=emdb.camera_lengthType,
                            units=const.U_MM,
                        )

                    cry_mic_dict_in = make_dict(const.EM_DIFFRACTION, const.K_IMAGING_ID)
                    if mic_id in cry_mic_dict_in:
                        cry_mic_in = cry_mic_dict_in[mic_id]
                        # element 1
                        set_el_camera_length(cryst_mic, cry_mic_in)
                        # choice 1
                        # CIF: _em_diffraction.tilt_angle_list 20,40,50,55
                        tilt_str_in = get_cif_value("tilt_angle_list", const.EM_DIFFRACTION, cry_mic_in)
                        if tilt_str_in is not None:
                            tilt_list_in = tilt_str_in.strip().split(",")
                            tilt_list = emdb.tilt_listType()
                            for tilt_in in tilt_list_in:
                                # XSD: <xs:element name="tilt_list" minOccurs="1">
                                tilt_list.add_angle(float(tilt_in))
                            if tilt_list.has__content():
                                cryst_mic.set_tilt_list(tilt_list)
                                # choice 2 (IS THIS NECESSARY ????????????)
                                # XSD: <xs:element name="tilt_series"
                                #      type="tilt_series_type" maxOccurs="unbounded" minOccurs="1">

                mic_list_in = self.cif.get(const.EM_IMAGING, None)
                if mic_list_in != []:
                    for mic_in in mic_list_in:
                        mic_id = get_cif_value(const.K_ID, const.EM_IMAGING, mic_in)
                        if em_method == const.EMM_SP:
                            sp_mic = emdb.single_particle_microscopy_type()
                            sp_mic.original_tagname_ = "single_particle_microscopy"
                            set_base_microscopy(mic_in, sp_mic, mic_id)
                            microscopy_list.add_microscopy(sp_mic)
                        elif em_method == const.EMM_HEL:
                            hel_mic = emdb.helical_microscopy_type()
                            hel_mic.original_tagname_ = "helical_microscopy"
                            set_base_microscopy(mic_in, hel_mic, mic_id)
                            microscopy_list.add_microscopy(hel_mic)
                        elif em_method == const.EMM_TOM:
                            tom_mic = emdb.tomography_microscopy_type()
                            tom_mic.original_tagname_ = "tomography_microscopy"
                            set_base_microscopy(mic_in, tom_mic, mic_id)
                            set_tilt_parameters(mic_id, tom_mic)
                            microscopy_list.add_microscopy(tom_mic)
                        elif em_method == const.EMM_STOM:
                            subtom_mic = emdb.tomography_microscopy_type()
                            subtom_mic.original_tagname_ = "subtomogram_averaging_microscopy"
                            set_base_microscopy(mic_in, subtom_mic, mic_id)
                            set_tilt_parameters(mic_id, subtom_mic)
                            microscopy_list.add_microscopy(subtom_mic)
                        elif em_method == const.EMM_EC:
                            cryst_mic = emdb.crystallography_microscopy_type()
                            cryst_mic.original_tagname_ = "crystallography_microscopy"
                            set_base_microscopy(mic_in, cryst_mic, mic_id)
                            set_tilt_parameters(mic_id, cryst_mic)
                            set_cryst_mic_specifics(mic_id, cryst_mic)
                            microscopy_list.add_microscopy(cryst_mic)
                    # all done now: set the list
                    struct_det.set_microscopy_list(microscopy_list)
                else:
                    txt = u"CIF category (%s) is missing." % const.EM_IMAGING
                    self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

            def set_el_image_processing(struct_det, microscopy_list):
                """
                XSD: <xs:element ref="image_processing" maxOccurs="unbounded">
                """

                def set_base_image_processing(ip_in, im_proc):
                    """
                    Set all elements for base image processing type

                    Parameters:
                    @param ip_in: cif dictionary item
                    @param im_proc: image processing object
                    XSD: <xs:complexType name="base_image_processing_type"> has
                    .. 1 attributes
                    .. 2 elements
                    """

                    def set_attr_id(ip_in, im_proc):
                        """
                        XSD: <xs:attribute name="image_processing_id" type="xs:positiveInteger" use="required"/>
                        CIF: _em_image_processing.id
                        """
                        set_cif_value(im_proc.set_image_processing_id, const.K_ID, const.EM_IMAGE_PROCESSING, cif_list=ip_in, fmt=int)

                    def set_el_image_recording_id(ip_in, im_proc):
                        """
                        XSD: <xs:element name="image_recording_id" type="xs:positiveInteger"/>
                        CIF: _em_image_processing.image_recording_id
                        """
                        set_cif_value(im_proc.set_image_recording_id, const.K_IMAGE_RECORDING_ID, const.EM_IMAGE_PROCESSING, cif_list=ip_in, fmt=int)

                    def set_el_details(ip_in, im_proc):
                        """
                        XSD: <xs:element name="details" type="xs:token" minOccurs="0"/>
                        CIF: _em_image_processing.details
                        """
                        set_cif_value(im_proc.set_details, "details", const.EM_IMAGE_PROCESSING, cif_list=ip_in)

                    # attribute 1
                    set_attr_id(ip_in, im_proc)
                    # element 1
                    set_el_image_recording_id(ip_in, im_proc)
                    # element 2
                    set_el_details(ip_in, im_proc)

                def set_particle_selection(ps_in, im_proc, cat_soft_dict_in):
                    """
                    Sets particle selection

                    @param ps_in: cif categories for single particle ip
                    @param im_proc: image processing object
                    @param cat_soft_dict_in
                    XSD: <xs:complexType name="particle_selection_type">
                    """

                    def set_particle_selection_type(part_sel):
                        """
                        XSD: <xs:complexType name="particle_selection_type"> has
                        ..5 elements
                        """

                        def set_el_num_particles_sel(part_sel, ps_in):
                            """
                            XSD: <xs:element name="number_selected" type="xs:positiveInteger"/>
                            CIF: _em_particle_selection.num_particles_selected 840
                            """
                            set_cif_value(part_sel.set_number_selected, "num_particles_selected", const.EM_PARTICLE_SELECTION, cif_list=ps_in, fmt=int)

                        def set_el_reference_model(part_sel, ps_in):
                            """
                            XSD: <xs:element name="reference_model" type="xs:token" minOccurs="0">
                            CIF: _em_particle_selection.reference_model
                            """
                            set_cif_value(part_sel.set_reference_model, "reference_model", const.EM_PARTICLE_SELECTION, cif_list=ps_in)

                        def set_el_method(part_sel, ps_in):
                            """
                            XSD: <xs:element name="method" type="xs:string" minOccurs="0" >
                            CIF: _em_particle_selection.method
                            """
                            set_cif_value(part_sel.set_method, "method", const.EM_PARTICLE_SELECTION, cif_list=ps_in)

                        def set_el_software_list(part_sel):
                            """
                            XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                            """
                            set_software_list(const.SOFT_PARTICLE_SELECTION, cat_soft_dict_in, part_sel.set_software_list)

                        def set_el_details(part_sel, ps_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _em_particle_selection.details .
                            """
                            set_cif_value(part_sel.set_details, "details", const.EM_PARTICLE_SELECTION, cif_list=ps_in)

                        # element 1
                        set_el_num_particles_sel(part_sel, ps_in)
                        # element 2
                        set_el_reference_model(part_sel, ps_in)
                        # element 3
                        set_el_method(part_sel, ps_in)
                        # element 4
                        set_el_software_list(part_sel)
                        # element 5
                        set_el_details(part_sel, ps_in)

                    part_sel = emdb.particle_selection_type()
                    set_particle_selection_type(part_sel)
                    if part_sel.has__content():
                        im_proc.add_particle_selection(part_sel)

                def set_ctfcorrection(ip_id_in, im_proc, ctf_corr_dict_in):
                    """
                    Sets CTF correction

                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param ctf_corr_dict_in: dictionary for CTF correction
                    XSD: <xs:complexType name="ctf_correction_type">
                    """

                    def set_ctf_correction_type(ctf_corr):
                        """
                        XSD: <xs:complexType name="ctf_correction_type"> has
                        .. 5 elements
                        """

                        def set_el_phase_reversal(ctf_corr, ctf_corr_in):
                            """
                            XSD: <xs:element name="phase_reversal" minOccurs="0">
                            CIF: _em_ctf_correction.phase_reversal YES/NO
                            """

                            def set_phase_reversal_type(ph_rev, ctf_corr_in):
                                """
                                XSD: <xs:element name="phase_reversal" minOccurs="0"> has
                                .. 2 elements
                                """

                                def set_el_anisotropic(ph_rev):
                                    """
                                    XSD: <xs:element name="anisotropic" type="xs:boolean" minOccurs="0"/>
                                    CIF: _em_ctf_correction.phase_reversal_anisotropic YES/NO
                                    """
                                    set_cif_value(ph_rev.set_anisotropic, "phase_reversal_anisotropic", const.EM_CTF_CORRECTION, cif_list=ctf_corr_in)

                                def set_el_correction_space(ph_rev, ctf_corr_in):
                                    """
                                    XSD:: <xs:element name="correction_space" type="correction_space_type"  minOccurs="0"/>
                                    CIF: _em_ctf_correction.phase_reversal_correction_space
                                        REAL/RECIPROCAL
                                    """
                                    pra = get_cif_value("phase_reversal_anisotropic", const.EM_CTF_CORRECTION, ctf_corr_in)
                                    if pra == "YES":
                                        # Yes if Anisotropic phase reversal (flipping) was performed
                                        set_cif_value(ph_rev.set_correction_space, "phase_reversal_correction_space", const.EM_CTF_CORRECTION, cif_list=ctf_corr_in)

                                # element 1
                                set_el_anisotropic(ph_rev)
                                # element 2
                                set_el_correction_space(ph_rev, ctf_corr_in)

                            ph_rev = get_cif_value("phase_reversal", const.EM_CTF_CORRECTION, ctf_corr_in)
                            if ph_rev == "YES":
                                # Yes if Phase reversal (flipping) was performed
                                ph_rev = emdb.phase_reversalType()
                                set_phase_reversal_type(ph_rev, ctf_corr_in)
                                ctf_corr.set_phase_reversal(ph_rev)

                        def set_el_amplitude_correction(ctf_corr, ctf_corr_in):
                            """
                            XSD: <xs:element name="amplitude_correction" minOccurs="0">
                            CIF: _em_ctf_correction.amplitude_correction YES/NO
                            """

                            def set_amplitude_correction_type(am_corr, ctf_corr_in):
                                """
                                XSD: <xs:element name="amplitude_correction" minOccurs="0"> has
                                .. 2 elements
                                """

                                def set_el_factor(am_corr, ctf_corr_in):
                                    """
                                    XSD: <xs:element name="factor" type="xs:float" minOccurs="0">
                                    CIF: _em_ctf_correction.amplitude_correction_factor
                                    """
                                    set_cif_value(am_corr.set_factor, "amplitude_correction", const.EM_CTF_CORRECTION, cif_list=ctf_corr_in)

                                def set_el_correction_space(am_corr, ctf_corr_in):
                                    """
                                    XSD: <xs:element name="correction_space" type="correction_space_type" minOccurs="0"/>
                                    CIF: _em_ctf_correction.amplitude_correction_space
                                        REAL/RECIPROCAL
                                    """
                                    set_cif_value(am_corr.set_correction_space, "amplitude_correction_space", const.EM_CTF_CORRECTION, cif_list=ctf_corr_in)

                                # element 1
                                set_el_factor(am_corr, ctf_corr_in)
                                # element 2
                                set_el_correction_space(am_corr, ctf_corr_in)

                            amp_corr = get_cif_value("amplitude_correction", const.EM_CTF_CORRECTION, ctf_corr_in)
                            if amp_corr is not None:
                                am_corr = emdb.amplitude_correctionType()
                                set_amplitude_correction_type(am_corr, ctf_corr_in)
                                ctf_corr.set_amplitude_correction(am_corr)
                            else:
                                # ADD LOGGER MESSAGE
                                pass

                        def set_el_correction_operation(ctf_corr, ctf_corr_in):
                            """
                            XSD: <xs:element name="correction_operation" minOccurs="0">
                            CIF: _em_ctf_correction.correction_operation
                            """
                            set_cif_value(ctf_corr.set_correction_operation, "correction_operation", const.EM_CTF_CORRECTION, cif_list=ctf_corr_in)

                        def set_el_software_list(ctf_corr, cat_soft_dict_in):
                            """
                            XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                            """
                            set_software_list(const.SOFT_CTF_CORRECTION, cat_soft_dict_in, ctf_corr.set_software_list)

                        def set_el_details(ctf_corr, ctf_corr_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _em_ctf_correction.details
                            """
                            set_cif_value(ctf_corr.set_details, "details", const.EM_CTF_CORRECTION, cif_list=ctf_corr_in)

                        ctf_corr_in = ctf_corr_dict_in[ip_id_in]
                        # element 1
                        set_el_phase_reversal(ctf_corr, ctf_corr_in)
                        # element 2
                        set_el_amplitude_correction(ctf_corr, ctf_corr_in)
                        # element 3
                        set_el_correction_operation(ctf_corr, ctf_corr_in)
                        # element 4
                        set_el_software_list(ctf_corr, cat_soft_dict_in)
                        # element 5
                        set_el_details(ctf_corr, ctf_corr_in)

                    ctf_corr = emdb.ctf_correction_type()
                    set_ctf_correction_type(ctf_corr)
                    if ctf_corr.has__content():
                        im_proc.set_ctf_correction(ctf_corr)

                def set_startup_model(ip_id_in, im_proc, st_mod_dict_in):
                    """
                    Sets startup model

                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param st_mod_dict_in: the dictionary for startup model
                    XSD: <xs:complexType name="starting_map_type">
                    """

                    def set_starting_map_type(st_map, sm_in):
                        """
                        XSD: <xs:complexType name="starting_map_type"> has
                        .. 1 attribute
                        .. 1 element of {1 choice of 6 elements + 1 element}
                        """

                        def set_attr_type_of_model(st_map, sm_in):
                            """
                            XSD: <xs:attribute name="type_of_model" type="xs:token"/>
                            CIF: _em_start_model.type
                            {EMDB MAP,INSILICO MODEL,NONE,ORTHOGONAL TILT,OTHER,
                            PDB ENTRY,RANDOM CONICAL TILT}
                            """
                            set_cif_value(st_map.set_type_of_model, "type", const.EM_START_MODEL, cif_list=sm_in)

                        def set_el_choices(st_map, sm_in):
                            """
                            XSD: <xs:complexType name="starting_map_type"> element 1 has
                            .. 1  choice of 6 elements and
                            .. 1 element
                            """

                            def set_choice_random_conical_tilt(rct, sm_in):
                                """
                                XSD: <xs:element name="random_conical_tilt" minOccurs="0"> has
                                .. 4 elements
                                """

                                def set_el_number_images(rct, sm_in):
                                    """
                                    XSD: <xs:element name="number_images" type="xs:positiveInteger" minOccurs="0"/>
                                    CIF: _em_start_model.random_conical_tilt_num_images 40
                                    """
                                    set_cif_value(rct.set_number_images, "random_conical_tilt_num_images", const.EM_START_MODEL, cif_list=sm_in, fmt=int)

                                def set_el_tilt_angle(rct, sm_in):
                                    """
                                    XSD: <xs:element name="tilt_angle" minOccurs="0">
                                    CIF: _em_start_model.random_conical_tilt_angle 60
                                    """
                                    set_cif_value(
                                        rct.set_tilt_angle,
                                        "random_conical_tilt_angle",
                                        const.EM_START_MODEL,
                                        cif_list=sm_in,
                                        constructor=emdb.tilt_angleType,
                                        units=const.U_DEGF,
                                    )

                                def set_el_software_list():
                                    """
                                    XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                                    CIF: ???
                                    setSoftwareList(const.SOFT_?????,
                                        cat_soft_dict_in, rct.set_software_list)
                                    """

                                def set_el_details():
                                    """
                                    XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                    CIF: not in CIF yet
                                    """

                                # element 1
                                set_el_number_images(rct, sm_in)
                                # element 2
                                set_el_tilt_angle(rct, sm_in)
                                # element 3
                                set_el_software_list()
                                # element 4
                                set_el_details()

                            def set_choice_orthogonal_tilt(orth_tilt, sm_in):
                                """
                                XSD: <xs:element name="orthogonal_tilt" minOccurs="0"> has
                                .. 5 elements
                                """

                                def set_el_software_list():
                                    """
                                    XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                                    CIF: not in CIF yet
                                    """

                                def set_el_number_images(orth_tilt, sm_in):
                                    """
                                    XSD: <xs:element name="number_images" type="xs:positiveInteger"/ minOccurs="0">
                                    CIF: _em_start_model.orthogonal_tilt_num_images
                                    """
                                    set_cif_value(orth_tilt.set_number_images, "orthogonal_tilt_num_images", const.EM_START_MODEL, cif_list=sm_in, fmt=int)

                                def set_el_tilt_angle1(orth_tilt, sm_in):
                                    """
                                    XSD: <xs:element name="tilt_angle1" minOccurs="0">
                                    CIF: _em_start_model.orthogonal_tilt_angle1 (-180 to 180)
                                    """
                                    set_cif_value(
                                        orth_tilt.set_tilt_angle1,
                                        "orthogonal_tilt_angle1",
                                        const.EM_START_MODEL,
                                        cif_list=sm_in,
                                        constructor=emdb.tilt_angle1Type,
                                        units=const.U_DEGF,
                                    )

                                def set_el_tilt_angle2(orth_tilt, sm_in):
                                    """
                                    XSD: <xs:element name="tilt_angle2"  minOccurs="0">
                                    CIF: _em_start_model.orthogonal_tilt_angle2 (-180 to 180)
                                    """
                                    set_cif_value(
                                        orth_tilt.set_tilt_angle2,
                                        "orthogonal_tilt_angle2",
                                        const.EM_START_MODEL,
                                        cif_list=sm_in,
                                        constructor=emdb.tilt_angle2Type,
                                        units=const.U_DEGF,
                                    )

                                def set_el_details():
                                    """
                                    XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                    CIF: not in CIF
                                    """

                                # element 1
                                set_el_software_list()
                                # element 2
                                set_el_number_images(orth_tilt, sm_in)
                                # element 3
                                set_el_tilt_angle1(orth_tilt, sm_in)
                                # element 4
                                set_el_tilt_angle2(orth_tilt, sm_in)
                                # element 5
                                set_el_details()

                            def set_choice_emdb_id(st_map, sm_in):
                                """
                                XSD: <xs:element name="emdb_id" type="emdb_id_type" minOccurs="0"/>
                                CIF: _em_start_model.emdb_id
                                """
                                set_cif_value(st_map.set_emdb_id, "emdb_id", const.EM_START_MODEL, cif_list=sm_in)

                            def set_choice_pdb_model(st_map, sm_in):
                                """
                                XSD: <xs:element name="pdb_model" type="pdb_model_type" maxOccurs="unbounded" minOccurs="0">
                                """

                                def set_pdb_model_type(pdb_mt):
                                    """
                                    XSD: <xs:complexType name="pdb_model_type"> has
                                    .. 2 elements
                                    """

                                    def set_el_pdb_id(pdb_mt):
                                        """
                                        XSD: <xs:element name="pdb_id" type="pdb_code_type"/>
                                        CIF: _em_start_model.pdb_id
                                        """
                                        set_cif_value(pdb_mt.set_pdb_id, "pdb_id", const.EM_START_MODEL, cif_list=sm_in, parent_el_req=False)

                                    def set_el_chain_id_list():
                                        """
                                        XSD: <xs:element name="chain_id_list" type="chain_type" minOccurs="0"/>
                                        CIF:???
                                        """

                                    # element 1
                                    set_el_pdb_id(pdb_mt)
                                    # element 1
                                    set_el_chain_id_list()

                                el_pdb_id = get_cif_value("pdb_id", const.EM_START_MODEL, cif_list=sm_in)
                                if el_pdb_id is not None:
                                    pdb_mt = emdb.pdb_model_type()
                                    set_pdb_model_type(pdb_mt)
                                    if pdb_mt.has__content():
                                        st_map.set_pdb_model(pdb_mt)

                            def set_choice_insilico_model(st_map, sm_in):
                                """
                                XSD: <xs:element name="insilico_model" type="xs:token" minOccurs="0"/>
                                CIF: _em_start_model.insilico_model
                                """
                                set_cif_value(st_map.set_insilico_model, "insilico_model", const.EM_START_MODEL, cif_list=sm_in)

                            def set_choice_other(st_map, sm_in):
                                """
                                XSD: <xs:element name="other" type="xs:string" minOccurs="0"/>
                                CIF: _em_start_model.type
                                """
                                set_cif_value(st_map.set_other, "other", const.EM_START_MODEL, cif_list=sm_in)

                            t_of_m = get_cif_value("type", const.EM_START_MODEL, sm_in)
                            if t_of_m == "RANDOM CONICAL TILT":
                                rct = emdb.random_conical_tiltType()
                                set_choice_random_conical_tilt(rct, sm_in)
                                if rct.has__content():
                                    st_map.set_random_conical_tilt(rct)
                            elif t_of_m == "ORTHOGONAL TILT":
                                orth_tilt = emdb.orthogonal_tiltType()
                                set_choice_orthogonal_tilt(orth_tilt, sm_in)
                                if orth_tilt.has__content():
                                    st_map.set_orthogonal_tilt(orth_tilt)
                            elif t_of_m == "EMDB MAP":
                                set_choice_emdb_id(st_map, sm_in)
                            elif t_of_m == "PDB ENTRY":
                                set_choice_pdb_model(st_map, sm_in)
                            elif t_of_m == "INSILICO MODEL":
                                set_choice_insilico_model(st_map, sm_in)
                            elif t_of_m in ["OTHER" "NONE"]:  # pylint: disable=implicit-str-concat
                                set_choice_other(st_map, sm_in)

                        def set_el_details(st_map, sm_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _em_start_model.details
                            """
                            set_cif_value(st_map.set_details, "details", const.EM_START_MODEL, cif_list=sm_in)

                        # attribute 1
                        set_attr_type_of_model(st_map, sm_in)
                        # choice element
                        set_el_choices(st_map, sm_in)
                        # element 1
                        set_el_details(st_map, sm_in)

                    for sm_in in st_mod_dict_in[ip_id_in]:
                        st_map = emdb.starting_map_type()
                        set_starting_map_type(st_map, sm_in)
                        im_proc.add_startup_model(st_map)

                def set_final_reconstruction(final_rec, ip_id_in, final_dicts):
                    """
                    Final reconstruction for every method but non-subtomographs

                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param final_dicts: a dictionary of the following dictionaries:
                    final_rec_dict_in, cat_soft_dict_in,
                    p_sym_dict_in, h_sym_dict_in
                    XSD: <xs:complexType name="final_reconstruction_type"> has
                    .. 8 elements
                    """

                    def set_el_number_classes_used(final_rec, final_rec_in):
                        """
                        XSD: <xs:element name="number_classes_used" type="xs:positiveInteger" minOccurs="0"/>
                        CIF: _em_3d_reconstruction.num_class_averages 300
                        """
                        set_cif_value(final_rec.set_number_classes_used, "num_class_averages", const.EM_3D_RECONSTRUCTION, cif_list=final_rec_in, fmt=int)

                    def set_el_applied_symmetry(final_rec, final_dicts):
                        """
                        XSD: <xs:element name="applied_symmetry" type="applied_symmetry_type">
                        CIF: _em_3d_reconstruction.symmetry_type
                            {2D CRYSTAL,3D CRYSTAL,HELICAL,POINT}
                        """

                        def set_applied_symmetry_type(app_sym, final_dicts):
                            """
                            XSD: <xs:complexType name="applied_symmetry_type"> has
                            .. 1 element (3 choices)
                            """

                            def set_choice_space_group():
                                """
                                This is being dealt in def set_applied_symmetry_type(app_sym):
                                XSD: <xs:element name="space_group" type="xs:int"/>
                                """

                            def set_choice_point_group(app_sym, p_sym_in):
                                """
                                XSD: <xs:element name="point_group">
                                CIF: _em_single_particle_entity.point_symmetry C1
                                """
                                set_cif_value(app_sym.set_point_group, "point_symmetry", const.EM_SINGLE_PARTICLE_ENTITY, cif_list=p_sym_in)

                            def set_choice_helical_parameters(h_sym_in):
                                """
                                XSD: <xs:element name="helical_parameters" type="helical_parameters_type">
                                """

                                def set_helical_parameters_type(h_sym, h_sym_in):
                                    """
                                    XSD: <xs:complexType name="helical_parameters_type"> has
                                    .. 4 elements
                                    """

                                    def set_el_delta_z(h_sym, h_sym_in):
                                        """
                                        XSD: <xs:element name="delta_z">
                                        CIF: _em_helical_entity.axial_rise_per_subunit 17.400000
                                        """
                                        set_cif_value(
                                            h_sym.set_delta_z,
                                            "axial_rise_per_subunit",
                                            const.EM_HELICAL_ENTITY,
                                            cif_list=h_sym_in,
                                            constructor=emdb.delta_zType,
                                            fmt=float,
                                            units=const.U_ANG,
                                        )

                                    def set_el_delta_phi(h_sym):
                                        """
                                        XSD: <xs:element name="delta_phi">
                                        CIF: _em_helical_entity.angular_rotation_per_subunit -34.616000
                                        """
                                        set_cif_value(
                                            h_sym.set_delta_phi, "angular_rotation_per_subunit", const.EM_HELICAL_ENTITY, cif_list=h_sym_in, constructor=emdb.delta_phiType, units=const.U_DEG
                                        )

                                    def set_el_axial_symmetry(h_sym, h_sym_in):
                                        """
                                        XSD: <xs:element name="axial_symmetry">
                                        CIF: _em_helical_entity.axial_symmetry C1
                                        """
                                        set_cif_value(h_sym.set_axial_symmetry, "axial_symmetry", const.EM_HELICAL_ENTITY, cif_list=h_sym_in)

                                    # element 1
                                    set_el_delta_z(h_sym, h_sym_in)
                                    # element 2
                                    set_el_delta_phi(h_sym)
                                    # element 3
                                    set_el_axial_symmetry(h_sym, h_sym_in)
                                    # element 4 - REMOVED
                                    # set_el_hand(h_sym)

                                h_sym = emdb.helical_parameters_type()
                                set_helical_parameters_type(h_sym, h_sym_in)
                                if h_sym.has__content():
                                    app_sym.set_helical_parameters(h_sym)

                            p_sym_dict_in = final_dicts["p_sym_dict_in"]
                            h_sym_dict_in = final_dicts["h_sym_dict_in"]
                            if sym_type_in in ["2D CRYSTAL" "3D CRYSTAL"]:  # pylint: disable=implicit-str-concat
                                # choice 1
                                set_choice_space_group()
                            elif sym_type_in == "POINT" and ip_id_in in p_sym_dict_in:
                                p_sym_in = p_sym_dict_in[ip_id_in]
                                # choice 2
                                set_choice_point_group(app_sym, p_sym_in)
                            elif sym_type_in == "HELICAL" and ip_id_in in h_sym_dict_in:
                                h_sym_in = h_sym_dict_in[ip_id_in]
                                if h_sym_in is not None:
                                    # choice 3
                                    set_choice_helical_parameters(h_sym_in)

                        sym_type_in = get_cif_value("symmetry_type", const.EM_3D_RECONSTRUCTION, final_rec_in)
                        if sym_type_in is not None:
                            app_sym = emdb.applied_symmetry_type()
                            set_applied_symmetry_type(app_sym, final_dicts)
                            if app_sym.has__content():
                                final_rec.set_applied_symmetry(app_sym)

                    def set_el_algorithm(final_rec, final_rec_in):
                        """
                        XSD: <xs:element name="algorithm" type="reconstruction_algorithm_type"  minOccurs="0">
                        CIF: _em_3d_reconstruction.algorithm
                            {ALGEBRAIC (ARTS),BACK PROJECTION,EXACT BACK PROJECTION,
                            FOURIER SPACE,SIMULTANEOUS ITERATIVE (SIRT)}
                        """
                        set_cif_value(final_rec.set_algorithm, "algorithm", const.EM_3D_RECONSTRUCTION, cif_list=final_rec_in)

                    def set_el_resolution(final_rec, final_rec_in):
                        """
                        XSD: <xs:element name="resolution" minOccurs="0">
                        CIF: _em_3d_reconstruction.resolution 8.9
                        """
                        set_cif_value(
                            final_rec.set_resolution,
                            "resolution",
                            const.EM_3D_RECONSTRUCTION,
                            cif_list=final_rec_in,
                            constructor=emdb.resolutionType,
                            fmt=float,
                            units=const.U_ANG,
                            res_type="BY AUTHOR",
                        )

                    def set_el_resolution_method(final_rec, final_rec_in):
                        """
                        XSD: <xs:element name="resolution_method" minOccurs="0">
                        CIF: _em_3d_reconstruction.resolution_method
                        DIFFRACTION PATTERN/LAYERLINES
                        """
                        set_cif_value(final_rec.set_resolution_method, "resolution_method", const.EM_3D_RECONSTRUCTION, cif_list=final_rec_in)

                    def set_el_reconstruction_filtering():
                        """
                        XSD: <xs:element name="reconstruction_filtering" type="reconstruction_filtering_type" minOccurs="0">
                        CIF:???
                        XSD: <xs:complexType name="reconstruction_filtering_type"> has
                        .. 5 elements
                        XSD: <xs:complexType name="reconstruction_filtering_type"> element 1: <xs:element name="background_masked" type="background_masked_type" minOccurs="0"/>
                        XSD: <xs:complexType name="reconstruction_filtering_type"> element 2: <xs:element name="spatial_filtering" minOccurs="0">
                        XSD: <xs:complexType name="reconstruction_filtering_type"> element 3: <xs:element name="sharpening" minOccurs="0">
                        XSD: <xs:complexType name="reconstruction_filtering_type"> element 4: <xs:element name="b-factorSharpening" minOccurs="0">
                        XSD: <xs:complexType name="reconstruction_filtering_type"> element 5: <xs:element name="other" minOccurs="0">
                        """

                    def set_el_software_list(final_rec, final_dicts):
                        """
                        XSD: <xs:element name="software_list" type="software_list_type"/>
                        """
                        cat_soft_dict_in = final_dicts["cat_soft_dict_in"]
                        set_software_list(const.SOFT_RECONSTRUCTION, cat_soft_dict_in, final_rec.set_software_list)

                    def set_el_details(final_rec, final_rec_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _em_3d_reconstruction.details
                        """
                        set_cif_value(final_rec.set_details, "details", const.EM_3D_RECONSTRUCTION, cif_list=final_rec_in)

                    final_rec_dict_in = final_dicts["final_rec_dict_in"]
                    final_rec_in = final_rec_dict_in[ip_id_in]
                    # element 1
                    set_el_number_classes_used(final_rec, final_rec_in)
                    # element 2
                    set_el_applied_symmetry(final_rec, final_dicts)
                    # element 3
                    set_el_algorithm(final_rec, final_rec_in)
                    # element 4
                    set_el_resolution(final_rec, final_rec_in)
                    # element 5
                    set_el_resolution_method(final_rec, final_rec_in)
                    # element 6
                    set_el_reconstruction_filtering()
                    # element 7
                    set_el_software_list(final_rec, final_dicts)
                    # element 8
                    set_el_details(final_rec, final_rec_in)

                def set_angle_assignments(im_proc, ip_id_in, ang_dict_in, cat_soft_dict_in, em_method, parent_req):
                    """
                    Initial or final Euler angles assignment in image processing

                    Params:
                    @param im_proc: image processing object
                    @param ip_id_in: image processing id
                    @param ang_dict_in: the dictionary containing angle assignments
                    XSD: <xs:complexType name="angle_assignment_type">
                    """

                    def set_angle_assignment_type(im_proc, ang, ang_in, cat_soft_dict_in, em_method, parent_req):
                        """
                        XSD: <xs:complexType name="angle_assignment_type"> has
                        .. 4 elements
                        """

                        def set_el_type(ang, ang_in, parent_req):
                            """
                            XSD: <xs:element name="type">
                            CIF: _em_euler_angle_assignment.type {ANGULAR RECONSTITUTION,COMMON LINE, NOT APPLICABLE,OTHER,PROJECTION MATCHING,RANDOM ASSIGNMENT,MAXIMUM LIKELIHOOD}
                            """
                            angle_type = get_cif_value("type", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)
                            if angle_type is not None:
                                if angle_type == "" or angle_type.isspace():
                                    txt = u"Cif item (_em_euler_angle_assignment.type) is not set. The value given should be one of: ANGULAR RECONSTITUTION, COMMON LINE, NOT APPLICABLE, " \
                                          u"OTHER, PROJECTION MATCHING, RANDOM ASSIGNMENT, MAXIMUM LIKELIHOOD."
                                    self.current_entry_log.error_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                    )
                                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                            set_cif_value(ang.set_type, "type", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in, parent_el_req=parent_req)

                        def set_el_proj_match_processing(ang, ang_in):
                            """
                            XSD: <xs:element name="projection_matching_processing" minOccurs="0">
                            """

                            def set_projection_match_proc_type(prj):
                                """
                                XSD: <xs:element name="projection_matching_processing" minOccurs="0"> has
                                .. 3 elements
                                """

                                def set_el_num_ref_projections(prj, ang_in):
                                    """
                                    XSD: <xs:element name="number_reference_projections" type="xs:positiveInteger" minOccurs="0"/>
                                    CIF: _em_euler_angle_assignment.proj_matching_num_projections
                                    """
                                    set_cif_value(
                                        prj.set_number_reference_projections,
                                        "proj_matching_num_projections",
                                        const.EM_EULER_ANGLE_ASSIGNMENT,
                                        cif_list=ang_in,
                                        fmt=int,
                                    )

                                def set_el_merit_function(prj, ang_in):
                                    """
                                    XSD: <xs:element name="merit_function" type="xs:token"  minOccurs="0"/>
                                    CIF: _em_euler_angle_assignment.proj_matching_merit_function 'Correlation coeficient (CC)'
                                    """
                                    set_cif_value(prj.set_merit_function, "proj_matching_merit_function", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)

                                def set_el_angular_sampling(prj, ang_in):
                                    """
                                    XSD: <xs:element name="angular_sampling" minOccurs="0">
                                    CIF: _em_euler_angle_assignment.proj_matching_angular_sampling
                                    """
                                    set_cif_value(
                                        prj.set_angular_sampling,
                                        "proj_matching_angular_sampling",
                                        const.EM_EULER_ANGLE_ASSIGNMENT,
                                        cif_list=ang_in,
                                        constructor=emdb.angular_samplingType,
                                        fmt=float,
                                        units=const.U_DEGF,
                                    )

                                # element 1
                                set_el_num_ref_projections(prj, ang_in)
                                # element 2
                                set_el_merit_function(prj, ang_in)
                                # element 3
                                set_el_angular_sampling(prj, ang_in)

                            ang_type = get_cif_value("type", const.EM_EULER_ANGLE_ASSIGNMENT, ang_in)
                            if ang_type == "PROJECTION MATCHING":
                                prj = emdb.projection_matching_processingType()
                                set_projection_match_proc_type(prj)
                                ang.set_projection_matching_processing(prj)

                        def set_el_software_list(im_proc, ang, ang_in, cat_soft_dict_in, em_method):
                            """
                            XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                            CIF: _em_euler_angle_assignment.order {FINAL,INITIAL}
                            """
                            order = get_cif_value("order", const.EM_EULER_ANGLE_ASSIGNMENT, ang_in)
                            if order == "INITIAL" and em_method == const.EMM_SP:
                                set_software_list(const.SOFT_INITIAL_EULER_ASSIGNMENT, cat_soft_dict_in, ang.set_software_list)
                                if ang.has__content():
                                    im_proc.set_initial_angle_assignment(ang)
                            elif order == "FINAL":
                                set_software_list(const.SOFT_FINAL_EULER_ASSIGNMENT, cat_soft_dict_in, ang.set_software_list)
                                if ang.has__content():
                                    im_proc.set_final_angle_assignment(ang)

                        def set_el_details(ang, ang_in):
                            """
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            CIF: _em_euler_angle_assignment.details
                            """
                            set_cif_value(ang.set_details, "details", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)

                        # element 1
                        set_el_type(ang, ang_in, parent_req)
                        # element 2
                        set_el_proj_match_processing(ang, ang_in)
                        # element 3
                        set_el_software_list(im_proc, ang, ang_in, cat_soft_dict_in, em_method)
                        # element 4
                        set_el_details(ang, ang_in)

                    for ang_in in ang_dict_in[ip_id_in]:
                        el_type = get_cif_value("type", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)
                        el_num_ref_projections = get_cif_value("proj_matching_number_reference_projections", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)
                        el_merit_function = get_cif_value("proj_matching_merit_function", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)
                        el_angular_sampling = get_cif_value("proj_matching_angular_sampling", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)
                        _el_software_list = get_cif_value("order", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)  # noqa: F841
                        el_details = get_cif_value("details", const.EM_EULER_ANGLE_ASSIGNMENT, cif_list=ang_in)
                        angle_assignment_type_list = [el_type, el_num_ref_projections, el_merit_function, el_angular_sampling, el_details]
                        if any(x is not None for x in angle_assignment_type_list):
                            ang = emdb.angle_assignment_type()
                            set_angle_assignment_type(im_proc, ang, ang_in, cat_soft_dict_in, em_method, parent_req)

                def set_crystal_parameters(ip_id_in, im_proc, cryst_dict_in, dict_category, parent_req):
                    """
                    Parameters for crystalography image processing

                    @param ip_id_in: image processing id in cif file
                    @param im_proc: image processing object
                    @param cryst_dict_in: the dictionary with crystalography parameters
                    @param dict_category: category in cifused for the dictionary
                    """

                    def set_crystal_parameters_type(cryst, cryst_in, dict_category, parent_req):
                        """
                        XSD: <xs:complexType name="crystal_parameters_type"> has
                        .. 1 element and
                        .. 1 choice of 2 elements
                        """

                        def set_el_unit_cell(cryst, cryst_in, dict_category, parent_req):
                            """
                            XSD: <xs:element name="unit_cell" type="unit_cell_type"/>
                            """

                            def set_unit_cell_type(u_cell, cryst_in, dict_category, parent_req):
                                """
                                XSD: <xs:complexType name="unit_cell_type"> has
                                .. 7 elements
                                """

                                def set_el_a(u_cell, cryst_in, dict_category):
                                    """
                                    XSD: <xs:element name="a" type="cell_type">
                                    CIF: _em_2d_crystal_entity.length_a 62.3
                                    CIF: _em_3d_crystal_entity.length_a 62.4
                                    """
                                    set_cif_value(u_cell.set_a, "length_a", dict_category, cif_list=cryst_in, constructor=emdb.cell_type, units=const.U_ANG, parent_el_req=parent_req)

                                def set_el_b(u_cell, cryst_in, dict_category):
                                    """
                                    XSD: <xs:element name="b" type="cell_type">
                                    CIF: _em_2d_crystal_entity.length_b 62.3
                                    CIF: _em_3d_crystal_entity.length_b 62.4
                                    """
                                    set_cif_value(u_cell.set_b, "length_b", dict_category, cif_list=cryst_in, constructor=emdb.cell_type, units=const.U_ANG, parent_el_req=parent_req)

                                def set_el_c(u_cell, cryst_in, dict_category):
                                    """
                                    XSD: <xs:element name="c" type="cell_type">
                                    CIF: _em_2d_crystal_entity.length_c 62.4
                                    CIF: _em_3d_crystal_entity.length_c
                                    """
                                    set_cif_value(u_cell.set_c, "length_c", dict_category, cif_list=cryst_in, constructor=emdb.cell_type, units=const.U_ANG, parent_el_req=parent_req)

                                def set_el_c_sampling_length(u_cell, cryst_in, dict_category):
                                    """
                                    XSD: <xs:element name="c_sampling_length" type="cell_type">
                                    CIF: _em_2d_crystal_entity.c_sampling_length
                                    CIF: _em_3d_crystal_entity.c_sampling_length
                                    """
                                    set_cif_value(
                                        u_cell.set_c_sampling_length, "c_sampling_length", dict_category, cif_list=cryst_in, constructor=emdb.cell_type, units=const.U_ANG
                                    )

                                def set_el_gamma(u_cell, cryst_in, dict_category, parent_req):
                                    """
                                    XSD: <xs:element name="gamma" type="cell_angle_type">
                                    CIF: _em_2d_crystal_entity.angle_gamma 120.0
                                    CIF: _em_3d_crystal_entity.angle_gamma 120.0
                                    """
                                    set_cif_value(
                                        u_cell.set_gamma, "angle_gamma", dict_category, cif_list=cryst_in, constructor=emdb.cell_angle_type, units=const.U_DEG, parent_el_req=parent_req
                                    )

                                def set_el_alpha(u_cell, cryst_in, dict_category):
                                    """
                                    XSD: <xs:element name="alpha" type="cell_angle_type">
                                    CIF: _em_3d_crystal_entity.angle_alpha 120.0
                                    """
                                    set_cif_value(u_cell.set_alpha, "alpha", dict_category, cif_list=cryst_in, constructor=emdb.cell_angle_type, units=const.U_DEG)

                                def set_el_beta(u_cell, cryst_in, dict_category):
                                    """
                                    XSD: <xs:element name="beta" type="cell_angle_type">
                                    CIF: _em_3d_crystal_entity.angle_beta 120.0
                                    """
                                    set_cif_value(u_cell.set_beta, "beta", dict_category, cif_list=cryst_in, constructor=emdb.cell_angle_type, units=const.U_DEG)

                                # element 1
                                set_el_a(u_cell, cryst_in, dict_category)
                                # element 2
                                set_el_b(u_cell, cryst_in, dict_category)
                                # element 3
                                set_el_c(u_cell, cryst_in, dict_category)
                                # element 4
                                set_el_c_sampling_length(u_cell, cryst_in, dict_category)
                                # element 5
                                set_el_gamma(u_cell, cryst_in, dict_category, parent_req)
                                # element 6
                                set_el_alpha(u_cell, cryst_in, dict_category)
                                # element 7
                                set_el_beta(u_cell, cryst_in, dict_category)

                            u_cell = emdb.unit_cell_type()
                            set_unit_cell_type(u_cell, cryst_in, dict_category, parent_req)
                            if u_cell.has__content():
                                cryst.set_unit_cell(u_cell)

                        def set_choice_plane_group():
                            """
                            XSD: <xs:element name="plane_group">
                            CIF:????
                            """

                        def set_choice_space_group(cryst, cryst_in, dict_category):
                            """
                            XSD: <xs:element name="space_group" type="xs:token">
                            CIF: _em_3d_crystal_entity.space_group_name  P 1
                            CIF: _em_2d_crystal_entity.space_group_name_H-M
                            """
                            if dict_category == const.EM_3D_CRYSTAL_ENTITY:
                                set_cif_value(cryst.set_space_group, "space_group_name", const.EM_3D_CRYSTAL_ENTITY, cif_list=cryst_in)
                            elif dict_category == const.EM_2D_CRYSTAL_ENTITY:
                                set_cif_value(cryst.set_plane_group, "space_group_name_H-M", const.EM_2D_CRYSTAL_ENTITY, cif_list=cryst_in)

                        # element 1
                        set_el_unit_cell(cryst, cryst_in, dict_category, parent_req)
                        # choice 1
                        set_choice_plane_group()
                        # choice 2
                        set_choice_space_group(cryst, cryst_in, dict_category)

                    cryst_in = cryst_dict_in[ip_id_in]
                    el_a = get_cif_value("a", dict_category, cif_list=cryst_in)
                    el_b = get_cif_value("b", dict_category, cif_list=cryst_in)
                    el_c = get_cif_value("c", dict_category, cif_list=cryst_in)
                    el_c_sampling_length = get_cif_value("c_sampling_length", dict_category, cif_list=cryst_in)
                    el_gamma = get_cif_value("gamma", dict_category, cif_list=cryst_in)
                    el_alpha = get_cif_value("alpha", dict_category, cif_list=cryst_in)
                    el_beta = get_cif_value("beta", dict_category, cif_list=cryst_in)
                    el_choice_space_group = None
                    if dict_category == const.EM_3D_CRYSTAL_ENTITY:
                        el_choice_space_group = get_cif_value("space_group", const.EM_3D_CRYSTAL_ENTITY, cif_list=cryst_in)
                    elif dict_category == const.EM_2D_CRYSTAL_ENTITY:
                        el_choice_space_group = get_cif_value("plane_group", const.EM_2D_CRYSTAL_ENTITY, cif_list=cryst_in)
                    crystal_parameters_type_list = [el_a, el_b, el_c, el_c_sampling_length, el_gamma, el_alpha, el_beta, el_choice_space_group]
                    if any(x is not None for x in crystal_parameters_type_list):
                        cryst = emdb.crystal_parameters_type()
                        set_crystal_parameters_type(cryst, cryst_in, dict_category, parent_req)
                        if cryst.has__content():
                            im_proc.set_crystal_parameters(cryst)

                def set_non_sub_tom_final_recon(ip_id_in, im_proc, non_st_dicts):
                    """
                    Final reconstruction for non-subtomographs

                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param non_st_dicts: a dictionary of the following dictionaries:
                    final_rec_dict_in, cat_soft_dict_in,
                    p_sym_dict_in, h_sym_dict_in
                    XSD: <xs:complexType name="non_subtom_final_reconstruction_type">
                    """

                    def set_non_subtom_final_rec_type(spfr):
                        """
                        XSD: <xs:complexType name="single_particle_final_reconstruction_type"> has
                        .. a base (reconstruction_type) and
                        .. 1 element
                        """

                        def set_el_number_images_used(spfr, f_rec_in):
                            """
                            XSD: <xs:element name="num_particles" type="xs:positiveInteger" minOccurs="0">
                            CIF: _em_3d_reconstruction.num_particles 300
                            """
                            set_cif_value(spfr.set_number_images_used, "num_particles", const.EM_3D_RECONSTRUCTION, cif_list=f_rec_in, fmt=int)

                        # base
                        final_rec_dicts = {
                            "final_rec_dict_in": non_st_dicts["final_rec_dict_in"],
                            "cat_soft_dict_in": non_st_dicts["cat_soft_dict_in"],
                            "p_sym_dict_in": non_st_dicts["p_sym_dict_in"],
                            "h_sym_dict_in": non_st_dicts["h_sym_dict_in"],
                        }
                        set_final_reconstruction(spfr, ip_id_in, final_rec_dicts)
                        # element 1
                        set_el_number_images_used(spfr, f_rec_in)  # Will never happen but should be passed to function.  #  pylint: disable=possibly-used-before-assignment

                    final_rec_dict_in = non_st_dicts["final_rec_dict_in"]
                    if ip_id_in in final_rec_dict_in:
                        f_rec_in = final_rec_dict_in[ip_id_in]  # used in set_non_subtom_final_reconstruction_type
                        spfr = emdb.non_subtom_final_reconstruction_type()
                        set_non_subtom_final_rec_type(spfr)
                        im_proc.set_final_reconstruction(spfr)

                def set_single_part_im_proc_specs(ip_id_in, im_proc, sp_dict_list):
                    """
                    Set elements specific for single particle image processing add group

                    Parameters:
                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param sp_dict_list: a dictionary of the following dictionaries:
                    part_sel_dict_in, cat_soft_dict_in, ctf_corr_dict_in, st_mod_dict_in,
                    final_rec_dict_in, p_sym_dict_in, h_sym_dict_in, ang_dict_in,
                    final_2d_class_dict_in, final_class_dict_in
                    XSD: <xs:group name="single_particle_proc_add_group"> has
                    .. 9 elements
                    """

                    def set_el_particle_selection(ip_id_in, im_proc, sp_dict_list):
                        """
                        XSD: <xs:element name="particle_selection" type="particle_selection_type" maxOccurs="unbounded"/>
                        """
                        part_sel_dict_in = sp_dict_list["part_sel_dict_in"]
                        if ip_id_in in part_sel_dict_in:
                            ps_list_in = part_sel_dict_in[ip_id_in]
                            for ps_in in ps_list_in:
                                set_particle_selection(ps_in, im_proc, part_sel_dict_in)

                    def set_el_ctf_correction(ip_id_in, im_proc, sp_dict_list):
                        """
                        XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                        """
                        if ip_id_in in sp_dict_list["ctf_corr_dict_in"]:
                            set_ctfcorrection(ip_id_in, im_proc, sp_dict_list["ctf_corr_dict_in"])

                    def set_el_startup_model(ip_id_in, im_proc, sp_dict_list):
                        """
                        XSD: <xs:element name="startup_model" type="starting_map_type" minOccurs="0" >
                        """
                        if ip_id_in in sp_dict_list["st_mod_dict_in"]:
                            set_startup_model(ip_id_in, im_proc, sp_dict_list["st_mod_dict_in"])

                    def set_el_final_reconstruction(ip_id_in, im_proc, sp_dict_list):
                        """
                        XSD: <xs:element name="final_reconstruction" type="single_particle_final_reconstruction_type"/>
                        """
                        non_st_dicts = {
                            "final_rec_dict_in": sp_dict_list["final_rec_dict_in"],
                            "cat_soft_dict_in": sp_dict_list["cat_soft_dict_in"],
                            "p_sym_dict_in": sp_dict_list["p_sym_dict_in"],
                            "h_sym_dict_in": sp_dict_list["h_sym_dict_in"],
                        }
                        set_non_sub_tom_final_recon(ip_id_in, im_proc, non_st_dicts)

                    def set_el_initial_angle_assignment():
                        """
                        XSD: <xs:element name="initial_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                        """

                    def set_el_final_angle_assignment(ip_id_in, im_proc, sp_dict_list):
                        """
                        XSD: <xs:element name="final_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                        """
                        if ip_id_in in sp_dict_list["ang_dict_in"]:
                            set_angle_assignments(im_proc, ip_id_in, sp_dict_list["ang_dict_in"], sp_dict_list["cat_soft_dict_in"], em_method, parent_req=False)

                    def set_el_final_multi_ref_align():
                        """
                        XSD: <xs:element name="final_multi_reference_alignment" minOccurs="0">
                        """

                    def set_el_final_two_d_class(ip_id_in, im_proc, sp_dict_list):
                        """
                        XSD: <xs:element name="final_two_d_classification" type="classification_type" minOccurs="0"/>
                        """
                        final_2d_class_dict_in = sp_dict_list["final_2d_class_dict_in"]
                        cat_soft_dict_in = sp_dict_list["cat_soft_dict_in"]
                        if ip_id_in in final_2d_class_dict_in:
                            f2dc = set_classification(ip_id_in, final_2d_class_dict_in, cat_soft_dict_in)
                            if f2dc.has__content():
                                im_proc.set_final_two_d_classification(f2dc)

                    def set_el_final_three_d_class(ip_id_in, im_proc, sp_dict_list):
                        """
                        XSD: <xs:element name="final_three_d_classification" type="classification_type" minOccurs="0"/>
                        """
                        final_class_dict_in = sp_dict_list["final_class_dict_in"]
                        cat_soft_dict_in = sp_dict_list["cat_soft_dict_in"]
                        if ip_id_in in final_class_dict_in:
                            f3dc = set_classification(ip_id_in, final_class_dict_in, cat_soft_dict_in)
                            if f3dc.has__content():
                                im_proc.set_final_three_d_classification(f3dc)

                    # element 1
                    set_el_particle_selection(ip_id_in, im_proc, sp_dict_list)
                    # element 2
                    set_el_ctf_correction(ip_id_in, im_proc, sp_dict_list)
                    # element 3
                    set_el_startup_model(ip_id_in, im_proc, sp_dict_list)
                    # element 4
                    set_el_final_reconstruction(ip_id_in, im_proc, sp_dict_list)
                    # element 5
                    set_el_initial_angle_assignment()
                    # element 6
                    set_el_final_angle_assignment(ip_id_in, im_proc, sp_dict_list)
                    # element 7
                    set_el_final_multi_ref_align()
                    # element 8
                    set_el_final_two_d_class(ip_id_in, im_proc, sp_dict_list)
                    # element 9
                    set_el_final_three_d_class(ip_id_in, im_proc, sp_dict_list)

                def set_subtom_av_im_proc_specifics(ip_id_in, im_proc, subtom_dicts, em_method):
                    """
                    Set elements specific for tomography image processing add group

                    Parameters:
                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param subtom_dicts: a dictionary of the following dictionaries:
                    final_rec_dict_in, cat_soft_dict_in, p_sym_dict_in, h_sym_dict_in,
                    vol_sel_dict_in, ctf_corr_dict_in, final_class_dict_in, ang_dict_in,
                    two_d_cryst_dict_in, three_d_cryst_dict_in
                    @param em_method:
                    XSD: <xs:group name="subtomogram_averaging_proc_add_group"> has
                    .. 7 elements
                    """

                    def set_el_final_reconstruction(im_proc, subtom_dicts):
                        """
                        XSD: <xs:element name="final_reconstruction" type="subtomogram_final_reconstruction_type" minOccurs="0"/>
                        """

                        def set_subtom_final_rec_type(stfr, f_rec_in, subtom_dicts):
                            """
                            XSD: <xs:complexType name="subtomogram_final_reconstruction_type"> has
                            .. a base (final_reconstruction_type) and
                            .. 1 element
                            """

                            def set_el_number_subtomograms_used(stfr, f_rec_in):
                                """
                                XSD: <xs:element name="number_subtomograms_used" type="xs:positiveInteger" minOccurs="0">
                                CIF: _em_3d_reconstruction.num_particles 300
                                """
                                set_cif_value(stfr.set_number_subtomograms_used, "num_particles", const.EM_3D_RECONSTRUCTION, cif_list=f_rec_in, fmt=int)

                            # base
                            final_rec_dicts = {
                                "final_rec_dict_in": subtom_dicts["final_rec_dict_in"],
                                "cat_soft_dict_in": subtom_dicts["cat_soft_dict_in"],
                                "p_sym_dict_in": subtom_dicts["p_sym_dict_in"],
                                "h_sym_dict_in": subtom_dicts["h_sym_dict_in"],
                            }
                            set_final_reconstruction(stfr, ip_id_in, final_rec_dicts)
                            # element 1
                            set_el_number_subtomograms_used(stfr, f_rec_in)

                        final_rec_dict_in = subtom_dicts["final_rec_dict_in"]
                        if ip_id_in in final_rec_dict_in:
                            f_rec_in = final_rec_dict_in[ip_id_in]
                            stfr = emdb.subtomogram_final_reconstruction_type()
                            set_subtom_final_rec_type(stfr, f_rec_in, subtom_dicts)
                            im_proc.set_final_reconstruction(stfr)

                    def set_el_extraction(im_proc, ip_id_in, subtom_dicts):
                        """
                        XSD:  <xs:element name="extraction">
                        """

                        def set_extraction_type(extraction, vs_in, subtom_dicts):
                            """
                            XSD:  <xs:element name="extraction"> has
                            .. 6 elements
                            """

                            def set_el_number_tomograms(extraction, vs_in):
                                """
                                XSD: <xs:element name="number_tomograms" type="xs:positiveInteger"/>
                                CIF: _em_volume_selection.num_tomograms 20
                                """
                                set_cif_value(extraction.set_number_tomograms, "num_tomograms", const.EM_VOLUME_SELECTION, cif_list=vs_in, fmt=int)

                            def set_el_number_images_used(extraction, vs_in):
                                """
                                XSD: <xs:element name="number_images_used" type="xs:positiveInteger"/>
                                CIF: _em_volume_selection.num_volumes_extracted 840
                                """
                                set_cif_value(extraction.set_number_images_used, "num_volumes_extracted", const.EM_VOLUME_SELECTION, cif_list=vs_in, fmt=int)

                            def set_el_reference_model(extraction, vs_in):
                                """
                                XSD: <xs:element name="reference_model" type="xs:token" minOccurs="0">
                                CIF: _em_volume_selection.reference_model
                                """
                                set_cif_value(extraction.set_reference_model, "reference_model", const.EM_VOLUME_SELECTION, cif_list=vs_in)

                            def set_el_method(extraction, vs_in):
                                """
                                XSD: <xs:element name="method" type="xs:string" minOccurs="0">
                                CIF: _em_volume_selection.method 'volumes picked interactively'
                                """
                                set_cif_value(extraction.set_method, "method", const.EM_VOLUME_SELECTION, cif_list=vs_in)

                            def set_el_software_list(extraction, subtom_dicts):
                                """
                                XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                                """
                                cat_soft_dict_in = subtom_dicts["cat_soft_dict_in"]
                                set_software_list("VOLUME SELECTION", cat_soft_dict_in, extraction.set_software_list)

                            def set_el_details(extraction, vs_in):
                                """
                                XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                CIF: _em_volume_selection.details
                                """
                                set_cif_value(extraction.set_details, "details", const.EM_VOLUME_SELECTION, cif_list=vs_in)

                            # element 1
                            set_el_number_tomograms(extraction, vs_in)
                            # element 2
                            set_el_number_images_used(extraction, vs_in)
                            # element 3
                            set_el_reference_model(extraction, vs_in)
                            # element 4
                            set_el_method(extraction, vs_in)
                            # element 5
                            set_el_software_list(extraction, subtom_dicts)
                            # element 6
                            set_el_details(extraction, vs_in)

                        vol_sel_dict_in = subtom_dicts["vol_sel_dict_in"]
                        if ip_id_in in vol_sel_dict_in:
                            vs_in = vol_sel_dict_in[ip_id_in]
                            extraction = emdb.extractionType()
                            set_extraction_type(extraction, vs_in, subtom_dicts)
                            if extraction.has__content():
                                im_proc.set_extraction(extraction)

                    def set_el_ctf_correction(im_proc, ip_id_in, subtom_dicts):
                        """
                        XSD: <xs:element name="ctf_correction" type="ctf_correction_type"/>
                        """
                        ctf_corr_dict_in = subtom_dicts["ctf_corr_dict_in"]
                        if ip_id_in in ctf_corr_dict_in:
                            set_ctfcorrection(ip_id_in, im_proc, ctf_corr_dict_in)

                    def set_el_final_multi_ref_align():
                        """
                        XSD: <xs:element name="final_multi_reference_alignment" minOccurs="0">
                        CIF:??
                        """

                    def set_el_final_three_d_class(im_proc, ip_id_in, subtom_dicts):
                        """
                        XSD: <xs:element name="final_three_d_classification" type="classification_type" minOccurs="0"/>
                        """
                        final_class_dict_in = subtom_dicts["final_class_dict_in"]
                        cat_soft_dict_in = subtom_dicts["cat_soft_dict_in"]
                        if ip_id_in in final_class_dict_in:
                            f3dc = set_classification(ip_id_in, final_class_dict_in, cat_soft_dict_in)
                            if f3dc.has__content():
                                im_proc.set_final_three_d_classification(f3dc)

                    def set_el_final_angle_assignment(im_proc, subtom_dicts, em_method):
                        """
                        XSD: <xs:element name="final_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                        """
                        ang_dict_in = subtom_dicts["ang_dict_in"]
                        cat_soft_dict_in = subtom_dicts["cat_soft_dict_in"]
                        if ip_id_in in ang_dict_in:
                            set_angle_assignments(im_proc, ip_id_in, ang_dict_in, cat_soft_dict_in, em_method, parent_req=False)

                    def set_el_crystal_parameters(im_proc, subtom_dicts):
                        """
                        XSD: <xs:element name="crystal_parameters" type="crystal_parameters_type" minOccurs="0"/>
                        """
                        if ip_id_in in subtom_dicts["two_d_cryst_dict_in"]:
                            set_crystal_parameters(ip_id_in, im_proc, subtom_dicts["two_d_cryst_dict_in"], const.EM_2D_CRYSTAL_ENTITY, parent_req=False)
                        if ip_id_in in subtom_dicts["three_d_cryst_dict_in"]:
                            set_crystal_parameters(ip_id_in, im_proc, subtom_dicts["three_d_cryst_dict_in"], const.EM_3D_CRYSTAL_ENTITY, parent_req=False)

                    # element 1
                    set_el_final_reconstruction(im_proc, subtom_dicts)
                    # element 2
                    set_el_extraction(im_proc, ip_id_in, subtom_dicts)
                    # element 3
                    set_el_ctf_correction(im_proc, ip_id_in, subtom_dicts)
                    # element 4
                    set_el_final_multi_ref_align()
                    # element 5
                    set_el_final_three_d_class(im_proc, ip_id_in, subtom_dicts)
                    # element 6
                    set_el_final_angle_assignment(im_proc, subtom_dicts, em_method)
                    # element 7
                    set_el_crystal_parameters(im_proc, subtom_dicts)

                def set_helical_proc_specifics(ip_id_in, im_proc, hel_dict_list):
                    """
                    Set elements specific for helical image processing add group

                    Parameters:
                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param hel_dict_list: a dictionary of the following dictionaries:
                    final_rec_dict_in, part_sel_dict_in, cat_soft_dict_in,
                    p_sym_dict_in, h_sym_dict_in, ctf_corr_dict_in, ang_dict_in
                    st_mod_dict_in, two_d_cryst_dict_in, three_d_cryst_dict_in
                    XSD: <xs:group name="helical_processing_add_group"> has
                    .. 9 elements
                    """

                    def set_el_final_reconstruction(im_proc, ip_id_in, hel_dict_list):
                        """
                        XSD: <xs:element name="final_reconstruction" type="non_subtom_final_reconstruction_type" minOccurs="0"/>
                        """

                        def set_non_subtom_final_rec_type(nstfr, f_rec_in):
                            """
                            XSD: <xs:complexType name="non_subtom_final_reconstruction_type"> has
                            .. a base (reconstruction_type) and
                            .. 1 element
                            """

                            def set_el_number_images_used(nstfr, f_rec_in):
                                """
                                XSD: <xs:element name="number_images_used" type="xs:positiveInteger" minOccurs="0">
                                CIF: _em_3d_reconstruction.num_particles 300
                                """
                                set_cif_value(nstfr.set_number_images_used, "num_particles", const.EM_3D_RECONSTRUCTION, cif_list=f_rec_in, fmt=int)

                            # base
                            final_rec_dicts = {
                                "final_rec_dict_in": hel_dict_list["final_rec_dict_in"],
                                "cat_soft_dict_in": hel_dict_list["cat_soft_dict_in"],
                                "p_sym_dict_in": hel_dict_list["p_sym_dict_in"],
                                "h_sym_dict_in": hel_dict_list["h_sym_dict_in"],
                            }
                            set_final_reconstruction(nstfr, ip_id_in, final_rec_dicts)
                            # element 1
                            set_el_number_images_used(nstfr, f_rec_in)

                        final_rec_dict_in = hel_dict_list["final_rec_dict_in"]
                        if ip_id_in in final_rec_dict_in:
                            f_rec_in = final_rec_dict_in[ip_id_in]
                            nstfr = emdb.non_subtom_final_reconstruction_type()
                            set_non_subtom_final_rec_type(nstfr, f_rec_in)
                            im_proc.set_final_reconstruction(nstfr)

                    def set_el_ctf_correction(im_proc, ip_id_in, hel_dict_list):
                        """
                        XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                        """
                        ctf_corr_dict_in = hel_dict_list["ctf_corr_dict_in"]
                        if ip_id_in in ctf_corr_dict_in:
                            set_ctfcorrection(ip_id_in, im_proc, ctf_corr_dict_in)

                    def set_el_segment_selection(im_proc, hel_dict_list):
                        """
                        XSD: <xs:element name="segment_selection" type="segment_selection_type" maxOccurs="unbounded" minOccurs="0">
                        """

                        def set_segment_selection_type(seg_sel, ps_in):
                            """
                            XSD: <xs:complexType name="segment_selection_type"> has
                            .. 6 elements
                            """

                            def set_el_number_segments(seg_sel, ps_in):
                                """
                                XSD: <xs:element name="number_selected" type="xs:positiveInteger" minOccurs="0"/>
                                CIF: _em_particle_selection.num_particles_selected
                                """
                                set_cif_value(seg_sel.set_number_selected, "num_particles_selected", const.EM_PARTICLE_SELECTION, cif_list=ps_in, fmt=int)

                            def set_el_segment_length():
                                """
                                XSD: <xs:element name="segment_length"  minOccurs="0">
                                CIF: ??
                                """

                            def set_el_segment_overlap():
                                """
                                XSD: <xs:element name="segment_overlap" minOccurs="0">
                                CIF: ??
                                """

                            def set_el_total_filament_length():
                                """
                                XSD: <xs:element name="total_filament_length" minOccurs="0">
                                CIF: ??
                                """

                            def set_el_software_list(seg_sel, hel_dict_list):
                                """
                                XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                                """
                                set_software_list(const.SOFT_PARTICLE_SELECTION, hel_dict_list["cat_soft_dict_in"], seg_sel.set_software_list)

                            def set_el_details(seg_sel, ps_in):
                                """
                                XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                CIF: _em_particle_selection.details
                                """
                                set_cif_value(seg_sel.set_details, "details", const.EM_PARTICLE_SELECTION, cif_list=ps_in)

                            # element 1
                            set_el_number_segments(seg_sel, ps_in)
                            # element 2
                            set_el_segment_length()
                            # element 3
                            set_el_segment_overlap()
                            # element 4
                            set_el_total_filament_length()
                            # element 5
                            set_el_software_list(seg_sel, hel_dict_list)
                            # element 6
                            set_el_details(seg_sel, ps_in)

                        part_sel_dict_in = hel_dict_list["part_sel_dict_in"]
                        if ip_id_in in part_sel_dict_in:
                            ps_list_in = part_sel_dict_in[ip_id_in]
                            for ps_in in ps_list_in:
                                seg_sel = emdb.segment_selection_type()
                                set_segment_selection_type(seg_sel, ps_in)
                                if seg_sel.has__content():
                                    im_proc.add_segment_selection(seg_sel)

                    def set_el_refinement(im_proc):
                        """
                        XSD: <xs:element name="refinement" type="refinement_type"/>
                        """

                        def set_refinement_type():
                            """
                            XSD: <xs:complexType name="refinement_type"> has
                            .. 4 elements

                            XSD: <xs:element name="startup_model" type="starting_map_type" maxOccurs="unbounded" minOccurs="0">
                            XSD: <xs:element name="starting_symmetry" maxOccurs="unbounded">
                            XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                            XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                            """

                        refinement = emdb.refinement_type()
                        set_refinement_type()
                        if refinement.has__content():
                            im_proc.set_refinement(refinement)

                    def set_el_startup_model(im_proc, ip_id_in, hel_dict_list):
                        """
                        XSD: <xs:element name="startup_model" type="starting_map_type" maxOccurs="unbounded" minOccurs="0">
                        """
                        st_mod_dict_in = hel_dict_list["st_mod_dict_in"]
                        if ip_id_in in st_mod_dict_in:
                            set_startup_model(ip_id_in, im_proc, st_mod_dict_in)

                    def set_el_helical_layer_lines(im_proc):
                        """
                        XSD: <xs:element name="helical_layer_lines" type="layer_lines_type" minOccurs="0"/>
                        XSD: <xs:complexType name="layer_lines_type"> has
                        .. 4 elements
                        """

                        def set_layer_lines_type():
                            """
                            XSD: <xs:element name="number_helices"/>
                            XSD: <xs:element name="helix_length">
                            XSD: <xs:element name="straightening" type="xs:boolean" minOccurs="0"/>
                            XSD: <xs:element name="indexing">
                            """

                        lay_lines = emdb.layer_lines_type()
                        set_layer_lines_type()
                        if lay_lines.has__content():
                            im_proc.set_helical_layer_lines(lay_lines)

                    def set_el_initial_angle_assignment():
                        """
                        XSD: <xs:element name="initial_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                        """

                    def set_el_final_angle_assignment(im_proc, ip_id_in, hel_dict_list, em_method):
                        """
                        XSD: <xs:element name="final_angle_assignment" type="angle_assignment_type" minOccurs="0"/>
                        """
                        ang_dict_in = hel_dict_list["ang_dict_in"]
                        cat_soft_dict_in = hel_dict_list["cat_soft_dict_in"]
                        if ip_id_in in ang_dict_in:
                            set_angle_assignments(im_proc, ip_id_in, ang_dict_in, cat_soft_dict_in, em_method, parent_req=False)

                    def set_el_crystal_parameters(ip_id_in, im_proc, hel_dict_list):
                        """
                        XSD: <xs:element name="crystal_parameters" type="crystal_parameters_type" minOccurs="0"/>
                        """
                        two_d_cryst_dict_in = hel_dict_list["two_d_cryst_dict_in"]
                        if ip_id_in in two_d_cryst_dict_in:
                            set_crystal_parameters(ip_id_in, im_proc, two_d_cryst_dict_in, const.EM_2D_CRYSTAL_ENTITY, parent_req=False)
                        three_d_cryst_dict_in = hel_dict_list["three_d_cryst_dict_in"]
                        if ip_id_in in three_d_cryst_dict_in:
                            set_crystal_parameters(ip_id_in, im_proc, three_d_cryst_dict_in, const.EM_3D_CRYSTAL_ENTITY, parent_req=False)

                    # element 1
                    set_el_final_reconstruction(im_proc, ip_id_in, hel_dict_list)
                    # element 2
                    set_el_ctf_correction(im_proc, ip_id_in, hel_dict_list)
                    # element 3
                    set_el_segment_selection(im_proc, hel_dict_list)
                    # element 4
                    set_el_refinement(im_proc)
                    # element 5
                    set_el_startup_model(im_proc, ip_id_in, hel_dict_list)
                    # element 6
                    set_el_helical_layer_lines(im_proc)
                    # element 7
                    set_el_initial_angle_assignment()
                    # element 8
                    set_el_final_angle_assignment(im_proc, ip_id_in, hel_dict_list, em_method)
                    # element 9
                    set_el_crystal_parameters(ip_id_in, im_proc, hel_dict_list)

                def set_cryst_proc_specifics(ip_id_in, im_proc, cryst_dicts):
                    """
                    Set elements specific for crystallography image processing add group

                    Parameters:
                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param cryst_dicts: a dictionary of the following dictionaries:
                    final_rec_dict_in, cat_soft_dict_in, p_sym_dict_in, h_sym_dict_in,
                    two_d_cryst_dict_in, three_d_cryst_dict_in, cry_shell_dict_in
                    st_mod_dict_in, ctf_corr_dict_in, cry_stats_dict_in
                    XSD: <xs:group name="crystallography_proc_add_group"> has
                    .. 9 elements
                    """

                    def set_el_final_reconstruction(im_proc, ip_id_in):
                        """
                        XSD: <xs:element name="final_reconstruction" type="non_subtom_final_reconstruction_type" minOccurs="0"/>
                        """
                        non_st_dicts = {
                            "final_rec_dict_in": cryst_dicts["final_rec_dict_in"],
                            "cat_soft_dict_in": cryst_dicts["cat_soft_dict_in"],
                            "p_sym_dict_in": cryst_dicts["p_sym_dict_in"],
                            "h_sym_dict_in": cryst_dicts["h_sym_dict_in"],
                        }
                        set_non_sub_tom_final_recon(ip_id_in, im_proc, non_st_dicts)

                    def set_el_crystal_parameters(im_proc, ip_id_in, cryst_dicts):
                        """
                        XSD: <xs:element name="crystal_parameters" type="crystal_parameters_type" minOccurs="0"/>
                        """
                        two_d_cryst_dict_in = cryst_dicts["two_d_cryst_dict_in"]
                        if ip_id_in in two_d_cryst_dict_in:
                            set_crystal_parameters(ip_id_in, im_proc, two_d_cryst_dict_in, const.EM_2D_CRYSTAL_ENTITY, parent_req=False)
                        three_d_cryst_dict_in = cryst_dicts["three_d_cryst_dict_in"]
                        if ip_id_in in three_d_cryst_dict_in:
                            set_crystal_parameters(ip_id_in, im_proc, three_d_cryst_dict_in, const.EM_3D_CRYSTAL_ENTITY, parent_req=False)

                    def set_el_startup_model(im_proc, ip_id_in):
                        """
                        XSD: <xs:element name="startup_model" type="starting_map_type" maxOccurs="unbounded" minOccurs="0">
                        """
                        st_mod_dict_in = cryst_dicts["st_mod_dict_in"]
                        if ip_id_in in st_mod_dict_in:
                            set_startup_model(ip_id_in, im_proc, st_mod_dict_in)

                    def set_el_ctf_correction(im_proc, ip_id_in, cryst_dicts):
                        """
                        XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                        """
                        ctf_corr_dict_in = cryst_dicts["ctf_corr_dict_in"]
                        if ip_id_in in ctf_corr_dict_in:
                            set_ctfcorrection(ip_id_in, im_proc, ctf_corr_dict_in)

                    def set_el_molecular_replacement(im_proc, cryst_dicts):
                        """
                        XSD: <xs:element name="molecular_replacement">
                        """

                        def set_molecular_replacement_type(mol_repl, cryst_dicts):
                            """
                            XSD: <xs:element name="molecular_replacement"> has
                            .. 3 elements
                            """

                            def set_el_starting_model():
                                """
                                XSD: <xs:element name="starting_model" maxOccurs="unbounded"> has
                                .. 2 elements
                                element 1: <xs:element name="access_code" type="pdb_code_type">
                                CIF:???
                                element 2: <xs:element name="chain" type="chain_type" maxOccurs="unbounded" minOccurs="0"/>
                                CIF:???
                                """

                            def set_el_resolution_range():
                                """
                                <xs:element name="resolution_range"> has
                                .. 2 elements
                                element 1: <xs:element name="high_resolution"> has
                                .. a base <xs:extension base="resolution_type"> and
                                .. 1 attribute <xs:attribute fixed="A" name="units" type="xs:token" use="required"/>
                                element 2: <xs:element name="low_resolution"> has
                                .. a base <xs:extension base="resolution_type"> and
                                .. 1 attribute <xs:attribute fixed="A" name="units" type="xs:token" use="required"/>
                                """

                            def set_el_software_list(mol_repl, cryst_dicts):
                                """
                                XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                                """
                                set_software_list(const.SOFT_MOLECULAR_REPLACEMENT, cryst_dicts["cat_soft_dict_in"], mol_repl.set_software_list)

                            # element 1
                            set_el_starting_model()
                            # element 2
                            set_el_resolution_range()
                            # element 3
                            set_el_software_list(mol_repl, cryst_dicts)

                        mol_repl = emdb.molecular_replacement_type()
                        set_molecular_replacement_type(mol_repl, cryst_dicts)
                        if mol_repl.has__content():
                            im_proc.set_molecular_replacement(mol_repl)

                    def set_el_lat_dist_corr_soft_list(im_proc, cryst_dicts):
                        """
                        XSD: <xs:element name="lattice_distortion_correction_software_list" type="software_list_type" minOccurs="0"/>
                        """
                        set_software_list(const.SOFT_LATTICE_DISTORTION_CORRECTION, cryst_dicts["cat_soft_dict_in"], im_proc.set_lattice_distortion_correction_software_list)

                    def set_el_sym_det_software_list(im_proc, cryst_dicts):
                        """
                        XSD: <xs:element name="symmetry_determination_software_list" type="software_list_type" minOccurs="0"/>
                        """
                        set_software_list(const.SOFT_SYMMETRY_DETERMINATION, cryst_dicts["cat_soft_dict_in"], im_proc.set_symmetry_determination_software_list)

                    def set_el_merging_software_list(im_proc, cryst_dicts):
                        """
                        XSD: <xs:element name="merging_software_list" type="software_list_type" minOccurs="0"/>
                        """
                        set_software_list(const.SOFT_CRYSTALLOGRAPHY_MERGING, cryst_dicts["cat_soft_dict_in"], im_proc.set_merging_software_list)

                    def set_el_cryst_statistics(im_proc, cryst_dicts):
                        """
                        XSD: <xs:element name="crystallography_statistics" type="crystallography_statistics_type" minOccurs="0"/>
                        """

                        def set_cryst_statistics_type(cry_stats, cry_stats_in):
                            """
                            XSD: <xs:complexType name="crystallography_statistics_type"> has
                            .. 11 elements
                            """

                            def set_el_num_int_measured(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="number_intensities_measured" type="xs:positiveInteger"/>
                                CIF: _em_diffraction_stats.num_intensities_measured 1590
                                """
                                set_cif_value(
                                    cry_stats.set_number_intensities_measured, "num_intensities_measured", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in, fmt=int
                                )

                            def set_el_number_structure_factors(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="number_structure_factors" type="xs:positiveInteger"/>
                                CIF: _em_diffraction_stats.num_structure_factors 1590
                                """
                                set_cif_value(cry_stats.set_number_structure_factors, "num_structure_factors", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in, fmt=int)

                            def set_el_fourier_space_coverage(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="fourier_space_coverage" type="xs:float">
                                CIF: _em_diffraction_stats.fourier_space_coverage 89.3
                                """
                                set_cif_value(cry_stats.set_fourier_space_coverage, "fourier_space_coverage", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in, fmt=float)

                            def set_el_r_sym(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="r_sym" type="xs:float"/>
                                CIF: _em_diffraction_stats.r_sym 0.244
                                """
                                set_cif_value(cry_stats.set_r_sym, "r_sym", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in, fmt=float)

                            def set_el_r_merge(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="r_merge" type="xs:float"/>
                                CIF: _em_diffraction_stats.r_merge 0.198
                                """
                                set_cif_value(cry_stats.set_r_merge, "r_merge", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in, fmt=float)

                            def set_el_overall_phase_error(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="overall_phase_error" type="xs:token" minOccurs="0"/>
                                CIF: _em_diffraction_stats.overall_phase_error "17.5"
                                """
                                set_cif_value(cry_stats.set_overall_phase_error, "overall_phase_error", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in)

                            def set_el_overall_phase_residual(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="overall_phase_residual" type="xs:float"/>
                                CIF: _em_diffraction_stats.overall_phase_residual 17.5
                                """
                                phase_residual = get_cif_value("overall_phase_residual", const.EM_DIFFRACTION_STATS, cry_stats_in)
                                if phase_residual is None:
                                    phase_residual = 0.0  # chosen default value
                                    txt = u"(_em_diffraction_stats.overall_phase_residual) is set to (%s) as no value is given and it is required." % phase_residual
                                    self.current_entry_log.warn_logs.append(
                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.change_title + txt)
                                    )
                                    self.log_formatted(self.warn_log_string, const.CHANGE_MADE + txt)
                                    set_cif_value(
                                        cry_stats.set_overall_phase_residual,
                                        "overall_phase_residual",
                                        const.EM_DIFFRACTION_STATS,
                                        cif_list=cry_stats_in,
                                        cif_value=phase_residual,
                                    )
                                else:
                                    set_cif_value(cry_stats.set_overall_phase_residual, "overall_phase_residual", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in, fmt=float)

                            def set_el_phase_err_rej_criteria(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="phase_error_rejection_criteria" type="xs:token"/>
                                CIF: _em_diffraction_stats.phase_error_rejection_criteria
                                """
                                set_cif_value(
                                    cry_stats.set_phase_error_rejection_criteria, "phase_error_rejection_criteria", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in
                                )

                            def set_el_high_resolution(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="high_resolution">
                                CIF: _em_diffraction_stats.high_resolution
                                """
                                set_cif_value(
                                    cry_stats.set_high_resolution,
                                    "high_resolution",
                                    const.EM_DIFFRACTION_STATS,
                                    cif_list=cry_stats_in,
                                    constructor=emdb.high_resolutionType,
                                    fmt=float,
                                    units=const.U_ANG,
                                )

                            def set_el_shell_list(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="shell_list" minOccurs="0"> has
                                .. 1 element
                                """

                                def set_shell_type(shell, cs_in):
                                    """
                                    XSD: <xs:element name="shell" maxOccurs="unbounded"> has
                                    .. 1 attribute and
                                    .. 6 elements
                                    """

                                    def set_attr_id(shell, cs_in):
                                        """
                                        XSD: <xs:attribute name="shell_id" type="xs:positiveInteger"/>
                                        CIF: _em_diffraction_shell.id
                                        """
                                        set_cif_value(shell.set_shell_id, "id", const.EM_DIFFRACTION_SHELL, cif_list=cs_in, fmt=int, parent_el_req=False)

                                    def set_el_high_resolution(shell, cs_in):
                                        """
                                        XSD: <xs:element name="high_resolution">
                                        CIF: _em_diffraction_shell.high_resolution 3.0
                                        """
                                        set_cif_value(
                                            shell.set_high_resolution,
                                            "high_resolution",
                                            const.EM_DIFFRACTION_SHELL,
                                            cif_list=cs_in,
                                            constructor=emdb.high_resolutionType,
                                            fmt=float,
                                            units=const.U_ANG,
                                            parent_el_req=False,
                                        )

                                    def set_el_low_resolution(shell, cs_in):
                                        """
                                        XSD: <xs:element name="low_resolution">
                                        CIF: _em_diffraction_shell.low_resolution 5.5
                                        """
                                        set_cif_value(
                                            shell.set_low_resolution,
                                            "low_resolution",
                                            const.EM_DIFFRACTION_SHELL,
                                            cif_list=cs_in,
                                            constructor=emdb.low_resolutionType,
                                            fmt=float,
                                            units=const.U_ANG,
                                            parent_el_req=False,
                                        )

                                    def set_el_number_structure_factors(shell, cs_in):
                                        """
                                        XSD: <xs:element name="number_structure_factors" type="xs:positiveInteger"/>
                                        CIF: _em_diffraction_shell.num_structure_factors 244
                                        """
                                        set_cif_value(
                                            shell.set_number_structure_factors,
                                            "num_structure_factors",
                                            const.EM_DIFFRACTION_SHELL,
                                            cif_list=cs_in,
                                            fmt=int,
                                            parent_el_req=False,
                                        )

                                    def set_el_phase_residual(shell, cs_in):
                                        """
                                        XSD: <xs:element name="phase_residual" type="xs:float"/>
                                        CIF: _em_diffraction_shell.phase_residual 13.5
                                        """
                                        set_cif_value(shell.set_phase_residual, "phase_residual", const.EM_DIFFRACTION_SHELL, cif_list=cs_in, fmt=float, parent_el_req=False)

                                    def set_el_fourier_space_coverage(shell, cs_in):
                                        """
                                        XSD: <xs:element name="fourier_space_coverage" type="xs:float">
                                        CIF: _em_diffraction_shell.fourier_space_coverage 93.2
                                        """
                                        set_cif_value(
                                            shell.set_fourier_space_coverage,
                                            "fourier_space_coverage",
                                            const.EM_DIFFRACTION_SHELL,
                                            cif_list=cs_in,
                                            fmt=float,
                                            parent_el_req=False,
                                        )

                                    def set_el_multiplicity(shell, cs_in):
                                        """
                                        XSD: <xs:element name="multiplicity" type="xs:float"/>
                                        CIF: _em_diffraction_shell.multiplicity 2.5
                                        """
                                        set_cif_value(shell.set_multiplicity, "multiplicity", const.EM_DIFFRACTION_SHELL, cif_list=cs_in, fmt=float, parent_el_req=False)

                                    # attribute 1
                                    set_attr_id(shell, cs_in)
                                    # element 1
                                    set_el_high_resolution(shell, cs_in)
                                    # element 2
                                    set_el_low_resolution(shell, cs_in)
                                    # element 3
                                    set_el_number_structure_factors(shell, cs_in)
                                    # element 4
                                    set_el_phase_residual(shell, cs_in)
                                    # element 5
                                    set_el_fourier_space_coverage(shell, cs_in)
                                    # element 6
                                    set_el_multiplicity(shell, cs_in)

                                cry_id_in = get_cif_value(const.K_ID, const.EM_DIFFRACTION_STATS, cry_stats_in)
                                if cry_id_in is not None:
                                    cry_shell_dict_in = cryst_dicts["cry_shell_dict_in"]
                                    if cry_id_in in cry_shell_dict_in:
                                        cry_shell_list_in = cry_shell_dict_in[cry_id_in]
                                        cry_shell_list = emdb.shell_listType()
                                        for cs_in in cry_shell_list_in:
                                            attr_id = get_cif_value("id", const.EM_DIFFRACTION_SHELL, cif_list=cs_in)
                                            el_high_resolution = get_cif_value("high_resolution", const.EM_DIFFRACTION_SHELL, cif_list=cs_in)
                                            el_low_resolution = get_cif_value("low_resolution", const.EM_DIFFRACTION_SHELL, cif_list=cs_in)
                                            el_number_structure_factors = get_cif_value("num_structure_factors", const.EM_DIFFRACTION_SHELL, cif_list=cs_in)
                                            el_phase_residual = get_cif_value("phase_residual", const.EM_DIFFRACTION_SHELL, cif_list=cs_in)
                                            el_fourier_space_coverage = get_cif_value("fourier_space_coverage", const.EM_DIFFRACTION_SHELL, cif_list=cs_in)
                                            el_multiplicity = get_cif_value("multiplicity", const.EM_DIFFRACTION_SHELL, cif_list=cs_in)
                                            shell_type_list = [
                                                attr_id,
                                                el_high_resolution,
                                                el_low_resolution,
                                                el_number_structure_factors,
                                                el_phase_residual,
                                                el_fourier_space_coverage,
                                                el_multiplicity,
                                            ]
                                            if any(x is not None for x in shell_type_list):
                                                shell = emdb.shellType()
                                                set_shell_type(shell, cs_in)
                                                if shell.has__content():
                                                    cry_shell_list.add_shell(shell)
                                        if cry_shell_list.has__content():
                                            cry_stats.set_shell_list(cry_shell_list)

                            def set_el_details(cry_stats, cry_stats_in):
                                """
                                XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                CIF: _em_diffraction_stats.details
                                """
                                set_cif_value(cry_stats.set_details, "details", const.EM_DIFFRACTION_STATS, cif_list=cry_stats_in)

                            # element 1
                            set_el_num_int_measured(cry_stats, cry_stats_in)
                            # element 2
                            set_el_number_structure_factors(cry_stats, cry_stats_in)
                            # element 3
                            set_el_fourier_space_coverage(cry_stats, cry_stats_in)
                            # element 4
                            set_el_r_sym(cry_stats, cry_stats_in)
                            # element 5
                            set_el_r_merge(cry_stats, cry_stats_in)
                            # element 6
                            set_el_overall_phase_error(cry_stats, cry_stats_in)
                            # element 7
                            set_el_overall_phase_residual(cry_stats, cry_stats_in)
                            # element 8
                            set_el_phase_err_rej_criteria(cry_stats, cry_stats_in)
                            # element 9
                            set_el_high_resolution(cry_stats, cry_stats_in)
                            # element 10
                            set_el_shell_list(cry_stats, cry_stats_in)
                            # element 11
                            set_el_details(cry_stats, cry_stats_in)

                        cry_stats_dict_in = cryst_dicts["cry_stats_dict_in"]
                        if ip_id_in in cry_stats_dict_in:
                            cry_stats_in = cry_stats_dict_in[ip_id_in]
                            cry_stats = emdb.crystallography_statistics_type()
                            set_cryst_statistics_type(cry_stats, cry_stats_in)
                            if cry_stats.has__content():
                                im_proc.set_crystallography_statistics(cry_stats)

                    # element 1
                    set_el_final_reconstruction(im_proc, ip_id_in)
                    # element 2
                    set_el_crystal_parameters(im_proc, ip_id_in, cryst_dicts)
                    # element 3
                    set_el_startup_model(im_proc, ip_id_in)
                    # element 4
                    set_el_ctf_correction(im_proc, ip_id_in, cryst_dicts)
                    # element 5
                    set_el_molecular_replacement(im_proc, cryst_dicts)
                    # element 6
                    set_el_lat_dist_corr_soft_list(im_proc, cryst_dicts)
                    # element 7
                    set_el_sym_det_software_list(im_proc, cryst_dicts)
                    # element 8
                    set_el_merging_software_list(im_proc, cryst_dicts)
                    # element 9
                    set_el_cryst_statistics(im_proc, cryst_dicts)

                def set_tom_proc_specifics(ip_id_in, im_proc, tomo_dicts):
                    """
                    Set elements specific for tomography image processing add group

                    Parameters:
                    @param ip_id_in: image processing id
                    @param im_proc: image processing object
                    @param tomo_dicts: a dictionary of the following dictionaries:
                    final_rec_dict_in, cat_soft_dict_in,
                    p_sym_dict_in, h_sym_dict_in, ctf_corr_dict_in,
                    two_d_cryst_dict_in, three_d_cryst_dict_in
                    XSD: <xs:group name="tomography_proc_add_group"> has
                    .. 4 elements
                    """

                    def set_el_final_reconstruction(im_proc, ip_id_in):
                        """
                        XSD: <xs:element name="final_reconstruction" type="non_subtom_final_reconstruction_type" minOccurs="0"/>
                        """
                        non_st_dicts = {
                            "final_rec_dict_in": tomo_dicts["final_rec_dict_in"],
                            "cat_soft_dict_in": tomo_dicts["cat_soft_dict_in"],
                            "p_sym_dict_in": tomo_dicts["p_sym_dict_in"],
                            "h_sym_dict_in": tomo_dicts["h_sym_dict_in"],
                        }
                        set_non_sub_tom_final_recon(ip_id_in, im_proc, non_st_dicts)

                    def set_el_alig_soft_list(im_proc, tomo_dicts):
                        """
                        XSD: <xs:element name="series_aligment_software_list" type="software_list_type" minOccurs="0"/>
                        """
                        cat_soft_dict_in = tomo_dicts["cat_soft_dict_in"]
                        set_software_list(const.SOFT_LATTICE_DISTORTION_CORRECTION, cat_soft_dict_in, im_proc.set_series_aligment_software_list)

                    def set_el_ctf_correction(im_proc, ip_id_in, tomo_dicts):
                        """
                        XSD: <xs:element name="ctf_correction" type="ctf_correction_type" minOccurs="0"/>
                        """
                        ctf_corr_dict_in = tomo_dicts["ctf_corr_dict_in"]
                        if ip_id_in in ctf_corr_dict_in:
                            set_ctfcorrection(ip_id_in, im_proc, ctf_corr_dict_in)

                    def set_el_crystal_parameters(im_proc, ip_id_in, tomo_dicts):
                        """
                        XSD: <xs:element name="crystal_parameters" type="crystal_parameters_type" minOccurs="0"/>
                        """
                        two_d_cryst_dict_in = tomo_dicts["two_d_cryst_dict_in"]
                        if ip_id_in in two_d_cryst_dict_in:
                            set_crystal_parameters(ip_id_in, im_proc, two_d_cryst_dict_in, const.EM_2D_CRYSTAL_ENTITY, parent_req=False)
                        three_d_cryst_dict_in = tomo_dicts["three_d_cryst_dict_in"]
                        if ip_id_in in three_d_cryst_dict_in:
                            set_crystal_parameters(ip_id_in, im_proc, three_d_cryst_dict_in, const.EM_3D_CRYSTAL_ENTITY, parent_req=False)

                    # element 1
                    set_el_final_reconstruction(im_proc, ip_id_in)
                    # element 2
                    set_el_alig_soft_list(im_proc, tomo_dicts)
                    # element 3
                    set_el_ctf_correction(im_proc, ip_id_in, tomo_dicts)
                    # element 4
                    set_el_crystal_parameters(im_proc, ip_id_in, tomo_dicts)

                # Create dictionaries indexed by em_image_processing_id
                soft_dict_in = make_list_of_dicts(const.EM_SOFTWARE, const.K_IMAGE_PROCESSING_ID)
                part_sel_dict_in = make_list_of_dicts(const.EM_PARTICLE_SELECTION, const.K_IMAGE_PROCESSING_ID, min_length=2)
                vol_sel_dict_in = make_dict(const.EM_VOLUME_SELECTION, const.K_IMAGE_PROCESSING_ID, min_length=2)
                ctf_corr_dict_in = make_dict(const.EM_CTF_CORRECTION, const.K_IMAGE_PROCESSING_ID)
                st_mod_dict_in = make_list_of_dicts(const.EM_START_MODEL, const.K_IMAGE_PROCESSING_ID)
                ang_dict_in = make_list_of_dicts(const.EM_EULER_ANGLE_ASSIGNMENT, const.K_IMAGE_PROCESSING_ID)
                final_class_dict_in = make_dict(const.EM_FINAL_CLASSIFICATION, const.K_IMAGE_PROCESSING_ID)
                final_2d_class_dict_in = make_dict(const.EM_FINAL_2D_CLASSIFICATION, const.K_IMAGE_PROCESSING_ID)
                final_rec_dict_in = make_dict(const.EM_3D_RECONSTRUCTION, const.K_IMAGE_PROCESSING_ID)
                p_sym_dict_in = make_dict(const.EM_SINGLE_PARTICLE_ENTITY, const.K_IMAGE_PROCESSING_ID)
                h_sym_dict_in = make_dict(const.EM_HELICAL_ENTITY, const.K_IMAGE_PROCESSING_ID)
                two_d_cryst_dict_in = make_dict(const.EM_2D_CRYSTAL_ENTITY, const.K_IMAGE_PROCESSING_ID)
                three_d_cryst_dict_in = make_dict(const.EM_3D_CRYSTAL_ENTITY, const.K_IMAGE_PROCESSING_ID)
                cry_stats_dict_in = make_dict(const.EM_DIFFRACTION_STATS, const.K_IMAGE_PROCESSING_ID)
                cry_shell_dict_in = make_list_of_dicts(const.EM_DIFFRACTION_SHELL, const.K_EM_DIFFRACTION_STATS_ID)
                im_rec_dict_in = make_dict(const.EM_IMAGE_RECORDING, const.K_ID)

                ip_list_in = assert_get_value(const.EM_IMAGE_PROCESSING, self.cif)
                for ip_in in ip_list_in:
                    # Create a dictionary of software keyed by category
                    cat_soft_dict_in = {}
                    ip_id_in = get_cif_value(const.K_ID, const.EM_IMAGE_PROCESSING, ip_in)
                    if ip_id_in in soft_dict_in:
                        cat_soft_dict_in = make_list_of_dicts(const.EM_SOFTWARE, "category", soft_dict_in[ip_id_in])
                    if em_method == const.EMM_SP:
                        sp_im_proc = emdb.singleparticle_processing_type()
                        sp_im_proc.original_tagname_ = "singleparticle_processing"
                        set_base_image_processing(ip_in, sp_im_proc)
                        sp_dicts = {
                            "cat_soft_dict_in": cat_soft_dict_in,
                            "ctf_corr_dict_in": ctf_corr_dict_in,
                            "part_sel_dict_in": part_sel_dict_in,
                            "st_mod_dict_in": st_mod_dict_in,
                            "ang_dict_in": ang_dict_in,
                            "p_sym_dict_in": p_sym_dict_in,
                            "h_sym_dict_in": h_sym_dict_in,
                            "final_class_dict_in": final_class_dict_in,
                            "final_2d_class_dict_in": final_2d_class_dict_in,
                            "final_rec_dict_in": final_rec_dict_in,
                        }
                        set_single_part_im_proc_specs(ip_id_in, sp_im_proc, sp_dicts)
                        struct_det.add_image_processing(sp_im_proc)
                    elif em_method == const.EMM_HEL:
                        hel_im_proc = emdb.helical_processing_type()
                        hel_im_proc.original_tagname_ = "helical_processing"
                        set_base_image_processing(ip_in, hel_im_proc)
                        helical_dicts = {
                            "part_sel_dict_in": part_sel_dict_in,
                            "final_rec_dict_in": final_rec_dict_in,
                            "cat_soft_dict_in": cat_soft_dict_in,
                            "p_sym_dict_in": p_sym_dict_in,
                            "h_sym_dict_in": h_sym_dict_in,
                            "ctf_corr_dict_in": ctf_corr_dict_in,
                            "ang_dict_in": ang_dict_in,
                            "st_mod_dict_in": st_mod_dict_in,
                            "two_d_cryst_dict_in": two_d_cryst_dict_in,
                            "three_d_cryst_dict_in": three_d_cryst_dict_in,
                        }
                        set_helical_proc_specifics(ip_id_in, hel_im_proc, helical_dicts)
                        struct_det.add_image_processing(hel_im_proc)
                    elif em_method == const.EMM_TOM:
                        tom_im_proc = emdb.tomography_processing_type()
                        tom_im_proc.original_tagname_ = "tomography_processing"
                        set_base_image_processing(ip_in, tom_im_proc)
                        tomo_dicts = {
                            "final_rec_dict_in": final_rec_dict_in,
                            "cat_soft_dict_in": cat_soft_dict_in,
                            "p_sym_dict_in": p_sym_dict_in,
                            "h_sym_dict_in": h_sym_dict_in,
                            "ctf_corr_dict_in": ctf_corr_dict_in,
                            "two_d_cryst_dict_in": two_d_cryst_dict_in,
                            "three_d_cryst_dict_in": three_d_cryst_dict_in,
                        }
                        set_tom_proc_specifics(ip_id_in, tom_im_proc, tomo_dicts)
                        struct_det.add_image_processing(tom_im_proc)
                    elif em_method == const.EMM_STOM:
                        subtom_im_proc = emdb.subtomogram_averaging_processing_type()
                        subtom_im_proc.original_tagname_ = "subtomogram_averaging_processing"
                        set_base_image_processing(ip_in, subtom_im_proc)
                        subtom_dicts = {
                            "final_rec_dict_in": final_rec_dict_in,
                            "cat_soft_dict_in": cat_soft_dict_in,
                            "p_sym_dict_in": p_sym_dict_in,
                            "h_sym_dict_in": h_sym_dict_in,
                            "vol_sel_dict_in": vol_sel_dict_in,
                            "ctf_corr_dict_in": ctf_corr_dict_in,
                            "final_class_dict_in": final_class_dict_in,
                            "ang_dict_in": ang_dict_in,
                            "two_d_cryst_dict_in": two_d_cryst_dict_in,
                            "three_d_cryst_dict_in": three_d_cryst_dict_in,
                        }
                        set_subtom_av_im_proc_specifics(ip_id_in, subtom_im_proc, subtom_dicts, em_method)
                        struct_det.add_image_processing(subtom_im_proc)
                    elif em_method == const.EMM_EC:
                        cryst_im_proc = emdb.crystallography_processing_type()
                        cryst_im_proc.original_tagname_ = "crystallography_processing"
                        set_base_image_processing(ip_in, cryst_im_proc)
                        cryst_dicts = {
                            "final_rec_dict_in": final_rec_dict_in,
                            "cat_soft_dict_in": cat_soft_dict_in,
                            "p_sym_dict_in": p_sym_dict_in,
                            "h_sym_dict_in": h_sym_dict_in,
                            "two_d_cryst_dict_in": two_d_cryst_dict_in,
                            "three_d_cryst_dict_in": three_d_cryst_dict_in,
                            "cry_shell_dict_in": cry_shell_dict_in,
                            "st_mod_dict_in": st_mod_dict_in,
                            "ctf_corr_dict_in": ctf_corr_dict_in,
                            "cry_stats_dict_in": cry_stats_dict_in,
                        }
                        set_cryst_proc_specifics(ip_id_in, cryst_im_proc, cryst_dicts)
                        struct_det.add_image_processing(cryst_im_proc)

                    # Check if microscopy software is there in the list
                    # XSD: in <xs:complexType name="base_microscopy_type">:
                    # XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                    if const.SOFT_IMAGE_ACQUISITION in cat_soft_dict_in:
                        ip_im_rec_id_in = get_cif_value(const.K_IMAGE_RECORDING_ID, const.EM_IMAGE_PROCESSING, ip_in)
                        if ip_im_rec_id_in is not None:
                            if ip_im_rec_id_in in im_rec_dict_in:
                                im_rec_in = im_rec_dict_in[ip_im_rec_id_in]
                                if im_rec_in is not None:
                                    mic_id_in = get_cif_value(const.K_IMAGING_ID, const.EM_IMAGE_RECORDING, im_rec_in)
                                    if mic_id_in is not None:
                                        mic_id = int(mic_id_in)
                                        mic_list = microscopy_list.get_microscopy()
                                        for mic in mic_list:
                                            if mic.get_microscopy_id() == mic_id:
                                                set_software_list(const.SOFT_IMAGE_ACQUISITION, cat_soft_dict_in, mic.set_software_list)

            # attribute 1
            set_attr_id(struct_det)
            # element 1
            em_method = set_el_method(struct_det)
            # element 2
            set_el_aggregation_state(struct_det)
            # element 3
            set_el_mol_and_complexes()
            # element 4
            set_el_specimen_prep_list(struct_det, em_method)
            # element 5
            microscopy_list = emdb.microscopy_listType()
            set_el_microscopy_list(struct_det, microscopy_list)
            # element 6
            set_el_image_processing(struct_det, microscopy_list)

        def make_map(map_in, map_name=None):
            """
            Make a map element from the cif EM_MAP category

            Parameters:
            @param map_in: cif dictionary item containing one row of EM_MAP
            @param map_name: name to give file -
                             if none just take the name found in map file
            @return emdb.map_type() element
            XSD: <xs:element name="map" type="map_type">
            """

            def set_map_type(em_map, map_in, map_name):
                """
                XSD: <xs:element name="map" type="map_type"> has
                .. 2 attributes
                .. 14 elements
                """

                def set_attr_format(em_map, map_in):
                    """
                    XSD: <xs:attribute name="format" fixed="CCP4" use="required"/>
                    CIF: _em_map.format CCP4
                    """
                    set_cif_value(em_map.set_format, "format", const.EM_MAP, cif_list=map_in)

                def set_attr_size_kbytes(em_map, map_in):
                    """
                    XSD: <xs:attribute name="size_kbytes" type="xs:positiveInteger" use="required"/>
                    CIF: _em_map.size_kb
                    """
                    set_cif_value(em_map.set_size_kbytes, "size_kb", const.EM_MAP, cif_list=map_in, fmt=lambda x: int(float(x) / 1000.0))

                def set_el_file(em_map, map_in, map_name_in):
                    """
                    XSD: <xs:element name="file">
                    CIF: _em_map.file emd_5470.map.gz
                    """
                    map_name = map_name_in
                    map_in_dict = dict(map_in)
                    map_type = map_in_dict.get("_em_map.type")
                    if map_type == "primary map":
                        if map_name_in is None:
                            map_name = get_cif_value("file", const.EM_MAP, map_in)
                        if map_name is None or map_name == "" or map_name.isspace():
                            map_name = self.emdb_id_u.lower() + ".map.gz"
                            txt = u"Map file name is not given for (_em_map.file). Map name is set to (%s)." % map_name
                            self.current_entry_log.warn_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.warn_title + txt))
                            self.log_formatted(self.warn_log_string, const.NOT_REQUIRED_ALERT + txt)
                    else:
                        if map_name == "" or map_name.isspace():
                            txt = u"The value for (_em_map.file) is not given in (%s)." % map_in
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                    set_cif_value(em_map.set_file, "file", const.EM_MAP, cif_list=map_in, cif_value=map_name, fmt=str.lower)

                def set_el_symmetry(em_map, map_in):
                    """
                    XSD: <xs:element name="symmetry" type="applied_symmetry_type"> has
                    .. 1 element out of 3 choices
                    """

                    def set_applied_symmetry_type(app_sym):
                        """
                        XSD: <xs:complexType name="applied_symmetry_type">
                        Only space_group choice has mmCIF value
                        CIF: _em_map.symmetry_space_group 1
                        """
                        set_cif_value(app_sym.set_space_group, "symmetry_space_group", const.EM_MAP, cif_list=map_in, fmt=int)

                    app_sym = emdb.applied_symmetry_type()
                    set_applied_symmetry_type(app_sym)
                    if app_sym.has__content():
                        em_map.set_symmetry(app_sym)

                def set_el_data_type(em_map, map_in):
                    """
                    XSD: <xs:element name="data_type" type="map_data_type"/>
                    CIF: _em_map.data_type 'Image stored as signed byte'
                    """
                    set_cif_value(em_map.set_data_type, "data_type", const.EM_MAP, cif_list=map_in, fmt=const.MAP_DATA_TYPE_CIF2XML)

                def set_el_dimensions(em_map, map_in):
                    """
                    XSD: <xs:element name="dimensions" type="integer_vector_map_type"/>
                    """

                    def set_integer_vector_map_type(dim, map_in):
                        """
                        XSD: <xs:complexType name="integer_vector_map_type"> has
                        .. 3 elements
                        """

                        def set_el_col(dim, map_in):
                            """
                            XSD: <xs:element name="col" type="xs:positiveInteger"/>
                            CIF: _em_map.dimensions_col
                            """
                            set_cif_value(dim.set_col, "dimensions_col", const.EM_MAP, cif_list=map_in, fmt=int)

                        def set_el_row(dim, map_in):
                            """
                            XSD: <xs:element name="row" type="xs:positiveInteger"/>
                            CIF: _em_map.dimensions_row
                            """
                            set_cif_value(dim.set_row, "dimensions_row", const.EM_MAP, cif_list=map_in, fmt=int)

                        def set_el_sec(dim, map_in):
                            """
                            XSD: <xs:element name="sec" type="xs:positiveInteger"/>
                            CIF: _em_map.dimensions_sec
                            """
                            set_cif_value(dim.set_sec, "dimensions_sec", const.EM_MAP, cif_list=map_in, fmt=int)

                        # element 1
                        set_el_col(dim, map_in)
                        # element 2
                        set_el_row(dim, map_in)
                        # element 3
                        set_el_sec(dim, map_in)

                    dim = emdb.integer_vector_map_type()
                    set_integer_vector_map_type(dim, map_in)
                    if dim.has__content():
                        em_map.set_dimensions(dim)

                def set_el_origin(em_map, map_in):
                    """
                    XSD: <xs:element name="origin">
                    """

                    def set_origin_type(orig, map_in):
                        """
                        XSD: <xs:element name="origin"> has
                        .. 3 elements
                        """

                        def set_el_col(orig, map_in):
                            """
                            XSD: <xs:element name="col" type="xs:integer"/>
                            CIF: _em_map.origin_col
                            """
                            set_cif_value(orig.set_col, "origin_col", const.EM_MAP, cif_list=map_in, fmt=int)

                        def set_el_row(orig, map_in):
                            """
                            XSD: <xs:element name="row" type="xs:integer"/>
                            CIF: _em_map.origin_row
                            """
                            set_cif_value(orig.set_row, "origin_row", const.EM_MAP, cif_list=map_in, fmt=int)

                        def set_el_sec(orig, map_in):
                            """
                            XSD: <xs:element name="sec" type="xs:integer"/>
                            CIF: _em_map.origin_sec
                            """
                            set_cif_value(orig.set_sec, "origin_sec", const.EM_MAP, cif_list=map_in, fmt=int)

                        # element 1
                        set_el_col(orig, map_in)
                        # element 2
                        set_el_row(orig, map_in)
                        # element 3
                        set_el_sec(orig, map_in)

                    orig = emdb.originType()
                    set_origin_type(orig, map_in)
                    em_map.set_origin(orig)

                def set_el_spacing(em_map, map_in):
                    """
                    XSD: <xs:element name="spacing">
                    """

                    def set_spacing_type(spc, map_in):
                        """
                        XSD: <xs:element name="spacing"> has
                        .. 3 elements
                        """

                        def set_el_x(spc, map_in):
                            """
                            XSD: <xs:element name="x" type="xs:positiveInteger"/>
                            CIF: _em_map.spacing_x
                            """
                            set_cif_value(spc.set_x, "spacing_x", const.EM_MAP, cif_list=map_in, fmt=int)

                        def set_el_y(spc, map_in):
                            """
                            XSD: <xs:element name="y" type="xs:nonNegativeInteger"/>
                            CIF: _em_map.spacing_y
                            """
                            set_cif_value(spc.set_y, "spacing_y", const.EM_MAP, cif_list=map_in, fmt=int)

                        def set_el_z(spc, map_in):
                            """
                            XSD: <xs:element name="z" type="xs:nonNegativeInteger"/>
                            CIF: _em_map.spacing_z
                            """
                            set_cif_value(spc.set_z, "spacing_z", const.EM_MAP, cif_list=map_in, fmt=int)

                        # element 1
                        set_el_x(spc, map_in)
                        # element 2
                        set_el_y(spc, map_in)
                        # element 3
                        set_el_z(spc, map_in)

                    spc = emdb.spacingType()
                    set_spacing_type(spc, map_in)
                    em_map.set_spacing(spc)

                def set_el_cell(em_map, map_in):
                    """
                    XSD: <xs:element name="cell">
                    """

                    def set_cell_type(cell, map_in):
                        """
                        XSD: <xs:element name="cell"> has
                        .. 6 elements
                        """

                        def set_el_a(cell, map_in):
                            """
                            XSD: <xs:element name="a" type="cell_type"/>
                            CIF: _em_map.cell_a
                            """
                            set_cif_value(cell.set_a, "cell_a", const.EM_MAP, cif_list=map_in, constructor=emdb.cell_type, fmt=float, units=const.U_ANG)

                        def set_el_b(cell, map_in):
                            """
                            XSD: <xs:element name="b" type="cell_type"/>
                            CIF: _em_map.cell_b
                            """
                            set_cif_value(cell.set_b, "cell_b", const.EM_MAP, cif_list=map_in, constructor=emdb.cell_type, fmt=float, units=const.U_ANG)

                        def set_el_c(cell, map_in):
                            """
                            XSD: <xs:element name="c" type="cell_type"/>
                            CIF: _em_map.cell_c
                            """
                            set_cif_value(cell.set_c, "cell_c", const.EM_MAP, cif_list=map_in, constructor=emdb.cell_type, fmt=float, units=const.U_ANG)

                        def set_el_alpha(cell, map_in):
                            """
                            XSD: <xs:element name="alpha" type="cell_angle_type"/>
                            CIF: _em_map.cell_alpha
                            """
                            set_cif_value(cell.set_alpha, "cell_alpha", const.EM_MAP, cif_list=map_in, constructor=emdb.cell_type, fmt=float, units=const.U_DEG)

                        def set_el_beta(cell, map_in):
                            """
                            XSD: <xs:element name="beta" type="cell_angle_type"/>
                            CIF: _em_map.cell_beta
                            """
                            set_cif_value(cell.set_beta, "cell_beta", const.EM_MAP, cif_list=map_in, constructor=emdb.cell_type, fmt=float, units=const.U_DEG)

                        def set_el_gamma(cell, map_in):
                            """
                            XSD: <xs:element name="gamma" type="cell_angle_type"/>
                            CIF: _em_map.cell_gamma
                            """
                            set_cif_value(cell.set_gamma, "cell_gamma", const.EM_MAP, cif_list=map_in, constructor=emdb.cell_type, fmt=float, units=const.U_DEG)

                        # element 1
                        set_el_a(cell, map_in)
                        # element 2
                        set_el_b(cell, map_in)
                        # element 3
                        set_el_c(cell, map_in)
                        # element 4
                        set_el_alpha(cell, map_in)
                        # element 5
                        set_el_beta(cell, map_in)
                        # element 6
                        set_el_gamma(cell, map_in)

                    cell = emdb.cellType()
                    set_cell_type(cell, map_in)
                    em_map.set_cell(cell)

                def set_el_axis_order(em_map, map_in):
                    """
                    XSD: <xs:element name="axis_order">
                    """

                    def set_axis_order_type(axis, map_in):
                        """
                        XSD: <xs:element name="axis_order"> has
                        .. 3 elements
                        """

                        def set_el_fast(axis, map_in):
                            """
                            XSD: <xs:element name="fast">
                            CIF: _em_map.axis_order_fast
                            """
                            set_cif_value(axis.set_fast, "axis_order_fast", const.EM_MAP, cif_list=map_in)

                        def set_el_medium(axis, map_in):
                            """
                            XSD: <xs:element name="medium">
                            CIF: _em_map.axis_order_medium
                            """
                            set_cif_value(axis.set_medium, "axis_order_medium", const.EM_MAP, cif_list=map_in)

                        def set_el_slow(axis, map_in):
                            """
                            XSD: <xs:element name="slow">
                            CIF: _em_map.axis_order_slow
                            """
                            set_cif_value(axis.set_slow, "axis_order_slow", const.EM_MAP, cif_list=map_in)

                        # element 1
                        set_el_fast(axis, map_in)
                        # element 2
                        set_el_medium(axis, map_in)
                        # element 3
                        set_el_slow(axis, map_in)

                    axis = emdb.axis_orderType()
                    set_axis_order_type(axis, map_in)
                    em_map.set_axis_order(axis)

                def set_el_statistics(em_map, map_in):
                    """
                    XSD: <xs:element name="statistics" type="map_statistics_type"/>
                    """

                    def set_map_statistics_type(stat, map_in):
                        """
                        XSD: <xs:element name="statistics" type="map_statistics_type"/> has
                        .. 4 elements
                        """

                        def set_el_minimum(stat, map_in):
                            """
                            XSD: <xs:element name="minimum" type="xs:float"/>
                            CIF: _em_map.statistics_minimum
                            """
                            set_cif_value(stat.set_minimum, "statistics_minimum", const.EM_MAP, cif_list=map_in, fmt=float)

                        def set_el_maximum(stat, map_in):
                            """
                            XSD: <xs:element name="maximum" type="xs:float"/>
                            CIF: _em_map.statistics_maximum
                            """
                            set_cif_value(stat.set_maximum, "statistics_maximum", const.EM_MAP, cif_list=map_in, fmt=float)

                        def set_el_average(stat, map_in):
                            """
                            XSD: <xs:element name="average" type="xs:float"/>
                            CIF: _em_map.statistics_average
                            """
                            set_cif_value(stat.set_average, "statistics_average", const.EM_MAP, cif_list=map_in, fmt=float)

                        def set_el_std(stat, map_in):
                            """
                            XSD: <xs:element name="std" type="xs:float"/>
                            CIF: _em_map.statistics_std
                            """
                            set_cif_value(stat.set_std, "statistics_std", const.EM_MAP, cif_list=map_in, fmt=float)

                        # element 1
                        set_el_minimum(stat, map_in)
                        # element 2
                        set_el_maximum(stat, map_in)
                        # element 3
                        set_el_average(stat, map_in)
                        # element 4
                        set_el_std(stat, map_in)

                    stat = emdb.map_statistics_type()
                    set_map_statistics_type(stat, map_in)
                    em_map.set_statistics(stat)

                def set_el_pixel_spacing(em_map, map_in):
                    """
                    XSD: <xs:element name="pixel_spacing">
                    """

                    def set_pixel_spacing_type(pix, map_in):
                        """
                        XSD: <xs:element name="pixel_spacing"> has
                        .. 3 elements
                        """

                        def set_el_x(pix, map_in):
                            """
                            XSD: <xs:element name="x" type="pixel_spacing_type"/>
                            CIF: _em_map.pixel_spacing_x
                            """
                            set_cif_value(pix.set_x, "pixel_spacing_x", const.EM_MAP, cif_list=map_in, constructor=emdb.pixel_spacing_type, fmt=float, units=const.U_ANG)

                        def set_el_y(pix, map_in):
                            """
                            XSD: <xs:element name="y" type="pixel_spacing_type"/>
                            CIF: _em_map.pixel_spacing_y
                            """
                            set_cif_value(pix.set_y, "pixel_spacing_y", const.EM_MAP, cif_list=map_in, constructor=emdb.pixel_spacing_type, fmt=float, units=const.U_ANG)

                        def set_el_z(pix, map_in):
                            """
                            XSD: <xs:element name="z" type="pixel_spacing_type"/>
                            CIF: _em_map.pixel_spacing_z
                            """
                            set_cif_value(pix.set_z, "pixel_spacing_z", const.EM_MAP, cif_list=map_in, constructor=emdb.pixel_spacing_type, fmt=float, units=const.U_ANG)

                        # element 1
                        set_el_x(pix, map_in)
                        # element 2
                        set_el_y(pix, map_in)
                        # element 3
                        set_el_z(pix, map_in)

                    pix = emdb.pixel_spacingType()
                    set_pixel_spacing_type(pix, map_in)
                    em_map.set_pixel_spacing(pix)

                def set_el_contour_list(em_map, map_in):
                    """
                    XSD: <xs:element name="contour_list" minOccurs="0">
                    """

                    def set_contour_type(cntr, map_in):
                        """
                        <xs:element name="contour" maxOccurs="unbounded"> has
                        .. 1 attribute and
                        .. 2 elements
                        """

                        def set_attr_primary(cntr):
                            """
                            XSD: <xs:attribute name="primary" type="xs:boolean" use="required"/>
                            IS THERE NO OTHER THAN PRIMARY?
                            """
                            set_cif_value(cntr.set_primary, cif_value=True)

                        def set_el_contour_level(cntr, map_in):
                            """
                            XSD: <xs:element name="level" type="xs:float"  minOccurs="0">
                            CIF: _em_map.contour_level {author,emdb,software}
                            Contour level had to be made non-mandatory as it's not given for tomograms
                            _em_experiment.reconstruction_method != "TOMOGRAPHY"
                            """
                            # check if the map contour level is not None
                            cntr_level = get_cif_value("contour_level", const.EM_MAP, cif_list=map_in)
                            if cntr_level == "None":
                                self.create_xml = False
                                txt = u'Contour level is "%s". This is not correct. The XML file is not going to be created now.' % cntr_level
                                self.current_entry_log.error_logs.append(
                                    self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                )
                                self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                            else:
                                # contour level can be set for the primary map
                                # get structure determination method
                                struct_det_method = get_cif_value("reconstruction_method", const.EM_EXPERIMENT)
                                # determine map type
                                map_type = get_cif_value("type", const.EM_MAP, cif_list=map_in)
                                if map_type == "primary map":
                                    if struct_det_method != "TOMOGRAPHY":
                                        cntr_level = get_cif_value("contour_level", const.EM_MAP, cif_list=map_in)
                                        if cntr_level is not None:
                                            if not isinstance(cntr_level, str):
                                                set_cif_value(cntr.set_level, "contour_level", const.EM_MAP, cif_list=map_in, fmt=float)
                                            else:
                                                # contour level is a string; check if the string can be converted
                                                if is_number(cntr_level.lstrip("+-")):
                                                    cl_float = float(cntr_level.lstrip("+-"))
                                                    set_cif_value(cntr.set_level, "contour_level", const.EM_MAP, cif_list=map_in, cif_value=cl_float)
                                                else:
                                                    txt = u"Contour level is given as a text value of %s. This is not correct. It should be a number." % cntr_level
                                                    self.current_entry_log.error_logs.append(
                                                        self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                                    )
                                                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                                        else:
                                            txt = u"Contour level is missing for %s." % struct_det_method
                                            self.current_entry_log.error_logs.append(
                                                self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                            )
                                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                                    else:
                                        txt = u"Contour level is not set for TOMOGRAPHY primary map."
                                        self.current_entry_log.info_logs.append(
                                            self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.info_title + txt)
                                        )
                                        self.log_formatted(self.info_log_string, const.INFO_ALERT + txt)
                            if map_type == "mask":
                                cntr_level = get_cif_value("contour_level", const.EM_MAP, cif_list=map_in)
                                if cntr_level is not None:
                                    if not isinstance(cntr_level, str):
                                        set_cif_value(cntr.set_level, "contour_level", const.EM_MAP, cif_list=map_in, fmt=float)
                                    else:
                                        # contour level is a string; check if the string can be converted
                                        if is_number(cntr_level.lstrip("+-")):
                                            cl_float = float(cntr_level.lstrip("+-"))
                                            set_cif_value(cntr.set_level, "contour_level", const.EM_MAP, cif_list=map_in, cif_value=cl_float)
                                        else:
                                            txt = u"Contour level is given as a text value of %s. This is not correct. It should be a number." % cntr_level
                                            self.current_entry_log.error_logs.append(
                                                self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt)
                                            )
                                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                        def set_el_source(cntr, map_in):
                            """
                            XSD: <xs:element name="source" minOccurs="0">
                            CIF: _em_map.contour_level_source
                            """
                            set_cif_value(cntr.set_source, "contour_level_source", const.EM_MAP, cif_list=map_in, fmt=str.upper)

                        # attribute 1
                        set_attr_primary(cntr)
                        # element 1
                        set_el_contour_level(cntr, map_in)
                        # element 2
                        set_el_source(cntr, map_in)

                    cntr_list = emdb.contour_listType()
                    cntr = emdb.contourType()
                    set_contour_type(cntr, map_in)
                    if cntr.has__content():
                        cntr_list.add_contour(cntr)
                    if cntr_list.has__content():
                        em_map.set_contour_list(cntr_list)

                def set_el_label(em_map, map_in):
                    """
                    XSD: <xs:element name="label" type="xs:token" minOccurs="0"/>
                    CIF: _em_map.label
                    """
                    set_cif_value(em_map.set_label, "label", const.EM_MAP, cif_list=map_in)

                def set_el_annotation_details(em_map, map_in):
                    """
                    XSD: <xs:element name="annotation_details" type="xs:string" minOccurs="0"/>
                    CIF: _em_map.annotation_details
                    """
                    set_cif_value(em_map.set_annotation_details, "annotation_details", const.EM_MAP, cif_list=map_in)

                def set_el_details():
                    """
                    XSD: <xs:element name="details" type="xs:string" minOccurs="0">
                    Deprecated (2014-10-21)
                    """

                # attribute 1
                set_attr_format(em_map, map_in)
                # attribute 2
                set_attr_size_kbytes(em_map, map_in)
                # element 1
                set_el_file(em_map, map_in, map_name)
                # element 2
                set_el_symmetry(em_map, map_in)
                # element 3
                set_el_data_type(em_map, map_in)
                # element 4
                set_el_dimensions(em_map, map_in)
                # element 5
                set_el_origin(em_map, map_in)
                # element 6
                set_el_spacing(em_map, map_in)
                # element 7
                set_el_cell(em_map, map_in)
                # element 8
                set_el_axis_order(em_map, map_in)
                # element 9
                set_el_statistics(em_map, map_in)
                # element 10
                set_el_pixel_spacing(em_map, map_in)
                # element 11
                set_el_contour_list(em_map, map_in)
                # element 12
                set_el_label(em_map, map_in)
                # element 13
                set_el_annotation_details(em_map, map_in)
                # element 14
                set_el_details()

            em_map = emdb.map_type()
            set_map_type(em_map, map_in, map_name)
            return em_map

        def get_canonical_map_name(map_in, map_type):
            """
            The map names will often be internal names of files in the D&A system.
            This routine will return the intended map name in archive

            Parameters:
            @param map_in: cif dictionary item containing one row of EM_MAP
            @param map_type: map types (PRIMARY, HALFMAP, ADDMAP, MASK)
            @return: canonical map name
            """
            map_name_in = get_cif_value("file", const.EM_MAP, map_in)
            map_name_out = ""
            if map_name_in is not None:
                if map_type == "PRIMARY":
                    map_name_out = "%s.map.gz" % self.emdb_id_u
                elif map_type == "HALFMAP":
                    auth_match = re.match(const.CIF_HALF_MAP_RE, map_name_in)
                    if auth_match is not None:
                        match_groups = auth_match.groups()
                        if match_groups is not None:
                            map_name_out = "%s_half_map_%s.map.gz" % (self.emdb_id_u, match_groups[0])
                        else:
                            map_name_out = map_name_in
                elif map_type == "ADDMAP":
                    auth_match = re.match(const.CIF_ADD_MAP_RE, map_name_in)
                    if auth_match is not None:
                        match_groups = auth_match.groups()
                        if match_groups is not None:
                            map_name_out = "%s_additional_%s.map.gz" % (self.emdb_id_u, match_groups[0])
                        else:
                            map_name_out = map_name_in
                elif map_type == "MASK":
                    auth_match = re.match(const.CIF_MSK_MAP_RE, map_name_in)
                    if auth_match is not None:
                        match_groups = auth_match.groups()
                        if match_groups is not None:
                            map_name_out = "%s_msk_%s.map" % (self.emdb_id_u, match_groups[0])
                        else:
                            map_name_out = map_name_in
                else:
                    map_name_out = map_name_in
            return map_name_out

        def set_entry_type():
            """
            Sets <xs:complexType name="entry_type"> that has:
            ..2 attributes
            ..a sequence of 7 elements
            """

            def set_attr_emdb_id():
                """
                XSD: <xs:attribute name="emdb_id" type="emdb_id_type" use="required">
                Note: emdb_id_type expects value="EMD-d{4,}
                CIF: _database_2.database_id _database_2.database_code
                CIF example: PDB   5LZZ
                CIF example: WWPDB D_1200001653
                CIF example: EMDB  EMD-4137
                """
                db2_in = assert_get_value(const.DATABASE_2, self.cif)
                # There are sometimes two emdb lines - why?
                db_id_dict = make_list_of_dicts(const.DATABASE_2, "database_id", db2_in, 2)
                if "EMDB" in db_id_dict:
                    emdb_db_id = db_id_dict["EMDB"][0]
                    if emdb_db_id is not None:
                        emdb_id = get_cif_value("database_code", const.DATABASE_2, emdb_db_id)
                        if emdb_id is not None:
                            # The following is as a prefix for file names in the archive
                            self.emdb_id_u = emdb_id.replace("-", "_")
                            self.xml_out.set_emdb_id(emdb_id)
                            if self.entry_in_translation_log.id == "":
                                self.entry_in_translation_log.id = emdb_id

                            if self.__show_log_id:
                                # Add entry ID into the warning logger messages
                                self.log_formatted(self.error_log_string, self.emdb_id_u)
                        else:
                            txt = u"The value for EMDB is not given in CIF for (_database_2.database_id _database_2.database_code). SOLUTION:  have (EMDB EMD-xxxx) instead of (EMDB .) or (EMDB ?)."
                            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                            self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                else:
                    txt = u"EMDB id is not given in CIF for (_database_2.database_id _database_2.database_code). SOLUTION: add (EMDB EMD-xxxx) in CIF."
                    self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

            def set_attr_version():
                """
                XSD: <xs:attribute name="version" type="xs:token" default="3.0.0.1"/>
                NOT IN CIF: The value is added in here -
                """
                # no need to set the version as the dafualt value is given
                # self.xml_out.set_version(const.XML_OUT_VERSION)
                pass  # pylint: disable=unnecessary-pass

            def set_el_admin():
                """
                XSD: <xs:element name="admin" type="admin_type">
                CIF: _em_admin
                """
                admin = emdb.admin_type()
                set_admin_type(admin)
                self.xml_out.set_admin(admin)

            def set_el_crossreferences():
                """
                XSD: <xs:element name="crossreferences" type="crossreferences_type"/>
                """
                cross_references = emdb.crossreferences_type()
                set_crossreferences_type(cross_references)
                self.xml_out.set_crossreferences(cross_references)

            def set_el_sample():
                """
                XSD: <xs:element name="sample" type="sample_type">
                """
                sample = emdb.sample_type()
                set_sample_type(sample)
                self.xml_out.set_sample(sample)

            def set_el_struct_det_list():
                """
                XSD: <xs:element name="structure_determination_list">
                """

                def set_struct_det_list(sd_list):
                    """
                    XSD: <xs:element name="structure_determination_list"> is
                    .. a sequence of 1 element
                    """
                    struct_det = emdb.structure_determination_type()
                    set_struct_determination_type(struct_det)
                    sd_list.add_structure_determination(struct_det)

                sd_list = emdb.structure_determination_listType()
                set_struct_det_list(sd_list)
                self.xml_out.set_structure_determination_list(sd_list)

            def set_el_map():
                """
                XSD: <xs:element name="map" type="map_type">
                """
                map_dict_in = make_list_of_dicts(const.EM_MAP, "type")
                pr_map_list_in = map_dict_in[const.MAP_PRIMARY] if const.MAP_PRIMARY in map_dict_in else []
                len_pr_map_list_in = len(pr_map_list_in)
                if len_pr_map_list_in != 1:
                    txt = u"There should be one and only one primary map. (%d) map(s) found." % len_pr_map_list_in
                    self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                else:
                    em_map = make_map(pr_map_list_in[0], get_canonical_map_name(pr_map_list_in[0], "PRIMARY"))
                    if em_map.has__content():
                        self.xml_out.set_map(em_map)
                    else:
                        txt = u"No information given for the primary map."
                        self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                        self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

            def set_el_interpretation():
                """
                XSD: <xs:element name="interpretation" type="interpretation_type"  minOccurs="0"/>
                """

                def set_interpretation_type(intrp, map_dict_in):
                    """
                    XSD: <xs:complexType name="interpretation_type"> has
                    .. 6 elements
                    """

                    def set_el_modelling_list(intrp):
                        """
                        XSD: <xs:element name="modelling_list" minOccurs="0">
                        Modelling - each modelling can have multiple models
                        """

                        def set_modelling_type(modelling, modelling_in, cat_soft_dict_in):
                            """
                            XSD: <xs:complexType name="modelling_type"> has
                            .. 8 elements
                            """

                            def set_el_initial_model(modelling, modelling_id_in, model_dict_in):
                                """
                                XSD: <xs:element name="initial_model" maxOccurs="unbounded">
                                """

                                def set_initial_model_type(model, model_in):
                                    """
                                    XSD: <xs:element name="initial_model" maxOccurs="unbounded"> has
                                    .. 3 elements
                                    """

                                    def set_el_access_code(model, model_in):
                                        """
                                        XSD: <xs:element name="access_code">
                                        CIF: _em_3d_fitting_list.pdb_entry_id  1EHZ
                                        CIF: _em_3d_fitting_list.entry_id  1EHZ
                                        pattern "d[dA-Za-z]{3}"
                                        """
                                        db_name = get_cif_value("source_name", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                        accession_code = get_cif_value("accession_code", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                        access_code = get_cif_value("pdb_entry_id", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                        pdb_pattern = re.compile(r"\d[\dA-Za-z]{3}|pdb_\d{5}[\dA-Za-z]{3}")
                                        if db_name == "PDB":
                                            if access_code is not None:
                                                if not pdb_pattern.match(str(access_code)):
                                                    txt = u"(%s) PDB id is not in the correct format" % access_code
                                                    self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                                                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)
                                        if accession_code is None and access_code is not None:
                                            set_cif_value(model.set_access_code, "pdb_entry_id", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                        if accession_code is not None and access_code is None:
                                            set_cif_value(model.set_access_code, "accession_code", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                        if accession_code is not None and access_code is not None:
                                            set_cif_value(model.set_access_code, "pdb_entry_id", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            if accession_code != access_code:
                                                txt = u"Cannot be two access_codes. If both pdb_entry_id and accession_code are populated, then both should be same."
                                                self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                                                self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                                    def set_el_chain(model, model_in):
                                        """
                                        XSD: <xs:element name="chain" maxOccurs="unbounded" minOccurs="0">
                                        .. a base chain_type and
                                        .. 1 element
                                        """
                                        # def set_chain_type(chain, model_in):
                                        #     """
                                        #     XSD: <xs:complexType name="chain_type"> has
                                        #     .. 2 elements
                                        #     """
                                        def set_el_id(chain, model_in):
                                            """
                                            XSD: <xs:element name="id" type="chain_pdb_id" minOccurs="0" maxOccurs="unbounded"/>
                                            CIF: _em_3d_fitting_list.pdb_chain_id A
                                            """
                                            ids_in = get_cif_value("pdb_chain_id", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            ch_ids = get_cif_value("chain_id", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            if ids_in is not None and ch_ids is None:
                                                ids = ids_in.split(",")
                                                for a_id in ids:
                                                    set_cif_value(chain.set_chain_id, "pdb_chain_id", const.EM_3D_FITTING_LIST, cif_list=model_in, cif_value=a_id)
                                            if ids_in is None and ch_ids is not None:
                                                c_ids = ch_ids.split(",")
                                                for b_id in c_ids:
                                                    set_cif_value(chain.set_chain_id, "chain_id", const.EM_3D_FITTING_LIST, cif_list=model_in, cif_value=b_id)
                                            if ids_in is not None and ch_ids is not None:
                                                ids = ids_in.split(",")
                                                for a_id in ids:
                                                    set_cif_value(chain.set_chain_id, "pdb_chain_id", const.EM_3D_FITTING_LIST, cif_list=model_in, cif_value=a_id)
                                                if ids_in != ch_ids:
                                                    txt = u"Cannot be two chain_ids. If both pdb_chain_id and chain_id are populated, then both should be same."
                                                    self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                                                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                                        def set_el_residue_range(chain, model_in):
                                            """
                                            XSD: <xs:element name="residue_range" maxOccurs="unbounded" minOccurs="0">
                                            CIF: _em_3d_fitting_list.pdb_chain_residue_range 5-545
                                            """
                                            p_res_range = get_cif_value("pdb_chain_residue_range", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            res_range = get_cif_value("chain_residue_range", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            if p_res_range is not None and res_range is None:
                                                set_cif_value(chain.set_residue_range, "pdb_chain_residue_range", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            if p_res_range is None and res_range is not None:
                                                set_cif_value(chain.set_residue_range, "chain_residue_range", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            if p_res_range is not None and res_range is not None:
                                                set_cif_value(chain.set_residue_range, "pdb_chain_residue_range", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                                if p_res_range != res_range:
                                                    txt = u"Cannot be two chain_residue_range. If both pdb_chain_residue_range and chain_residue_range are populated, then both should be same."
                                                    self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                                                    self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                                        def set_el_numofcps_in_final_model(chain, model_in):
                                            """
                                            XSD: <xs:element name="number_of_copies_in_final_model" type="xs:positiveInteger"  minOccurs="0">
                                            CIF: ?????? How do you calculate this? It's not in mmCIf
                                            """
                                            set_cif_value(chain.set_number_of_copies_in_final_model, "number_of_copies_in_final_model", const.EM_3D_FITTING_LIST, cif_list=model_in)

                                        def set_el_source_name(chain, model_in):
                                            """
                                            XSD: <xs:element name="source_name" type="xs:string"  minOccurs="0" maxOccurs="1">
                                            CIF: _em_3d_fitting_list.source_name SwissModel
                                            """
                                            access_code = get_cif_value("pdb_entry_id", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            accession_code = get_cif_value("accession_code", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            sourcename = get_cif_value("source_name", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            set_cif_value(chain.set_source_name, "source_name", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                            if sourcename == "PDB" and access_code is None and accession_code is None:
                                                txt = u"Error! Missing PDB ID. If initial model is from PDB, then access code is mandatory."
                                                self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
                                                self.log_formatted(self.error_log_string, const.REQUIRED_ALERT + txt)

                                        def set_el_initial_model_type(chain, model_in):
                                            """
                                            XSD: <xs:element name="initial_model_type" type="xs:string"  minOccurs="0" maxOccurs="1">
                                            CIF: _em_3d_fitting_list.type experimental model
                                            """
                                            set_cif_value(chain.set_initial_model_type, "type", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                        # element 2
                                        chain = emdb.chainType()
                                        set_el_id(chain, model_in)
                                        set_el_residue_range(chain, model_in)
                                        set_el_numofcps_in_final_model(chain, model_in)
                                        set_el_source_name(chain, model_in)
                                        set_el_initial_model_type(chain, model_in)
                                        if chain.has__content():
                                            model.set_chain(chain)

                                    def set_el_details(model, model_in):
                                        """
                                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                        CIF: _em_3d_fitting_list.details
                                        """
                                        set_cif_value(model.set_details, "details", const.EM_3D_FITTING_LIST, cif_list=model_in)
                                    # element 1
                                    set_el_access_code(model, model_in)
                                    # element 2
                                    set_el_chain(model, model_in)
                                    # element 3
                                    set_el_details(model, model_in)
                                if modelling_id_in in model_dict_in:
                                    for model_in in model_dict_in[modelling_id_in]:
                                        model = emdb.initial_modelType()
                                        set_initial_model_type(model, model_in)
                                        if model.has__content():
                                            modelling.add_initial_model(model)

                            def set_el_final_model():
                                """
                                XSD: <xs:element name="final_model" minOccurs="0"> has
                                .. 3 elements
                                CIF: there is no _emd_modelling_final_model !!!
                                """

                            def set_el_refinement_protocol(modelling, modelling_in):
                                """
                                XSD: <xs:element name="refinement_protocol" minOccurs="0">
                                CIF: _em_3d_fitting.ref_protocol {AB INITIO MODEL,...}
                                """
                                set_cif_value(modelling.set_refinement_protocol, "ref_protocol", const.EM_3D_FITTING, cif_list=modelling_in)

                            def set_el_software_list(modelling, cat_soft_dict_in):
                                """
                                XSD: <xs:element name="software_list" type="software_list_type" minOccurs="0"/>
                                """
                                set_software_list(const.SOFT_MODEL_FITTING, cat_soft_dict_in, modelling.set_software_list)

                            def set_el_details(modelling, modelling_in):
                                """
                                XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                CIF: _em_3d_fitting.details
                                """
                                set_cif_value(modelling.set_details, "details", const.EM_3D_FITTING, cif_list=modelling_in)

                            def set_el_target_criteria(modelling, modelling_in):
                                """
                                XSD: <xs:element name="target_criteria" type="xs:token" minOccurs="0">
                                CIF: _em_3d_fitting.target_criteria 'Correlation coefficient'
                                """
                                set_cif_value(modelling.set_target_criteria, "target_criteria", const.EM_3D_FITTING, cif_list=modelling_in)

                            def set_el_refinement_space(modelling, modelling_in):
                                """
                                XSD: <xs:element name="refinement_space" type="xs:token" minOccurs="0">
                                CIF: _em_3d_fitting.ref_space {REAL,RECIPROCAL}
                                """
                                set_cif_value(modelling.set_refinement_space, "ref_space", const.EM_3D_FITTING, cif_list=modelling_in)

                            def set_el_overall_bvalue(modelling, modelling_in):
                                """
                                XSD: <xs:element name="overall_bvalue" type="xs:float" minOccurs="0">
                                CIF: _em_3d_fitting.overall_b_value 200
                                """
                                set_cif_value(modelling.set_overall_bvalue, "overall_b_value", const.EM_3D_FITTING, cif_list=modelling_in, fmt=float)

                            # element 1
                            set_el_initial_model(modelling, modelling_id_in, model_dict_in)
                            # element 2
                            set_el_final_model()
                            # element 3
                            set_el_refinement_protocol(modelling, modelling_in)
                            # element 4
                            set_el_software_list(modelling, cat_soft_dict_in)
                            # element 5
                            set_el_details(modelling, modelling_in)
                            # element 6
                            set_el_target_criteria(modelling, modelling_in)
                            # element 7
                            set_el_refinement_space(modelling, modelling_in)
                            # element 8
                            set_el_overall_bvalue(modelling, modelling_in)

                        # Create a dictionary of software keyed by category
                        cat_soft_dict_in = {}
                        ip_list_in = assert_get_value(const.EM_IMAGE_PROCESSING, self.cif)
                        for ip_in in ip_list_in:
                            ip_id_in = get_cif_value(const.K_ID, const.EM_IMAGE_PROCESSING, ip_in)
                            soft_dict_in = make_list_of_dicts(const.EM_SOFTWARE, const.K_IMAGE_PROCESSING_ID)
                            if ip_id_in in soft_dict_in:
                                cat_soft_dict_in = make_list_of_dicts(const.EM_SOFTWARE, "category", soft_dict_in[ip_id_in])

                        modelling_list = emdb.modelling_listType()
                        modelling_list_in = self.cif.get(const.EM_3D_FITTING)
                        model_dict_in = make_list_of_dicts(const.EM_3D_FITTING_LIST, const.K_3D_FITTING_ID)
                        for modelling_in in modelling_list_in:
                            modelling_id_in = get_cif_value(const.K_ID, const.EM_3D_FITTING, modelling_in)
                            modelling = emdb.modelling_type()
                            set_modelling_type(modelling, modelling_in, cat_soft_dict_in)
                            if modelling.has__content():
                                modelling_list.add_modelling(modelling)
                        if modelling_list.has__content():
                            intrp.set_modelling_list(modelling_list)

                    def set_el_figure_list():
                        """
                        XSD: <xs:element name="figure_list" minOccurs="0">
                        CIF: ??
                        """

                    def set_el_segmentation_list(intrp):
                        """
                        XSD: <xs:element name="segmentation_list" minOccurs="0"> has
                        .. 1 element
                        """

                        def set_segmentation_type(seg, msk):
                            """
                            XSD: <xs:element name="segmentation" maxOccurs="unbounded"> has
                            .. 3 elements
                            """

                            def set_el_file(seg, msk):
                                """
                                XSD: <xs:element name="file">
                                pattern [emd_d{4,}]+.*
                                """
                                msk_file = msk.get_file()
                                set_cif_value(seg.set_file, cif_value=msk_file)

                            def set_el_details():
                                """
                                XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                                CIF: ????
                                """

                            def set_el_mask_details():
                                """
                                XSD: <xs:element name="mask_details" type="map_type" minOccurs="0">
                                Deprecated (2014-10-21)
                                """

                            # element 1
                            set_el_file(seg, msk)
                            # element 2
                            set_el_details()
                            # element 3
                            set_el_mask_details()

                        mask_list_in = map_dict_in[const.MAP_MASK] if const.MAP_MASK in map_dict_in else []
                        seg_list = emdb.segmentation_listType()
                        for msk_in in mask_list_in:
                            msk = make_map(msk_in, get_canonical_map_name(msk_in, "MASK"))
                            if msk.has__content():
                                seg = emdb.segmentationType()
                                set_segmentation_type(seg, msk)
                                seg_list.add_segmentation(seg)
                        if seg_list.has__content():
                            intrp.set_segmentation_list(seg_list)

                    def set_el_slices_list():
                        """
                        XSD: <xs:element name="slices_list" minOccurs="0">
                        CIF: Deprecated (2014-10-21)
                        """

                    def set_el_additional_map_list(intrp, map_dict_in):
                        """
                        XSD: <xs:element name="additional_map_list" minOccurs="0">
                        """
                        add_map_list_in = map_dict_in[const.MAP_ADD] if const.MAP_ADD in map_dict_in else []
                        add_map_list = emdb.additional_map_listType()
                        for add_map_in in add_map_list_in:
                            add_map = make_map(add_map_in, get_canonical_map_name(add_map_in, "ADDMAP"))
                            if add_map.has__content():
                                add_map_list.add_additional_map(add_map)
                        if add_map_list.has__content():
                            intrp.set_additional_map_list(add_map_list)

                    def set_el_half_map_list(intrp, map_dict_in):
                        """
                        XSD: <xs:element name="half_map_list" minOccurs="0">
                        """
                        hf_map_list_in = map_dict_in[const.MAP_HALF] if const.MAP_HALF in map_dict_in else []
                        hf_map_list = emdb.half_map_listType()
                        for hf_map_in in hf_map_list_in:
                            hf_map = make_map(hf_map_in, get_canonical_map_name(hf_map_in, "HALFMAP"))
                            if hf_map.has__content():
                                hf_map_list.add_half_map(hf_map)
                        if hf_map_list.has__content():
                            intrp.set_half_map_list(hf_map_list)

                    # element 1
                    set_el_modelling_list(intrp)
                    # element 2
                    set_el_figure_list()
                    # element 3
                    set_el_segmentation_list(intrp)
                    # element 4
                    set_el_slices_list()
                    # element 5
                    set_el_additional_map_list(intrp, map_dict_in)
                    # element 6
                    set_el_half_map_list(intrp, map_dict_in)

                map_dict_in = make_list_of_dicts(const.EM_MAP, "type")
                intrp = emdb.interpretation_type()
                set_interpretation_type(intrp, map_dict_in)
                if intrp.has__content():
                    self.xml_out.set_interpretation(intrp)

            def set_el_validation():
                """
                XSD: <xs:element name="validation" minOccurs="0">
                There should be 4 validation methods:
                1. crystallography validation
                2. fsc curve
                3. layer lines
                4. structure factors
                Only fsc curve is implemented here???
                """
                def set_fsc_curve(fsc, fsc_in):
                    """
                    XSD: <xs:complexType name="fsc_curve_validation_type"> has
                    .. a base validation_type only
                    XSD: <xs:complexType name="validation_type"> has
                    .. 2 elements with pattern for a fsc curve only!!!???
                    """
                    def set_el_file(fsc, fsc_in):
                        """
                        XSD: <xs:element name="file">
                        CIF: _em_fsc_curve.file
                        """
                        set_cif_value(fsc.set_file, "file", const.EM_FSC_CURVE, cif_list=fsc_in)

                    def set_el_details(fsc, fsc_in):
                        """
                        XSD: <xs:element name="details" type="xs:string" minOccurs="0"/>
                        CIF: _em_fsc_curve.details
                        """
                        set_cif_value(fsc.set_details, "details", const.EM_FSC_CURVE, cif_list=fsc_in)

                    # element 1
                    set_el_file(fsc, fsc_in)
                    # element 2
                    set_el_details(fsc, fsc_in)

                validation = emdb.validationType()
                fsc_list_in = self.cif.get(const.EM_FSC_CURVE, None)
                for fsc_in in fsc_list_in:
                    fsc = emdb.fsc_curve_validation_type()
                    fsc.original_tagname_ = "fsc_curve"
                    set_fsc_curve(fsc, fsc_in)
                    if fsc.has__content():
                        validation.add_validation_method(fsc)
                if validation.has__content():
                    self.xml_out.set_validation(validation)

            # attribute 1
            set_attr_emdb_id()
            # attribute 2
            set_attr_version()
            # element 1
            set_el_admin()
            # element 2
            set_el_crossreferences()
            # element 3
            set_el_sample()
            # element 4
            set_el_struct_det_list()
            # element 5
            set_el_map()
            # element 6
            set_el_interpretation()
            # element 7
            set_el_validation()

        def get_entry_id_from_input_file():
            cif_re = re.compile(r"EMD\-([0-9]){4,}")
            search_result = re.search(cif_re, self.cif_file_name)
            if search_result is not None:
                return search_result.group()
            else:
                return ""

        if self.cif_file_read:
            const = self.Constants
            self.entry_in_translation_log = self.EntryLog(get_entry_id_from_input_file())
            self.create_xml = True
            self.xml_out = emdb.entry_type()
            set_entry_type()
            self.translation_log.logs.append(self.entry_in_translation_log)
        else:
            print("cif file NOT read")
            # TODO  # pylint: disable=fixme
            # txt = u'Translation cannot be performed. The cif file (%s) cannot be read.' % self.cif_file_name
            # self.current_entry_log.error_logs.append(self.ALog(log_text='(' + self.entry_in_translation_log.id + ')' +txt))
            # self.log_formatted(self.error_log_string, self.Constants.REQUIRED_ALERT + txt)

    def translate(self, in_cif, out_xml):
        """
        Reads input cif file, translates it and writes out XML
        """
        self.read_cif_in_file(in_cif)
        self.translate_cif_to_xml()
        if self.create_xml:
            self.write_xml_out_file(out_xml)
        else:
            txt = u"The XML output file (%s) has not been created. See previous messages for the reason." % out_xml
            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.error_title + txt))
            self.log_formatted(self.error_log_string, self.Constants.REQUIRED_ALERT + txt)

    def validate_file(self, the_parser, xml_filename):
        """
        Method to validate any schema against any file
        """
        try:
            # Python3 reqires a byte string for lxml as encoding in file
            xml_file = open(xml_filename, "rb")
            try:
                etree.fromstring(xml_file.read(), the_parser)
            except etree.XMLSyntaxError:
                return False
            except etree.XMLSchemaError:
                return False
            return True
        except IOError as exp:
            txt = u"Error (%s) occured. Arguments (%s)." % (str(exp), exp.args)
            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.validation_title + txt))
            self.log_formatted(self.error_log_string, self.Constants.VALIDATION_ERROR + txt)
            return False
        finally:
            xml_file.close()

    def show_validation_errors(self, in_xml, in_schema_filename):
        """
        Called if the validation of the in_xml file against
        the schema in_schema_filename fails. Shows the list of validation errors.
        """
        try:
            xml_doc = etree.parse(in_xml)
            xsd = etree.parse(in_schema_filename)
            xml_schema = etree.XMLSchema(xsd)
            xml_schema.assertValid(xml_doc)
            txt = u"File %s validates." % in_xml
            self.current_entry_log.info_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + txt))
            self.log_formatted(self.info_log_string, txt)
        except etree.XMLSyntaxError as exp:
            txt = u"PARSING ERROR: %s." % exp
            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + txt))
            self.log_formatted(self.error_log_string, txt)
        except etree.DocumentInvalid as exp:
            i = 1
            for err in exp.error_log:  # pylint: disable=no-member
                txt = u"%d: %s" % (i, err)
                self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.validation_title + txt))
                self.log_formatted(self.error_log_string, self.Constants.VALIDATION_ERROR + txt)
                i = i + 1

    def validation_logger_header(self, log_str, in_schema_filename):  # pylint: disable=unused-argument
        """
        Makes the validation logger pretty
        """
        if self.logger is not None and log_str is not None:
            if log_str.closed:
                if log_str == self.info_log_string:
                    self.info_log_string = self.open_log_stream(self.Constants.INFO_LOG_STRING)
                    log_str = self.info_log_string
                if log_str == self.warn_log_string:
                    self.warn_log_string = self.open_log_stream(self.Constants.WARN_LOG_STRING)
                    log_str = self.warn_log_string
                if log_str == self.error_log_string:
                    self.error_log_string = self.open_log_stream(self.Constants.ERROR_LOG_STRING)
                    log_str = self.error_log_string
            log_str.write(u"\n")

    def validate(self, in_xml, in_schema_filename):
        """
        Validate in_xml against in_schema
        """
        try:
            # For python3, as encoding in file, lxml requires bytes string
            in_schema = open(in_schema_filename, "rb")
        except IOError as exp:
            txt = u"Error (%s) occured. Arguments (%s)." % (str(exp), exp.args)
            self.current_entry_log.error_logs.append(self.ALog(log_text="(" + self.entry_in_translation_log.id + ")" + self.current_entry_log.validation_title + txt))
            self.log_formatted(self.error_log_string, self.Constants.VALIDATION_ERROR + txt)
            return False
        else:
            schema_doc = in_schema.read()
            schema_root = etree.XML(schema_doc)
            the_schema = etree.XMLSchema(schema_root)

            xml_parser = etree.XMLParser(schema=the_schema)
            validates = self.validate_file(xml_parser, in_xml)
            if not validates:
                self.validation_logger_header(self.error_log_string, in_schema_filename)
                self.show_validation_errors(in_xml, in_schema_filename)
            return validates
        finally:
            in_schema.close()

    def translate_and_validate(self, in_cif, out_xml, in_schema=None):
        """
        Validates XML file after reading an input cif file, translating it and creating XML
        """
        if in_schema is None:
            es = EmdbSchema()
            in_schema = es.getSchemaPath()

        self.translate(in_cif, out_xml)
        if self.create_xml:
            self.validate(out_xml, in_schema)


def main():
    """Translate a cif file to EMDB xml 3.x.x.x"""
    # Handle command line options
    usage = """
            cif_emdb_translator.py [options]
            Convert a cif file to EMDB xml 3.x.x.x or vice versa

            Examples:
            python cif_emdb_translator.py --xml -i input_file -o outputFile

            Typical run:
            python cif_emdb_translator.py --xml -i in.cif -o out.xml
                in.cif is assumed to be a cif file in emd space
                an XML file following EMDB XML schema 3.0 is produced and written out to out.xml

            Enabling console logging information:
            python cif_emdb_translator.py --xml -l -i in.cif -o out.xml

            Turn off logging files:
            python cif_emdb_translator.py --xml --info-log-off --warn-log-off --error-log-off -i in.cif -o out.xml
            or
            python cif_emdb_translator.py --xml -f -w -e -i in.cif -o out.xml
                -f or --info-log-off: no information written into INFO_cifEMDBTranslation.log
                -w or --warn-log-off: no information written into WARN_cifEMDBTranslation.log
                -e or --error-log-off: no information written into ERROR_cifEMDBTranslation.log

            Combine turning off log files and enabling console log:
            python cif_emdb_translator.py --xml -w -e -i in.cif -o out.xml
            """

    logging.info("*** cif_emdb_translator.py version %s ***", __version__)
    parser = OptionParser(usage=usage, version=__version__)
    parser.add_option("-o", "--out-file", action="store", type="string", metavar="FILE", dest="outputFile", help="Write output to FILE")
    parser.add_option("-i", "--in-file", action="store", type="string", metavar="FILE", dest="inputFile", help="Write input to FILE")
    parser.add_option("-x", "--xml", action="store_true", dest="translate2Xml", help="Translate cif file to XML", default=True)
    parser.add_option("-c", "--cif", action="store_false", dest="translate2Xml", help="Translate XML file to cif")
    parser.add_option("-f", "--info-log", action="store_false", dest="info_log", help="Logging to INFO file flag")
    parser.add_option("-w", "--warn-log", action="store_false", dest="warn_log", help="Logging to WARN file flag")
    parser.add_option("-e", "--error-log", action="store_false", dest="err_log", help="Logging to ERROR file flag")
    parser.add_option("-l", "--console-log", action="store_false", dest="console_log", help="Logging to console turned on")
    parser.add_option("-p", "--private-include", action="store_false", dest="private_inc", help="Private elements included in the xml output file")
    parser.add_option("-s", "--show-log-id", action="store_false", dest="show_log_ids", help="Entry ids are printed into logs")

    (options, args) = parser.parse_args()

    # Check for sensible/supported options
    if len(args) == 0:
        print(usage)
        sys.exit("No input options given!")
    #     else:
    #         input_file = args[0]

    if options.translate2Xml is True:
        translator = CifEMDBTranslator()
        translator.set_logger_logging(options.info_log, options.warn_log, options.err_log, options.console_log)
        translator.set_show_private(options.private_inc)
        translator.set_show_log_id(options.show_log_ids)
        translator.read_cif_in_file(options.inputFile)
        translator.translate_cif_to_xml()
        translator.write_xml_out_file(options.outputFile)
    else:
        logging.info("XML to CIF translation not implemented yet. Coming sooooon :)")


if __name__ == "__main__":
    main()
