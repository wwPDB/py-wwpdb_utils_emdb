
import os
import re
import emdb_xml2cif_translator.input_files


class Mappings(object):
    """
    This class contains tools and containers for serving mappings described in the input file MAPPINGS_FILENAME.

    Its tools manage the load of the mappings from MAPPINGS_FILENAME to populate a dictionary with the "mappings" logic

    The logic is used when preparing a container (dictionary) of values for the CIF output

    """
    class Const(object):
        """

        """
        # input text file containing the mapping logic from any EMDB v3.x XML header file into an _emd mmcif file
        MAPPINGS_LOGIC_FILENAME = 'emdb-xml2cif-mappings.txt'

        # mappings logic dictionary constants
        XML_MAPPING = 'xml_mapping'
        XML_VALUE = 'xml_value'
        LOGIC = 'logic'
        CIF_MAPPINGS = 'cif_mappings'
        CATEGORY_ID = 'category_id'
        ITEMS = 'items'
        DATA = 'data'

        # mappings logic mappings keys
        MAP_EMD_EMDB_ID = 'emd@emdb_id'

        # other constants
        XML_VALUE_UPPER = 'XML_VALUE'
        ID = 'ID'
        N = 'N'
        B = 'B'
        D = 'D'
        M = 'M'
        E = 'E'
        U = 'U'
        S = 'S'
        I = "I"

    def __init__(self):
        """
        Initialises (xml to cif) mappings as an empty dictionary in the following format:
        self.mappings =
        {
            xml mapping 1: {
                'xml_value': '',
                'logic': {
                    multiples: False,
                    cif category: {
                        'items_logic': [],
                        'data_logic': []
                    }
                },
                'cif_mappings': {
                    cif category: {
                        'items': [],
                        'data': []
                    }
                    ...
                }
            }

            ...
            xml mapping N: {

            }
        }
        where
            a. xml mapping (1..N) corresponds to
            b. cif category
        a. and b. are found in one mapping (a line of text) from MAPPINGS_LOGIC_FILENAME (one whole line of text)
        e.g. emd.admin.title emd_admin.title.XML_VALUE

        Populates the "mappings" logic dictionary with a call to the class method load_mappings()
        """
        # Mapping logic is a dictionary of dictionaries (one for each mapping)
        self.mappings = {}
        # read the input text file with the mapping logic and populate the "mappings" logic dictionary
        self.load_mappings()

    def load_mappings(self):
        """
        This method opens and reads the text file containing the xml to mmcif emd mappings
        located in the input_files project folder
        :param self:
        :return loaded: a boolean; True if the input text file can be opened and is read
        """
        loaded = False
        mappings_file = os.path.join(os.path.dirname(emdb_xml2cif_translator.input_files.__file__),
                                     self.Const.MAPPINGS_LOGIC_FILENAME)
        if os.path.isfile(mappings_file):
            f = open(mappings_file, 'r').read().split("\n")
            for line in f:
                if len(line) != 0:
                    if line[0] != '#':
                        self.read_one_mapping(line)

            loaded = True

        return loaded

    def read_one_mapping(self, line):
        """
        Helper method.
        Interprets one line from MAPPINGS_LOGIC_FILENAME and stores the logic and values into mappings logic dictionary
        :param line: a string; one line of text from the MAPPINGS_LOGIC_FILENAME file; it's not an empty string;
                    contains the logic of how one xml attribute or element is to be represented in the mmcif file;
                    xml mapping to cif mapping is a one-to-many relationship;
                    mappings are separated by a space;
                    e.g. 'emd.admin.title emd_admin.title.XML_VALUE'
                        where
                        'emd.admin.title' is the representation of the title sub-element within the XML input file
                        and
                        'emd_admin.title.XML_VALUE' represents a notation where
                        the item 'title' will be written within the _emd_admin category and
                        its value is what is read from 'emd.admin.title'
        :return read: a boolean; True if the mapping logic given in line is read
        """
        read = False
        # flag noting if the XML value is coming from a list
        is_list = False
        # flag noting if the XML value is coming from a dictionary
        is_dict = False
        # a flag for multiple instances
        is_multiple = False

        # check if XML mapping is coming from a list
        if line:
            is_multiple = None
            if line[0] in ['N']:
                is_list = True
            elif line[0] in ['D']:
                is_dict = True

        if is_list or is_dict or is_multiple:
            # remove the leading anchor (e.g. 'N:') from the xml mapping
            line = line[2:]

        # save the mappings into a list
        line_mappings = line.rstrip().split(' ')

        xml_maps= []
        # the first mapping is always a xml mapping
        if 'S$' in line_mappings[0]:
            substitution = re.split('S\$|\$S', line_mappings[0])
            count = len(substitution)
            sub_ele = substitution[1].split("|")
            for sub in sub_ele:
                if count == 2:
                    xml_map = substitution[0]+sub
                    xml_maps.append(xml_map)
                if count == 3:
                    xml_map = substitution[0]+sub+substitution[2]
                    xml_maps.append(xml_map)
        else:
            xml_maps.append(line_mappings[0])


        mappings = {}
        if xml_maps:
            for xml_mapp in xml_maps:
                # xml mapping is unique; use it as a dict key
                if "_supramolecule" in xml_mapp:
                    xml_mapping = re.sub(r'(cell_supramolecule|complex_supramolecule|organelle_or_cellular_component_supramolecule|sample_supramolecule|tissue_supramolecule|virus_supramolecule)', 'all_supramolecules', xml_mapp)
                    xml_mapping = re.sub(r'sci_species_name', 'natural_source', xml_mapping)
                elif "macromolecule_list" in xml_mapp:
                    xml_mapping = re.sub(r'(protein_or_peptide|em_label|ligand|other_macromolecule|dna|rna|saccharide)', 'all_macromolecules', xml_mapp)
                elif ("map" in xml_mapp or "mask_details" in xml_mapp) and not "map_release" in xml_mapp:
                    if "&" in xml_mapp:
                        rep_item = xml_mapp.rsplit(".", 1)[1].split("&", 1)[1]
                    elif "$I" in xml_mapp:
                        rep_item = xml_mapp.rsplit(".", 1)[1].split("$I", 1)[1]+".$I"
                    else:
                        rep_item = "_".join(xml_mapp.rsplit(".", 2)[1:])
                        if xml_mapp.rsplit('.', 2)[1] in ['map', 'mask_details', 'half_map', 'additional_map']:
                            rep_item = xml_mapp.rsplit(".", 1)[1]
                    xml_mapping = "emd.all_maps."+rep_item
                    if "H$" in xml_mapp:
                        xml_mapping = xml_mapp
                else:
                    xml_mapping = xml_mapp
                # create this line's mapping now
                if is_list:
                    mapping = {
                        xml_mapping: {
                            self.Const.XML_VALUE: [],
                            self.Const.LOGIC: {},
                            self.Const.CIF_MAPPINGS: {}}
                    }
                    if mapping:
                        mappings.update(mapping)
                elif is_dict:
                    mapping = {
                        xml_mapping: {
                            self.Const.XML_VALUE: {},
                            self.Const.LOGIC: {},
                            self.Const.CIF_MAPPINGS: {}
                        }
                    }
                    if mapping:
                        mappings.update(mapping)
                else:
                    mapping = {
                        xml_mapping: {
                            self.Const.XML_VALUE: '',
                            self.Const.LOGIC: {},
                            self.Const.CIF_MAPPINGS: {}
                        }
                    }
                    if mapping:
                        mappings.update(mapping)

        # the rest of the mappings are cif mappings
        if len(line_mappings) > 1:
            # there is at least one cif mapping
            logic = self.read_logic(line_mappings[1:], is_multiple)
            if logic:
                for map in mappings:
                    mappings.get(map).get(self.Const.LOGIC).update(logic)
                    # both xml mapping and cif mappings are now read; set the return flag
                    read = True
        # the line information is read and its mapping logic dictionary is created,
        # add it to the dictionary containing all mappings logic
        if mappings:
            self.mappings.update(mappings)

        return read

    def read_logic(self, input_all_logic, is_multiple):
        """
        Helper method.
        Reads all cif mappings given in one line within the text input file containing all xml to cif mappings logic

        :param is_multiple: A flag; If True one cif item will have multiple values
        :param input_all_logic: a string representing cif mapping as 'cif category'.'cif item 1'.'value to use'
        :return cif_mappings: a dictionary;
            contains dictionaries (one for each cif category found in the cif mappings)

                            logic = {
                                'cif category': {
                                    self.ITEMS: ['cif item 1'],
                                    self.DATA: ['value to use']
                                }
                            }
            e.g. for cases where one xml value writes to two different cif categories:
            'emd@emdb_id database_2.database_id.EMDB database_2.database_code.XML_VALUE' the dictionary is:

                    logic = {
                                database_2: {
                                    self.ITEMS: [database_id, database_code],
                                    self.DATA: [EMDB, XML_VALUE]
                                }
                            }
            where XML_VALUE is the value read from emd@emdb_id
            (the value given for 'emdb_id attribute in the 'emd' element)

            e.g. for case where the cif mappings are different cif categories
            'emd.admin.sites.deposition emd_admin.deposition_site.XML_VALUE pdbx_database_status.deposit_site.XML_VALUE'
                    logic = {
                                emd_admin: {
                                    self.ITEMS: [deposition_site],
                                    self.DATA: [XML_VALUE]
                                },
                                pdbx_database_status: {
                                    self.ITEMS: [deposit_site],
                                    self.DATA: [XML_VALUE]
                                }
                            }
                            TODO: Add example for a dictionary here


            For multiples:
                    logic = {
                                'emd.crossreferences.emdb_list.emdb_reference.emdb_id': {
                                    self.ITEMS: ['pdbx_database_related.db_name', 'pdbx_database_related.db_id'],
                                    self.DATA: [['EMDB'], [XML_VALUE]]
                                }
                            }

        """
        all_logic = {}
        for input_logic in input_all_logic:
            items_ids = []
            items_values = []
            logic = {}

            # get the cif category and its items; here referred as cif components
            cif_components = input_logic.split('.')
            if len(cif_components) > 0:
                # the first cif component is the cif category
                category_id = cif_components[0]
                if len(cif_components) > 1:
                    # the second cif component is the cif item
                    items_ids.append(cif_components[1])
                    if len(cif_components) > 2:
                        items_values.append(cif_components[2])

                # create dictionary to hold the cif mappings
                if category_id not in all_logic.keys():
                    logic = {
                        category_id: {
                            self.Const.ITEMS: [],
                            self.Const.DATA: []
                        }
                    }

                all_logic.update(logic)
                # add read items and values
                all_logic.get(category_id).get(self.Const.ITEMS).extend(items_ids)
                all_logic.get(category_id).get(self.Const.DATA).extend(items_values)

        return all_logic

    def get_mapping_logic_value(self, mapping_logic, mapping_logic_key):
        """
        Method to be called by other classes. Providing a key for a mapping logic, its value is returned
        :param mapping_logic: a constant defined in this class;
                              a value from the first column from the input text file containing log.
                              e.g. MAP_EMD_EMDB_ID
        :param mapping_logic_key: a constant defined in this class and contained in each mapping logic;
                                  can be {XML_MAPPING, XML_VALUE, EXTRAS, CIF_MAPPINGS}
        :return: The value for mapping_logic_key within mapping_logic; None if mapping_logic doesn't exist
        """
        # get the dictionary containing mapping for requested mapping logic
        mapping_logic_dict = self.mappings.get(mapping_logic)
        if mapping_logic_dict:
            # the mapping logic requested exists in the mappings dictionary;
            # return the value for the mapping logic key
            return mapping_logic_dict.get(mapping_logic_key)
        else:
            # mapping logic is not in the mappings logic dictionary; None to return
            return None

    def set_cif_mappings_values(self):
        """
        This method sets the values for each category item in the mappings logic
        :return values_set: a boolean; True when the xml values are set into all mappings
        """
        values_set = False
        for mapping_key, mapping in self.mappings.items():
            xml_value = mapping.get(self.Const.XML_VALUE)
            all_logic = mapping.get(self.Const.LOGIC)
            if all_logic:
                for a_logic_key, a_logic_value in all_logic.items():
                    logic_keys = a_logic_value.get(self.Const.ITEMS)
                    logic_values = a_logic_value.get(self.Const.DATA)
                    logic_k, logic_v = '', ''
                    list_values, source_val = [], []
                    list_item = []
                    for logic_value in logic_values:
                        if '$' in logic_value:
                            logic_k, logic_v = [str(i) for i in logic_value.split("$")]
                        if logic_value == self.Const.XML_VALUE_UPPER:
                            list_values = list(map(lambda x: x if x != logic_value else xml_value, logic_values))
                            if list_values[0]:
                                if any('FULLOVERLAP' in sublist or 'unknown' in sublist for sublist in list_values):
                                    replacement_mapping = {'FULLOVERLAP': 'IN FRAME', 'unknown': 'OTHER'}
                                    list_values = [[replacement_mapping.get(value, value) for value in sublist] for sublist in list_values]
                                if 'DOI' in [item for sublist in list_values for item in sublist]:
                                    doi_value = next((item['DOI'] for item in list_values if 'DOI' in item), None)
                                    new_doi = doi_value.replace('doi:', '')
                                    [entry.update({'DOI': new_doi}) for entry in list_values if 'DOI' in entry]
                                if any("IMAGE STORED AS" in element for sublist in list_values for element in sublist):
                                    list_values = [[original_string.capitalize() for original_string in sublist] for sublist in list_values]
                        if logic_value == self.Const.U:
                            if xml_value == "twoDArray":
                                xml_value = "2D ARRAY"
                            if xml_value == "threeDArray":
                                xml_value = "3D ARRAY"
                            if xml_value == "helicalArray":
                                xml_value = "HELICAL ARRAY"
                            if xml_value == "singleParticle":
                                xml_value = "SINGLE PARTICLE"
                            if xml_value == "subtomogramAveraging":
                                xml_value = "SUBTOMOGRAM AVERAGING"
                            if xml_value == "electronCrystallography":
                                xml_value = "CRYSTALLOGRAPHY"
                            else:
                                xml_value = xml_value
                            list_values.append(xml_value.upper())
                        if logic_value == self.Const.N:
                            n_values = [xml_value.index(i) + 1 for i in xml_value]
                            list_values.append(n_values)
                            break
                        if logic_value == self.Const.B:
                            if xml_value == "true":
                                list_values.append("N")
                            elif xml_value == "false":
                                list_values.append("Y")
                        if logic_value == self.Const.S:
                            for val in xml_value:
                                if val == '':
                                    source_val.append('')
                                if val == "nat":
                                    source_val.append("NATURAL")
                                if val == "man":
                                    source_val.append("RECOMBINANT")
                                if val == "syn":
                                    source_val.append("SYNTHETIC")
                            list_values = [source_val]
                        if logic_value == self.Const.E:
                            modified_names = [f"{n.split(' ', 1)[0]}, {'.'.join(n.split(' ', 1)[1])}." if ' ' in n else n for n in xml_value]
                            list_values.append(modified_names)
                        if logic_value == self.Const.I:
                            new_list = [re.search(r'\d+', item).group(0).zfill(4) for item in xml_value]
                            new_list = ['EMD-' + item for item in new_list]
                            list_values.remove("I")
                            list_values.append(new_list)
                        elif logic_k == self.Const.M:
                            if self.Const.XML_VALUE_UPPER in logic_value:
                                list_item.append(xml_value)
                            else:
                                list_item.append([logic_v]*(len(xml_value)))
                            list_values = list_item

                    cif_mapping = self.create_cif_mapping_dict(a_logic_key, logic_keys, list_values)
                    mapping.get(self.Const.CIF_MAPPINGS).update(cif_mapping)

                    for logic_key in logic_keys:
                        if logic_key == self.Const.D:
                            new_items, dict_data = [], []
                            for dict_key in xml_value:
                                dict_data.append(xml_value.get(dict_key))
                                if dict_key == "PUBMED":
                                    new_items.append("pdbx_database_id_PubMed")
                                if dict_key == "DOI":
                                    new_items.append("pdbx_database_id_DOI")
                                if dict_key == "PATENT":
                                    new_items.append("pdbx_database_id_patent")
                                if dict_key == "ISSN":
                                    new_items.append("journal_id_ISSN")
                                if dict_key == "CSD":
                                    new_items.append("journal_id_CSD")
                                if dict_key == "ASTM":
                                    new_items.append("journal_id_ASTM")
                            cif_mapping = self.create_cif_mapping_dict(a_logic_key, new_items, dict_data)
                            mapping.get(self.Const.CIF_MAPPINGS).update(cif_mapping)

                values_set = True

        return values_set

    def create_cif_mapping_dict(self, cif_cat_key, cif_items, cif_values):
        """

        :param cif_cat_key:
        :param cif_items:
        :param cif_values:
        :return:
        """
        one_mapping = {
            cif_cat_key:
            {
                self.Const.ITEMS: cif_items,
                self.Const.DATA: cif_values
            }
        }

        return one_mapping

    def map_xml_value_to_code(self, xml_value, xml_mapping_code, optional_value=None):
        """
        Stores the value for the XML element or attribute for the mapping that corresponds to the mapping code
        :param optional_value:
        :param xml_value: a string; a value from XML
        :param xml_mapping_code: a string; its format is
            a. 'element'.'child'.'grandchild' for elements; e.g. emd.admin.title
            b. 'element'.'child'.'grandchild'@'grandchild attribute' for attributes; e.g. emd@emdb_id
        :return map_successful: a boolean; True when value is set
        """
        map_successful = False
        for mapping, sub in self.mappings.items():
            if mapping == xml_mapping_code:
                self.set_mapping_xml_value(xml_value, mapping)
                map_successful = True
            elif optional_value:
                if xml_mapping_code in mapping and '&' in mapping:
                    self.set_mapping_xml_value({xml_value: optional_value}, mapping)
                    map_successful = True

        return map_successful

    def set_mapping_xml_value(self, value, mapping):
        """
        This method sets the value for XML_VALUE in a specific mapping
        e.g. from 'emd.admin.title emd_admin.title.XML_VALUE' to
        'emd_admin': {
                      self.ITEMS: ['title'],
                      self.DATA: [XML_VALUE]
        }
        The mapping is depends on weather a list or a single value is expected
        :param value: the value to be assigned to XML_VALUE in
        :param mapping: a dict key; a string defined
        :return value_set: a boolean; True if the mapping exists and the xml value is set
        """
        value_set = False
        if self.mappings.get(mapping):
            xml_value_object = self.mappings.get(mapping).get(self.Const.XML_VALUE)
            if isinstance(xml_value_object, list):
                xml_value_object.append(value)
            elif isinstance(xml_value_object, dict):
                xml_value_object.update(value)
            else:
                self.mappings.get(mapping)[self.Const.XML_VALUE] = value
            value_set = True

        return value_set
