
import os
import re
from emdb_xml2cif_translator.translator_classes.Mappings import Mappings
import xml.etree.ElementTree as ET
import emdb_xml2cif_translator.input_files


class EMDBMetadata(object):
    """
    This class stores data needed by both, the input EMDB XML and the mmCIF output file
    There are two dictionaries for collecting data:
    1. mappings_in: a dictionary containing xml-to-cif mappings as read from the mappings input file
    2. emd: a dictionary containing values parsed from the input XML file
    """

    def __init__(self):
        """
        Initialisation prompts loading mapping logic into the map_utils object
        """
        self.emdb_id_from_filename = None
        self.filename_in = None
        self.xml_tree_created = False
        self.mappings_in = Mappings()
        self.cif = None
        self.emd = None
        self.MAPPINGS_LOGIC_FILENAME = 'emdb-xml2cif-mappings.txt'

    def parse_into_xml_tree(self, emdb_header_file_in):
        """
        Method for parsing an XML file into a dictionary
        :param emdb_header_file_in: an XML file; an EMDB entry v3.x header file
        :return xml_tree_created: a boolean; True when the input XML file is parsed and
                                  stored into a dictionary (self.emd)
        """
        # store the input file name; the EMDB ID contained in it might be different to the one in the file
        self.filename_in = emdb_header_file_in
        if not self.xml_tree_created:
            tree = ET.parse(self.filename_in)
            self.emd = tree.getroot()
            if self.emd is not None:
                self.xml_tree_created = True

        return self.xml_tree_created

    def process(self):
        """
        This script processes XML values using the mapping logic and stores them ready for writing into cif
        :return processed: a boolean; True when processing is finished (a bit useless)
        """
        processed = False
        if self.xml_tree_created:
            ###use XML values to store them in the mappings logic
            mappings_file = os.path.join(os.path.dirname(emdb_xml2cif_translator.input_files.__file__),
                                         self.MAPPINGS_LOGIC_FILENAME)
            if self.find_xml_values_to_mappings(mappings_file):
                #by having the xml values, use the mappings logic to create the cif ready dictionary
                if self.mappings_in.set_cif_mappings_values():
                    # cif data is now ready;
                    # before starting with cif mappings prepare cif container using EMDB ID value from XML file
                    if self.cif.prepare_container(self.mappings_in.get_mapping_logic_value(
                            self.mappings_in.Const.MAP_EMD_EMDB_ID,
                            self.mappings_in.Const.XML_VALUE).replace('-', '_').lower()):
                        # insert data into the cif object container
                        self.add_data_into_cif_container()
                        processed = True

        return processed

    def find_xml_values_to_mappings(self, mappings_file):
        """
       This method reads through the XML (emd) object, finds an element or attribute with
       its counterpart in mappings (self.map_utils.mappings_logic)
       If the mapping requires XML_VALUE, the XML value from emd is written
       """
        root = self.emd
        if os.path.isfile(mappings_file):
            f = open(mappings_file, 'r').read().split("\n")
            for line in f:
                if len(line) != 0:
                    if line[0] != '#':
                        xml_part = line.split(" ")[0]

                        if ':' in xml_part:
                            xml_part = xml_part.split(":")[1]

                        if '@' in xml_part and any(keyword not in xml_part for keyword in ['S$', 'H$', '^']):
                            xml_elem = xml_part.rsplit('@', 1)
                            if not '.' in xml_elem[0]:
                                for attrib_key, attrib_val in root.attrib.items():
                                    attrib_tag = '@'.join((xml_elem[0], attrib_key))
                                    self.mappings_in.map_xml_value_to_code(attrib_val, attrib_tag, root.text)
                            else:
                                el = xml_elem[0].split(".", 1)[1].replace('.', '/')
                                attrib_key = xml_elem[1]
                                el = root.findall(el)
                                if el:
                                    for at in el:
                                        attrib_val = at.get(attrib_key)
                                        if attrib_key == "size_kbytes":
                                            self.mappings_in.map_xml_value_to_code(str(int(attrib_val)*10e2), xml_part, at.text)
                                        else:
                                            self.mappings_in.map_xml_value_to_code(attrib_val, xml_part, at.text)
                        elif '^' in xml_part:
                            self.multiple_same_element(root, xml_part)
                        elif 'S$' in xml_part:
                            self.substitution_groups(root, xml_part)
                        elif 'H$' in xml_part:
                            xml_elem, item = xml_part.replace("H$", '').rsplit(".", 1)
                            map_type = ''
                            rtags = ['map', 'interpretation/half_map_list/half_map', 'interpretation/additional_map_list/additional_map', 'interpretation/segmentation_list/segmentation/mask_details']
                            for re_tags in rtags:
                                for elem in root.findall(re_tags):
                                    sub_elem = elem.find(item)
                                    if sub_elem is not None:
                                        if re_tags == "map":
                                            map_type = "primary map"
                                        elif "half_map" in re_tags:
                                            map_type = "half map"
                                        elif "additional_map" in re_tags:
                                            map_type = "additional map"
                                        elif "mask" in re_tags:
                                            map_type = "mask"
                                        self.mappings_in.map_xml_value_to_code(map_type, xml_part)
                        else:
                            elem = xml_part.split(".", 1)[1].replace('.', '/')
                            tags, item = elem.rsplit('/', 1)
                            att = ''
                            if "U$" in tags:
                                tags = tags.replace("U$", '')
                            elif "$I" in item:
                                item = item.replace("$I", '')
                            elif "R$" in item:
                                att, item = item.rsplit('|', 1)
                                item = item.replace("R$", '')
                            if "$I" in elem or "R$" in elem:
                                self.primary_and_reference_ids(xml_part, root, tags, item, att)
                            else:
                                for elem in root.findall(tags):
                                    sub_elem = elem.find(item)
                                    sub_elements = '' if sub_elem is None else str(sub_elem.text)
                                    self.mappings_in.map_xml_value_to_code(sub_elements, xml_part)

        return True

    def multiple_same_element(self, root, xml_part):
        """
        This method writes the XML_VALUES if there are multiple child elements with same name in the parent element
        """
        elem_dict = {}
        tags = xml_part.split('^', 1)[0].split('.', 1)[1].replace('.', '/')
        item = xml_part.rsplit('^', 1)[1].split('&', 1)[0]
        for elem in root.findall(tags):
            sub_elem = elem.findall(item)
            for sub in sub_elem:
                if '&' in xml_part:
                    attrib_key = xml_part.rsplit('^', 1)[1].split('&', 1)[1]
                    attrib_val = sub.get(attrib_key)
                    elem_dict[attrib_val] = sub.text
                    self.mappings_in.map_xml_value_to_code(elem_dict, xml_part)
                else:
                    self.mappings_in.map_xml_value_to_code(sub.text, xml_part)

    def substitution_groups(self, root, xml_part):
        """
        This method writes the substitution group's XML_VALUES for all the special anchors used in the input file
        """
        tags, item, att, parent_elem, sub_elements, enantio = '', '', '', '', '', ''
        xml_slices, other_slices, list_id, macro_list_elem, attrib_index, macro_list_index, reference_ids = [], [], [], [], [], [], []
        substitution = re.split('S\$|\$S', xml_part)
        count = len(substitution)
        sub_ele = substitution[1].split("|")
        for sub in sub_ele:
            if count == 2:
                xml_slice = substitution[0]+sub+"!"
                xml_slices.append(xml_slice)
            if count == 3:
                xml_slice = substitution[0]+sub+"!"+substitution[2]
                xml_slices.append(xml_slice)
        for slice in xml_slices:
            parent_elem = slice.split("!", 1)[0].split(".", 1)[1].replace(".", "/")
            elem = slice.replace("!", '').split(".", 1)[1].replace(".", "/")
            xslice = slice.replace("!", '')
            find_type = root.findall(parent_elem)
            if "_supramolecule" in slice:
                slice = re.sub(r'(cell_supramolecule|complex_supramolecule|organelle_or_cellular_component_supramolecule|sample_supramolecule|tissue_supramolecule|virus_supramolecule)', 'all_supramolecules', xslice)
                slice = re.sub(r'sci_species_name', 'natural_source', slice)
            elif "macromolecule_list" in slice:
                slice = re.sub(r'(protein_or_peptide|em_label|ligand|other_macromolecule|dna|rna|saccharide)', 'all_macromolecules', xslice)
            elif "map" in slice or "mask_details" in slice:
                if "&" in slice:
                    rep_item = xslice.rsplit(".", 1)[1].split("&", 1)[1]
                elif "$I" in slice:
                    rep_item = xslice.rsplit(".", 1)[1].split("$I", 1)[1]+".$I"
                else:
                    rep_item = xslice.rsplit(".", 1)[1]
                slice = "emd.all_maps."+rep_item
            else:
                slice = xslice

            if '>' in slice and find_type:
                for mol_type in find_type:
                    molecule_type = (mol_type.tag).split("_", 1)[0]
                    if molecule_type == "complex":
                        ribo = mol_type.find("ribosome-details")
                        if ribo is not None:
                            molecule_type = "ribosome"
                    if "_supramolecule" in mol_type.tag:
                        if molecule_type == "organelle":
                            molecule_type = "ORGANELLE OR CELLULAR COMPONENT"
                        self.mappings_in.map_xml_value_to_code(str(molecule_type).upper(), slice)
                    else:
                        if molecule_type in ["protein", "em", "other", "dna", "rna", "saccharide"]:
                            self.mappings_in.map_xml_value_to_code("polymer", slice)
                        if molecule_type == "ligand":
                            self.mappings_in.map_xml_value_to_code("non-polymer", slice)

            if '<' in slice and find_type:
                tags, item = xslice.split(".", 1)[1].replace(".", "/").replace("<", '').rsplit("/", 1)
                for element in root.findall(tags):
                    mol_type = (element.tag).split("_", 1)[0]
                    if mol_type == "protein" or mol_type == "other":
                        selem = element.find(item)
                        if selem.text == "LEVO":
                            enantio = "polypeptide(L)"
                        elif selem.text == "DEXTRO":
                            enantio = "polypeptide(D)"
                    elif mol_type == "saccharide":
                        selem = element.find(item)
                        if selem.text == "DEXTRO":
                            enantio = "polysaccharide(D)"
                        elif selem.text == "LEVO":
                            enantio = "polysaccharide(L)"
                    elif mol_type == "rna":
                        se = element.find("classification")
                        if se is not None:
                            enantio = "polyribonucleotide"
                    elif mol_type == "dna":
                        se = element.find("classification")
                        if se is not None:
                            enantio = "polydeoxyribonucleotide"
                    self.mappings_in.map_xml_value_to_code(enantio, slice)

            if ';' in slice and find_type:
                svalue = ''
                tags = xslice.split(".", 1)[1].replace(".", "/").replace(";", '')
                items = ['space_group', 'point_group', 'helical_parameters']
                for element in root.findall(tags):
                    for item in items:
                        sub_elem = element.find(item)
                        if sub_elem is not None:
                            if "space_group" in sub_elem.tag:
                                svalue = "2D CRYSTAL"
                            elif "point_group" in sub_elem.tag:
                                svalue = "POINT"
                            elif "helical_paramenters" in sub_elem.tag:
                                svalue = "HELICAL"
                            self.mappings_in.map_xml_value_to_code(str(svalue), slice)

            if '@' in slice:
                tags, attrib_key = elem.split('@', 1)
                el = root.find(tags)
                if el is not None:
                    attrib_val = el.get(attrib_key)
                    self.mappings_in.map_xml_value_to_code(attrib_val, slice, el.text)

            if not any(char in xml_part for char in ['@', '$I', 'R$', '>', '<', 'E$', 'A$', '%', '+', 'T$', 'H$']):
                if '&' in xslice:
                    tags, item = elem.rsplit('&', 1)
                else:
                    tags, item = elem.rsplit('/', 1)
                find_elem = root.findall(tags)
                find_parent_elem = root.findall(parent_elem)
                if find_elem:
                    for element in root.findall(tags):
                        if '&' in xslice:
                            sub_elements = element.get(item)
                        else:
                            sub_elem = element.find(item)
                            sub_elements = '' if sub_elem is None or str(sub_elem.text) == "None" else str(sub_elem.text)
                        self.mappings_in.map_xml_value_to_code(sub_elements, slice)
                else:
                    if find_parent_elem:
                        for l in range(len(find_parent_elem)):
                            self.mappings_in.map_xml_value_to_code('', slice)

            if any(char in slice for char in ['$I', 'R$']) and not any(char in slice for char in ['E$', 'A$', 'T$', 'H$']):
                if '$I' in slice:
                    tags, item = elem.split("$I", 1)
                elif 'R$' in slice:
                    tags, items = elem.split("R$", 1)
                    att, item = items.rsplit('|', 1)
                self.primary_and_reference_ids(slice, root, tags, item, att)

            if '%' in slice and '+' in slice:
                tags, items = elem.rsplit('%', 1)
                item, att = items.rsplit('+', 1)
                parent_element = root.findall(parent_elem)
                if parent_element:
                    for par_attrib in parent_element:
                        target_id = par_attrib.get(att)
                        xpath_expression = './/complex_supramolecule[@supramolecule_id="{}"]/macromolecule_list/macromolecule'.format(target_id)
                        matching_elements = root.findall(xpath_expression)
                        index = len(matching_elements)
                        slice_ind = slice + "&S&" + str(index)
                        macro_list_index.append(slice_ind)

                for element in root.findall(tags):
                    selem = element.find(item)
                    macro_list_elem.append(selem.text)
                    index = len(root.findall(tags))

            if 'E$' in slice:
                e_tags, e_item = [], []
                either_one = re.split('E\$|\$E', xslice)
                other_slices = self.spliting_anchors(either_one)
                for sl in other_slices:
                    if '&' in sl:
                        tags, item = sl.split(".", 1)[1].replace(".", "/").rsplit("&", 1)
                    elif '$I' in sl:
                        tags, item = sl.split(".", 1)[1].replace(".", "/").rsplit('$I', 1)
                        e_item.append(item)
                    elif 'R$' in sl:
                        tags, items = sl.split(".", 1)[1].replace(".", "/").rsplit('R$', 1)
                        att, it = items.split('|', 1)
                        e_item.append(it)
                    else:
                        tags, item = sl.split(".", 1)[1].replace(".", "/").rsplit("/", 1)

                if "$I" in slice or 'R$' in slice:
                    self.primary_and_reference_ids(slice, root, tags, e_item, att)
                else:
                    for element in root.findall(tags):
                        selem = element.find(item)
                        if '&' in slice:
                            sub_elements = element.get(item)
                            self.mappings_in.map_xml_value_to_code(str(sub_elements), slice)
                        elif '?' in slice:
                            if item == "theoretical?":
                                self.mappings_in.map_xml_value_to_code('NO', slice)
                            elif item == "experimental?":
                                self.mappings_in.map_xml_value_to_code('YES', slice)
                            else:
                                self.mappings_in.map_xml_value_to_code('', slice)
                        else:
                            if selem is None:
                                self.mappings_in.map_xml_value_to_code('', slice)
                            else:
                                matches = ["experimental", "theoretical"]
                                if any(x in slice for x in matches):
                                    if "supramolecule_list" in tags:
                                        self.mappings_in.map_xml_value_to_code(str(selem.text), slice)
                                    else:
                                        svalue = round(float(selem.text) * 1e6, 3)
                                        self.mappings_in.map_xml_value_to_code(str(svalue), slice)


            if 'A$' in slice:
                one_sl = re.split('A\$|\$A', xslice)
                all_sl = self.spliting_anchors(one_sl)
                for sl in all_sl:
                    tags, items = sl.split(".", 1)[1].replace(".", "/").rsplit("/", 1)
                    item, att = items.rsplit('+', 1)
                    parent_element = root.findall(parent_elem)
                    xpath_expression, nelem, relem = '', '', ''
                    if parent_element:
                        for par_attrib in parent_element:
                            target_id = par_attrib.get(att)
                            if "supramolecule" in parent_elem:
                                npath_expression = parent_elem+'[@supramolecule_id="{}"]/natural_source/organism'.format(target_id)
                                nelem = root.find(npath_expression)
                                rpath_expression = parent_elem+'[@supramolecule_id="{}"]/recombinant_expression/recombinant_organism'.format(target_id)
                                relem = root.find(rpath_expression)
                            elif "macromolecule" in parent_elem:
                                npath_expression = parent_elem+'[@macromolecule_id="{}"]/natural_source/organism'.format(target_id)
                                nelem = root.find(npath_expression)
                                rpath_expression = parent_elem+'[@macromolecule_id="{}"]/recombinant_expression/recombinant_organism'.format(target_id)
                                relem = root.find(rpath_expression)

                            if nelem is None and relem is None:
                                self.mappings_in.map_xml_value_to_code('', slice)
                            else:
                                if relem is not None:
                                    svalue = relem.text
                                    if svalue == "synthetic construct":
                                        svalue = "syn"
                                    else:
                                        svalue = "man"
                                    self.mappings_in.map_xml_value_to_code(str(svalue), slice)
                                elif nelem is not None and relem is None:
                                    self.mappings_in.map_xml_value_to_code("nat", slice)

            if 'T$' in slice:
                one_sl = re.split('T\$|\$T', xslice)
                all_sl = self.spliting_anchors(one_sl)
                for sl in all_sl:
                    if "$I" in sl:
                        tags, items = sl.split(".", 1)[1].replace(".", "/").rsplit('$I', 1)
                    elif 'R$' in sl:
                        tags, item = sl.split(".", 1)[1].replace(".", "/").rsplit('R$', 1)
                        att, items = item.split('|', 1)
                    else:
                        tags, items = sl.split(".", 1)[1].replace(".", "/").rsplit("/", 1)
                    chk_parent = root.findall(tags.rsplit("/", 1)[0])
                    chk_element = root.findall(tags)
                    if chk_parent and chk_element:
                        for elem in root.findall(tags):
                            selem = elem.find(items)
                            if '*' in items:
                                if "intial" in items:
                                    self.mappings_in.map_xml_value_to_code("INITIAL", slice)
                                elif "final" in items:
                                    self.mappings_in.map_xml_value_to_code("FINAL", slice)
                            elif "$I" in slice or 'R$' in slice:
                                self.primary_and_reference_ids(slice, root, tags, items, att)
                            else:
                                if selem is not None:
                                   self.mappings_in.map_xml_value_to_code(str(selem.text), slice)
                    else:
                        if not "R$" in slice and not "*" in slice:
                            if chk_parent and not chk_element:
                                self.mappings_in.map_xml_value_to_code('', slice)

            # elif '+' in xml_part:
            #     xml_elem, item = xml_part.replace("+", '').rsplit(".", 1)
            #     map_type = ''
            #     rtags = ['map', 'interpretation/half_map_list/half_map', 'interpretation/additional_map_list/additional_map', 'interpretation/segmentation_list/segmentation/mask_details']
            #     for re_tags in rtags:
            #         for elem in root.findall(re_tags):
            #             sub_elem = elem.find(item)
            #             if sub_elem is not None:
            #                 if re_tags == "map":
            #                     map_type = "primary map"
            #                 elif "half_map" in re_tags:
            #                     map_type = "half map"
            #                 elif "additional_map" in re_tags:
            #                     map_type = "additional map"
            #                 elif "mask" in re_tags:
            #                     map_type = "mask"
            #                 self.mappings_in.map_xml_value_to_code(map_type, xml_part)

        if "%" in xml_part:
            index = 0
            for macro_list in macro_list_index:
                macro_slice, size = macro_list.split("&S&", 1)
                sublist = macro_list_elem[index:index + int(size)]
                macro_sublist = ', '.join(sublist)
                index += int(size)
                self.mappings_in.map_xml_value_to_code(macro_sublist, macro_slice)


    def primary_and_reference_ids(self, slice, root, tags, e_item, att=None):
        index = 0
        if not isinstance(e_item, list):
            list_item = [e_item]
        else:
            list_item = e_item
        chk_parent = root.findall(tags.rsplit("/", 1)[0])
        chk_element = root.findall(tags)
        if chk_parent and chk_element:
            for item in list_item:
                for element in root.findall(tags):
                    sub_elem = element.findall(item)
                    if '$I' in slice:
                        for ind in range(1, len(sub_elem)+1):
                            index += 1
                            self.mappings_in.map_xml_value_to_code(str(index), slice)
                    elif 'R$' in slice:
                        if sub_elem:
                            attrib = element.get(att)
                            for ind in range(1, len(sub_elem)+1):
                                self.mappings_in.map_xml_value_to_code(str(attrib), slice)
                        else:
                            self.mappings_in.map_xml_value_to_code("", slice)
        else:
            if "$I" in slice:
                if chk_parent and not chk_element:
                    self.mappings_in.map_xml_value_to_code('', slice)


    def spliting_anchors(self, either_one):
        other_slices = []
        combined = either_one[1].split("|")
        for each in combined:
            other_slice = either_one[0]+each+either_one[2]
            other_slices.append(other_slice)
        return other_slices

    def add_data_into_cif_container(self):
        """
        This method adds all data from all cif mappings to the cif container
        :return:
        """
        cif_mappings = self.collate_cif_categories()
        for cif_category_id, cif_category_data in cif_mappings.items():
            self.cif.insert_data_into_category(cif_category_id,
                                               cif_category_data.get(self.mappings_in.Const.ITEMS),
                                               cif_category_data.get(self.mappings_in.Const.DATA))

    def collate_cif_categories(self):
        """
        This method reorganises mappings dictionary to another dictionary where
        data is better organised for sending to cif writer

        :return cif_ready_mappings: a dictionary; collates all items and data within one cif category; key: category id
        """
        cif_ready_mappings = {}
        for mapping_code, mapping in self.mappings_in.mappings.items():
            cif_mappings = mapping.get(self.mappings_in.Const.CIF_MAPPINGS)
            for cif_mapping, cif_values in cif_mappings.items():
                cif_v = list(dict.fromkeys(cif_values['items']))
                cif_values['items'] = cif_v
                if cif_ready_mappings.get(cif_mapping) is None:
                    cif_ready_mappings.update({cif_mapping: cif_values})
                else:
                    # the category already exists; append lists of items and data to those in the category
                    items_append = cif_values.get(self.mappings_in.Const.ITEMS)
                    data_append = cif_values.get(self.mappings_in.Const.DATA)
                    cif_ready_mappings.get(cif_mapping).get(self.mappings_in.Const.ITEMS).extend(items_append)
                    cif_ready_mappings.get(cif_mapping).get(self.mappings_in.Const.DATA).extend(data_append)

        return cif_ready_mappings