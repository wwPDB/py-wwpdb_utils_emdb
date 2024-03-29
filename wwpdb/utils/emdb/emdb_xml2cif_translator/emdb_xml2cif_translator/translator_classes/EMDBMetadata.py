
from emdb_xml2cif_translator.translator_classes.Mappings import Mappings
from lxml import objectify
import xml.etree.ElementTree as ET


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
            with open(emdb_header_file_in) as f:
                xml = f.read().encode()
            self.emd = objectify.fromstring(xml)
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
            # use XML values to store them in the mappings logic
            if self.store_xml_values_to_mappings_recursively(self.emd):
                # by having the xml values, use the mappings logic to create the cif ready dictionaryœ
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

    def store_xml_values_to_mappings_recursively(self, el, parent_tag=None):
        """
        This method iterates through the XML (emd) object, matches an element or attribute with
        its counterpart in mappings (self.map_utils.mappings_logic)
        If the mapping requires XML_VALUE, the XML value from emd is written

        :param parent_tag: a string; generated from appending elements names for the use in recursive calls
        :param el: a key in self.emd; represents an element of an XML file
        :return finished_successfully: a boolean; True if the recursion is finished successfully
        """
        finished_successfully = False

        if el is not None:
            acc_tag = '.'.join(filter(None, (parent_tag, el.tag)))
            # map element's attributes values
            if el.attrib:
                for attrib_key, attrib_val in el.attrib.items():
                    attrib_tag = '@'.join((acc_tag, attrib_key))
                    self.mappings_in.map_xml_value_to_code(attrib_val, attrib_tag, el.text)
            # map element's values
            self.mappings_in.map_xml_value_to_code(el.text, acc_tag)
            # check is el has children (sub-elements)
            children = el.getchildren()
            if children:
                # el has children; call itself for each child to find their values, attributes and children
                for child in children:
                    self.store_xml_values_to_mappings_recursively(child, acc_tag)
            # all done - success?
            finished_successfully = True

        return finished_successfully

    def list_mapping(self, xml_mapping):
        sub_elements = []
        tree = ET.parse(self.filename_in)
        root = tree.getroot()
        xml_elem = xml_mapping.rsplit(".", 1)
        for elem in root.findall(xml_elem[0]):
            sub_elem = elem.find(xml_elem[1])
            if sub_elem is None:
                sub_elements.append('')
            else:
                sub_elements.append(sub_elem.text)
        return sub_elements

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
                cif_v = list(set(cif_values['items']))
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
