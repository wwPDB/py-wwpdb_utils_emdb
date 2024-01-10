
from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from mmcif.io.PdbxWriter import PdbxWriter


class CIF(object):
    """
    This class uses the mmcif library to create an mmCIF-like object
    Each object has one container, with a container ID and a list of DataCategory objects
    """
    DUMMY_CONTAINER_ID = "emd_0000"

    def __init__(self, cif_name_name):
        self.filename = cif_name_name
        # self.__dataList needed for PDBxWriter
        self.__dataList = []
        self.__container_id = None
        self.__container = None
        self.__dataMap = {}

    def write(self):
        """
        Given a file name, a pdbx writer is used to write data stored in self.__dataList
        :return written: a boolean; True when pdf writer is finished
        """
        written = False
        if self.filename:
            ofh = open(self.filename, "w")
            pdbx_writer = PdbxWriter(ofh)
            pdbx_writer.write(self.__dataList)
            ofh.close()
            written = True

        return written

    def add_container(self, container_id):
        """
        This method provides the basic functionality to set up a container
        :param container_id: a string; mmcif category e.g. 'emd_admin'
        :return:
        """
        added = False
        self.__container_id = container_id
        self.__container = DataContainer(container_id)
        self.__dataMap[container_id] = len(self.__dataList)
        self.__dataList.append(self.__container)
        if self.__container is not None:
            added = True

        return added

    def prepare_container(self, container_id):
        """
        Creates a container is it doesn't exist using either provided value or the dummy value
        :param container_id: a string; mmcif category e.g. 'emd_admin'
        :return:
        """
        if not self.__container:
            if container_id is None:
                container_id = self.DUMMY_CONTAINER_ID

            return self.add_container(container_id)

    def add_category(self, category_id, items):
        """
        This method creates a data category object, adds all items to it and appends it to the container
        :param category_id: a string; mmcif category e.g. 'emd_admin'
        :param items: a list of strings; each element in the list is an item of mmcif category as defined by category_id
        :return: a list of strings; each element represents a value for the corresponding element in data_items
        """
        category = DataCategory(category_id)
        for item in items:
            category.appendAttribute(item)
        self.__container.append(category)

    def insert_data(self, category_id, data_list):
        """
        This method appends the data in data_list to the container labeled category_id
        :param category_id: a string; mmcif category e.g. 'emd_admin'
        :param data_list:
        :return:
        """
        cat_obj = self.__container.getObj(category_id)
        if cat_obj is None:
            return

        if all(isinstance(i, list) for i in data_list):
            list_values = [list(t) for t in zip(*data_list)]
            # if all(item == "" for sublist in list_values for item in sublist):
            #     list_val = list_values[0]
            #     cat_obj.append(list_val)
            # else:
            cat_obj.extend(list_values)
        else:
            cat_obj.append(data_list)
        # print("END-MAPPINGS", cat_obj)

    def insert_data_into_category(self, category_id, data_items, data_list):
        """
        Helper method: calls two other methods, one to add a category and its items into a container and
        another to insert the data for the category items
        :param category_id: a string; mmcif category e.g. 'emd_admin'
        :param data_items: a list of strings; each element in the list is an item of mmcif category as
                           defined by category_id
        :param data_list: a list of strings; each element represents a value for the corresponding element in data_items
        :return:
        """
        data_ids, data_values = [], []
        if len(data_items) == len(set(data_items)):
            data_ids = data_items
            data_values = data_list

        elif len(data_items) != len(set(data_items)):
            for j, i in enumerate(data_items):
                if i not in data_ids:
                    dict_items = {i:data_items.count(i) for i in data_items}
                    nsub = max(dict_items.values())
                    data_ids.append(i)
                    data_i = data_items[j: j+2]
                    data_nsub = data_items[j: j+nsub]
                    data_v = data_list[j: j+nsub]
                    if data_i[0] != data_i[1]:
                        data_values.append(data_list[j])
                    if all(x==data_nsub[0] for x in data_nsub):
                        data_v = ['' if i is None else i for i in data_v]
                        if all(ele == '' for ele in data_v):
                            data_values.append('')
                        else:
                            data_value = [i for i in data_v if i]
                            for d in data_value:
                                data_values.append(d)

        deduped_data_values = self.deduping_empty_values(data_ids, data_values)
        self.add_category(category_id, data_ids)
        self.insert_data(category_id, deduped_data_values)

    def deduping_empty_values(self, data_ids, data_values):
        """
        This method keeps only one empty list is all the sublists are empty in the data_values
        :return:
        """
        deduped_data_values = []
        empty_data_values = all(all(not element for element in sub_list) for sub_list in data_values)
        if empty_data_values is True:
            if len(data_values) > 1:
                subs = iter(data_values)
                subs_len = len(next(subs))
                equal_subs = all(len(sub) == subs_len for sub in subs)
                if equal_subs:
                    if len(data_ids) == len(data_values[0]):
                        deduped_data_values = data_values[0]
                    else:
                        deduped_data_values = [None] * len(data_ids)
            else:
                deduped_data_values = [None] * len(data_ids)
        else:
            deduped_data_values = data_values
        return deduped_data_values

