
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

    #
    # def update_single_row_value(self, category_id, item_name, row, value):
    #     """Update value in single row
    #     """
    #     catObj = self.__container.getObj(category_id)
    #     if catObj is None:
    #         return
    #     #
    #     catObj.setValue(value, item_name, row)
    #
    # def update_multiple_rows_value(self, category_id, item_name, value):
    #     """Update value in multiple rows
    #     """
    #     cat_obj = self.__container.getObj(category_id)
    #     if cat_obj is None:
    #         return
    #     #
    #     row_no = cat_obj.getRowCount()
    #     for row in range(0, row_no):
    #         cat_obj.setValue(value, item_name, row)

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
        if any(isinstance(el, list) for el in data_list):
            if all(isinstance(i, int) for i in data_list[0]):
                cat_id = [x for x in data_list[0]]
                new_list = [list(t) for t in zip(cat_id, data_list[1])]
                cat_obj.extend(new_list)
            else:
                updated_list = [v for i in data_list for v in (i if isinstance(i, list) else [i])]
                cat_obj.append(updated_list)
        else:
            cat_obj.append(data_list)
        print("ENDMAP", cat_obj)

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
        # print('INSERT DATA INTO CATEGORY:', category_id, data_items, data_list)
        self.add_category(category_id, data_items)
        self.insert_data(category_id, data_list)
