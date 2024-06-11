# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:59:54 2022

@author: Dustin Garner
"""


import os
import subprocess
import copy
import time
import datetime
from pathlib import Path
import numpy as np
import pandas as pd
import pydantic
import matplotlib.pyplot as plt
import seaborn as sns
from fafbseg import flywire


def time_elapsed(func):
    """A wrapper function to print the elapsed time of a function.

    Parameters
    ----------
    func : callable
        The function of which to measure the duration.
    """
    def wrap(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f"Elapsed Time: {end-start}")
        return result
    return wrap


def attempt_func(func, repeat_times=15):
    """
    A wrapper function to repeat a function a number of times if it fails.
    
    Parameters
    ----------
    func : callable
        The function to repeat if needed.
    repeat_times : int
        The number of times to repeat function if it doesn't work.
    """
    def wrap(*args, **kwargs):
        attempts = 0
        while(True):
            try:
                result = func(*args, **kwargs)
            except KeyboardInterrupt:
                raise
            except:
                print("Trying")
                attempts += 1
                if attempts>repeat_times:
                    print("It didn't work")
                    break
            else:
                return result
        return
    return wrap

def notify_finished(func):
    """A wrapper function that makes a noise when the function is done.

    Parameters
    ----------
    func : callable
        The function to carry out.
    """
    def wrap(*args, **kwargs):
        result = func(*args, **kwargs)
        print('\a')
        return result
    return wrap
    

def copy_ids(neur_ids, make_unique=False, output=True):
    """A function that copies neuron IDs to the clipboard.
    
    Parameters
    ----------
    neur_ids : int or str or list-like
        The neurons to copy.
    make_unique : bool, optional
        Whether you want the list of neurons to be unique. The default
        is False.
    output : bool, optional
        Whether to print the number of neurons copied. The default is True.
    """
    if isinstance(neur_ids, int) or isinstance(neur_ids, str):
        neur_ids = np.array([neur_ids], dtype="U50")
    neur_ids = np.array(neur_ids, dtype="U50")
    if make_unique:
        neur_ids = np.unique(neur_ids)
    id_string = ", ".join(neur_ids)
    subprocess.run("pbcopy", text=True, input=id_string)
    if output:
        print(f"Copied {len(neur_ids)} IDs to the clipboard.")


def row_to_str(id_row):
    """Takes in a row of IDs that are int and converts to str for excel.
    
    Parameters: 
    ----------
    id_row : pd.Series
        A pd.Series that contains neuron IDs.

    Returns
    -------
    new_row: np.array
        A np.array of neuron IDs that are str.
    """
    new_row = np.array(id_row, dtype="U50")
    return new_row


def syn_df_to_str(syn_df):
    """Takes in a synapses dataframe from flywire.get_synapses and makes
        the pre and post columns into strings.

    Parameters
    ----------
    syn_df : pd.DataFrame
        Synapse dataframe from flywire.get_synapses

    Returns
    -------
    syn_df : pd.DataFrame
        The same df with pre and post columns being strings.
    """
    syn_df["pre"] = row_to_str(syn_df["pre"])
    syn_df["post"] = row_to_str(syn_df["post"])
    return syn_df


def coord_str_to_list(coord):
    """Takes in a string of coordinates and returns them as a list.
    
    Parameters
    ----------
    coords: str
        "(x, y, z)"
    
    Returns
    -------
    coord_list : list
        [x, y, z]
    """
    coord = coord.replace("(", "")
    coord = coord.replace(")", "")
    coord_list = coord.split(",")
    coord_list = [int(x.strip()) for x in coord_list]
    return coord_list


def coord_str_to_arr(coord):
    """Takes in a string of coordinates and returns them as a numpy array,
        ready for flywire.locs_to_segments.
    
    Parameters
    ----------
    coords: str
        "(x, y, z)"
    
    Returns
    -------
    coord_arr : arr
        np.array([[x, y, z]])
    """
    coord = coord_str_to_list(coord)
    coord_arr = np.array([coord])
    return coord_arr


def coord_list_to_str(coord):
    """Takes in a coordinate as a list and returns it as a string.
    
    Parameters
    ----------
    coord : np.array
        The coordinate.

    Returns
    -------
    coord_str : str
        The coordinate.
    """
    assert coord.size==3, "Coordinate should be of size 3."
    return f"{coord[0]!s}, {coord[1]!s}, {coord[2]!s}"


def coord_column_to_array(column):
    """Takes in an excel column of coordinates and returns array ready for Flywire.
    
    Parameters
    ----------
    column : pd.Series
        A series of coords stored as strings.
    
    Returns
    -------
    coord_arr : np.array
        The coords stored as an array, ready for flywire.locs_to_segments().
    """
    coord_arr = np.array([[]], dtype=np.uint64)
    for i in column:
        current_coord = np.array([coord_str_to_list(i)], dtype=np.uint64)
        if coord_arr.size==0:
            coord_arr = current_coord
        else:
            coord_arr = np.append(coord_arr, current_coord, axis=0)
    return coord_arr


def basemodel_list_to_df(basemodels):
    """Takes in a list of pydantic BaseModels and returns a pd.DataFrame of them.
    
    Parameters
    ----------
    basemodels : list[pydantic.BaseModel]
        A list of any dataclass inheriting from pydantic.BaseModel.
    
    Returns
    -------
        A pd.DataFrame in which every attribute of the BaseModel is a column,
        and each individual member of basemodels is a row.
    """
    if len(basemodels)==0:
        return pd.DataFrame()
    data = {x: [] for x in basemodels[0].__fields__.keys()}
    for i in basemodels:
        for j in data:
            data[j].append(getattr(i, j))
    for i in data:
        data[i] = np.array(data[i])
    return pd.DataFrame(data)
    

def write_excel(df, file_name):
    """Writes a dataframe to excel in the Excel-Plots folder.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to convert to excel.
    file_name : str
        The name of the file.
    """
    absolute_path = os.path.dirname(__file__)
    relative_path = os.path.join("Excel-Plots", f"{file_name}.xlsx")
    file_path = os.path.join(absolute_path, relative_path)
    df.to_excel(file_path)
    

def usable_coords(x_coords, y_coords, z_coords, divide):
    """Takes in three lists of x, y, and z coordinate components and returns
        them as usable coordinates for importing into Flywire.
    
    Parameters
    ----------
    x_coords : list-like
        Array of x components.
    y_coords : list-like
        Array of y components.
    z_coords : list-like
        Array of z components.
    divide : bool
        if True, will divide components by (4, 4, 40).

    Returns
    -------
    coord_array: np.array
        Array of coordinates.
    """
    coord_array = np.array([], dtype="U50")
    for x, y, z in zip(x_coords, y_coords, z_coords):
        if divide:
            temp_x, temp_y, temp_z = x//4, y//4, z//40
        else:
            temp_x, temp_y, temp_z = int(x), int(y), int(z)
        
        coord = f"({temp_x!s}, {temp_y!s}, {temp_z!s})"
        coord_array = np.append(coord_array, coord)
    return coord_array


def make_importable_syn_coords(syn_df, pre_or_post="pre", 
                         file_name="Importable Coordinates",
                         save_file=True):
    """Take a synapse dataframe and turns it into an importable coordinate
    spreadsheet to Flywire.

    Parameters
    ----------
    syn_df : pd.DataFrame
        A synapse dataframe from flywire.get_synapses().
    pre_or_post : str, options
        "pre" or "post". Whether you want the importable coordinates to be the
        presynaptic or postsynaptic sites. The default is "pre".
    file_name : str, optional
        The name of the file. The default is "Importable Coordinates".
    save_file : bool
        Whether to save the file.
    """
    coords = usable_coords(syn_df[f"{pre_or_post}_x"],
                           syn_df[f"{pre_or_post}_y"],
                           syn_df[f"{pre_or_post}_z"], True)
    return importable_coords(coords, file_name = file_name)


def importable_coords(coords, file_name="Importable Coordinates",
                      save_file=True):
    """Makes a csv file that can be imported into Flywire as annotations.
        Stores it in the Importable-Coords folder.

    Parameters
    ----------
    coords : list-like
        List of coordinates.
    file_name : str
        Name of file to export.
    save_file : bool
        Whether to save the file.
    """
    length = len(coords)
    empty_arr = np.empty(length, dtype="U50")
    coord_dict = {
        "Coordinate 1": coords,
        "Coordinate 2": empty_arr,
        "Ellipsoid Dimensions": empty_arr,
        "Tags": empty_arr,
        "Description": empty_arr,
        "Segment IDs": empty_arr,
        "Parent ID": empty_arr,
        "Type": np.full(length, "Point", dtype="U50"),
        "ID": empty_arr}
    coord_df = pd.DataFrame(coord_dict)
    absolute_path = os.path.dirname(__file__)
    relative_path = os.path.join("Importable-Coords", f"{file_name}.csv")
    file_path = os.path.join(absolute_path, relative_path)
    coord_df.to_csv(file_path, index=False)
    return coord_df


def importable_coord_lines(coords1, coords2, file_name="Importable Lines"):
    """Makes a csv file that can be imported into Flywire as annotation lines.
        Stores it in the Importable-Coords folder.

    Parameters
    ----------
    coords1 : list-like
        Array of coordinates for the first line.
    coords2 : list-like
        Array of coordinates for the second line.
    file_name : str
        Name of file to export.
    """
    length = len(coords1)
    empty_arr = np.empty(length, dtype="U50")
    coord_dict = {
        "Coordinate 1": coords1,
        "Coordinate 2": coords2,
        "Ellipsoid Dimensions": empty_arr,
        "Tags": empty_arr,
        "Description": empty_arr,
        "Segment IDs": empty_arr,
        "Parent ID": empty_arr,
        "Type": np.full(length, "Line", dtype="U50"),
        "ID": empty_arr}
    coord_df = pd.DataFrame(coord_dict)
    absolute_path = os.path.dirname(__file__)
    relative_path = f"Importable-Coords/{file_name}.csv"
    file_path = os.path.join(absolute_path, relative_path)
    coord_df.to_csv(file_path, index=False)


def count_instances(neurons, min_neurs=5):
    """Returns a dict of the count of neurons in the input array, and returns a 
       unique array in the order of number of instances.
    
    Parameters
    ----------
    neurons : np.array
        An array of neuron IDs.
    min_neurs : int, Optional
        If the number of instances is lower than this, the IDs will not be 
        returned.
    
    Returns
    -------
    neur_dict : dict
        Keys (int): Neuron ID.
        Values (int): The number of instances of ID.
    """
    neurons, counts = np.unique(neurons, return_counts=True)
    neur_dict = dict(zip(neurons, counts))
    neurons = list(neur_dict.keys())
    for i in neurons:
        if neur_dict[i]<min_neurs:
            del neur_dict[i]
    neur_dict = dict(sorted(neur_dict.items(), key=lambda x: x[1], 
                            reverse=True))
    return neur_dict


def build_tree(values):
    """Builds a BST.

    Parameters
    ----------
    values : list-like
        List of values (likely neuron IDs) to include in the BST.

    Returns
    -------
    tree : SearchTree
        A binary search tree containing the values.
    """
    tree = SearchTree()
    for i in values:
        tree.insert(i)
    return tree

class TreeNode:
    """
    Node of a binary search tree.
    """
    def __init__(self, value):
        self.value = value
        self.parent = None
        self.left = None
        self.right = None
        self.size = 1
        self.height = 1
    
    def insert_left(self, node):
        self.left = node
        if node!= None:
            node.parent = self
        self.update_all()
    
    def insert_right(self, node):
        self.right = node
        if node!=None:
            node.parent = self
        self.update_all()
    
    def update_all(self):
        self.update_height()
        self.update_size()
    
    def update_size(self):
        left = 0 if self.left==None else self.left.size
        right = 0 if self.right==None else self.right.size
        self.size = left+right+1
    
    def update_height(self):
        left = 0 if self.left==None else self.left.height
        right = 0 if self.right==None else self.right.height
        self.height = max(left, right) + 1
            
    def get_skew(self):
        left = 0 if self.left==None else self.left.height
        right = 0 if self.right==None else self.right.height
        return right-left

    def rotate_right(self):
        if self.left==None:
            return
        grandparent = self.parent
        new_parent = self.left
        new_child = new_parent.right
        
        if grandparent==None:
            new_parent.parent = None
        elif grandparent.left==self:
            grandparent.insert_left(new_parent)
        else:
            grandparent.insert_right(new_parent)
        new_parent.insert_right(self)
        self.insert_left(new_child)
        
        for i in [self, new_parent, grandparent]:
            if i==None:
                continue
            i.update_all()
        
    def rotate_left(self):
        if self.right==None:
            return
        grandparent = self.parent
        new_parent = self.right
        new_child = new_parent.left
        
        if grandparent==None:
            new_parent.parent = None
        elif grandparent.right==self:
            grandparent.insert_right(new_parent)
        else:
            grandparent.insert_left(new_parent)
        new_parent.insert_left(self)
        self.insert_right(new_child)
        
        for i in [self, new_parent, grandparent]:
            if i==None:
                continue
            i.update_all()

    def make_balanced(self, tree):
        self.update_all()
        skew = self.get_skew()
        if skew>=-1 and self.get_skew()<=1:
            return
        if skew==2 and self.right.get_skew()==1:
            self.rotate_left()
        elif skew==2:
            self.right.rotate_right()
            self.rotate_left()
        elif skew==-2 and self.left.get_skew()==-1:
            self.rotate_right()
        else:
            self.left.rotate_left()
            self.rotate_right()
        if tree.root==self:
            tree.root = self.parent
    
    def __str__(self):
        return f"Value: {self.value}\n"\
               f"Parent: {'None' if self.parent==None else self.parent.value}\n"\
               f"Left: {'None' if self.left==None else self.left.value}\n"\
               f"Right: {'None' if self.right==None else self.right.value}\n"
        

class SearchTree:
    """
    A binary search tree.
    """
    def __init__(self):
        self.root = None
    
    def __str__(self):
        return f"[{', '.join([str(x) for x in self.get_values()])}]"
    
    def print_all(self):
        def recurse_print_all(node, depth=0):
            if node.left!=None:
                recurse_print_all(node.left, depth=depth+1)
            print(f"{node}Depth: {depth}\n")
            if node.right!=None:
                recurse_print_all(node.right, depth=depth+1)
        recurse_print_all(self.root)
    
    def get_values(self):
        def get_values(node):
            values = []
            if node==None:
                return values
            values += get_values(node.left)
            values.append(node.value)
            values += get_values(node.right)
            return values
        return get_values(self.root)
    
    def has_value(self, value):
        def does_node_exist(value, node):
            if node==None:
                return False
            if node.value==value:
                return True
            if node.value<value:
                return does_node_exist(value, node=node.right)
            return does_node_exist(value, node=node.left)
        return does_node_exist(value, self.root)
    
    def insert(self, value):
        new_node = TreeNode(value)
        if self.root==None:
            self.root = new_node
            return
        def do_insertion(node):
            if node.value==value:
                return
            elif node.value<value and node.right!=None:
                do_insertion(node.right)
            elif node.value<value:
                node.insert_right(new_node)
            elif node.value>value and node.left!=None:
                do_insertion(node.left)
            elif node.value>value:
                node.insert_left(new_node)
            node.make_balanced(self)
        do_insertion(self.root)
            








