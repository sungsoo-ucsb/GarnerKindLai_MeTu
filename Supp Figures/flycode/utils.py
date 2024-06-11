# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:59:54 2022

@author: dusti
"""


import os
import copy
import time
import datetime
from pathlib import Path
import numpy as np
import pandas as pd
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


def write_excel(df, file_name):
    """Writes a dataframe to excel in the Excel Plots folder.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to convert to excel.
    file_name : str
        The name of the file.
    """
    absolute_path = os.path.dirname(__file__)
    relative_path = os.path.join("Excel Plots", f"{file_name}.xlsx")
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
            temp_x, temp_y, temp_z = int(x/4), int(y/4), int(z/40)
        else:
            temp_x, temp_y, temp_z = int(x), int(y), int(z)
        
        coord = f"({temp_x!s}, {temp_y!s}, {temp_z!s})"
        coord_array = np.append(coord_array, coord)
    return coord_array


def make_importable_syn_coords(syn_df, pre_or_post="pre", 
                         file_name="Importable Coordinates"):
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
    """
    coords = usable_coords(syn_df[f"{pre_or_post}_x"],
                           syn_df[f"{pre_or_post}_y"],
                           syn_df[f"{pre_or_post}_z"], True)
    importable_coords(coords, file_name = file_name)


def importable_coords(coords, file_name="Importable Coordinates"):
    """Makes a csv file that can be imported into Flywire as annotations.
        Stores it in the Importable Coords folder.

    Parameters
    ----------
    coords : list-like
        List of coordinates.
    file_name : str
        Name of file to export.
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
    relative_path = os.path.join("Importable Coords", f"{file_name}.csv")
    file_path = os.path.join(absolute_path, relative_path)
    coord_df.to_csv(file_path, index=False)


def importable_coord_lines(coords1, coords2, file_name="Importable Lines"):
    """Makes a csv file that can be imported into Flywire as annotation lines.
        Stores it in the Importable Coords folder.

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
    relative_path = f"Importable Coords/{file_name}.csv"
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









