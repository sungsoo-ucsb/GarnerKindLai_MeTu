<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 20:46:46 2022

@author: dusti
"""


import os
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from fafbseg import flywire
import flycode.utils as utils


def import_file(file_name, sheet_name=0, usecols=None, 
                file_type="xlsx"):
    """Imports a file from the Readable folder as a pandas DataFrame.

    Parameters
    ----------
    file_name : str
        The name of the file in Readable (excluding .xlsx)
    sheet_name : str or int
        Optional, the sheet name.
    usecols : str, list-like
        Optional, the columns to be parsed.
    file_type : str
        Optional, "xlsx" if you are reading and Excel file. "csv" if csv.

    Returns
    -------
    file : pd.DataFrame
        The file as a pandas DataFrame.
    """
    absolute_path = os.path.dirname(__file__)
    relative_path = f"Readable/{file_name}.{file_type}"
    file_path = os.path.join(absolute_path, relative_path)
    
    if file_type == "xlsx":
        file = pd.read_excel(file_path, sheet_name = sheet_name, 
                             usecols = usecols)
    else:
        file = pd.read_csv(file_path)
    return file


def import_input():
    """Imports the MeTu_input file and renames its columns to align with other
        versions that had different column names.

    Returns
    -------
    input_df : pd.DataFrame
        The MeTu_input file as a dataframe.
    """
    input_df = import_file("MeTu_input")
    input_df = input_df.rename(columns = {"type": "pre_type", \
                        "hemisphere": "pre_hemisphere", "weight": "syn_count"})
    return input_df


def import_proofreading(just_connectivity=False):
    """Imports the MeTu_proofreading file, depending on which sheet you want.
    
    Parameters
    ----------
    just_connectivity : bool, optional
        Whether you want sheet 0 (False) or sheet 1 (True). The default is False.
    
    Returns
    -------
    pd.DataFrame
        The MeTu_proofreading file as a dataframe.
    """
    sheet_name = 0 if not just_connectivity else 1
    return import_file("MeTu_proofreading", sheet_name=sheet_name)


def import_joined():
    """Imports joined_avp_table.

    Returns
    -------
    pd.DataFrame
        The joined_avp_table as a dataframe.

    """
    return import_file("joined_avp_table")
    

def import_regions():
    """Imports the regions with their include and exclude points.
    
    Returns
    -------
    regions_dict : Dictionary
        Makes a dictionary of the coordinates areas included and excluded
        from each region.
    """
    include_coords = import_file("Regions", sheet_name="Include")
    exclude_coords = import_file("Regions", sheet_name="Exclude")
    regions_dict = {}
    
    for i in include_coords:
        include_arr = np.array(include_coords[i].dropna())
        include_arr = np.array([utils.coord_str_to_list(x) for x in include_arr])
        
        regions_dict[i] = {"include": include_arr}
        if i in exclude_coords:
            exclude_arr = np.array(exclude_coords[i].dropna())
            exclude_arr = np.array([utils.coord_str_to_list(x) for x in exclude_arr])
            regions_dict[i]["exclude"] = exclude_arr
        else:
            regions_dict[i]["exclude"] = np.array([])
    return regions_dict
    
    
def import_coords():
    """Imports the coordinates from Neuron Spreadsheet.
    
    Returns
    -------
    neur_types : Dictionary
        Keys are neuron types from Neuron Spreadsheet and values are lists of
        coordinates of those types.
    """
    neur_file = import_file("Neuron Spreadsheet", sheet_name="Neurons")
    neur_types = {}
    prev_ids = {}
    
    for i in range(len(neur_file["Coordinates"])):
        if not isinstance(neur_file["Coordinates"][i], str):
            continue
        current_coord = utils.coord_str_to_list(neur_file["Coordinates"][i])
        current_id = neur_file["Current ID"][i]
        if current_id == current_id:
            current_id = np.int64(current_id)
        current_type = neur_file["Type"][i]
        
        if current_type in neur_types:
            neur_types[current_type] = np.concatenate((neur_types[current_type], 
                                np.array([current_coord], dtype=int)))
            
        else:
            neur_types[current_type] = np.array([current_coord], dtype=int)
        if current_type in prev_ids:
            prev_ids[current_type] = np.concatenate((prev_ids[current_type],
                            np.array([current_id], dtype=np.int64)))
        elif current_id == current_id:
            prev_ids[current_type] = np.array([current_id], dtype=np.int64)
    return neur_types, prev_ids


def import_neuprint_data():
    """Imports the Neuprint comparison data from the Comparison sheet of
        Neuron Spreadsheet.
    
    Returns
    -------
    compare_file: pd.DataFrame
        A dataframe from the "Comparison" sheet of the Neuron Spreadsheet, 
        limited to the Hemibrain dataset.
    """
    compare_file = import_file("Neuron Spreadsheet", sheet_name="Comparison")
    return compare_file[compare_file["Dataset"]=="Hemibrain"]
    

def import_colors():
    """Imports the colors from the Colors spreadsheet as pyplot.ListedColormap.
    
    Returns
    -------
    new_colors: dict
        Contains various extra colors as keys for colormaps.
    """
    new_colors = {}
    colors = ["Purple", "Orange", "Green", "TPurple", "TOrange1", "TOrange2", 
              "TGreen", "TBlue", "TRed", "TPink", "TCyan", "RegionTRed",
              "RegionTBlue", "RegionTGreen", "RegionTYellow", "RegionSolidRed",
              "RegionSolidBlue", "RegionSolidGreen", "RegionSolidYellow",
              "RegionSolidGray"]
    for i in colors:
        new_colors[i] = ListedColormap(np.array(import_file("Colors", \
                sheet_name = i, usecols = ["R", "G", "B", "A"])))
    return new_colors



=======
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 20:46:46 2022

@author: Dustin Garner
"""

import os
import math
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from fafbseg import flywire
import flycode.utils as utils


def import_file(file_name, sheet_name=0, usecols=None, 
                file_type="xlsx", dtype=None):
    """Imports a file from the Readable folder as a pandas DataFrame.

    Parameters
    ----------
    file_name : str
        The name of the file in Readable (excluding .xlsx)
    sheet_name : str or int
        Optional, the sheet name.
    usecols : str, list-like
        Optional, the columns to be parsed.
    file_type : str
        Optional, "xlsx" if you are reading and Excel file. "csv" if csv.
    dtype
        Optional, the datatype you would like the datasheet to be imported as.

    Returns
    -------
    file : pd.DataFrame
        The file as a pandas DataFrame.
    """
    absolute_path = os.path.dirname(__file__)
    relative_path = f"Readable/{file_name}.{file_type}"
    file_path = os.path.join(absolute_path, relative_path)
    
    if file_type == "xlsx":
        file = pd.read_excel(file_path, sheet_name = sheet_name, 
                             usecols = usecols, dtype=dtype)
    else:
        file = pd.read_csv(file_path, dtype=dtype)
    return file


def import_input():
    """Imports the MeTu_input file and renames its columns to align with other
        versions that had different column names.

    Returns
    -------
    input_df : pd.DataFrame
        The MeTu_input file as a dataframe.
    """
    input_df = import_file("MeTu_input")
    input_df = input_df.rename(columns = {"type": "pre_type", \
                        "hemisphere": "pre_hemisphere", "weight": "syn_count"})
    return input_df


def import_proofreading(just_connectivity=False):
    """Imports the MeTu_proofreading file, depending on which sheet you want.
    
    Parameters
    ----------
    just_connectivity : bool, optional
        Whether you want sheet 0 (False) or sheet 1 (True). The default is False.
    
    Returns
    -------
    pd.DataFrame
        The MeTu_proofreading file as a dataframe.
    """
    sheet_name = 0 if not just_connectivity else 1
    return import_file("MeTu_proofreading", sheet_name=sheet_name)


def import_joined():
    """Imports joined_avp_table.

    Returns
    -------
    pd.DataFrame
        The joined_avp_table as a dataframe.

    """
    return import_file("joined_avp_table")
    

def import_regions():
    """Imports the regions with their include and exclude points.
    
    Returns
    -------
    regions_dict : Dictionary
        Makes a dictionary of the coordinates areas included and excluded
        from each region.
    """
    include_coords = import_file("Regions", sheet_name="Include")
    exclude_coords = import_file("Regions", sheet_name="Exclude")
    regions_dict = {}
    
    for i in include_coords:
        include_arr = np.array(include_coords[i].dropna())
        include_arr = np.array([utils.coord_str_to_list(x) for x in include_arr])
        
        regions_dict[i] = {"include": include_arr}
        if i in exclude_coords:
            exclude_arr = np.array(exclude_coords[i].dropna())
            exclude_arr = np.array([utils.coord_str_to_list(x) for x in exclude_arr])
            regions_dict[i]["exclude"] = exclude_arr
        else:
            regions_dict[i]["exclude"] = np.array([])
    return regions_dict
    
    
def import_coords():
    """Imports the coordinates from Neuron Spreadsheet.
    
    Returns
    -------
    neur_types : Dictionary
        Keys are neuron types from Neuron Spreadsheet and values are lists of
        coordinates of those types.
    """
    neur_file = import_file("Neuron Spreadsheet", sheet_name="Neurons")
    neur_types = {}
    prev_ids = {}
    
    for i in range(len(neur_file["Coordinates"])):
        if not isinstance(neur_file["Coordinates"][i], str):
            continue
        current_coord = utils.coord_str_to_list(neur_file["Coordinates"][i])
        current_id = neur_file["Current ID"][i]
        if current_id == current_id:
            current_id = np.int64(current_id)
        current_type = neur_file["Type"][i]
        
        if current_type in neur_types:
            neur_types[current_type] = np.concatenate((neur_types[current_type], 
                                np.array([current_coord], dtype=int)))
            
        else:
            neur_types[current_type] = np.array([current_coord], dtype=int)
        if current_type in prev_ids:
            prev_ids[current_type] = np.concatenate((prev_ids[current_type],
                            np.array([current_id], dtype=np.int64)))
        elif current_id == current_id:
            prev_ids[current_type] = np.array([current_id], dtype=np.int64)
    return neur_types, prev_ids


def import_neuprint_data():
    """Imports the Neuprint comparison data from the Comparison sheet of
        Neuron Spreadsheet.
    
    Returns
    -------
    compare_file: pd.DataFrame
        A dataframe from the "Comparison" sheet of the Neuron Spreadsheet, 
        limited to the Hemibrain dataset.
    """
    compare_file = import_file("Neuron Spreadsheet", sheet_name="Comparison")
    return compare_file[compare_file["Dataset"]=="Hemibrain"]
    

def import_colors():
    """Imports the colors from the Colors spreadsheet as pyplot.ListedColormap.
    
    Returns
    -------
    new_colors: dict
        Contains various extra colors as keys for colormaps.
    """
    new_colors = {}
    colors = ["Purple", "Orange", "Green", "TPurple", "TOrange1", "TOrange2", 
              "TGreen", "TBlue", "TRed", "TPink", "TCyan", "RegionTRed",
              "RegionTBlue", "RegionTGreen", "RegionTYellow", "RegionSolidRed",
              "RegionSolidBlue", "RegionSolidGreen", "RegionSolidYellow",
              "RegionSolidGray"]
    for i in colors:
        new_colors[i] = ListedColormap(np.array(import_file("Colors", \
                sheet_name = i, usecols = ["R", "G", "B", "A"])))
    return new_colors


def import_synapses():
    """Imports the synapse spreadsheet from Codex.

    Returns
    -------
    A dataframe with all of the synapses, limited to 5 connections.
    """
    df = import_file("Codex Synapse Coordinates", file_type="csv", dtype=str)
    def convert_series(series):
        series = np.asarray(series)
        temp_value = series[0]
        for indi, i in enumerate(series):
            if math.isnan(float(i)):
                series[indi] = temp_value
            else:
                temp_value = i
        series = np.asarray(series, dtype=np.int64)
        return pd.Series(series)
    
    for i in [f"{x}_root_id" for x in ["pre", "post"]]:
        df[i] = convert_series(df[i])
    df = df.rename(columns = {f"{x}_root_id": x for x in ["pre", "post"]})    
    for i in "xyz":
        df[i] = df[i].astype("int")
    return df















>>>>>>> Stashed changes
