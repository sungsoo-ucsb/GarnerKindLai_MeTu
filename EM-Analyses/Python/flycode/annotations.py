#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 12:26:37 2023

@author: Dustin Garner
"""
import math
import datetime
import numpy as np
import pandas as pd
from fafbseg import flywire
import flycode.mapping as mapping
import flycode.utils as utils
import flycode.flywire_functions as fw


client = fw.client
version = client.materialize.version

user_ids = {
    2: ["Sven Dorkenwald", "Sebastian Seung Lab"],
    4: ["Claire McKellar", "Mala Murthy Lab, Sebastian Seung Lab"],
    8: ["Szi-chieh Yu", "Mala Murthy Lab, Sebastian Seung Lab"],
    13: ["Jay Gager", "Mala Murthy Lab, Sebastian Seung Lab"],
    14: ["James Hebditch", "Mala Murthy Lab, Sebastian Seung Lab"],
    17: ["Austin T Burke", "Mala Murthy Lab, Sebastian Seung Lab"],   
    27: ["Kyle Patrick Willie", "Mala Murthy Lab, Sebastian Seung Lab"],   
    28: ["Celia D", "Mala Murthy Lab, Sebastian Seung Lab"],   
    29: ["M Sorek", "Mala Murthy Lab, Sebastian Seung Lab"],
    60: ["Greg Jefferis", "Greg Jefferis Lab"],
    62: ["Philipp Schlegel", "Greg Jefferis Lab"],
    94: ["Gizem Sancer", "Mathias Wernet Lab"],
    95: ["Emil Kind", "Mathias Wernet Lab"],
    96: ["Doug Bland", "Mala Murthy Lab, Sebastian Seung Lab"],
    100: ["Dustin Garner", "Sung Soo Kim Lab"],
    103: ["Ben Silverman", "Mala Murthy Lab, Sebastian Seung Lab"],
    119: ["Zhihao Zheng", "Sebastian Seung Lab"],
    129: ["Gerit Linneweber", "Gerit Linneweber Lab"],
    153: ["Amalia Braun", "Alexander Borst Lab"],
    355: ["Alexander Bates", "Greg Jefferis Lab, Rachel Wilson Lab"],
    372: ["Yijie Yin", "Greg Jefferis Lab"],
    392: ["Krzysztof Kruk", "Eyewire"],
    700: ["a5hm0r", "Eyewire"],
    957: ["Varun Sane", "Greg Jefferis Lab"],
    1012: ["AzureJay", "Eyewire"],
    1063: ["Lab Members", "Volker Hartenstein Lab"],
    1230: ["Ben Gorko", "Sung Soo Kim Lab"],
    1321: ["J. Dolorosa", "Mala Murthy Lab, Sebastian Seung Lab"],
    1323: ["J. Anthony Ocho", "Mala Murthy Lab, Sebastian Seung Lab"],
    1324: ["Zairene Lenizo", "Mala Murthy Lab, Sebastian Seung Lab"],
    1329: ["Nash Hadjerol", "Mala Murthy Lab, Sebastian Seung Lab"],
    1491: ["Lucy Houghton", "Sung Soo Kim Lab"],
    1927: ["Ariel Dagohoy", "Mala Murthy Lab, Sebastian Seung Lab"],
    1928: ["remer tancontian", "Mala Murthy Lab, Sebastian Seung Lab"],
    1929: ["regine salem", "Mala Murthy Lab, Sebastian Seung Lab"],
    1931: ["Rey Adrian Candilada", "Mala Murthy Lab, Sebastian Seung Lab"],
    1932: ["Kendrick Joules Vinson", "Mala Murthy Lab, Sebastian Seung Lab"],
    1933: ["Joshua Bañez", "Mala Murthy Lab, Sebastian Seung Lab"],
    2102: ["Maria Ioannidou", "Marion Silies"],
    2395: ["Sebastian Mauricio Molina Obando", "Marion Silies Lab"],
    2455: ["annkri", "Eyewire"],
    2654: ["Burak Gur", "Marion Silies Lab"],
    2815: ["Nseraf", "Eyewire"],
    2843: ["TR77", "Eyewire"],
    2935: ["Lena Lörsch", "Marion Silies Lab"],
    3431: ["juan felipe vargas fique", "Marion Silies Lab"],
    3504: ["Jenna Joroff", "Wei-Chung Lee lab"]
}


def test_types(neur_types):
    """Tests whether the counts of neurons in Neuron Spreadsheet match up with
        the annotated neurons.

    Parameters
    ----------
    neur_types : list-like
        Neuron types in Neuron Spreadsheet.
    """
    
    df = client.materialize.query_table("neuron_information_v2", \
                        split_positions=True, materialization_version=version)
    df['created'] = df['created'].dt.tz_localize(None) #Allows to be exported to excel
    
    counter = 0
    fixed_types = []
    for i in neur_types:
        fixed_types.append(i[:-2])
    fixed_types = list(set(fixed_types))
    for i in fixed_types:
        mapping.add_types([f"{i}_L", f"{i}_R"], output=False)
        counter += len(df[df["tag"] == i])
        l = len(mapping.neur_ids[f"{i}_L"])
        r = len(mapping.neur_ids[f"{i}_R"])
        length = len(df[(df["tag"]==i) & (df["user_id"]==100)])
        print(f'Type: {i} has {length} / {l+r} neurons with the same name as my scheme in the annotations spreadsheet.')
    

def compare_annotations(neur_types):
    """Makes a spreadsheet to compare all annotations of the neuron types.

    Parameters
    ----------
    neur_types : list-like
        Types to check in spreadsheet.

    Returns
    -------
    label_df : pd.DataFrame
        The spreadsheet that is exported to excel.
    """
    neur_ids = {}
    for i in neur_types:
        neur_ids[i] = fw.locs_to_segments(mapping.neur_coords[i])
    
    anno_df = client.materialize.query_table("neuron_information_v2", \
            split_positions=True, 
            materialization_version=version)
    anno_df["created"] = anno_df["created"].dt.tz_localize(None)
    label_dict = {"type": np.array([], dtype=str),
                  "hemisphere": np.array([], dtype=str),
                  "x_coord": np.array([], dtype=float),
                  "y_coord": np.array([], dtype=float),
                  "z_coord": np.array([], dtype=float),
                 "id": np.array([], dtype=str),
                 "annotations": np.array([], dtype=str),
                 "timestamp": np.array([], dtype="datetime64[ns]"),
                 "user_id": np.array([], dtype=str),
                 "user_name": np.array([], dtype=str),
                 "user_affiliation": np.array([], dtype=str)}
    
    type_dict = {}
    for i in neur_types:
        for j in neur_ids[i]:
            type_dict[j] = i
        
    coord_axes = ["x", "y", "z"]
    for neur_id, neur_type in type_dict.items():
        temp_df = anno_df[anno_df["pt_root_id"]==neur_id]
        
        if len(temp_df) == 0:
            print(f"{neur_id} has no annotation")
            continue
        
        temp_df.loc[(temp_df["user_id"].isin([100,95,1491])) & 
                    (temp_df["created"]>np.datetime64('2022-12-14')), "created"]\
            = np.datetime64('2022-12-13T00:00:00.000000')
        temp_df = temp_df.sort_values("created")
        
        root_added = False
        for indk, k in temp_df.iterrows():
            anno_type = ""
            hemis = ""
            anno_id = ""
            coords = {w: math.nan for w in coord_axes}
            if not root_added:
                anno_type = neur_type[:-2]
                hemis = neur_type[-1]
                anno_id = str(neur_id)
                root_added = True
                for l in coord_axes:
                    divisor = 40 if l == "z" else 4
                    coords[l] = k[f"pt_position_{l}"] // divisor
                
            label_dict["type"] = np.append(label_dict["type"], anno_type)
            label_dict["hemisphere"] = np.append(label_dict["hemisphere"], hemis)
            for l in coord_axes:
                label_dict[f"{l}_coord"] = np.append(label_dict[f"{l}_coord"],
                                            coords[l])
            label_dict["id"] = np.append(label_dict["id"], anno_id)
            label_dict["annotations"] = np.append(label_dict["annotations"], 
                                                  k["tag"])
            label_dict["timestamp"] = np.append(label_dict["timestamp"], 
                                                 k["created"])
            label_dict["user_id"] = np.append(label_dict["user_id"], 
                                                k["user_id"])
            temp_id = k["user_id"]
            if not temp_id in user_ids:
                print(f"{temp_id} is not in user_ids dict")
            else:
                label_dict["user_name"] = np.append(label_dict["user_name"],
                                                    user_ids[temp_id][0])
                label_dict["user_affiliation"] = \
                        np.append(label_dict["user_affiliation"],
                                                    user_ids[temp_id][1])
    label_df = pd.DataFrame(label_dict)
    utils.write_excel(label_df, "Codex Annotations History")
    return label_df

        
        
    
    
    


