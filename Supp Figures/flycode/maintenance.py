<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 18:52:51 2023

@author: dusti
"""
import numpy as np
import pandas as pd
from fafbseg import flywire
import flycode.readfiles as readfiles
import flycode.mapping as mapping
import flycode.utils as utils
import flycode.reduction as reduction
import flycode.flywire_functions as fw


mt1_l = ["MeTu1_L"]
mt1_r = ["MeTu1_R"]
mt2_l = ["MeTu2a_L", "MeTu2b_L"]
mt2_r = ["MeTu2a_R", "MeTu2b_R"]
mt3_l = ["MeTu3a_L", "MeTu3b_L", "MeTu3c_L"]
mt3_r = ["MeTu3a_R", "MeTu3b_R", "MeTu3c_R"]
mt4_l = ["MeTu4a_L", "MeTu4b_L", "MeTu4c_L", "MeTu4d_L"]
mt4_r = ["MeTu4a_R", "MeTu4b_R", "MeTu4c_R", "MeTu4d_R"]
mt_l_general = ["MeTu1_L", "MeTu2_L", "MeTu3_L", "MeTu4_L"]
mt_r_general = ["MeTu1_R", "MeTu2_R", "MeTu3_R", "MeTu4_R"]
mt_l = mt1_l + mt2_l + mt3_l + mt4_l
mt_r = mt1_r + mt2_r + mt3_r + mt4_r


def find_undone_input(min_syns=4):
    """Uses MeTu input table and finds which neurons have not been classified based
        on the minimum number of synapses.

    Parameters
    ----------
    min_syn : int, optional
        Minimum number of synapses you want to find which neurons have not been
        classified.

    Returns
    -------
    np.array
        The neurons which have not been classified.

    """
    input_df = readfiles.import_input()
    unknown_input = input_df[(input_df["pre_type"]!=input_df["pre_type"]) &
                             (input_df["syn_count"]>=min_syns)]    
    return np.array(unknown_input["pre_id"])


def check_joined_concurrence(types):
    """Checks concurrence with joined_avp_table. Input list of types.

    Parameters
    ----------
    types : list-like
        Neuron types in Neuron Spreadsheet.
    """
    full_joined = readfiles.import_joined()
    no_issues = np.array([], dtype=str)
    for i in types:
        current_type = i[:-2]
        current_hemis = i[-1]
        if current_type == "MeTu3c":
            extra_types = ["MeTu3c_ventral", "MeTu3c_dorsal"]
            temp_df = full_joined[(full_joined["type"].isin(extra_types)) & \
                        (full_joined["post_switch_hemisphere"]==current_hemis)]
        else:
            temp_df = full_joined[(full_joined["type"]==current_type) & \
                        (full_joined["post_switch_hemisphere"]==current_hemis)]
        temp_coords = utils.coord_column_to_array(temp_df["XYZ"])
        all_good = True
        for j in temp_coords:
            if not j in mapping.neur_coords[i]:
                print(f"{utils.coord_list_to_str(j)}        {i} \
                      is in joined but not in my list.")
                all_good = False
        for j in mapping.neur_coords[i]:
            if not j in temp_coords:
                print(f"{utils.coord_list_to_str(j)}        {i} \
                      is in my list but not in joined.")
                all_good = False
        if all_good:
            no_issues = np.append(no_issues, current_type)
    #print(f"\nThese types had no issues: {no_issues})                


def check_proof_concurrence():
    """
    Checks if proofreading spreadsheet matches with neuron spreadsheet.
    """
    full_proof = readfiles.import_proofreading()
    types = [f"MeTu{x}_R" for x in 
             ["1","2a","2b","3a","3b","3c","4a","4b","4c","4d"]]
    for i in types:
        current_type = i[4:-2]
        temp_df = full_proof[full_proof["subtype"] == current_type]
        temp_coords = utils.coord_column_to_array(temp_df["XYZ"])
        for j in temp_coords:
            if j in mapping.neur_coords[i]:
                continue
            print(f"{list(j)} {current_type} is in proofreading but not in list.")
        for j in mapping.neur_coords[i]:
            if j in temp_coords:
                continue
            print(f"{list(j)} {current_type} is in list but not in proofreading.")


def check_neur_connections(types, dend_region, axon_region):
    """Checks if neuron types have inputs in dend_region and outputs in axon_region.
    
    Parameters
    ----------
    types : list-like
        Neuron types in Neuron Spreadsheet.
    dend_region : str
        The region in which these types have dendrites.
    axon_region: str
        The region in which these types have axons.

    Returns
    -------
    issues : np.array
        The neurons which do not have inputs or outputs in the respective
        regions.
    """
    
    mapping.add_types(types)
    print("Done adding")
    all_neurs = np.array([], dtype = np.int64)
    left_neurs, right_neurs = np.copy(all_neurs), np.copy(all_neurs)
    for i in types:
        all_neurs = np.append(all_neurs, mapping.neur_ids[i])
        if i[-1]=="L":
            left_neurs = np.append(left_neurs, mapping.neur_ids[i])
        else:
            right_neurs = np.append(right_neurs, mapping.neur_ids[i])
    counter = len(all_neurs)
    issues = np.array([], dtype = np.int64)
    region_hemis = lambda x,y : "EB" if x=="EB" else f"{x}_{y}"
    for i in all_neurs:
        syn_df = fw.fetch_synapses(i)
        hemis = "L" if i in left_neurs else "R"
        temp_dend_region = region_hemis(dend_region, hemis)
        temp_axon_region = region_hemis(axon_region, hemis)
        dend_df = reduction.lim_region(syn_df, temp_dend_region)
        dend_df = dend_df[dend_df["post"]==i]
        axon_df = reduction.lim_region(syn_df, temp_axon_region)
        axon_df = axon_df[axon_df["pre"]==i]
        if dend_df.size==0:
            issues = np.append(issues, np.array([], dtype=np.int64))
            print(f"{i} does not have dendrites in the {dend_region}.")
        if axon_df.size==0:
            issues = np.append(issues, np.array([], dtype=np.int64))
            print(f"{i} does not have axons in the {axon_region}.")
        counter -= 1
        if not counter%10:
            print(counter)
    return issues


def check_duplicates():
    """
    Checks whether any types in Neuron Spreadsheet have duplicate neurons.
    """
    neur_file = readfiles.import_file("Neuron Spreadsheet")
    types = np.unique(np.array(neur_file["Type"], dtype=str))
    mapping.add_types(types)
    for i in types:
        if not i in mapping.neur_ids:
            continue
        temp_ids = mapping.neur_ids[i]
        unique = np.unique(temp_ids)
        if len(temp_ids)==len(unique):
            continue
        for j in temp_ids:
            count = np.where(temp_ids==j)
            if len(count)==1:
                continue
            print(f"Duplicate: {i} {count}")
    





    
=======
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 18:52:51 2023

@author: Dustin Garner
"""

"""
-----------
Maintenance
-----------

maintenance.find_undone_input()
================
#Last Run 5/8/23, identified 9 neurons. 


maintenance.check_joined_concurrence(all_mt_l + all_mt_r + tb_l + tb_r + \
                                     er_all + aotu46 + tutu + upstream_neurons + \
                                     ["Mi1_L", "Mi1_R"])
================
#Last Run 8/23/23
No issues

maintenance.check_proof_concurrence()
================
#Last Run 6/10/24
No issues


maintenance.check_duplicates()
"""


import numpy as np
import pandas as pd
from fafbseg import flywire
import flycode.readfiles as readfiles
import flycode.mapping as mapping
import flycode.utils as utils
import flycode.reduction as reduction
import flycode.flywire_functions as fw


mt1_l = ["MeTu1_L"]
mt1_r = ["MeTu1_R"]
mt2_l = ["MeTu2a_L", "MeTu2b_L"]
mt2_r = ["MeTu2a_R", "MeTu2b_R"]
mt3_l = ["MeTu3a_L", "MeTu3b_L", "MeTu3c_L"]
mt3_r = ["MeTu3a_R", "MeTu3b_R", "MeTu3c_R"]
mt4_l = ["MeTu4a_L", "MeTu4b_L", "MeTu4c_L", "MeTu4d_L"]
mt4_r = ["MeTu4a_R", "MeTu4b_R", "MeTu4c_R", "MeTu4d_R"]
mt_l_general = ["MeTu1_L", "MeTu2_L", "MeTu3_L", "MeTu4_L"]
mt_r_general = ["MeTu1_R", "MeTu2_R", "MeTu3_R", "MeTu4_R"]
mt_l = mt1_l + mt2_l + mt3_l + mt4_l
mt_r = mt1_r + mt2_r + mt3_r + mt4_r


def find_undone_input(min_syns=4):
    """Uses MeTu input table and finds which neurons have not been classified based
        on the minimum number of synapses.

    Parameters
    ----------
    min_syn : int, optional
        Minimum number of synapses you want to find which neurons have not been
        classified.

    Returns
    -------
    np.array
        The neurons which have not been classified.

    """
    input_df = readfiles.import_input()
    unknown_input = input_df[(input_df["pre_type"]!=input_df["pre_type"]) &
                             (input_df["syn_count"]>=min_syns)]    
    return np.array(unknown_input["pre_id"])


def check_joined_concurrence(types):
    """Checks concurrence with joined_avp_table. Input list of types.

    Parameters
    ----------
    types : list-like
        Neuron types in Neuron Spreadsheet.
    """
    full_joined = readfiles.import_joined()
    no_issues = np.array([], dtype=str)
    for i in types:
        current_type = i[:-2]
        current_hemis = i[-1]
        if current_type == "MeTu3c":
            extra_types = ["MeTu3c_ventral", "MeTu3c_dorsal"]
            temp_df = full_joined[(full_joined["type"].isin(extra_types)) & \
                        (full_joined["post_switch_hemisphere"]==current_hemis)]
        else:
            temp_df = full_joined[(full_joined["type"]==current_type) & \
                        (full_joined["post_switch_hemisphere"]==current_hemis)]
        temp_coords = utils.coord_column_to_array(temp_df["XYZ"])
        all_good = True
        for j in temp_coords:
            if not j in mapping.neur_coords[i]:
                print(f"{utils.coord_list_to_str(j)}        {i} \
                      is in joined but not in my list.")
                all_good = False
        for j in mapping.neur_coords[i]:
            if not j in temp_coords:
                print(f"{utils.coord_list_to_str(j)}        {i} \
                      is in my list but not in joined.")
                all_good = False
        if all_good:
            no_issues = np.append(no_issues, current_type)
    #print(f"\nThese types had no issues: {no_issues})                


def check_proof_concurrence():
    """
    Checks if proofreading spreadsheet matches with neuron spreadsheet.
    """
    full_proof = readfiles.import_proofreading()
    types = [f"MeTu{x}_R" for x in 
             ["1","2a","2b","3a","3b","3c","4a","4b","4c","4d"]]
    for i in types:
        current_type = i[4:-2]
        temp_df = full_proof[full_proof["subtype"] == current_type]
        temp_coords = utils.coord_column_to_array(temp_df["XYZ"])
        for j in temp_coords:
            if j in mapping.neur_coords[i]:
                continue
            print(f"{list(j)} {current_type} is in proofreading but not in list.")
        for j in mapping.neur_coords[i]:
            if j in temp_coords:
                continue
            print(f"{list(j)} {current_type} is in list but not in proofreading.")


def check_neur_connections(types, dend_region, axon_region):
    """Checks if neuron types have inputs in dend_region and outputs in axon_region.
    
    Parameters
    ----------
    types : list-like
        Neuron types in Neuron Spreadsheet.
    dend_region : str
        The region in which these types have dendrites.
    axon_region: str
        The region in which these types have axons.

    Returns
    -------
    issues : np.array
        The neurons which do not have inputs or outputs in the respective
        regions.
    """
    
    mapping.add_types(types)
    print("Done adding")
    all_neurs = np.array([], dtype = np.int64)
    left_neurs, right_neurs = np.copy(all_neurs), np.copy(all_neurs)
    for i in types:
        all_neurs = np.append(all_neurs, mapping.neur_ids[i])
        if i[-1]=="L":
            left_neurs = np.append(left_neurs, mapping.neur_ids[i])
        else:
            right_neurs = np.append(right_neurs, mapping.neur_ids[i])
    counter = len(all_neurs)
    issues = np.array([], dtype = np.int64)
    region_hemis = lambda x,y : "EB" if x=="EB" else f"{x}_{y}"
    for i in all_neurs:
        syn_df = fw.fetch_synapses(i)
        hemis = "L" if i in left_neurs else "R"
        temp_dend_region = region_hemis(dend_region, hemis)
        temp_axon_region = region_hemis(axon_region, hemis)
        dend_df = reduction.lim_region(syn_df, temp_dend_region)
        dend_df = dend_df[dend_df["post"]==i]
        axon_df = reduction.lim_region(syn_df, temp_axon_region)
        axon_df = axon_df[axon_df["pre"]==i]
        if dend_df.size==0:
            issues = np.append(issues, np.array([], dtype=np.int64))
            print(f"{i} does not have dendrites in the {dend_region}.")
        if axon_df.size==0:
            issues = np.append(issues, np.array([], dtype=np.int64))
            print(f"{i} does not have axons in the {axon_region}.")
        counter -= 1
        if not counter%10:
            print(counter)
    return issues


def check_duplicates():
    """
    Checks whether any types in Neuron Spreadsheet have duplicate neurons.
    """
    neur_file = readfiles.import_file("Neuron Spreadsheet")
    types = np.unique(np.array(neur_file["Type"], dtype=str))
    mapping.add_types(types)
    for i in types:
        if not i in mapping.neur_ids:
            continue
        temp_ids = mapping.neur_ids[i]
        unique = np.unique(temp_ids)
        if len(temp_ids)==len(unique):
            continue
        for j in temp_ids:
            count = np.where(temp_ids==j)
            if len(count)==1:
                continue
            print(f"Duplicate: {i} {count}")
    





    
>>>>>>> Stashed changes
