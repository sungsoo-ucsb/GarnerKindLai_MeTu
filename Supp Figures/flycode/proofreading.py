<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 08:10:52 2023

@author: dusti
"""

import os
import numpy as np
import pandas as pd
import caveclient
from fafbseg import flywire
import flycode.mapping as mapping
import flycode.flywire_functions as fw
import flycode.utils as utils


client = fw.client


#1575500097250 is timestamp before all MeTu edits

def empty_change_df():
    """Makes an empty df that can be combined with the edits of other neurons,
        specifically for neurons that have issues with accessing their data.
    
    Returns
    -------
    pd.DataFrame
        A dataframe in the style of client.chunkedgraph.get_tabular_change_log()
    """
    return pd.DataFrame(columns = ['operation_id', 'timestamp', 'user_id', 
                    'before_root_ids', 'after_root_ids', 'is_merge', 
                    'user_name', 'user_affiliation'])


def get_change_log(ids):
    """Retrieves the change log of the IDs of the neurons.
    
    Parameters
    ----------
    ids : list-like
        Neuron IDs from which to retrieve the change log.

    Returns
    -------
    proof_dict : dict
        Each key is a neuron ID, each value is a df of its edit history.
    """
    try:
        proof_dict = client.chunkedgraph.get_tabular_change_log(ids)
    except:
        proof_dict = individual_change_logs(ids)
    return proof_dict


def individual_change_logs(ids):
    """Similar to get_change_log(), but queries the IDs individual so if there
        is a problem neuron, it gives the neuron a blank df instead.
    
    Parameters
    ----------
    ids : list-like
        Neuron IDs from which to retrieve the change log.

    Returns
    -------
    proof_dict : dict
        Each key is a neuron ID, each value is a df of its edit history.
    """
    proof_dict = {}
    for i in ids:
        try:
            temp_dict = client.chunkedgraph.get_tabular_change_log(ids)
        except:
            temp_dict = {i: empty_change_df()}
        proof_dict[i] = temp_dict[i]
    return proof_dict


def condense_change_log(change_log):
    """Condenses a dictionary of change dfs into a single change df.
    
    Parameters
    ----------
    change_log : dict
        A log of changes retrieved from client.chunkedgraph.get_tabular_change_log()

    Returns
    -------
    change_df : pd.DataFrame
        The change_log condensed into a single dataframe.

    """
    change_df = pd.concat([change_log[x] for x in change_log])
    return change_df


def full_proof_table(ids):
    """Returns a condensed change log based on IDs.

    Parameters
    ----------
    ids : list-like
        A list of neuron IDs.

    Returns
    -------
    pd.DataFrame
        df of all edits of the neurons.
    """
    proof_dict = get_change_log(ids)
    return condense_change_log(proof_dict)


def proof_types(types):
    """Returned a change log table based on neuron types.

    Parameters
    ----------
    types : list-like
        Neuron types in mapping.neur_coords.

    Returns
    -------
    full_df : pd.DataFrame
        df of all edits of the neurons.
    """
    mapping.add_types(types)
    full_df = empty_change_df()
    for i in types:
        full_df = pd.concat([full_df, full_proof_table(mapping.neur_ids[i])])
    return full_df


def proof_counts(full_df):
    """Gives the edit counts of each user.
    
    Parameters
    ----------
    full_df : pd.DataFrame
        Dataframe given by proof_types().
        
    Returns
    -------
    proof_df : pd.DataFrame
        Dataframe with the edit counts of each person.
    """
    names = np.unique(np.array(full_df["user_name"]))
    proof_counts = {"user_id": np.array([], dtype=int),
                    "name": np.array([], dtype=str),
                    "lab": np.array([], dtype=str),
                    "edit_count": np.array([], dtype=int)}
    for i in names:
        temp_df = full_df[full_df["user_name"]==i]
        user_id = temp_df["user_id"][0]
        lab = temp_df["user_affiliation"][0]
        proof_counts["user_id"] = np.append(proof_counts["user_id"], user_id)
        proof_counts["name"] = np.append(proof_counts["name"], i)
        proof_counts["lab"] = np.append(proof_counts["lab"], lab)
        proof_counts["edit_count"] = np.append(proof_counts["edit_count"], 
                                               len(temp_df))
    proof_df = pd.DataFrame(proof_counts)
    proof_df = proof_df.sort_values("edit_count", ascending=False)
    return proof_df


def proof_counts_types(types):
    """Runs proof_counts() based on neuron types.

    Parameters
    ----------
    types : list-like
        A list of neuron types in Neuron Spreadsheet.

    Returns
    -------
    pd.DataFrame
        Dataframe with the edit counts of each person.
    """
    full_df = proof_types(types)
    return proof_counts(full_df)
        

def full_percent_spread(neur_types):
    """Shows how much each lab contributed to certain neurons.
    
    Parameters
    ----------
    neur_types : list-like
        A list of neuron types that are in mapping.neur_coords

    Returns
    -------
    edit_df : pd.DataFrame
        Breaks down how much each lab contributed to editing each neuron out
        of neur_types, assuming they contributed >=0.1 of the edits.
    """
    mapping.add_types(neur_types)
    all_ids = mapping.ids_from_types(neur_types)
    proof_dict = get_change_log(all_ids)
    edit_dict = {"id": np.array([], dtype=str),
                 "neuron_type": np.array([], dtype=str),
                 "total_edit_count": np.array([], dtype=int),
                 "lab": np.array([], dtype=str),
                 "lab_edit_count": np.array([], dtype=int),
                 "lab_percent_edited": np.array([], dtype=float)}
    for i in all_ids:
        temp_df = proof_dict[i]
        if len(temp_df)==0:
            continue
        current_type = mapping.find_neur_type(i)
        temp_df.replace("Mathias Wernet lab", "Mathias Wernet Lab")
        temp_df.replace("Greg Jefferis lab", "Greg Jefferis Lab")
        labs = np.unique(temp_df["user_affiliation"])
        total_edits = len(temp_df)
        for j in labs:
            temp_df2 = temp_df[temp_df["user_affiliation"]==j]
            if len(temp_df2) < (total_edits/10):
                continue
            edit_dict["id"] = np.append(edit_dict["id"], str(i))
            edit_dict["neuron_type"] = np.append(\
                    edit_dict["neuron_type"], current_type)
            edit_dict["total_edit_count"] = np.append(\
                    edit_dict["total_edit_count"], total_edits)
            edit_dict["lab"] = np.append(edit_dict["lab"], j)
            edit_dict["lab_edit_count"] = np.append(\
                    edit_dict["lab_edit_count"], len(temp_df2))
            edit_dict["lab_percent_edited"] = np.append(\
                    edit_dict["lab_percent_edited"], len(temp_df2)/total_edits)
    edit_df = pd.DataFrame(edit_dict)
    return edit_df
    

def neuron_counts(neur_types):
    """Shows the total edits each person performed on each neuron.
    
    Parameters
    ----------
    neur_types : list-like
        A list of neuron types that are in mapping.neur_coords

    Returns
    -------
    user_df : pd.DataFrame
        A dataframe showing the total number of edits that each person performed
        among those neurons (Only if they did >=0.1 of the edits of that neuron).
    """
    mapping.add_types(neur_types, output=False)
    lab_dict = {}
    edit_count_dict = {}
    for i in neur_types:
        ids = np.copy(mapping.neur_ids[i])
        change_log = get_change_log(ids)
        for j in ids:
            temp_df = change_log[j]
            if len(temp_df)==0:
                continue
            editors = np.unique(temp_df["user_name"])
            for k in editors:
                temp_df2 = temp_df[temp_df["user_name"]==k]
                if len(temp_df2)/len(temp_df) < 0.1:
                    continue
                if not k in lab_dict:
                    user_lab = np.array(temp_df2["user_affiliation"])[0]
                    user_lab = "Mathias Wernet Lab" if user_lab==\
                            "Mathias Wernet lab" else user_lab
                    user_lab = "Greg Jefferis Lab" if user_lab==\
                            "Greg Jefferis lab" else user_lab
                    lab_dict[k] = user_lab
                if not k in edit_count_dict:
                    edit_count_dict[k] = 0
                edit_count_dict[k] += len(temp_df2)
        print(i, " done")
    user_dict = {"Name": np.array([], dtype=str), 
                 "Lab": np.array([], dtype=str), 
                 "Total Edits": np.array([], dtype=int)}
    for i in lab_dict:
        user_dict["Name"] = np.append(user_dict["Name"], i)
        user_dict["Lab"] = np.append(user_dict["Lab"], lab_dict[i])
        user_dict["Total Edits"] = np.append(user_dict["Total Edits"], 
                                             edit_count_dict[i])
    user_df = pd.DataFrame(user_dict)
    return user_df


def all_codex_names(neur_types, keep_our_names = True):
    mapping.add_types(neur_types)
    absolute_path = os.path.dirname(__file__)
    relative_path = "Readable/Codex Labels.csv"
    file_path = os.path.join(absolute_path, relative_path)
    codex_csv = pd.read_csv(file_path)
    label_dict = {"type": np.array([], dtype = str),
                 "id": np.array([], dtype = str),
                 "neuron_name": np.array([], dtype = str),
                 "date_named": np.array([], dtype = str),
                 "user_name": np.array([], dtype = str),
                 "user_lab": np.array([], dtype = str)}
    our_names = ["Dustin Garner", "Emil Kind", "Ben Gorko", "Lucy Houghton"]
    for i in neur_types:
        for j in mapping.neur_ids[i]:
            temp_df = codex_csv[codex_csv["root_id"] == j]
            root_added = False
            if len(temp_df) == 0:
                label_dict["type"] = np.append(label_dict["type"], i)
                label_dict["id"] = np.append(label_dict["id"], j)
                label_dict["neuron_name"] = np.append(label_dict["neuron_name"], 
                                                      "No label in Codex")
                label_dict["date_named"] = np.append(label_dict["date_named"], "")
                label_dict["user_name"] = np.append(label_dict["user_name"], "")
                label_dict["user_lab"] = np.append(label_dict["user_lab"], "")
                
            
            for indk, k in temp_df.iterrows():
                if not keep_our_names and k["user_name"] in our_names:
                    continue
                type_name = ""
                id_name = ""
                if not root_added:
                    type_name = i
                    id_name = str(j)
                    root_added = True
                label_dict["type"] = np.append(label_dict["type"], type_name)
                label_dict["id"] = np.append(label_dict["id"], id_name)
                label_dict["neuron_name"] = np.append(label_dict["neuron_name"], 
                                                      k["label"])
                label_dict["date_named"] = np.append(label_dict["date_named"], 
                                                     k["date_created"])
                label_dict["user_name"] = np.append(label_dict["user_name"], 
                                                    k["user_name"])
                label_dict["user_lab"] = np.append(label_dict["user_lab"], 
                                                   k["user_affiliation"])
    label_df = pd.DataFrame(label_dict)
    utils.write_excel(label_df, "Codex Naming History")




=======
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 08:10:52 2023

@author: Dustin Garner
"""

import os
import numpy as np
import pandas as pd
import caveclient
from fafbseg import flywire
import flycode.mapping as mapping
import flycode.flywire_functions as fw
import flycode.utils as utils
import flycode.readfiles as readfiles


client = fw.client
#1575500097250 is timestamp before all MeTu edits


def empty_change_df():
    """Makes an empty df that can be combined with the edits of other neurons,
        specifically for neurons that have issues with accessing their data.
    
    Returns
    -------
    pd.DataFrame
        A dataframe in the style of client.chunkedgraph.get_tabular_change_log()
    """
    return pd.DataFrame(columns = ['operation_id', 'timestamp', 'user_id', 
                    'before_root_ids', 'after_root_ids', 'is_merge', 
                    'user_name', 'user_affiliation'])


def get_change_log(ids):
    """Retrieves the change log of the IDs of the neurons.
    
    Parameters
    ----------
    ids : list-like
        Neuron IDs from which to retrieve the change log.

    Returns
    -------
    proof_dict : dict
        Each key is a neuron ID, each value is a df of its edit history.
    """
    try:
        proof_dict = client.chunkedgraph.get_tabular_change_log(ids)
    except:
        proof_dict = individual_change_logs(ids)
    return proof_dict


def individual_change_logs(ids):
    """Similar to get_change_log(), but queries the IDs individual so if there
        is a problem neuron, it gives the neuron a blank df instead.
    
    Parameters
    ----------
    ids : list-like
        Neuron IDs from which to retrieve the change log.

    Returns
    -------
    proof_dict : dict
        Each key is a neuron ID, each value is a df of its edit history.
    """
    proof_dict = {}
    for i in ids:        
        try:
            temp_dict = client.chunkedgraph.get_tabular_change_log(ids)
        except:
            temp_dict = {i: empty_change_df()}
        proof_dict[i] = temp_dict[i]
    return proof_dict


def condense_change_log(change_log):
    """Condenses a dictionary of change dfs into a single change df.
    
    Parameters
    ----------
    change_log : dict
        A log of changes retrieved from client.chunkedgraph.get_tabular_change_log()

    Returns
    -------
    change_df : pd.DataFrame
        The change_log condensed into a single dataframe.

    """
    change_df = pd.concat([change_log[x] for x in change_log])
    return change_df


def full_proof_table(ids):
    """Returns a condensed change log based on IDs.

    Parameters
    ----------
    ids : list-like
        A list of neuron IDs.

    Returns
    -------
    pd.DataFrame
        df of all edits of the neurons.
    """
    proof_dict = get_change_log(ids)
    return condense_change_log(proof_dict)


def proof_types(types):
    """Returned a change log table based on neuron types.

    Parameters
    ----------
    types : list-like
        Neuron types in mapping.neur_coords.

    Returns
    -------
    full_df : pd.DataFrame
        df of all edits of the neurons.
    """
    mapping.add_types(types)
    full_df = empty_change_df()
    for i in types:
        full_df = pd.concat([full_df, full_proof_table(mapping.neur_ids[i])])
    return full_df


def proof_counts(full_df):
    """Gives the edit counts of each user.
    
    Parameters
    ----------
    full_df : pd.DataFrame
        Dataframe given by proof_types().
        
    Returns
    -------
    proof_df : pd.DataFrame
        Dataframe with the edit counts of each person.
    """
    names = np.unique(np.array(full_df["user_name"]))
    proof_counts = {"user_id": np.array([], dtype=int),
                    "name": np.array([], dtype=str),
                    "lab": np.array([], dtype=str),
                    "edit_count": np.array([], dtype=int)}
    for i in names:
        temp_df = full_df[full_df["user_name"]==i]
        user_id = temp_df["user_id"][0]
        lab = temp_df["user_affiliation"][0]
        proof_counts["user_id"] = np.append(proof_counts["user_id"], user_id)
        proof_counts["name"] = np.append(proof_counts["name"], i)
        proof_counts["lab"] = np.append(proof_counts["lab"], lab)
        proof_counts["edit_count"] = np.append(proof_counts["edit_count"], 
                                               len(temp_df))
    proof_df = pd.DataFrame(proof_counts)
    proof_df = proof_df.sort_values("edit_count", ascending=False)
    return proof_df


def proof_counts_types(types):
    """Runs proof_counts() based on neuron types.

    Parameters
    ----------
    types : list-like
        A list of neuron types in Neuron Spreadsheet.

    Returns
    -------
    pd.DataFrame
        Dataframe with the edit counts of each person.
    """
    full_df = proof_types(types)
    return proof_counts(full_df)
        

def full_percent_spread(neur_types):
    """Shows how much each lab contributed to certain neurons.
    
    Parameters
    ----------
    neur_types : list-like
        A list of neuron types that are in mapping.neur_coords

    Returns
    -------
    edit_df : pd.DataFrame
        Breaks down how much each lab contributed to editing each neuron out
        of neur_types, assuming they contributed >=0.1 of the edits.
    """
    all_ids = mapping.ids_from_types(neur_types)
    edit_dict = {"id": np.array([], dtype=str),
                 "neuron_type": np.array([], dtype=str),
                 "total_edit_count": np.array([], dtype=int),
                 "lab": np.array([], dtype=str),
                 "lab_edit_count": np.array([], dtype=int),
                 "lab_percent_edited": np.array([], dtype=float)}
    counter = 0
    for i in all_ids:
        if counter % 40 == 0:
            print(f"{counter} completed")
        counter += 1
        temp_df = get_change_log([i])[i]
        if len(temp_df)==0:
            continue
        current_type = mapping.find_neur_type(i)
        temp_df.replace("Mathias Wernet lab", "Mathias Wernet Lab")
        temp_df.replace("Greg Jefferis lab", "Greg Jefferis Lab")
        labs = np.unique(temp_df["user_affiliation"])
        total_edits = len(temp_df)
        for j in labs:
            temp_df2 = temp_df[temp_df["user_affiliation"]==j]
            if len(temp_df2) < (total_edits/10):
                continue
            edit_dict["id"] = np.append(edit_dict["id"], str(i))
            edit_dict["neuron_type"] = np.append(\
                    edit_dict["neuron_type"], current_type)
            edit_dict["total_edit_count"] = np.append(\
                    edit_dict["total_edit_count"], total_edits)
            edit_dict["lab"] = np.append(edit_dict["lab"], j)
            edit_dict["lab_edit_count"] = np.append(\
                    edit_dict["lab_edit_count"], len(temp_df2))
            edit_dict["lab_percent_edited"] = np.append(\
                    edit_dict["lab_percent_edited"], len(temp_df2)/total_edits)
    edit_df = pd.DataFrame(edit_dict)
    return edit_df
    

def neuron_counts(neur_types):
    """Shows the total edits each person performed on each neuron.
    
    Parameters
    ----------
    neur_types : list-like
        A list of neuron types that are in mapping.neur_coords

    Returns
    -------
    user_df : pd.DataFrame
        A dataframe showing the total number of edits that each person performed
        among those neurons (Only if they did >=0.1 of the edits of that neuron).
    """
    mapping.add_types(neur_types, output=False)
    lab_dict = {}
    edit_count_dict = {}
    for i in neur_types:
        if not i in mapping.neur_ids:
            continue
        ids = np.copy(mapping.neur_ids[i])
        change_log = get_change_log(ids)
        for j in ids:
            temp_df = change_log[j]
            if len(temp_df)==0:
                continue
            editors = np.unique(temp_df["user_name"])
            for k in editors:
                temp_df2 = temp_df[temp_df["user_name"]==k]
                if len(temp_df2)/len(temp_df) < 0.1:
                    continue
                if not k in lab_dict:
                    user_lab = np.array(temp_df2["user_affiliation"])[0]
                    user_lab = "Mathias Wernet Lab" if user_lab==\
                            "Mathias Wernet lab" else user_lab
                    user_lab = "Greg Jefferis Lab" if user_lab==\
                            "Greg Jefferis lab" else user_lab
                    lab_dict[k] = user_lab
                if not k in edit_count_dict:
                    edit_count_dict[k] = 0
                edit_count_dict[k] += len(temp_df2)
        print(i, " done")
    user_dict = {"Name": np.array([], dtype=str), 
                 "Lab": np.array([], dtype=str), 
                 "Total Edits": np.array([], dtype=int)}
    for i in lab_dict:
        user_dict["Name"] = np.append(user_dict["Name"], i)
        user_dict["Lab"] = np.append(user_dict["Lab"], lab_dict[i])
        user_dict["Total Edits"] = np.append(user_dict["Total Edits"], 
                                             edit_count_dict[i])
    user_df = pd.DataFrame(user_dict)
    return user_df


def all_codex_names(neur_types, keep_our_names = True):
    """Makes a spreadsheet of all relevant named neurons in Codex.

    Parameters
    ----------
    neur_types : list-like
        The types to include in the spreadsheet..
    keep_our_names : bool, optional
        Whether to keep the names of people in Kim and Wernet Labs. The default 
        is True.

    Returns
    -------
    label_df : pd.DataFrame
        The dataframe containing name information.
    """
    mapping.add_types(neur_types)
    codex_csv = readfiles.import_file("Codex Labels 783", file_type = "csv")
    columns = ["type",
               "hemisphere",
               "x_coord",
               "y_coord",
               "z_coord",
               "neur_id",
               "annotations",
               "timestamp",
               "user_id",
               "user_name",
               "user_affiliation"]
    label_dict = {x: [] for x in columns}
    
    our_names = ["Dustin Garner", "Emil Kind", "Ben Gorko", "Lucy Houghton"]
    original_time = "2022-12-13 00:00:00"
    user_ids = {100: ["Dustin Garner", "Sung Soo Kim Lab"],
                27: ["Kyle Patrick Willie", "Mala Murthy Lab, Sebastian Seung Lab"],
                392: ["Krzysztof Kruk", "Eyewire"],
                1929: ["regine salem", "Mala Murthy Lab, Sebastian Seung Lab"],
                28: ["Celia D", "Mala Murthy Lab, Sebastian Seung Lab"],
                96: ["Doug Bland", "Mala Murthy Lab, Sebastian Seung Lab"],
                13: ["Jay Gager", "Mala Murthy Lab, Sebastian Seung Lab"]}
    
    #exceptions are neurons that for some reason were not in the Labels spreadsheet,
    #despite having annotations in Flywire. Here I added them all in manually.
    exceptions = {720575940616396191: [["MeTu1", 100, original_time]],
                  720575940637556581: [["MeTu1", 100, original_time]],
                  
                  720575940631367877: [["Distal medulla 2a; Dm2a", 27, "2023-11-06 13:00:39"],
                                       ["Dm2", 100, original_time]],
                  720575940625284872: [["Distal medulla 2a; Dm2a", 27, "2023-11-06 12:50:34"],
                                       ["Dm2", 100, original_time]],
                  720575940638653254: [["Mi15", 100, original_time]],
                  720575940627943957: [["Mti_unknown_3", 100, original_time]],
                  720575940635179164: [["Mti_unknown_3", 100, original_time]],
                  720575940638867290: [["Retinula axon; R7", 1929, "2024-02-22 23:48:12"],
                                       ["R7", 100, original_time]],
                  720575940609212427: [["R7", 100, original_time],
                                       ["Retinula axon 8; R8", 28, "2023-02-01 15:07:05"]],
                  720575940629362929: [["Retinula axon; putative R7-DRA", 96, "2023-10-26 04:40:20"],
                                       ["R7", 100, original_time]],
                  720575940621194889: [["Retinula axon; R7", 96, "2023-11-10 22:13:29"],
                                       ["R7", 100, original_time]],
                  720575940638085232: [["Retinula axon; R7", 96, "2023-09-08 22:14:17"],
                                       ["R7", 100, original_time]],
                  720575940617613562: [["(correction) Unknown cell type", 96, "2023-12-24 15:35:47"],
                                       ["R7", 100, original_time],
                                       ["Putative Medulla intrinsic 15; Mi15", 
                                        13, "2023-08-24 13:17:45"]]}
    for i in [720575940621041217, 720575940631084800, 720575940640751960,
              720575940613967853]:
        exceptions[i] = [["Dm2", 100, original_time]]
    for i in [720575940618706136, 720575940624707767,720575940631170876,
              720575940618296374]:
        exceptions[i] = [["Mi15", 100, original_time]]
    for i, j in zip([720575940614414767, 720575940645359639, 720575940632269419,
                     720575940627537210, 720575940653528566, 720575940645582884,
                     720575940624363783],
                    [[41, 25], [41, 20], [42, 27],
                     [41, 15], [41, 15], [41, 11],
                     [41, 25]]):
        exceptions[i] = [["Mi15", 100, original_time],
                         ["Mi15; Medullary intrinsic neuron 15; FBbt_00111285",
                          392, f"2023-03-25 17:{j[0]}:{j[1]}"]]
    for i in [720575940629044042, 720575940616562754, 720575940624102282]:
        exceptions[i] = [["R7", 100, original_time]]
    
    for i in neur_types:
        if not i in mapping.neur_ids:
            continue
        for neur_id, coord in zip(mapping.neur_ids[i], mapping.neur_coords[i]):
            temp_type = i[:-2]
            temp_hemis = i[-1]
            x, y, z = np.array(coord, dtype=str)
            temp_id = str(neur_id)
            
            temp_df = codex_csv[codex_csv.root_id==neur_id]
            if len(temp_df) == 0:
                if not neur_id in exceptions:
                    print("No annotation: ", neur_id)
                    temp_df.loc[0] = [temp_id, "No annotation in Codex",
                                      "", "", "", "", "", "", ""]
                for indj, j in enumerate(exceptions.get(neur_id, [])):
                    user_name, user_affiliation = user_ids[j[1]]
                    temp_df.loc[indj] = [temp_id, j[0], j[1], "", "", "",
                                         j[2], user_name, user_affiliation]
            temp_df.loc[(temp_df.user_name.isin(our_names)) &
                        (temp_df.date_created>original_time), 
                        "date_created"] = original_time
            temp_df = temp_df.sort_values(by=["date_created"])
            
            for indk, k in temp_df.iterrows():
                if not keep_our_names and k.user_name in our_names:
                    continue
                label_dict["type"].append(temp_type)
                label_dict["hemisphere"].append(temp_hemis)
                label_dict["x_coord"].append(x)
                label_dict["y_coord"].append(y)
                label_dict["z_coord"].append(z)
                label_dict["neur_id"].append(temp_id)
                label_dict["annotations"].append(k.label)
                label_dict["timestamp"].append(np.datetime64(k.date_created))
                label_dict["user_id"].append(str(k.user_id))
                label_dict["user_name"].append(k.user_name)
                label_dict["user_affiliation"].append(k.user_affiliation)
                
                temp_type = ""
                temp_hemis = ""
                x, y, z = "", "", ""
                temp_id = ""
            
    label_df = pd.DataFrame(label_dict)
    return label_df




>>>>>>> Stashed changes
