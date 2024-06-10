#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 09:24:49 2023

@author: dustin
"""

from collections import defaultdict 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fafbseg import flywire
import flycode.utils as utils
import flycode.readfiles as readfiles
import flycode.flywire_functions as fw


def make_change_log():
    """
    Returns the MeTu_Change_Log df ready to be exported as excel
    """
    proof = readfiles.import_proofreading()
    proof = proof[proof["proofreading"]=="yes"]
    change_dict={"coords": np.array([], dtype=str),
                "current_id": np.array([], dtype=str),
                "id_at_timestamp": np.array([], dtype=str),
                "timing": np.array([], dtype=str),
                "round": np.array([], dtype=str),
                "timestamp": np.array([], dtype=np.int64),
                }
    for j in proof.aotu_XYZ:
        coord = utils.coord_str_to_arr(j)
        current_id = fw.locs_to_segments(coord)
        for i in range(6):
            timing = "start" if i%2==0 else "end"
            edit_round = str(i//2 + 1)
            temp_df = proof[proof["aotu_XYZ"]==j]
            timestamp = \
                    np.array(temp_df[f"Round_{edit_round}_{timing}_Unix"])[0]
            id_at_timestamp = \
                    flywire.locs_to_segments(coord, timestamp=timestamp)
            
            timing = "before" if timing=="start" else "after"
            change_dict["coords"] = \
                    np.append(change_dict["coords"], j)
            change_dict["current_id"] = \
                    np.append(change_dict["current_id"], current_id)
            change_dict["id_at_timestamp"] = \
                    np.append(change_dict["id_at_timestamp"], id_at_timestamp)
            change_dict["timing"] = \
                    np.append(change_dict["timing"], timing)
            change_dict["round"] = \
                    np.append(change_dict["round"], edit_round)
            change_dict["timestamp"] = \
                    np.append(change_dict["timestamp"], timestamp)
    change_df = pd.DataFrame(change_dict)
    return change_df


def l2_info(save_file=True):
    """Gets the l2_info of neurons in MeTu_Change_Log.xlsx.
    
    Parameters
    ----------
    save_file : bool
        Whether to save the dataframe as an excel file.
        
    Returns
    -------
    l2 : pd.DataFrame
        The retrieved l2 information.
    """
    changes = readfiles.import_file("MeTu_Change_Log")
    neur_ids = np.asarray(changes["id_at_timestamp"])
    l2 = flywire.get_l2_info(neur_ids)
    l2["root_id"] = utils.row_to_str(l2["root_id"])
    if save_file:
        utils.write_excel(l2, "MeTu Proofreading Volumes")
    return l2
    

def get_f1_score(true_positive, false_negative, false_positive):
    """Returns and F1 score based on the given values.

    Parameters
    ----------
    true_positive : int
        The number of items in both the pre-edited and post-edited neurons.
    false_negative : int
        The number of items in the post-edited but not pre-edited neuron.
    false_positive : int
        The number of items in the pre-edited but not post-edited neuron.

    Returns
    -------
    f1_score : float
        The F1 score based on the given values.
    """
    f1_score = (true_positive) / \
            (true_positive + (0.5*(false_positive+false_negative)))
    return f1_score
    
        
def get_f1_nodes(before_id, after_id):
    """Get the F1 score comparing the nodes of a neuron before and after 
    editing.
    
    Parameters
    ----------
    before_id : np.int64
        The neuron ID from before editing.
    after_id : np.int64
        The neuron ID from after editing.

    Returns
    -------
    float
        The f1 score of the skeletal nodes comparing the before and after ID.
    """
    try:
        before_skeleton = flywire.get_l2_skeleton(before_id)
    except:
        return 0
    before_nodes = before_skeleton.nodes
    after_skeleton = flywire.get_l2_skeleton(after_id)
    after_nodes = after_skeleton.nodes
    
    true_positive = 0
    false_negative = 0
    false_positive = 0
    
    for x,y,z in zip(np.asarray(before_nodes["x"]), 
            np.asarray(before_nodes["y"]), np.asarray(before_nodes["z"])):
        temp_df = after_nodes[(after_nodes["x"]==x) & 
            (after_nodes["y"]==y) & (after_nodes["z"]==z)]
        if temp_df.empty:
            false_positive+=1
        else:
            true_positive+=1
    for x,y,z in zip(np.asarray(after_nodes["x"]), 
            np.asarray(after_nodes["y"]), np.asarray(after_nodes["z"])):
        temp_df = before_nodes[(before_nodes["x"]==x) & 
            (before_nodes["y"]==y) & (before_nodes["z"]==z)]
        if temp_df.empty:
            false_negative+=1
    
    return get_f1_score(true_positive, false_negative, false_positive)
    

def get_f1_synapses(before_id, after_id):
    """Get the F1 Score comparing the synapses of a neuron before and after 
    editing.
    
    Parameters
    ----------
    before_id : np.int64
        The neuron ID from before editing.
    after_id : np.int64
        The neuron ID from after editing.

    Returns
    -------
    float
        The f1 score of the synapses comparing the before and after ID.
    """
    before_voxels = flywire.roots_to_supervoxels(before_id)[before_id]
    after_voxels = flywire.roots_to_supervoxels(after_id)[after_id]
    
    true_positive_voxels = np.array([], dtype=np.int64)
    false_positive_voxels = np.array([], dtype=np.int64)
    false_negative_voxels = np.array([], dtype=np.int64)

    for i in before_voxels:
        if i in after_voxels:
            true_positive_voxels = np.append(true_positive_voxels, np.int64(i))
        else:
            false_positive_voxels = \
                    np.append(false_positive_voxels, np.int64(i))
    for i in after_voxels:
        if i in true_positive_voxels:
            continue
        false_negative_voxels = np.append(false_negative_voxels, np.int64(i))
    
    all_voxels = np.concatenate((true_positive_voxels,
                                 false_positive_voxels,
                                 false_negative_voxels))
    
    all_ids = np.unique(flywire.supervoxels_to_roots(all_voxels))
    syn_df = fw.fetch_synapses(all_ids)
    
    true_positive = len(\
            syn_df[(syn_df["pre_supervoxel"].isin(true_positive_voxels))
                 | (syn_df["post_supervoxel"].isin(true_positive_voxels))])
    false_positive = len(\
            syn_df[(syn_df["pre_supervoxel"].isin(false_positive_voxels))
                 | (syn_df["post_supervoxel"].isin(false_positive_voxels))])
    false_negative = len(\
            syn_df[(syn_df["pre_supervoxel"].isin(false_negative_voxels))
                 | (syn_df["post_supervoxel"].isin(false_negative_voxels))])
    
    return get_f1_score(true_positive, false_negative, false_positive)


@utils.time_elapsed
def compare_f1_scores():
    """Get the F1 cores of the nodes and synapses of MeTu_Change_Log neurons.

    Returns
    -------
    node_f1_scores : dict
        Arrays of F1 scores of nodes between three rounds of proofreading.
    syn_f1_scores : dict
        Arrays of F1 scores of synapses between three rounds of proofreading.
    """
    changes = readfiles.import_file("MeTu_Change_Log")
    
    ids = np.unique(changes["current_id"])
    rounds = ("before_first", "after_first", "after_second")
    node_f1_scores = {x: np.array([], dtype=float) for x in rounds}
    syn_f1_scores = {x: np.array([], dtype=float) for x in rounds}
    for i in ids:
        temp_df = changes[changes["current_id"]==i]
        temp_df = temp_df.set_index(["timing", "round"])
        last_id = temp_df.loc[("after",3), "id_at_timestamp"]
        
        for j, k in zip(["before","after","after"], [1,1,2]):
            old_id = temp_df.loc[(j,k),"id_at_timestamp"]
            node_f1 = get_f1_nodes(old_id, last_id)
            syn_f1 = get_f1_synapses(old_id, last_id)
            
            temp_round = f"{j}_{'first' if k==1 else 'second'}"
            node_f1_scores[temp_round] = \
                    np.append(node_f1_scores[temp_round], node_f1)
            syn_f1_scores[temp_round] = \
                    np.append(syn_f1_scores[temp_round], syn_f1)
    
    return node_f1_scores, syn_f1_scores


def plot_f1(f1_scores):
    """Plots the node and synapse F1 score charts.

    Parameters
    ----------
    f1_scores : tuple
        Items given by compare_f1_scores().
    """
    for i,j in zip(f1_scores, ["Skeletal Node", "Synaptic"]):
        temp_scores = [i[x] for x in i]
        x_labels = np.array(['Round 1', 'Round 2', 'Round 3'])
        
        fig, ax = plt.subplots(figsize=(7, 4))
        bplot = ax.boxplot(temp_scores,
                             vert=True,
                             patch_artist=True)        
        colors = ["#D62728", "#1F77B4", "#2CA02C"]
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
            
        for median in bplot['medians']:
            median.set_color('black')
        
        ax.yaxis.grid(True)
        ax.set_xlabel('Proofreading Round')
        ax.set_ylabel(f'{j} F1 Score')
        plt.xticks(ticks=ax.get_xticks(), labels = x_labels, fontsize = 6)
        #utils.create_figure(fig, ax, 
        #                    plot_name=f'{j} F1 Score Per Proofreading Round',
        #                    rotation=0, 
        #                    fix_x_labels=False)
        plt.show()










