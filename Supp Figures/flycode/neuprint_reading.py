# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 10:20:49 2022

@author: dusti
"""


import os
import collections
import itertools
from enum import Enum
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import bokeh
import holoviews as hv
import hvplot.pandas
import bokeh.palettes
from bokeh.plotting import figure, show, output_notebook
from neuprint import Client
from neuprint import fetch_neurons, fetch_adjacencies, fetch_synapses, \
    NeuronCriteria as NC
import flycode.readfiles as readfiles
import flycode.utils as utils
import flycode.figures as figures


output_notebook()
c = Client('neuprint.janelia.org',
           dataset='hemibrain:v1.2.1', \
           token="eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6ImR1c3Rpbm"\
                 "dhcm5lcjY0QGdtYWlsLmNvbSIsImxldmVsIjoibm9hdXRoIiwiaW1hZ2Utd"\
                 "XJsIjoiaHR0cHM6Ly9saDMuZ29vZ2xldXNlcmNvbnRlbnQuY29tL2EvQUFU"\
                 "WEFKdzdLWHBiblkyb0VyNTJuVHhael9QWkFDZndZa3EzbzZwekVmbm09czk"\
                 "2LWM_c3o9NTA_c3o9NTAiLCJleHAiOjE4MjM4NDIyNTh9.klIdAPp7VE4nU"\
                 "X0QoYQrcisGYOEwdGnsA3ZoONDLk-4")
c.fetch_version()

neuprint_neuron_dict = {}
comparison_data = readfiles.import_file("Neuron Spreadsheet", 
                                        sheet_name = "Comparison")

class Ratio(Enum):
    METU_TO_TUBU = 1
    TUBU_TO_METU = 2
    RING_TO_TUBU = 3
    TUBU_TO_RING = 4
    
ratio_pairs = {}
ratio_pairs[Ratio.RING_TO_TUBU] = [["TuBu01", "ER4m"], ["TuBu02", "ER3a_ad"], 
                                    ["TuBu02", "ER3m"], ["TuBu03", "ER3d_a"],
                                    ["TuBu03", "ER3d_c"], ["TuBu03", "ER3d_d"],
                                    ["TuBu04", "ER3d_b"], ["TuBu05", "ER3p_ab"], 
                                    ["TuBu06", "ER5"], ["TuBu07", "ER3w_ab"], 
                                    ["TuBu08", "ER4d"], ["TuBu09", "ER2_ad"], 
                                    ["TuBu09", "ER2_b"], ["TuBu10", "ER2_c"]]
ratio_pairs[Ratio.TUBU_TO_RING] = [[x[1], x[0]] for x in 
                                   ratio_pairs[Ratio.RING_TO_TUBU]]
ratio_pairs[Ratio.TUBU_TO_METU] = [["MeTu1", "TuBu08"], ["MeTu2a", "TuBu01"], 
                                   ["MeTu2a", "TuBu01"],["MeTu2b", "TuBu06"], 
                                   ["MeTu3ab", "TuBu07"], ["MeTu3c", "TuBu09"],
                                   ["MeTu3c", "TuBu10"], ["MeTu4a", "TuBu03"], 
                                   ["MeTu4a", "TuBu04"], ["MeTu4b", "TuBu02"], 
                                   ["MeTu4c", "TuBu05"], ["MeTu4d", "TuBu02"],
                                   ["MeTu4d", "TuBu05"]]
ratio_pairs[Ratio.METU_TO_TUBU] = [[x[1], x[0]] for x in 
                                   ratio_pairs[Ratio.TUBU_TO_METU]]


"""
add_types(["MC61", "MC64", "TuBu01", "TuBu02", "TuBu03", "TuBu04", "TuBu05", 
           "TuBu06", "TuBu07", "TuBu08", "TuBu09", "TuBu10"])
neuread.neur_to_type(["MC61", "MC64"], 
                     ["TuBu08", "TuBu01", "TuBu06", "TuBu07", "TuBu09", 
                      "TuBu10", "TuBu02", "TuBu03", "TuBu04", "TuBu05"], 
                     file_name="MeTu to TuBu Neuprint")
"""


def add_types(types):
    """Adds types to neuprint_neuron_dict
    
    Parameters
    ----------
    types : list-like
        List of Neuprint types to add to neuron_dict.
    """
    for i in types:
        if i in neuprint_neuron_dict:
            continue
        neuron_df, roi_counts_df = fetch_neurons(NC(type=f"{i}.*", regex=True))
        neuprint_neuron_dict[i] = np.array(neuron_df["bodyId"], dtype=np.uint64)


def neur_to_type(pre_types, post_types, file_name="Neuprint Sum Array"):
    """Finds the sum of all post connections of given post types to individual
        pre neurons. Then exports to excel.
    
    Parameters
    ----------
    pre_types : list-like
        A list of pre_types.
    post_types : list-like
        A list of post_types.
    file_name : str, optional
        The name of the output excel file. The default is "Neuprint Sum Array".
    
    Returns
    -------
    syn_df : pandas.DataFrame
        Information of the neurons to the types.
    """
    add_types(pre_types)
    add_types(post_types)
    print("Done adding")
    
    pre_neurs = np.array([], dtype = np.uint64)
    for i in pre_types:
        pre_neurs = np.concatenate((pre_neurs, neuprint_neuron_dict[i]))
    
    syn_dict = {}
    syn_dict["Type"] = np.array([], dtype=str)
    for i in pre_types:
        neur_count = len(neuprint_neuron_dict[i])
        temp_arr = np.array([i] * neur_count, dtype=str)
        syn_dict["Type"] = np.concatenate((syn_dict["Type"], temp_arr))
    syn_dict["ID"] = utils.row_to_str(pre_neurs)
    
    neuron_df, conn_df = fetch_adjacencies(sources=pre_types, targets=post_types)
    for i in post_types:
        syn_dict[i] = np.array([], dtype=str)
        temp_df = conn_df[conn_df["bodyId_post"].isin(neuprint_neuron_dict[i])]
        for indj, j in enumerate(pre_neurs):
            temp_df2 = temp_df[temp_df["bodyId_pre"]==j]
            syn_count = str(temp_df2["weight"].sum())
            syn_dict[i] = np.append(syn_dict[i], str(syn_count))
    syn_df = pd.DataFrame(syn_dict)
    utils.write_excel(syn_df, file_name=file_name)
    return syn_df


def hemi_syn_counts(types, region, outliers=[], y_axis="Synapse Count"):
    """Gives the counts of each neuron's synapses within a given region.

    Parameters
    ----------
    types : TYPE
        DESCRIPTION.
    region : TYPE
        DESCRIPTION.
    outliers : TYPE, optional
        DESCRIPTION. The default is [].
    y_axis : TYPE, optional
        DESCRIPTION. The default is "Synapse Count".

    Returns
    -------
    neur_dict : dict
        A df-ready dict that contains each neuron in types and its number
        of synapses within region.
    """
    neuprint_data = readfiles.import_neuprint_data()
    neur_data = neuprint_data[neuprint_data["Generalized Type"].isin(types)]
    neur_data["Identifier"] = np.array(neur_data["Identifier"], dtype=np.int64)
    
    neur_ids = np.array(neur_data["Identifier"], dtype = np.int64)
    neur_dict = {
        "ID": np.array([], dtype="U50"),
        "Neuron Type": np.array([], dtype="U50"),
        "Synapse Type": np.array([], dtype="U50"),
        y_axis: np.array([], dtype=int)}
    unused, syn_df = fetch_neurons(neur_ids)
    for i in neur_ids:
        temp_row = neur_data[neur_data["Identifier"]==i]
        hemisphere = np.array(temp_row["Hemisphere"])[0][0]
        neur_type = f"{np.asarray(temp_row['Generalized Type'])[0]}_{hemisphere}"
        temp_region = f"{region}({hemisphere})"
        
        temp_df = syn_df[syn_df["bodyId"]==i]
        temp_df = temp_df[temp_df["roi"]==temp_region]
        
        pre_count = sum(np.array(temp_df["downstream"]))
        post_count = sum(np.array(temp_df["upstream"]))
        
        if i in outliers:
            pre_count, post_count = 0, 0
        
        neur_dict["ID"] = np.concatenate((neur_dict["ID"], \
                np.array([i, i], dtype="U50")))
        neur_dict["Neuron Type"] = np.concatenate((neur_dict["Neuron Type"],
                np.array([neur_type, neur_type], dtype="U50")))
        neur_dict["Synapse Type"] = np.concatenate((neur_dict["Synapse Type"],
                np.array(["Presynaptic", "Postsynaptic"], dtype="U50")))
        neur_dict[y_axis] = np.concatenate((neur_dict[y_axis],
                np.array([pre_count, post_count], dtype=int)))
    return neur_dict


def lobula_counts(just_metu4, plot_name="", save_figure=True):
    """Makes a dataframe with the lobula counts of all MeTu4 neurons or all
        MeTu neurons.
    
    Parameters
    ----------
    just_metu4 : bool
        Whether the figure will only conatin MeTu4 or all MeTu in Neuprint.
    plot_name : str, optional
        The name of the outputted plot. The default is "".
    save_figure : bool, optional
        Whether the figure will be saved. The default is True.
        
    Returns
    -------
    df : pd.DataFrame
    """
    neu_data = readfiles.import_neuprint_data()
    mt_data = neu_data[neu_data["Broad Type"]=="MeTu4"] if just_metu4\
                                       else neu_data[neu_data["Class"]=="MeTu"]
    mt_dict = {}
    mt_dict["ID"] = np.array(mt_data["Identifier"], dtype=np.uint64)
    mt_dict["Hemibrain Type"] = np.array(mt_data["Specific Type"], dtype = str)
    mt_dict["Generalized Type"] = np.array(mt_data["Generalized Type"], dtype = str)
    
    pre_lo = np.array([], dtype = np.uint64)
    post_lo = np.array([], dtype = np.uint64)
    for i in mt_dict["ID"]:
        unused, id_data = fetch_neurons([i])
        rois = np.array(id_data["roi"])
        if "LO(R)" in rois: 
            lo_row = id_data[id_data["roi"]=="LO(R)"]
            pre_count = np.array(lo_row["downstream"], dtype=np.uint64)
            post_count = np.array(lo_row["upstream"], dtype=np.uint64)
        else:
            pre_count = np.array([0], dtype=np.uint64)
            post_count = np.array([0], dtype=np.uint64)
        pre_lo = np.concatenate((pre_lo, pre_count))
        post_lo = np.concatenate((post_lo, post_count))
    mt_dict["Pre LO"] = pre_lo
    mt_dict["Post LO"] = post_lo
    
    double_dict = {}
    double_dict["ID"] = np.concatenate((mt_dict["ID"], mt_dict["ID"]))
    double_dict["Hemibrain Type"] = np.concatenate((mt_dict["Hemibrain Type"],
                                                    mt_dict["Hemibrain Type"]))
    double_dict["Generalized Type"] = np.concatenate((mt_dict["Generalized Type"],
                                                      mt_dict["Generalized Type"]))
    double_dict["Synapse Type"] = np.concatenate((
            np.array(["Presynaptic" for x in mt_dict["Pre LO"]], dtype=str),
            np.array(["Postsynaptic" for x in mt_dict["Post LO"]], dtype=str)))
    double_dict["Lobula Synapse Count"] = np.concatenate((mt_dict["Pre LO"], 
                                                   mt_dict["Post LO"]))
    df = pd.DataFrame(double_dict)
    neur_types = ["MC61", "MC64"]
    figures.strip_plot(data=df, 
                       x_label="Hemibrain Type", 
                       y_label="Lobula Synapse Count",
                       order=neur_types, 
                       hue="Synapse Type", 
                       plot_name=plot_name,
                       palette=["#FF0000", "#00FFFF"],
                       save_figure=save_figure)
    return df



def ratio_plot(ratio_type, plot_name, save_figure=True, fig_size=(2.5,1.25),
               y_ticks=np.arange(0,4,1)):
    """
    Plots ratios between two sets of neurons.

    Parameters
    ----------
    ratio_type : Ratio
        One of the ratios in the enum Ratio.
    plot_name : str
        Name of the plot.
    save_figure : bool, optional
        Whether to save the figure. The default is True.
    fig_size : tuple, optional
        The size of the figure. The default is (2.5, 1.25).

    Returns
    -------
    ratio_df : pd.DataFrame
        A dataframe of the ratios.
    """
    neuron_pairs = ratio_pairs[ratio_type]
    datasets = [("FAFB", "Left"), ("FAFB", "Right"), ("Hemibrain", "Right")]
    
    ratio_dict = {}
    for i in datasets:
        temp_df = comparison_data[(comparison_data["Dataset"]==i[0])&
                                 (comparison_data["Hemisphere"]==i[1])]
        temp_dataset = f"{i[0]} {i[1]}"
        temp_ratios = []
        for j in neuron_pairs:
            denominator_count = len(temp_df[temp_df["Generalized Type"]==j[0]])
            numerator_count = len(temp_df[temp_df["Generalized Type"]==j[1]])
            ratio = numerator_count / denominator_count
            temp_ratios.append(ratio)
        ratio_dict[temp_dataset] = np.array([temp_ratios], dtype=np.float64)
    
    types = ["", ""]
    if ratio_type == Ratio.METU_TO_TUBU:
        types = ["MeTu", "TuBu"]
    elif ratio_type == Ratio.TUBU_TO_METU:
        types = ["TuBu", "MeTu"]
    elif ratio_type == Ratio.RING_TO_TUBU:
        types = ["Ring", "TuBu"]
    elif ratio_type == Ratio.TUBU_TO_RING:
        types = ["TuBu", "Ring"]
    
    figures.create_bar_graph(data=ratio_dict, 
                           x_label="Neuron Ratio",
                           y_label=f"{types[0]} Counts to {types[1]} Counts",
                           x_ticks=np.array([f"{x[1]} / {x[0]}" for x in 
                                             neuron_pairs], dtype = "U50"),
                           y_ticks=y_ticks.round(1),
                           colors=["#1f77b4","#2ca02c","#d62728"],
                           fig_size=fig_size,
                           plot_name=plot_name,
                           save_figure=save_figure)
    
    return ratio_dict





def plot_comparison(compare_types, plot_name="Neuron Counts", 
                    y_ticks=np.arange(0,3,1),
                    save_figure=True, fig_size=(2.0,1.25)):
    """Compares a type between hemispheres and datasets.
    
    Parameters
    ----------
    compare_types : list-like
        The types that will be compared in the plot.
    plot_name : str, optional
        The name of the outputted plot. The default is "Neuron Counts".
    save_figure : bool, optional
        Whether the figure will be saved. The default is True.
    fig_size : tuple, optional
        The size of the figure. The default is (2.0, 1.25).
        
    Returns
    -------
    neuron_counts : dict
        The count of each neuron type, with their datasets and hemispheres.
    """
    neuron_counts = collections.defaultdict(lambda: np.array([], dtype=int))
    datasets = ["FAFB", "Hemibrain"]
    hemispheres = ["Left", "Right"]
    for i,j,k in itertools.product(compare_types,datasets,hemispheres):
        neur_count = len(comparison_data[
                    (comparison_data["Generalized Type"]==i) &
                    (comparison_data["Dataset"]==j) &
                    (comparison_data["Hemisphere"] == k)])
        if j=="Hemibrain" and k=="Left" and neur_count==0:
            continue
        
        type_label = f"{j} {k}"
        neuron_counts[type_label] = np.append(neuron_counts[type_label], 
                                               neur_count)
    for i in neuron_counts:
        neuron_counts[i] = np.expand_dims(neuron_counts[i], axis=0)

    dataset_hemispheres = [f"{x} {y}" for x in datasets for y in hemispheres]
    all_colors = {x:y for x,y in zip(dataset_hemispheres, 
                                 ["#1f77b4","#2ca02c","#ffa90e","#d62728"])}
    colors = [all_colors[x] for x in neuron_counts]
    
    figures.create_bar_graph(data=neuron_counts, 
                           x_label="Neuron Type",
                           y_label="Neuron Count",
                           x_ticks=compare_types,
                           y_ticks=y_ticks,
                           colors=colors,
                           fig_size=fig_size,
                           plot_name=plot_name,
                           save_figure=save_figure)
    
    return neuron_counts
    

