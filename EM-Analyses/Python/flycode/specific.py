# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 17:28:36 2023

@author: Dustin Garner
"""


import os
import copy
import math
import collections
import itertools
import numpy as np
import pandas as pd
from pydantic import BaseModel
import matplotlib.pyplot as plt
import seaborn as sns
from fafbseg import flywire
import flycode.utils as utils
import flycode.readfiles as readfiles
import flycode.mapping as mapping
import flycode.neuprint_reading as neuread
import flycode.reduction as reduction
import flycode.utils as utils
import flycode.figures as figures
import flycode.flywire_functions as fw


aotu_subregions = {
        "Posterior Lateral": ["TuBu08", "MeTu1"],
        "Posterior Central": ["TuBu01", "TuBu06", "MeTu2"],
        "Anterior": ["TuBu07", "TuBu09", "TuBu10", "MeTu3"],
        "Medial": ["TuBu02", "TuBu03", "TuBu04", "TuBu05", "MeTu4"],
        }
bulb_subregions = {
        "Superior": ["TuBu06", "TuBu07", "TuBu08", "TuBu09", "TuBu10",
                     "ER5", "ER3w_ab", "ER4d", "ER2_ad", "ER2_b", "ER2_c"],
        "Inferior": ["TuBu02", "TuBu03", "TuBu04", "TuBu05", "ER3a_ad",
                             "ER3m", "ER3d_a", "ER3d_b", "ER3d_c", "ER3d_d", "ER3p_ab"],
        "Anterior": ["TuBu01","ER4m"]
        }


def sm17_map_by_dv_axis(plot_name="", plot_folder="", save_figure=True): 
    """Plots MeTu3b and MeTu3c synaptic counts with respect to their offset 
        from the medulla centroid along the d-v axis.
    
    Parameters
    ----------
    plot_name : str, optional
        The name of the output plot. The default is "".
    plot_folder : str, optional
        The name of the folder to save the plots to. The default is "".
    save_figure : bool, optional
        Whether to save the figure. The default is True.
        
    Returns
    -------
    type_dict : dict
        A dictionary containing the neuron IDs and their positions.
    """
    mapping.add_types(["Sm17_R"], output=False)
    type_dict = {}
    for i in ["MeTu3b", "MeTu3c"]:
        type_dict[f"{i} y"] = np.array([], dtype=float)
        type_dict[f"{i} sm17"] = np.array([], dtype=int)
    
    dv_table = readfiles.import_file(file_name="MeTu_dorsal_ventral")
    dv_table = dv_table[(dv_table["subtype"]=="MeTu3b") |\
                        (dv_table["subtype"]=="MeTu3c")]
    
    coords = utils.coord_column_to_array(dv_table["XYZ"])
    metu_ids = fw.locs_to_segments(coords)
    syn_df = fw.fetch_synapses(metu_ids, pre=False)
    for i, j, k in zip(metu_ids, dv_table["subtype"], dv_table["dv"]):
        type_dict[f"{j} y"] = np.append(type_dict[f"{j} y"], k)
        temp_df = syn_df[(syn_df["post"]==i) & 
                    (syn_df["pre"].isin(mapping.neur_ids["Sm17_R"]))]
        type_dict[f"{j} sm17"] = np.append(type_dict[f"{j} sm17"], len(temp_df))
    
    fig, ax = plt.subplots(figsize = (2.0,1.25))
    for i, j in zip(["MeTu3b", "MeTu3c"], ["#D627B5", "#1F77B4"]):
        ax.scatter(type_dict[f"{i} y"],
                   type_dict[f"{i} sm17"],
                   marker="o",
                   c=j,
                   s=4.5)
    
    font_size, label_pad = 6, 2
    ax.set_xlabel("Relative Offset from Medulla Centroid", 
                  fontsize=font_size, 
                  labelpad=label_pad)
    ax.set_ylabel("Sm17 Presynapse Count", 
                  fontsize=font_size, 
                  labelpad=label_pad)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    
    ax.spines[['right', 'top']].set_visible(False)
    figures.add_legend(legend=["MeTu3b", "MeTu3c"], marker_scale=1.6)
    if save_figure:
        figures.save_fig(fig, plot_name=plot_name, folder_path=[plot_folder])
    return type_dict
        

     

def fafb_syn_counts(types, region, plot_name="Synapse Counts", plot_folder="",
                    y_axis="Synapse Count", save_figure=True):
    """Compares regional synaptic counts in the FAFB dataset.

    Parameters
    ----------
    types : list-like
        The MeTu neuron types to compare.
    region : str
        The region in which to compare synaptic counts.
    plot_name : str, optional
        The name of the plot. The default is "Synapse Counts".
    plot_folder : str, optional
        The name of the folder to save the plots to. The default is "".
    y_axis : str, optional
        The name of the y-axis. The default is "Synapse Count".
    save_figure : bool, optional
        Whether to save the figure. The default is True.

    Returns
    -------
    neur_dict : dict
        A dictionary that contains information about the neuron synaptic
        counts.
    """
    mapping.add_types(types, output=False)
    
    class NeurModel(BaseModel):
        root_id: str
        neur_type: str
        syn_type: str
        y_ax: int
    neur_models = []
    for i in types:
        current_neurons = mapping.neur_ids[i]
        hemisphere = i[-1]
        syn_df = fw.fetch_synapses(current_neurons, 
                                      region=f"{region}_{hemisphere}")
        for j, k in itertools.product(current_neurons, ["pre", "post"]):
            syn_count = len(syn_df[syn_df[k]==j])
            temp_model = NeurModel(root_id=str(j),
                                   neur_type=i,
                                   syn_type=f"{k}synaptic".capitalize(),
                                   y_ax=syn_count)
            neur_models.append(temp_model)
    
    full_df = utils.basemodel_list_to_df(neur_models)
    full_df.rename(columns={"root_id": "ID",
                            "neur_type": "Neuron Type",
                            "syn_type": "Synapse Type",
                            "y_ax": y_axis},
                   inplace=True)
    figures.strip_plot(data=full_df, 
                       x_label="Neuron Type", 
                       y_label=y_axis,
                       order=types, 
                       hue="Synapse Type", 
                       plot_name=plot_name,
                       folder_path=[plot_folder],
                       palette=["#FF0000","#00FFFF"], 
                       save_figure=save_figure)
    return full_df


def full_comparison(types, region, plot_name="Synapse Counts", plot_folder="",
                    neu_outliers=[], y_axis="Synapse Count", 
                    save_figure=True):
    """Compares regional synaptic counts between both datasets.

    Parameters
    ----------
    types : list-like
        The MeTu neuron types to compare.
    region : str
        The region in which to compare synaptic counts.
    plot_name : str, optional
        The name of the plot. The default is "Synapse Counts".
    plot_folder : str, optional
        The name of the folder to save the plots to. The default is "".
    neu_outliers : list-like, optional
        Neurons to exclude from the Hemibrain search. The default is [].
    y_axis : str, optional
        The name of the y-axis. The default is "Synapse Count".
    save_figure : bool, optional
        Whether to save the figure. The default is True.

    Returns
    -------
    full_df : pd.DataFrame
        Dataframe with all neurons and their synaptic counts.
    """
    fafb_types, hemi_types = [], []
    for i in types:
        fafb_types += [f"{i}_L", f"{i}_R"]
        general_type = "MeTu3ab" if (i=="MeTu3a" or i=="MeTu3b") else i
        hemi_types.append(general_type)
    fafb_dict = fafb_syn_counts(fafb_types, region, y_axis=y_axis, 
                                save_figure=False)
    hemi_dict = neuread.hemi_syn_counts(types, region, 
                                        outliers=neu_outliers, y_axis=y_axis)
    
    left_hemi = False
    for i, j in zip([fafb_dict, hemi_dict], ["FAFB", "Hemibrain"]):
        for indk, k in enumerate(i["Neuron Type"]):
            i["Neuron Type"][indk] = f"{j} {k}"
            if k[:-2] in ["MeTu3a", "MeTu3b"]:
                i["Neuron Type"][indk] = f"MeTu3ab{k[-2:]}"
            if not left_hemi and j=="Hemibrain" and k[-1]=="L":
                left_hemi = True
            
    full_dict = {}
    for i in fafb_dict:
        full_dict[i] = np.concatenate((fafb_dict[i], hemi_dict[i]))
    full_df = pd.DataFrame(full_dict)

    order = []
    hemispheres = ["FAFB _L", "FAFB _R", "Hemibrain _R"]
    if left_hemi:
        hemispheres.insert(2, "Hemibrain _L")
    for i in types:
        temp_type = "MeTu3ab" if (i=="MeTu3a" or i=="MeTu3b") else i
        for j in hemispheres:
            temp_type2 = f"{j[:-2]}{temp_type}{j[-2:]}"
            order.append(temp_type2)
    
    figures.strip_plot(data=full_df, 
                       x_label="Neuron Type", 
                       y_label=y_axis,
                       hue="Synapse Type", 
                       order=order, 
                       plot_name=plot_name,
                       folder_path=[plot_folder],
                       save_figure=save_figure, 
                       palette=["#FF0000", "#00FFFF"],
                       fig_size=(2.5,1.25))
    return full_df
    

def mt3_pre_connections(pre_types, plot_name="", plot_folder="",
                        save_figure=True):
    """Stripplot of MeTu3 subtypes and their presynaptic connections to chosen
        neurons

    Parameters
    ----------
    pre_types : list-like
        Presynaptic neurons to compare connections to MeTu3.
    plot_name : str, optional
        The name of the plot. The default is "".
    plot_folder : str, optional
        The name of the folder to save the plots to. The default is "".
    save_figure : bool, optional
        Whether to save the figure. The default is True.

    Returns
    -------
    neur_df : pd.DataFrame
        A dataframe with the presynapse counts to MeTu3 subtypes from pre_types.

    """
    metu3 = [f"MeTu3{x}_{y}" for x in "abc" for y in "LR"]
    all_pre_types = list(np.array([[f"{x}_L",f"{x}_R"] \
                                   for x in pre_types]).flatten())
    mapping.add_types(metu3+all_pre_types, output=False)
    
    class NeurModel(BaseModel):
        root_id: str
        neuron_type: str
        presyn_type: str
        syn_count: int
        syn_weight: float
    neur_models = []
    
    all_metu3 = mapping.ids_from_types(metu3)
    syn_df = fw.fetch_synapses(all_metu3, pre=False)
    for i in all_metu3:
        temp_type = mapping.find_neur_type(i, search=metu3)
        ipsi_hemis = temp_type[-1]
        temp_df = syn_df[syn_df.post==i]
        temp_df = reduction.lim_region(temp_df, f"ME_{ipsi_hemis}")
        total_me_syns = len(temp_df)
        contra_hemis = "R" if ipsi_hemis=="L" else "L"
        for j in pre_types:
            
            pre_type = f"{j}_{contra_hemis}" if j=="MeMeDRA" \
                        else f"{j}_{ipsi_hemis}"
            pre_count = len(temp_df[temp_df["pre"].isin(\
                                                mapping.neur_ids[pre_type])])
            
            temp_model = NeurModel(root_id=str(i),
                                   neuron_type=temp_type,
                                   presyn_type=pre_type[:-2],
                                   syn_count=pre_count,
                                   syn_weight=pre_count/total_me_syns)
            neur_models.append(temp_model)
    
    neur_df = utils.basemodel_list_to_df(neur_models)
    neur_df.rename(columns={"root_id": "ID",
                            "neuron_type": "Neuron Type",
                            "presyn_type": "Presynaptic Type",
                            "syn_count": "Synapse Count",
                            "syn_weight": "Synaptic Weight"},
                   inplace=True)
    
    for i in range(2):
        plot_type = "Count" if i==0 else "Weight"
        synapse_label = "Synapse" if i==0 else "Synaptic"
        temp_save_fig = False if i==0 else save_figure
        figures.strip_plot(data=neur_df,
                           x_label="Neuron Type",
                           y_label=f"{synapse_label} {plot_type}",
                           hue="Presynaptic Type",
                           order=metu3, 
                           palette=["#1F77B4","#D627B5","#2CA02C"],
                           plot_name=f"{plot_name} {plot_type}",
                           folder_path=[plot_folder],
                           save_figure=temp_save_fig)    
    return neur_df


    
def tutu_comparison(save_figure=True, plot_name="TuTu Synapse Counts",
                    plot_folder=""):
    """
    Makes a bar graph comparing the percentages of TuTu connections that are
        MeTu, TuBu, or other on the ipsilateral and contralateral sides.

    Parameters
    ----------
    save_figure : bool, optional
        Whether to save the figure. The default is True.
    plot_name : str, optional
        What the plot name would be upon saving it. The default is "TuTu Synapse
        Counts".
    plot_folder : str, optional
        The name of the folder to save the plots to. The default is "".
        
    Returns
    -------
    conn_dict : dict
        A dictionary containing the TuTu connection counts to MeTu, TuBu, and
        other.
    """
    
    tutu_types = [f"TuTuB_{x}_{y}" for x in "ab" for y in "LR"]
    metu_types = [f"MeTu{x}_{y}" for x in [1,2,3,4] for y in "LR"]
    tubu_types = [f"TuBu0{x!s}_{y}" for x in range(1,10) for y in "LR"] +\
                 [f"TuBu10_{x}" for x in "LR"]
    mapping.add_types(tutu_types + metu_types + tubu_types, output=False)
    
    tutu_ids = np.concatenate([mapping.neur_ids[x] for x in tutu_types])
    metu_ids = np.concatenate([mapping.neur_ids[x] for x in metu_types])
    tubu_ids = np.concatenate([mapping.neur_ids[x] for x in tubu_types])
    
    conn_dict = {}
    syn_df = fw.fetch_synapses(tutu_ids)
    for i in tutu_types:
        tutu_id = mapping.neur_ids[i][0]
        syn_counts = [[],[],[]]
        ipsi_contra = [i[-1], "R" if i[-1]=="L" else "L"]
        for j,k in itertools.product(ipsi_contra, ["pre", "post"]):
            temp_df = reduction.lim_region(syn_df, f"AOTU_{j}")
            temp_df = temp_df[temp_df[k]==tutu_id]
            
            other_site = "pre" if k=="post" else "post"
            syn_counts[0].append(\
                               len(temp_df[~(temp_df[other_site].isin(metu_ids)) & 
                                           ~(temp_df[other_site].isin(tubu_ids))]))
            syn_counts[1].append(len(temp_df[temp_df[other_site].isin(tubu_ids)]))
            syn_counts[2].append(len(temp_df[temp_df[other_site].isin(metu_ids)]))

        conn_dict[i] = np.array(syn_counts)
        
    figures.create_bar_graph(data=conn_dict,
                           x_label="Connection Type",
                           y_label="Synapse Count",
                           x_ticks=["Ipsilateral Pre", "Ipsilateral Post", 
                                    "Contralateral Pre", "Contralateral Post"],
                           y_ticks=np.arange(0, 4000, 1000),
                           colors=["#FF7F0E","#27B0D6","#D627B5"],
                           fig_size=(2.0, 1.25),
                           plot_name=plot_name,
                           folder_path=[plot_folder],
                           save_figure=save_figure,
                           color_axis=1)
    
    return conn_dict


def get_region_avg_counts(bihem_type, region, side, relevant_subregions={}):
    """Retrieves the average number of pre and post connections of  ipsilateral 
    or contralateral side from a given bihemispheric neuron type in a certain
    region.

    Parameters
    ----------
    bihem_type : str
        The type of bihemispheric neuron without a hemisphere. (For instance:
        'TuTuB_a').
    region : str
        The region to search without a hemisphere. (For instance: AOTU).
    side : str
        Either 'ipsi' or 'contra' for ipsilateral and contralateral connections
        respectively.
    relevant_subregions : dict, optional
        A dictionary that contains the names of subregions as keys with lists
        of neurons within those regions as values. The variables aotu_subregions
        and bulb_subregions can be passed for the AOTU and Bulb respecitvely.
        The default is {}.

    Returns
    -------
    syn_counts : dict
        The numnber of presynaptic and postsynaptic connections from the given
        neuron to the given subregions. Any that are outside the subregions
        (or if relevant_subregions was passed {}) are given the key 'other'.
    """
    assert side in ["ipsi", "contra"], "side must be 'ipsi' or 'contra'."
    ids_separated = {x: mapping.ids_from_types([f"{bihem_type}_{x}"]) for x in "LR"}
    bihem_ids = np.concatenate([ids_separated[x] for x in ids_separated])
    syn_dfs = {x: fw.fetch_synapses(ids_separated[x]) for x in "LR"}
    other_side = lambda x: "L" if x=="R" else "R"
    for i in syn_dfs:
        hemis = i if side=="ipsi" else other_side(i)
        assert hemis in "LR"
        syn_dfs[i] = reduction.lim_region(syn_dfs[i], f"{region}_{hemis}")
    full_df = pd.concat([syn_dfs[x] for x in syn_dfs])
    pre_df = full_df[full_df["pre"].isin(bihem_ids)]
    post_df = full_df[full_df["post"].isin(bihem_ids)]
    outside_pre, outside_post = len(pre_df), len(post_df)
    syn_counts = {}
    for i in relevant_subregions:
        partner_types = [f"{x}_{y}" for x in relevant_subregions[i] for y in "LR"]
        partner_ids = mapping.ids_from_types(partner_types)
        pre_count = len(pre_df[pre_df["post"].isin(partner_ids)])
        post_count = len(post_df[post_df["pre"].isin(partner_ids)])
        outside_pre -= pre_count
        outside_post -= post_count
        syn_counts[i] = {"pre": pre_count, "post": post_count}
    syn_counts["Other"] = {"pre": outside_pre, "post": outside_post}
    for i in syn_counts:
        for j in syn_counts[i]:
            syn_counts[i][j] /= len(bihem_ids)
    return syn_counts
        

#radius_factor = 1/2000
fig_size = tuple((34*4.8)/275 for _ in range(2))


def make_pie_chart(syn_counts, bihem_type, region, side, folder_path=[]):
    """Makes pie charts based on the return dict of get_region_avg_counts().
    
    Parameters
    ----------
    syn_counts : dict
        A dictionary of presynaptic and postsynaptic connections, as given by
        get_region_avg_counts().
    bihem_type : str
        The type of bihemispheric neuron without a hemisphere. (For instance:
        'TuTuB_a').
    region : str
        The region to search without a hemisphere. (For instance: AOTU).
    side : str
        Either 'ipsi' or 'contra' for ipsilateral and contralateral connections
        respectively.
    """
    colors = np.array(["#FF0000", "#00FFFF"])
    side = side.capitalize()
    radius_factor = 1/1000 if bihem_type=="AOTU046" else 1/2000
    folder_path += ["Pie Charts", region]
    for i in syn_counts:
        fig, ax = plt.subplots(figsize=fig_size)
        temp_dict = syn_counts[i]
        pie_arr = np.array([temp_dict["pre"], temp_dict["post"]])
        total = temp_dict["pre"] + temp_dict["post"]
        if total==0:
            return
        radius = total*radius_factor
        plt.pie(pie_arr, colors=colors, radius=radius)
        file_name = f"{bihem_type} {side} {i} {total}"
        figures.save_fig(fig, plot_name=file_name, folder_path=folder_path,
                         file_type=figures.FileType.SVG,
                         transparent=True)
    

def make_bihem_pie_charts(plot_folder=""):
    """Makes all relevant bihemispheric pie charts.
    
    Parameters
    ----------
    plot_folder : str
        The name of the folder to save the plots to.
    """
    bihem_types = ["AOTU046", "TuTuB_a", "TuTuB_b"]
    sides = ["ipsi", "contra"]
    regions = ["AOTU", "BU", "SPS"]
    
    highest_count = 0
    for i in bihem_types:
        relevant_regions = regions if i=="AOTU046" else regions[:1]
        for j, k in itertools.product(relevant_regions, sides):
            subregions = {}
            if j=="BU":
                subregions = bulb_subregions
            elif j=="AOTU":
                subregions = aotu_subregions
            avg_counts = get_region_avg_counts(i, j, k, 
                                               relevant_subregions=subregions)
            make_pie_chart(avg_counts, i, j, k, folder_path=[plot_folder])
            all_counts = [avg_counts[x][y] for x in avg_counts 
                          for y in ["pre","post"]]
            highest_count = max(highest_count, max(all_counts))
    for i,j in itertools.product(range(100, math.floor(highest_count)+500, 100),
                                 ("AOTU046","TuTu")):
        radius_factor = 1/1000 if j=="AOTU046" else 1/2000
        fig, ax = plt.subplots(figsize=fig_size)
        radius = i*radius_factor
        color = ["#FF0000"]
        plt.pie([i], colors=color, radius=radius)
        folder_path = [plot_folder, "Pie Charts", f"Legend Markers {j}"]
        figures.save_fig(fig, plot_name=i, folder_path=folder_path, 
                         file_type=figures.FileType.SVG,
                         transparent=True)
                                          
 
def get_bihem_weights(min_weight=0.05):
    """Gets the average ipsilateral and contralateral connection weight between
    bihemispheric neurons and MeTu and TuBu neurons in the AOTU.

    Parameters
    ----------
    min_weight : float, optional
        The lowest value that will be included in the returned dictionary. 
        The default is 0.05.

    Returns
    -------
    weight_dict : dict
        The weights of connections between neurons in the AOTU.
    """
    bihem_types = ["AOTU046", "TuTuB_a", "TuTuB_b"]
    
    mt_tb_types = [f"MeTu{x}" for x in ["1","2a","2b","3a","3b",
                                        "3c","4a","4b","4c","4d"]] + \
                  [f"TuBu{x}" for x in [f"0{y!s}" if y<10 else str(y) \
                                        for y in range(1,11)]]
    ipsi_dict = {}
    for i in "LR":
        temp_types = [f"{x}_{i}" for x in bihem_types] + [f"{x}_{i}" for x in mt_tb_types]
        ipsi_dict[i] = mapping.ConnectionMap(pre_types=temp_types,
                                             post_types=temp_types,
                                             region=f"AOTU_{i}")
    contra_dict = {}
    for i in "LR":
        opposite = "L" if i=="R" else "R"
        temp_types = [f"{x}_{opposite}" for x in bihem_types] + [f"{x}_{i}" for x in mt_tb_types]
        contra_dict[i] = mapping.ConnectionMap(pre_types=temp_types,
                                               post_types=temp_types,
                                               region=f"AOTU_{i}")
    ipsi_weights = (ipsi_dict["L"].type_weight_map + ipsi_dict["R"].type_weight_map) / 2
    contra_weights = (contra_dict["L"].type_weight_map + contra_dict["R"].type_weight_map) / 2
    all_types = bihem_types + mt_tb_types
    weight_dict = {}
    for indi, i in enumerate(all_types):
        for indj, j in enumerate(all_types):
            if indi > 2 and indj > 2:
                break
            temp_ipsi = ipsi_weights[indi][indj]
            temp_contra = contra_weights[indi][indj]
            if temp_ipsi >= min_weight:
                weight_dict[f"Ipsi {i} to {j}"] = temp_ipsi
            if temp_contra >= min_weight:
                weight_dict[f"Contra {i} to {j}"] = temp_contra
    return weight_dict
            

def get_bihem_weight_line_width(min_weight=0.05):
    """Gets the line width of arrows for the bihemispheric weight figure, using
    the formula y=2.5x+0.455.

    Parameters
    ----------
    min_weight : min_weight, optional
        The lowest value that will be queried in get_bihem_weights(). 
        The default is 0.05.

    Returns
    -------
    line_dict : dict
        The width of each line between neuron types.
    """
    weight_dict = get_bihem_weights(min_weight=min_weight)
    line_dict = {}
    for i in weight_dict:
        thickness = 2.5*weight_dict[i] + 0.455
        line_dict[i] = round(thickness, 1)
    return line_dict







