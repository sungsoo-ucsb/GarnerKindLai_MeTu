#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:05:50 2024

@author: Dustin Garner
"""


import os
import numpy as np
import pandas as pd
import flycode.mapping as mapping
import flycode.reduction as reduction
import flycode.readfiles as readfiles



#full_df = get_full_df()
#cx_df = get_cx_df(full_df)

all_tb_types = [f"0{x!s}" if x<10 else str(x) for x in range(1,11)]
tb_l = [f"TuBu{x}_L" for x in all_tb_types]
tb_r = [f"TuBu{x}_R" for x in all_tb_types]


def get_full_df():
    """Imports the spreadsheet with all neurons.

    Returns
    -------
    full_df : TYPE
        DESCRIPTION.
    """
    full_df = readfiles.import_file("Full783", file_type="feather")
    full_df.rename(columns={"pre_pt_root_id": "pre", "post_pt_root_id": "post"},
                   inplace=True)
    return full_df


def get_cx_df(full_df):
    """Gets the CX neurons from the full_df.

    Parameters
    ----------
    full_df : pd.DataFrame
        The dataframe retrieved from get_full_df().

    Returns
    -------
    cx_df : pd.DataFrame
        A dataframe with just the neurons in the CX.

    """
    info_df = readfiles.import_file("CX Neurons", sheet_name="CX with Outside")
    cx_neurs = np.asarray(info_df.neurs)
    cx_df = full_df[full_df.post.isin(cx_neurs)]
    cx_df = cx_df[(~(cx_df.neuropil.isin(["PB", "FB", "EB", "NO"])))| \
                  (cx_df.pre.isin(mapping.ids_from_types(tb_l+tb_r)))]
    cx_df = reduction.remove_min_synapses(cx_df, min_syns=5, pre_or_post="post")
    cx_df = reduction.remove_min_synapses(cx_df, min_syns=5, pre_or_post="pre")
    return cx_df







""" #Getting all the stuff


usable_cx_neurs = []
for i in cx_neurs:
    if len(cx_df[cx_df.post==i])==0:
        continue
    usable_cx_neurs.append(i)
#layer2 = pd.read_excel(desktop+"CX Neurons.xlsx", sheet_name="Layer2")
#layer2 = np.asarray(layer2.neurs)
layer2 = np.unique(cx_df.pre)

#XXX layer2_optic_l = np.asarray(pd.read_excel(desktop+"Layer 2 Left Optic Lobe Neurs.xlsx").neurs)
#XXX layer2_optic_r = np.asarray(pd.read_excel(desktop+"Layer 2 Right Optic Lobe Neurs.xlsx").neurs)
usable_layer2 = []
for i in layer2:
    if len(full_df[full_df.post==i])==0:
        continue
    if i in layer2_optic_l or i in layer2_optic_r:
        continue
    usable_layer2.append(i)
usable_layer2 = np.unique(usable_layer2)

layer3_df = full_df[full_df.post.isin(usable_layer2)]
layer3_df = reduction.remove_min_synapses(layer3_df, min_syns=5)
layer3 = np.unique(layer3_df.pre)

next_df = full_df[full_df.post.isin(layer3)]
next_df_l = next_df[next_df.neuropil.isin([f"{x}_L" for x in ["ME", "LO", "LOP"]])]
next_df_r = next_df[next_df.neuropil.isin([f"{x}_R" for x in ["ME", "LO", "LOP"]])]

counter = 0
left_neurs = []
right_neurs = []
both_side_neurs = []
for i in layer3:
    if counter%1000==0:
        print(counter//1000, time.time()-start)
        start = time.time()
    counter += 1
    
    left_df = next_df_l[next_df_l.post==i]
    right_df = next_df_r[next_df_r.post==i]

    if len(left_df)>0 and len(right_df)>0:
        both_side_neurs.append(i)
    elif len(left_df)>0:
        left_neurs.append(i)
    elif len(right_df)>0:
        right_neurs.append(i)


layer2_df_pre_optic_l = layer2_df[layer2_df.pre.isin(layer3_optic_l)]
layer2_df_pre_optic_r = layer2_df[layer2_df.pre.isin(layer3_optic_r)]
layer2_weights = {x: [] for x in ["neur", "right_weight", "left_weight", "total_weight"]}
for i in layer2:
    if len(layer2_df[layer2_df.post==i])==0:
        continue
    if i in layer2_optic_l or i in layer2_optic_r:
        continue
    temp_df = layer2_df[layer2_df.post==i]
    optic_l_df = layer2_df_pre_optic_l[layer2_df_pre_optic_l.post==i]
    optic_r_df = layer2_df_pre_optic_r[layer2_df_pre_optic_r.post==i]
    left_weight = len(optic_l_df) / len(temp_df)
    right_weight = len(optic_r_df) / len(temp_df)
    layer2_weights["neur"].append(i)
    layer2_weights["right_weight"].append(right_weight)
    layer2_weights["left_weight"].append(left_weight)
    layer2_weights["total_weight"].append(right_weight+left_weight)
layer2_weights["neur"] = np.asarray(layer2_weights["neur"], dtype=str)
for i in [f"{x}_weight" for x in ["left", "right", "total"]]:
    layer2_weights[i] = np.asarray(layer2_weights[i])
layer2_weights_df = pd.DataFrame(layer2_weights)



weight_dicts = {x: [] for x in ["neur", "left_weight", "right_weight", 
                                "Optic Lobe Weight", "neur_type", "CX Neuron Type"]}
all_rings = mapping.ids_from_types(er_all)
for neur in cx_neurs:
    temp_df = cx_df[cx_df.post==neur]
    weight_dict = {"L": 0, "R": 0}
    for i in np.unique(temp_df.pre):
        if not i in layer2:
            continue
        temp_df2 = temp_df[temp_df.pre==i]
        temp_weight = len(temp_df2) / len(temp_df)
        if i in layer2_optic_l:
            weight_dict["L"] += temp_weight
            continue
        elif i in layer2_optic_r:
            weight_dict["R"] += temp_weight
            continue
        for j in "LR":
            weight_dict[j] += temp_weight * layer2_weights[i][j]
    weight_dicts["neur"].append(str(neur))
    weight_dicts["left_weight"].append(weight_dict["L"])
    weight_dicts["right_weight"].append(weight_dict["R"])
    weight_dicts["Optic Lobe Weight"].append(weight_dict["L"] + weight_dict["R"])
    weight_dicts["neur_type"].append("ring" if neur in all_rings else "not ring")
    neur_type = mapping.find_neur_type(neur, 
                                       search=er_all_all+exr+
                                               [f"{x}_{y}" for x in ["FB8B", "IbSpsP", 
                                                                     "LNO2", "OA-AL2i1"]
                                                for y in "LR"], 
                                       unidentified_type="Unidentified")
    neur_type = neur_type if neur_type=="Unidentified" else neur_type[:-2]
    weight_dicts["CX Neuron Type"].append(neur_type)

weight_df = pd.DataFrame(weight_dicts)
weight_df = weight_df.sort_values("Optic Lobe Weight", ascending=False, ignore_index=True)
weight_df.to_excel(desktop+"CX Weights.xlsx")
"""

#####################################3

""" #This is the way to get the percentages
full_df = pd.read_csv(desktop+"Full783.csv", header=None, names=["id","pre_pt_root_id",
                                                                 "post_pt_root_id",
                                                                 "connection_score",
                                                                 "cleft_score",
                                                                 "gaba","ach","glut",
                                                                 "oct","ser","da",
                                                                 "pre_pt_supervoxel_id",
                                                                 "post_pt_supervoxel_id",
                                                                 "neuropil",
                                                                 "post_pt_position_x",
                                                                 "post_pt_position_y",
                                                                 "post_pt_position_z",
                                                                 "pre_pt_position_x",
                                                                 "pre_pt_position_y",
                                                                 "pre_pt_position_z"])
full_df = full_df.rename(columns = {"pre_pt_root_id": "pre", "post_pt_root_id": "post"})
cx_neurs = np.asarray(pd.read_excel(desktop+"CX Neurons.xlsx", 
                                    sheet_name="CX with Outside").neurs)

cx_df = full_df[full_df.post.isin(cx_neurs)]
cx_df = cx_df[(~(cx_df.neuropil.isin(["PB", "FB", "EB", "NO"])))| \
              (cx_df.pre.isin(mapping.ids_from_types(tb_l+tb_r)))]
cx_df = reduction.remove_min_synapses(cx_df, min_syns=5, pre_or_post="post")
cx_df = reduction.remove_min_synapses(cx_df, min_syns=5, pre_or_post="pre")

layer2 = np.asarray(pd.read_excel(desktop+"CX Neurons.xlsx", sheet_name="Layer2").neurs)
layer2_optic_l = np.asarray(pd.read_excel(desktop+"CX Neurons.xlsx", 
                                          sheet_name="Layer2 Optic Left").neurs)
layer2_optic_r = np.asarray(pd.read_excel(desktop+"CX Neurons.xlsx", 
                                          sheet_name="Layer2 Optic Right").neurs)
layer2_df = full_df[full_df.post.isin(layer2)]

layer3_optic_l = np.asarray(pd.read_excel(desktop+"CX Neurons.xlsx", 
                                          sheet_name="Layer3 Optic Left").neurs)
layer3_optic_r = np.asarray(pd.read_excel(desktop+"CX Neurons.xlsx", 
                                          sheet_name="Layer3 Optic Right").neurs)
layer2_df_pre_optic_l = layer2_df[layer2_df.pre.isin(layer3_optic_l)]
layer2_df_pre_optic_r = layer2_df[layer2_df.pre.isin(layer3_optic_r)]

layer2_weights_df = pd.read_excel(desktop+"CX Neurons.xlsx", sheet_name="Layer2 Weights")
layer2_weights = {}
for i, l, r in zip(np.asarray(layer2_weights_df.neur),
                   np.asarray(layer2_weights_df.left_weight),
                   np.asarray(layer2_weights_df.right_weight)):
    layer2_weights[i] = {"L": l, "R": r}


weight_df = pd.read_excel(desktop+"CX Neurons.xlsx", sheet_name="CX Weights")




def get_layer2(neur, min_weight=0.1):
    temp_df = cx_df[cx_df.post==neur]
    weight_dict = {x: 0 for x in np.unique(temp_df.pre)}
    for i in np.unique(temp_df.pre):
        if not i in layer2:
            continue
        temp_df2 = temp_df[temp_df.pre==i]
        temp_weight = len(temp_df2) / len(temp_df)
        if i in layer2_optic_l or i in layer2_optic_r:
            weight_dict[i] = temp_weight
            continue
        for j in "LR":
            weight_dict[i] += temp_weight * layer2_weights[i][j]
    weight_dict = dict(sorted(weight_dict.items(), key = lambda x: x[1], reverse = True))
    for i in weight_dict.copy():
        if weight_dict[i]<min_weight:
            del weight_dict[i]
    return weight_dict


def get_layer3(neur, min_weight=0.1):
    temp_df = layer2_df[layer2_df.post==neur]
    weight_dict = {x: 0 for x in np.unique(temp_df.pre)}
    
    for i in np.unique(temp_df.pre):
        hemis = ""
        if not (i in layer3_optic_l or i in layer3_optic_r):
            continue
        temp_df2 = temp_df[temp_df.pre==i]
        temp_weight = len(temp_df2) / len(temp_df)
        weight_dict[i] += temp_weight
    weight_dict = dict(sorted(weight_dict.items(), key = lambda x: x[1], reverse = True))
    for i in weight_dict.copy():
        if weight_dict[i]<min_weight:
            del weight_dict[i]
    return weight_dict


def get_type_weights(neur):
    weights = get_layer2(neur, 0)
    weight_dict = collections.defaultdict(int)
    for i in weights:
        neur_type = mapping.find_neur_type(i, search=avp_all, unidentified_type="Unidentified")
        weight_dict[neur_type] += weights[i]["L"] + weights[i]["R"]
    weight_dict = dict(sorted(weight_dict.items(), key = lambda x:x[1], reverse=True))
    return weight_dict




def strip_optic(neur_types, pre_search, palette=["#e6194B", "#3E7748", "#EBC400", "#4400dd"],
                include_unidentified=True, fig_size=(4.0,1.5), plot_name="Optic Plot"):
    weight_dicts = {}
    for i in mapping.ids_from_types(neur_types):
        temp_layer2 = get_layer2(i, 0)
        temp_weight_dict = collections.defaultdict(int)
        for j in temp_layer2:
            neur_type = mapping.find_neur_type(j, 
                                               search=pre_search, 
                                               unidentified_type="Unidentified")
            if not include_unidentified and neur_type=="Unidentified":
                continue
            temp_weight_dict[neur_type] += temp_layer2[j]
        for j in temp_weight_dict.copy():
            if temp_weight_dict[j] < 0.01:
                del temp_weight_dict[j]
        weight_dicts[i] = temp_weight_dict
    pre_types = []
    for i in weight_dicts:
        pre_types += list(weight_dicts[i].keys())
    pre_types = list(np.unique(pre_types))
    weight_df = {x: [] for x in ["neur", "Neur Type", "Presynaptic Type", "Optic Lobe Weight"]}
    for i in mapping.ids_from_types(neur_types):
        for j in pre_types:
            weight_df["neur"].append(i)
            weight_df["Neur Type"].append(mapping.find_neur_type(i, 
                                                            search=pre_search))
            weight_df["Presynaptic Type"].append(j)
            weight_df["Optic Lobe Weight"].append(weight_dicts[i][j])
    weight_df = pd.DataFrame(weight_df)
    
    figures.strip_plot(weight_df, 
                       "Presynaptic Type", 
                       "Optic Lobe Weight", 
                       pre_types, 
                       "neur", 
                       palette, 
                       plot_name, 
                       dodge=False,
                       fig_size=fig_size)
    return weight_df




def exr_strip_optic(palette=["#e6194B", "#3E7748", "#EBC400", "#4400dd"],
                include_unidentified=True, fig_size=(4.0,1.5), plot_name="Optic Plot"):
    weight_dicts = {}
    neur_types = ["ExR1_R", "ExR1_L"]
    for i in mapping.ids_from_types(neur_types):
        temp_layer2 = get_layer2(i, 0)
        temp_weight_dict = collections.defaultdict(int)
        for j in temp_layer2:
            neur_type = mapping.find_neur_type(j, 
                                               search=aotu46+tb_l+tb_r,
                                               unidentified_type="Unidentified")
            if neur_type in [f"TuBu{x}_{y}" for x in ["04","05","08","09","10"] for y in "LR"]:
                neur_type = "Other TuBu"
            if not include_unidentified and neur_type=="Unidentified":
                continue
            temp_weight_dict[neur_type] += temp_layer2[j]
        #for j in temp_weight_dict.copy():
        #    if temp_weight_dict[j] < 0.01:
        #        del temp_weight_dict[j]
        weight_dicts[i] = temp_weight_dict
    pre_types = []
    for i in weight_dicts:
        pre_types += list(weight_dicts[i].keys())
    pre_types = [f"AOTU046_{x}" for x in "LR"] + [f"TuBu0{x}_{y}" for x in "12367" for y in "LR"] + \
                ["Other TuBu", "Unidentified"]
    weight_df = {x: [] for x in ["neur", "Neur Type", "Presynaptic Type", "Optic Lobe Weight"]}
    for i in mapping.ids_from_types(neur_types):
        for j in pre_types:
            weight_df["neur"].append(i)
            weight_df["Neur Type"].append(mapping.find_neur_type(i, 
                                                            search=pre_search))
            weight_df["Presynaptic Type"].append(j)
            weight_df["Optic Lobe Weight"].append(weight_dicts[i][j])
    weight_df = pd.DataFrame(weight_df)
    
    figures.strip_plot(weight_df, 
                       "Presynaptic Type", 
                       "Optic Lobe Weight", 
                       pre_types, 
                       "neur", 
                       palette, 
                       plot_name, 
                       dodge=False,
                       fig_size=fig_size)
    return weight_df




exr_strip_optic()

pre_search = [f"{x}_{y}" for x in 
              ["H1","MeLoPlp_unknown_4","PS125","Olt","MeLoLo_unknown_1",
               "PlpSpsSpsGng_unknown_1"] for y in "LR"]
for i in [True, False]:
    plot_name = f"OA-AL2i1 Optic Weight{' No Unidentified' if not i else ''}"
    strip_optic(["OA-AL2i1_L", "OA-AL2i1_R"], 
                pre_search, 
                palette=["#4400dd","#e6194B",],
                include_unidentified=i,
                plot_name = plot_name)

strip_optic(["FB8B_L", "FB8B_R"],
            ["MePlpMb_unknown_1_R","MePlpMb_unknown_1_L",],
            palette=["#DD00C9","#EBC400", "#3E7748",  "#4400dd","#e6194B"],
            fig_size=(2.0,1.5),
            plot_name = "FB8B Optic Weight")





relevant_er = [x[:-2] for x in er_all_all if x[-1]=="R"]
relevant_exr = [x[:-2] for x in exr if x[-1]=="R"]
extra_types = [f"{x}" for x in ["FB8B", "OA-AL2i1", "IbSpsP", "LNO2",]]
figures.strip_plot(weight_df, "CX Neuron Type", "Optic Lobe Weight", 
                     relevant_er+relevant_exr+extra_types
                     +["Unidentified"],
                     None, 
                     #["#4400dd" for _ in relevant_er]+
                     #["#e6194B" for _ in relevant_exr]+
                     #["#3E7748" for _ in extra_types]+
                     #["#EBC400"], 
                     ["#e6194B"]+
                     ["#4400dd" for _ in range(4)]+
                     ["#e6194B" for _ in range(2)]+
                     ["#4400dd" for _ in range(10)]+
                     ["#e6194B"]+
                     ["#4400dd"]+
                     ["#e6194B" for _ in range(7)]+
                     ["#3E7748" for _ in extra_types]+
                     ["#EBC400"],
                     "CX Neuron Optic Lobe Weights", 
                     False, 
                     fig_size=(8.0, 2),
                     show_means=True,
                     mean_type="_")


"""

