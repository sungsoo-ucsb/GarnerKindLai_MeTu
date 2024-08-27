#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 11:42:24 2023

@author: Dustin Garner
"""

from enum import Enum
import numpy as np
import pandas as pd
from pydantic import BaseModel
import matplotlib.pyplot as plt
import seaborn as sns
from fafbseg import flywire
import flycode.readfiles as readfiles
import flycode.utils as utils
import flycode.mapping as mapping
import flycode.reduction as reduction
import flycode.figures as figures
import flycode.flywire_functions as fw


font_size = 6
dot_size =0.4
dv_table = readfiles.import_file(file_name = "MeTu_dorsal_ventral")
dv_table = dv_table.rename(columns = {
        "den_count_me": "Medulla Postsynapse Count",
        "den_count_aotu": "AOTU Postsynapse Count",
        "pre_syn_count_me": "Medulla Presynapse Count",
        "pre_syn_count_aotu": "AOTU Presynapse Count",
        "column_count": "Column Count",
        "subtype": "MeTu Subtype",
        "dv": "Relative Offset from Medulla Centroid (V-D)",
        "ap": "Relative Offset from Medulla Centroid (A-P)",
        "eli_ratio": "Ellipse Ratio",
        "a": "Ellipse Major Axis Length (nm)",
        "b": "Ellipse Minor Axis Length (nm)",
        "aotu_centroid_y": "AOTU Dorsal-Ventral",
        "aotu_centroid_x": "AOTU Medial-Lateral"})
dv_table["ME Dorsal-Ventral"] = dv_table["Relative Offset from Medulla Centroid (V-D)"] * -1
dv_table["ME Anterior-Posterior"] = dv_table["Relative Offset from Medulla Centroid (A-P)"]
spacing = {
    "Medulla Postsynapse Count": [0, 1200],
    "AOTU Postsynapse Count": [0, 180],
    "Medulla Presynapse Count": [0, 1400],
    "AOTU Presynapse Count": [0, 800],
    "Column Count": [0, 120],
    "Relative Offset from Medulla Centroid (V-D)": [-1, 1],
    "Relative Offset from Medulla Centroid (A-P)": [-1, 1],
    "AOTU Dorsal-Ventral": [min(dv_table["AOTU Dorsal-Ventral"]),
                            max(dv_table["AOTU Dorsal-Ventral"])],
    "AOTU Medial-Lateral": [min(dv_table["AOTU Medial-Lateral"]),
                            max(dv_table["AOTU Medial-Lateral"])],
    "ME Dorsal-Ventral": [min(dv_table["ME Dorsal-Ventral"]),
                          max(dv_table["ME Dorsal-Ventral"])],
    "ME Anterior-Posterior": [min(dv_table["ME Anterior-Posterior"]),
                              max(dv_table["ME Anterior-Posterior"])]
    }
colors = {
    "MeTu1": "#e6194B",
    "MeTu2a": "#f58231",
    "MeTu2b": "#EBC400",
    "MeTu3a": "#bfef45",
    "MeTu3b": "#3cb44b",
    "MeTu3c": "#3E7748",
    "MeTu4a": "#42d4f4",
    "MeTu4b": "#7587C9",
    "MeTu4c": "#4400dd",
    "MeTu4d": "#911eb4"
    }
metu_types = list(colors.keys())
nts = ['gaba', 'acetylcholine', 'glutamate',
       'octopamine', 'serotonin', 'dopamine']
adjusted_name = {
    "Position Along Medulla Ventral-Dorsal Axis": \
        "Position Along Medulla\nVentral-Dorsal Axis"
    }

class FormatType(Enum):
    ME_POSITION = 1
    RETINOTOPY = 2


def compare_metu(column, save_figure=True, plot_folder=""):
    """Makes a strip plot comparing all MeTu subtypes based on their values
        in a column within dv_table.

    Parameters
    ----------
    column : str
        The column in dv_table.
    save_fig : bool, optional
        Whether to save the figure. The default is True.
    plot_folder : str, optional
        The folder to save the plots to. The default is "".
    """
    figures.strip_plot(data=dv_table, 
                       x_label="MeTu Subtype", 
                       y_label=column,
                       order=colors.keys(), 
                       hue="MeTu Subtype",
                       palette=colors.values(),
                       plot_name=column, 
                       folder_path=[plot_folder, "Comparison"], 
                       save_figure=save_figure,
                       dodge=False, 
                       show_means=True)


def scatter_plots(x, y, type_palette={}, plot_name="", plot_folder="", 
                  save_figure=True, fig_size = (1.9,1.5), 
                  format_type=FormatType.ME_POSITION):
    """Makes scatter plots comparing two columns of dv_table for each MeTu
        subtype, as well as a plot with all subtypes together and their
        line of best fits together as two different subplots.

    Parameters
    ----------
    x : str
        The column in dv_table to be plotted on the x-axis.
    y : str
        The column in dv_table to be plotted on the y-axis.
    type_palette : dict, optional
        Includes types you want included and their plotted colors. It defaults
        to the variable colors.
    plot_name : str, optional
        The name of the saved figure. The default is "".
    plot_folder : str, optional
        The folder which the plots get save to. The default is "".
    save_figure : bool, optional
        Whether to save the figure as a file. The default is True.
    fig_size : tuple, optional
        The size of the figure. The default is (1.9, 1.25).
    format_type : int, optional
        The format type of the scatter plots. This is one of the enum values of
        FormatType. The default is FormatType.ME_POSITION.
    """
    fig, ax = plt.subplots(figsize=fig_size)
    type_palette = colors if type_palette=={} else type_palette
    for i in type_palette:
        temp_df = dv_table[dv_table["MeTu Subtype"]==i]
        a, b = np.polyfit(temp_df[x], temp_df[y], 1)
        ax.plot(temp_df[x], a*temp_df[x]+b, c=type_palette[i], linewidth=0.5, zorder=0)
        plt.xlim(spacing[x][0], spacing[x][1])
        if y in spacing:
            plt.ylim(spacing[y][0], spacing[y][1])
        ax.scatter(x=temp_df[x], y=temp_df[y],
                    s=dot_size, c=type_palette[i], zorder=1)
        
    if format_type==FormatType.ME_POSITION:
        medulla = "Medulla"
        end_x_pos = x.find(medulla) + len(medulla)
        xlabel = f"{x[:end_x_pos]}\n{x[end_x_pos+1:]}"
        fig.text(0.565, -0.02, xlabel, ha='center', fontsize=font_size)
        fig.text(0.05, 0.4, y, ha='center', rotation=90, fontsize=font_size)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
    elif format_type==FormatType.RETINOTOPY:
        end_y_pos = y.find(" ") + 1
        y = f"{y[:end_y_pos]}\n{y[end_y_pos:]}"
        fig.text(0.5, 0.032, x, ha='center', fontsize=font_size)
        fig.text(0.025, 0.24, y, ha='center', rotation=90, fontsize=font_size)
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.yaxis.set_tick_params(labelleft=False)
        plt.tick_params(axis='both', which='both', length=0)
        ax.spines[['right', 'top']].set_visible(False)
    fig.tight_layout()
    if not save_figure:
        return
    folder_path = [plot_folder, "Scatter"]
    figures.save_fig(fig, plot_name=plot_name, folder_path=folder_path)



def nt_by_types(neur_types, plot_names=["NTs"], plot_folder="",
                region="Connectome", palette=["#000000"], save_figure=True,
                fig_size=(1.6,1.25), separate_plots=False):
    """Makes a plot of the neurotransmitters for the neuron types.

    Parameters
    ----------
    neur_types : list-like
        Neuron types in Neuron Spreadsheet.
    plot_names : str, optional
        What to name the plot. The default is "NTs".
    plot_folder : str, optional
        The name of the folder to save the plots to. The default is "".
    region : str, optional
        The region by which to limit the NT probabilities. 
        The default is "Connectome".
    palette : list-like, optional
        The colors to make the neuron types. The default is ["#000000"].
    save_figure : bool, optional
        Whether to save the figure. The default is True.
    fig_size : tuple, optional
        The size of the figure. The default is (1.6, 1.25).
    separate_plots : bool, optional
        Whether the neuron types should be plotted separately. 
        The default is False.

    Returns
    -------
    df : pd.DataFrame
        Each neuron's probability of each NT.

    """
    mapping.add_types(neur_types, output=False)
    
    nt_dict = {
        "Subtype": np.array([], dtype=str),
        "ID": np.array([], dtype=np.int64),
        "Neurotransmitter": np.array([], dtype=str),
        "Average Prediction": np.array([], dtype=float)
        }
    for i in neur_types:
        temp_df = fw.fetch_synapses(mapping.neur_ids[i], region=region,
                                        post=False, transmitters=True)
        temp_df.insert(0, "Subtype", i)
        
        for j in mapping.neur_ids[i]:
            temp_df2 = temp_df[temp_df["pre"]==j]
            for k in nts:
                percent = np.average(np.asarray(temp_df2[k]))
                nt_dict["Subtype"] = np.append(nt_dict["Subtype"], i)
                nt_dict["ID"] = np.append(nt_dict["ID"], j)
                k = k.upper() if k=="gaba" else k.capitalize()
                nt_dict["Neurotransmitter"] = \
                        np.append(nt_dict["Neurotransmitter"], k)
                nt_dict["Average Prediction"] = \
                        np.append(nt_dict["Average Prediction"], percent)
    df = pd.DataFrame(nt_dict)
    all_dfs = [df[df["Subtype"]==x] for x in neur_types] if separate_plots \
                                                         else [df]
    for indi, i in enumerate(all_dfs):
        temp_palette = [palette[indi]] if separate_plots else palette
        figures.strip_plot(data=i, 
                           x_label="Neurotransmitter",
                           y_label="Average Prediction",
                           order=None,
                           hue="Subtype",
                           palette=temp_palette,
                           plot_name=f"{plot_names[indi]}",
                           folder_path=[plot_folder],
                           dodge=False,
                           fig_size=fig_size,
                           size=2,
                           lim_range=(0,1),
                           save_figure=save_figure)
    return df


def get_mt4e():
    """Makes a dictionary of Other MeTu4a and Putative MeTu4e.

    Returns
    -------
    types : dict
        Keys are "Other MeTu4a" and "Putative MeTu4e", and the values are arrays
        of neurons..
    """
    mt4e_file = readfiles.import_file("MeTu4a_4e_candidates", file_type="csv")
    root_ids = np.array(mt4e_file.root_id)
    mt4e_file = mt4e_file.set_index("root_id")
    types = {"Other MeTu4a": [], "Putative MeTu4e": []}
    for i in root_ids:
        temp_type = mt4e_file.loc[i, "candidate_type"]
        if temp_type=="MeTu4a":
            types["Other MeTu4a"].append(i)
        else:
            types["Putative MeTu4e"].append(i)
    for i in types:
        types[i] = np.array(types[i])
    return types


def mt4e_partner_comparison(plot_name, plot_folder="", save_figure=False):
    """Makes a graph comparing putative MeTu4e partners to other MeTu4a.

    Parameters
    ----------
    plot_name : str
        The name of the exported plot.
    plot_folder : str, optional
        The name of the folder to save the plots to. The default is "".
    save_figure : bool, optional
        Whether to save the figure. The default is False.

    Returns
    -------
    df : pd.DataFrame
        The resulting comparison data.
    """
    mt4e_file = readfiles.import_file("MeTu4a_4e_candidates", file_type="csv")
    root_ids = np.array(mt4e_file.root_id)
    mt4e_file = mt4e_file.set_index("root_id")
    class NeurModel(BaseModel):
        root_id: int
        neur_type: str
        syn_type: str
        partner_type: str
        syn_count: int
    models = []
    partner_types = ['Dm2-IN', 'MeTu4a-IN', 'MTe01a-IN', 'Sm08-OUT', 'Sm43-IN',
           'Sm14-IN', 'LTe22-OUT', 'Li12-OUT', 'MC65-IN', 'LT55-OUT']
    for i in root_ids:
        temp_type = mt4e_file.loc[i, "candidate_type"]
        temp_type = "Other MeTu4a" if temp_type=="MeTu4a" else "Putative MeTu4e"
        for j in partner_types:
            syn_type = "Postsynaptic" if j[-2:]=="IN" else "Presynaptic"
            syn_count = mt4e_file.loc[i, j]
            temp_model = NeurModel(root_id=i,
                                   neur_type=temp_type,
                                   syn_type=syn_type,
                                   partner_type=j.split("-")[0],
                                   syn_count=syn_count)
            models.append(temp_model)
    df = utils.basemodel_list_to_df(models)
    df.rename(columns={"neur_type": "Putative Neuron Type",
                       "partner_type": "Synaptic Partner",
                       "syn_type": "Synapse Type",
                       "syn_count": "Synapse Count"},
              inplace=True)
    df.sort_values(by="Putative Neuron Type",
                   inplace=True)
    for i in [f"{x}synaptic" for x in ["Pre", "Post"]]:
        temp_df = df[df["Synapse Type"]==i]
        opposite = "Pre" if i=="Postsynaptic" else "Post"
        temp_partner_type = f"{opposite}synaptic Partner"
        temp_synapse_type = f"{i[:-3]}se Count"
        temp_df.rename(columns={"Synapse Count": temp_synapse_type,
                                "Synaptic Partner": temp_partner_type}, 
                       inplace=True)
        end_letter = "N" if i=="Postsynaptic" else "T"
        fig_width = 3.5 if i=="Postsynaptic" else 2.5
        figures.strip_plot(data=temp_df, 
                           x_label=temp_partner_type, 
                           y_label=temp_synapse_type,
                           hue="Putative Neuron Type",
                           order=[x.split("-")[0] for x in partner_types \
                                          if x[-1]==end_letter], 
                           plot_name=f"{plot_name} MeTu4 are {i}",
                           folder_path=[plot_folder, "MeTu4e"],
                           palette=["#1F77B4", "#FF7F0E"],
                           dodge=True,
                           save_figure=save_figure,
                           fig_size=(fig_width, 1.25))
    return df


def mt4e_dorsal_comparison(plot_name, plot_folder="", save_figure=True):
    """Makes a graph comparing putative MeTu4e to other MeTu4a by medulla
    position along the D-V axis.

    Parameters
    ----------
    plot_name : str
        The name of the exported plot.
    plot_folder : str, optional
        The name of the folder to save the plots to. The default is "".
    save_figure : bool, optional
        Whether to save the figure. The default is False.

    Returns
    -------
    dorsal_locs : pd.DataFrame
        The resulting dataframe.
    """
    mt4e = get_mt4e()
    dorsal_locs = readfiles.import_file("MeTu_dorsal_ventral")
    dorsal_locs = dorsal_locs[dorsal_locs.subtype=="MeTu4a"]
    new_types = []
    root_ids = fw.locs_to_segments(utils.coord_column_to_array(dorsal_locs.XYZ))
    for i in root_ids:
        temp_type = "Other MeTu4a" if i in mt4e["Other MeTu4a"] else "Putative MeTu4e"
        new_types.append(temp_type)
    dorsal_locs["Putative Type"] = new_types
    rename =   {"den_count_me": "Medulla Postsynapse Count",
                #"den_count_aotu": "AOTU Postsynapse Count",
                "pre_syn_count_me": "Medulla Presynapse Count",
                #"pre_syn_count_aotu": "AOTU Presynapse Count",
                "column_count": "Medulla Column Count",
                #"den_centroid_x": "Medulla x Position",
                #"den_centroid_y": "Medulla y Position",
                #"den_centroid_z": "Medulla z Position",
                #"den_centroid_x_m6_xform": "Medulla x M6 Xform",
                #"den_centroid_y_m6_xform": "Medulla y M6 Xform",
                #"den_centroid_z_m6_xform": "Medulla z M6 Xform",
                "dv": "Relative Offset from Medulla Centroid (D-V)",
                #"ap": "Relative Offset from Medulla Centroid (A-P)",
                #"aotu_centroid_x": "AOTU x Position",
                #"aotu_centroid_y": "AOTU y Position",
                #"aotu_centroid_z": "AOTU z Position",
                #"angle": "Ellipse Angle",
                #"angle_deg": "Ellipse Angle Degrees",
                #"eli_ratio": "Ellipse Ratio",
                "a": "Ellipse Major Axis Length (nm)",
                "b": "Ellipse Minor Axis Length (nm)"}
    dorsal_locs.rename(columns=rename, inplace=True)
    for i in list(rename.values()):
        figures.strip_plot(data=dorsal_locs,
                         x_label="Putative Type",
                         y_label = i,
                         hue="Putative Type",
                         order=["Other MeTu4a", "Putative MeTu4e"],
                         plot_name=f"{plot_name} {i}",
                         folder_path=[plot_folder, "MeTu4e"],
                         palette=["#1F77B4", "#FF7F0E"],
                         save_figure=save_figure,
                         dodge=False,
                         fig_size=(0.5,0.75))
    return dorsal_locs













