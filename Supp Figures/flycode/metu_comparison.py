#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 11:42:24 2023

@author: dustin
"""

import numpy as np
import pandas as pd
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
dot_size = 0.75
dv_table = readfiles.import_file(file_name = "MeTu_dorsal_ventral")
dv_table = dv_table.rename(columns = {
        "den_count_me": "Medulla Postsynapse Counts",
        "den_count_aotu": "AOTU Postsynapse Counts",
        "pre_syn_count_me": "Medulla Presynapse Counts",
        "pre_syn_count_aotu": "AOTU Presynapse Counts",
        "column_count": "Column Counts",
        "close_col": "Close Column Counts",
        "subtype": "MeTu Subtype",
        "dv": "Position Along Medulla Ventral-Dorsal Axis",
        "ap": "Position Along Medulla Posterior-Anterior Axis",
        "eli_ratio": "Ellipse Ratio"})
spacing = {
    "Medulla Postsynapse Counts": [0, 1200],
    "AOTU Postsynapse Counts": [0, 180],
    "Medulla Presynapse Counts": [0, 1400],
    "AOTU Presynapse Counts": [0, 800],
    "Column Counts": [0, 120],
    "Position Along Medulla Ventral-Dorsal Axis": [-1, 1],
    "Position Along Medulla Posterior-Anterior Axis": [-1, 1]
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
nts = ['gaba', 'acetylcholine', 'glutamate',
       'octopamine', 'serotonin', 'dopamine']
adjusted_name = {
    "Position Along Medulla Ventral-Dorsal Axis": \
        "Position Along Medulla\nVentral-Dorsal Axis"
    }


def compare_metu(column, save_figure=True):
    """Makes a strip plot comparing all MeTu subtypes based on their values
        in a column within dv_table.

    Parameters
    ----------
    column : str
        The column in dv_table.
    save_fig : bool, optional
        Whether to save the figure. The default is True.
    """
    plot_name = f"Fig S6/Comparison/{column}"
    figures.strip_plot(data=dv_table, x_label="MeTu Subtype", y_label=column,
                    order=colors.keys(), hue="MeTu Subtype",
                    palette=colors.values(),
                    plot_name=plot_name, save_figure=save_figure,
                    dodge=False, show_means=True)
    
    
def single_fig(df, x, y, neur_type, plot_folder):
    """Makes an individual scatter plot based on dv_table limited to a specific
        MeTu subtype.

    Parameters
    ----------
    df : pd.DataFrame
        DESCRIPTION.
    x : str
        The column in dv_table to be plotted on the x-axis.
    y : str
        The column in dv_table to be plotted on the y-axis.
    neur_type : str
        The MeTu type the plot is regarding.
    plot_folder : str, optional
        The folder within MeTu Comparison in which to put the plots. 
        The default is "".
    """
    fig, ax = plt.subplots(figsize=(0.75,0.75))
    a, b = np.polyfit(df[x], df[y], 1)
    ax.scatter(x=df[x], y=df[y],
                s=dot_size, c=colors[neur_type])
    ax.plot(df[x], a*df[x]+b, c=colors[neur_type])
    plt.xlim(spacing[x][0], spacing[x][1])
    plt.ylim(spacing[y][0], spacing[y][1])
    
    x_label = adjusted_name[x] if x in adjusted_name else x
    plt.xlabel(x_label, fontsize=font_size)
    plt.ylabel(y, fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    figures.save_fig(fig, f"MeTu Comparison/{plot_folder}{neur_type}")


def scatter_plots(x, y, plot_folder="", save_figure=True, save_singles=False):
    """Makes scatter plots comparing two columns of dv_table for each MeTu
        subtype, as well as a plot with all subtypes together and their
        line of best fits together as two different subplots.

    Parameters
    ----------
    x : str
        The column in dv_table to be plotted on the x-axis.
    y : str
        The column in dv_table to be plotted on the y-axis.
    plot_folder : str, optional
        The folder within MeTu Comparison in which to put the plots. 
        The default is "".
    save_figure : bool, optional
        Whether to save the figure as a file. The default is True.
    save_singles : bool, optional
        Whether to save the individual scatter plots.
    """
    fig, (ax1, ax2) = plt.subplots(2, figsize=(1.5,2.4))
    for i in colors:
        temp_df = dv_table[dv_table["MeTu Subtype"]==i]
        a, b = np.polyfit(temp_df[x], temp_df[y], 1)
        ax1.scatter(x=temp_df[x], y=temp_df[y],
                    s=dot_size, c=colors[i])
        ax2.plot(temp_df[x], a*temp_df[x]+b, c=colors[i])
        plt.xlim(spacing[x][0], spacing[x][1])
        plt.ylim(spacing[y][0], spacing[y][1])
        if not save_singles:
            continue
        single_fig(temp_df, x, y, neur_type=i, plot_folder=plot_folder)
    
    plt.sca(ax1)
    fig.text(0.5, 0.02, x, ha='center', fontsize=font_size)
    fig.text(0.05, 0.3, y, ha='center', rotation=90, fontsize=font_size)
    for i in (ax1, ax2):
        plt.sca(i)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
    fig.tight_layout()
    if not save_figure:
        return
    figures.save_fig(fig, f"Fig S6/Scatter/{plot_folder.strip('/')}")


def nt_by_types(neur_types, plot_names=["NTs"], region="Connectome", 
                palette=["#000000"], save_figure=True,
                fig_size=(1.6,1.25), separate_plots=False):
    """Makes a plot of the neurotransmitters for the neuron types.

    Parameters
    ----------
    neur_types : list-like
        Neuron types in Neuron Spreadsheet.
    plot_names : str, optional
        What to name the plot. The default is "NTs".
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
                           dodge=False,
                           fig_size=fig_size,
                           size=2,
                           lim_range=(0,1),
                           save_figure=save_figure)
    return df


