#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:44:44 2024

@author: Dustin Garner
"""


import os
import numpy as np
import pandas as pd
import functools
from enum import Enum
import matplotlib.pyplot as plt
import seaborn as sns


border = 0
text_pad = 0
label_pad = 2
font_size = 6
xtick_rotation = 45
ha = "right"
edge_color = "white"
rotation_mode = "anchor"
loc = "lower left"
bbox = (1.04, 0)


class FileType(Enum):
    PDF = 1
    PNG = 2
    SVG = 3


def verify_folder_path(folder_path):
    """Creates a directory in the order of the passed list.
    
    Parameters
    ----------
    folder_path : list-like[str]
        A list of directory names. If they do not exist, they will be created.
        
    Returns
    -------
    directory : str
        The resulting directory based on the given folder_path.
    """
    directory = os.path.dirname(__file__)
    directory = os.path.join(directory, "Generated-Figures")
    if isinstance(folder_path, str):
        folder_path = [folder_path]
    for i in folder_path:
        if i == "":
            continue
        subfolders = [x.path for x in os.scandir(directory) if x.is_dir()]
        directory = os.path.join(directory, i)
        if not directory in subfolders:
            os.mkdir(directory)
    return directory


def get_figure_path(figure_name, folder_path=[], file_type=FileType.PDF):
    """Makes an exportable figure path, depending on the file type.
    
    Parameters
    ----------
    figure_name : str
        The name of the figure.
    folder_path : list-like, optional
        The folder path the figure will be saved to. The default is [].
    file_type : int, optional
        One of the values of the FileType enum. The default is FileType.PDF.
    
    Returns
    -------
    str
        file_path to save the figure to.
    """
    absolute_path = verify_folder_path(folder_path)
    if file_type==FileType.PDF:
        extension = "pdf"
    elif file_type==FileType.PNG:
        extension = "png"
    elif file_type==FileType.SVG:
        extension = "svg"
    figure_name = f"{figure_name}.{extension}"
    return os.path.join(absolute_path, figure_name)


def save_fig(fig, plot_name, folder_path=[], file_type=FileType.PDF, dpi=300, 
             transparent=False):
    """Saves a figure with the given file path.
    
    Parameters
    ----------
    fig : matplotlib.pyplot.figure
        The figure to be saved.
    plot_name : str
        The name of the file upon being exported.
    folder_path : list-like, optional
        The folder path the figure will be saved to. The default is [].
    file_type : int, optional
        One of the values of the FileType enum. The default is FileType.PDF.
    dpi : int, optional
        The dpi. The default is 300.
        folder_path : list-like, optional
            The folder path the figure will go into.
    transparent : bool, optional
        Whether to export the file with a transparent background. The default
        is False.
    """
    file_path = get_figure_path(plot_name, folder_path=folder_path, 
                                file_type=file_type)
    fig.savefig(file_path, 
                dpi=dpi, 
                bbox_inches='tight', 
                transparent=transparent)
    print(f"Figure saved at {file_path}")


def add_legend(legend=None, marker_scale=1.2):
    """Makes a legend for a figure.

    Parameters
    ----------
    legend : matplotlib.pyplot.Legend, optional
        The pre-existing legend, if there is one. The default is None.
    marker_scale : float, optional
        The size of legend marker scales. The default is 1.2.

    Returns
    -------
    leg : matplotlib.pyplot.Legend
        The resulting legend.

    """
    plt_legend = functools.partial(plt.legend,
                                   fontsize=font_size,
                                   bbox_to_anchor=bbox,
                                   loc=loc, 
                                   borderaxespad=border, 
                                   edgecolor=edge_color,
                                   handletextpad=text_pad, 
                                   markerscale=marker_scale)
    leg = plt_legend() if legend==None else plt_legend(legend)
    return leg


def create_strip_plot(data, x_label, y_label, order, hue, palette, plot_name, 
               folder_path=[], dodge=True, save_figure=True, fig_size=(2.0, 1.25),
              show_means=False, mean_type=".", size=2.5, lim_range=None):
    """Makes a strip plot.

    Parameters
    ----------
    data : pd.DataFrame
        The data.
    x_label : str
        The name of the data on the x-axis.
    y_label : str
        The name of the data on the y-axis.
    order : list-like
        The order of the data columns.
    hue : str
        What to base the different colors on.
    palette : list-like
        The colors.
    plot_name : str
        The name of the plot to be exported.
    folder_path : list-like, optional
        The folder path the figure will be saved to. The default is [].
    dodge : bool, optional
        Whether the hue strips are separated. The default is True.
    save_figure : bool, optional
        Whether the figure will be saved. The default is True.
    fig_size : tuple, optional
        The size of the figure. The default is (2.0, 1.25).
    show_means : bool, optional
        Whether to show the means. The default is False.
    mean_type : str, optional
        The type of mean line to have. The default is ".".
    size : float, optional
        Radius of the markers. The default is 2.5.
    lim_range : list-like, optional
        The range by which to limit the plot values. The default is None.
    """
    
    fig, ax = plt.subplots(figsize=fig_size)
    strip = sns.stripplot(data=data,
                          x=x_label, 
                          y=y_label, 
                          hue=hue, 
                          order=order,
                          palette=palette,
                          dodge=dodge, 
                          size=size)
    plt.xlabel(x_label,
               fontsize=font_size)
    plt.ylabel(y_label,
               fontsize=font_size)
    plt.xticks(rotation=xtick_rotation,
               ha=ha,
               rotation_mode=rotation_mode,
               fontsize=font_size)
    plt.yticks(fontsize=font_size)
    
    if not lim_range is None:
        plt.ylim(lim_range[0], lim_range[1])
    
    if show_means:
        sns.boxplot(showmeans=True,
                    meanprops={"markerfacecolor": "k",
                               "markeredgecolor": "k",
                               "marker": mean_type,
                               "markersize": 4},
                    medianprops={'visible': False},
                    whiskerprops={'visible': False},
                    zorder=10,
                    x=x_label,
                    y=y_label,
                    data=data,
                    showfliers=False,
                    showbox=False,
                    showcaps=False,
                    ax=strip)
    
    add_legend(marker_scale=1.6)
    
    ax.spines[['right', 'top']].set_visible(False)
    
    if save_figure:
        save_fig(fig, plot_name=plot_name, folder_path=folder_path)

strip_plot = create_strip_plot


def create_bar_graph(data, x_label, y_label, x_ticks, y_ticks, colors,
                     fig_size, plot_name, save_figure, folder_path=[],
                     color_axis=0):
    """Takes in bar graph data and makes a graph.

    Parameters
    ----------
    data : dict
        Keys are sets of sets of bars, one for each within a bar cluster.
        Values are 2D numpy arrays: one dimension for bars stacked on top of
        one another, the second for values of each bar.
    x_label : str
        The x-axis label of the graph.
    y_label : str
        The y-axis label of the graph.
    x_ticks : list-like
        The names of clusters of bars.
    y_ticks : np.array
        The spacing of number values on the y-axis.
    colors : list
        Values are the colors that each bar corresponds to.
    fig_size : tuple
        The size of the figure.
    plot_name : str
        The name of the outputted plot.
    save_figure : bool
        Whether to save the figure.
    folder_path : list-like, optional
        The folder path the figure will be saved to. The default is [].
    color_axis : int, optional
        Along what axis the color scheme occurs. The default is 0.
    """
    fig, ax = plt.subplots(figsize = fig_size)
    
    space_multiplier = 2
    bar_count = len(list(data.values())[0][0])
    bar_spacing = np.arange(0, bar_count*space_multiplier, 
                            space_multiplier)
    width = 0.35
    spacing = [x+(-0.5*(len(data)-1)) for x in range(len(data))]
    bar_placement = {x: y*width for x,y in zip(data.keys(), spacing)}
    bottom = {x: np.zeros(bar_count) for x in data.keys()}
    
    bars = {}
    color = "#808080"
    for indi, i in enumerate(data):
        if color_axis==0:
            color = colors[indi]
        for indj, j in enumerate(data[i]):
            if color_axis==1:
                color = colors[indj]
            bars[i] = ax.bar(bar_spacing+bar_placement[i],
                             j,
                             width,
                             label=i,
                             color=color,
                             bottom=bottom[i])
            bottom[i]+=j
    
    ax.set_xlabel(x_label, fontsize = font_size, labelpad = label_pad)
    ax.set_ylabel(y_label, fontsize = font_size, labelpad = label_pad)
    ax.set_xticks(bar_spacing, x_ticks, rotation = xtick_rotation, ha = ha, 
               rotation_mode = rotation_mode, fontsize = font_size)
    ax.set_yticks(y_ticks, y_ticks, fontsize = font_size)
    ax.spines[['right','top']].set_visible(False)
    add_legend()
    
    if save_figure:
        save_fig(fig, plot_name=plot_name, folder_path=folder_path)





