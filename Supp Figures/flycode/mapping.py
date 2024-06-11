<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 16:24:03 2022

@author: dusti
"""


import os
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fafbseg
from fafbseg import flywire
import flycode.reduction as reduction
import flycode.utils as utils
import flycode.readfiles as readfiles
import flycode.figures as figures
import flycode.flywire_functions as fw


neur_coords, prev_ids = readfiles.import_coords()
neur_ids = {}
new_colors = readfiles.import_colors()


def add_types(types, output=False):
    """Adds types to the dictionary neur_ids
    
    Parameters
    ----------
    types : list-like
        List of types that are in neur_coords. They will be added to neur_ids 
        in the same order.
    output : bool, Optional
        Whether to show if neurons in Neuron Spreadsheet have been updated
        since the IDs were last added.
    """
    not_updated = np.array([], dtype=np.int64)
    not_in_spread = np.array([], dtype=str)
    for i in types:
        if i not in neur_coords:
            print(f"{i} is not a valid neuron type.")
            continue
        elif i in neur_ids:
            pass
        else:
            neur_ids[i] = fw.locs_to_segments(neur_coords[i])
        if not i in prev_ids:
            not_in_spread = np.append(not_in_spread, i)
            continue
        for j, k in zip(neur_ids[i], prev_ids[i]):
            if not j == k:
                not_updated = np.append(not_updated, k)
    if output:
        print(f"No neurons in spreadsheet: {not_in_spread}\n")
        print(f"Not updated: {not_updated}\n")


def ids_from_types(types):
    """Returns a list of ids from the given types.
    
    Parameters
    ----------
    types: list-like 
        List of types that are in neur_coords.
    
    Returns
    -------
    ids: np.array
        Array of neur ids from types.
    """
    add_types(types, output=False)
    return np.concatenate(([neur_ids[x] for x in types]))


def find_neur_type(neur_id, unidentified_type=""):
    """Finds what neuron type the neuron ID is, assuming it is in neur_ids.

    Parameters
    ----------
    neur_id : int
        The neuron ID.
    unidentified_type : str, optional
        The type that is returned if the ID is not in neur_ids. The default is "".

    Returns
    -------
    str
        The type or unidentified type if ID is not in neur_ids.

    """
    for i in neur_ids:
        if int(neur_id) in neur_ids[i]:
            return i
    return unidentified_type


def check_existence_in_region(neur_types, regions, min_syns=10):
    """Checks whether any neurons from the given types do not have at least 
        min_syns in the given regions.
    
    Parameters
    ----------
    neur_types : list-like
        Neuron types in Neuron Spreadsheet.
    regions : list-like
        The regions to test.
    min_syns : int, Optional
        Below this number the neuron is displayed in output.
    """
    ids = ids_from_types(neur_types)
    syn_df = fw.fetch_synapses(ids)
    for i in ids:
        temp_df = syn_df[(syn_df["pre"]==i) | (syn_df["post"]==i)]
        for j in regions:
            temp_df2 = reduction.lim_region(temp_df, j)
            number_syns = len(temp_df2)
            if number_syns<min_syns:
                print(f"{i!s} has {number_syns!s} connections in {j}")


def find_partners(neurons, region, pre_or_post, min_neurs=5):
    """Finds the pre_or_post partners of all of the given neurons. Returns a 
        dictionary of the partner with number of connection and an array of 
        just the partners.
    
    Parameters
    ----------
    neurons : list-like
        Neuron IDs.
    region : str
        The region to limit the search.
    pre_or_post : str
        Whether pre or postsynaptic connections are wanted
    min_neurs : int, Optional
        The lowest number of connections to be shown in the output.
        
    Returns
    -------
    neur_dict : dict
        A dictionary. Keys are IDs and values are the number of connections from
        the input neurons.
    """
    pre = False if pre_or_post=="pre" else True
    post = not pre
    syn_df = fw.fetch_synapses(neurons, pre=pre, post=post)
    syn_df = reduction.lim_region(syn_df, region)
    all_partners = np.unique(np.array(syn_df[pre_or_post]))
    for i in all_partners:
        temp_df = syn_df[syn_df[pre_or_post]==i]
        if len(temp_df)<min_neurs:
            syn_df = syn_df[syn_df[pre_or_post]!=i]
    new_partners = np.array(syn_df[pre_or_post])
    neur_dict = utils.count_instances(new_partners, min_neurs=min_neurs)
    return neur_dict
    

def partner_df(ids, region="Connectome", pre=True, post=True, 
           autapses=False, min_synapses=1):
    """Retrieves the partner spreadsheet from Flywire.
    
    Parameters
    ----------
    ids : int or list of int
        Neurons that you want to find the partners of.
    region : str
        The region to limit the map to.
    pre : bool
        Whether you want the neurs to be pre.
    post : bool
        Whether you want the neurs to be post.
    autapses : bool
        If False, removes autapses.
    min_synapses : int
        Limits the df to this number of neurons.

    Returns
    -------
    syn_df: pd.DataFrame
        The resulting synapse dataframe.
    """
    syn_df = fw.fetch_synapses(ids, pre=pre, post=post)
    syn_df = reduction.lim_region(syn_df, region)
    if not autapses:
        syn_df = reduction.remove_autapses(syn_df)
    syn_df = reduction.remove_min_synapses(syn_df, minimum_syns=min_synapses)
    return syn_df


def make_syn_maps(df, pre_neurs, post_neurs, min_synapses=3, exclude=[]):
    """Makes connectivity and weight maps given a dataframe from Flywire.

    Parameters
    ----------
    df : pandas.DataFrame
        A neuron dataframe given by flywire.get_neurons.
    pre_neurs : list-like
        The pre-synaptic neurons.
    post_neurs : list-like
        The post-synaptic neurons.
    min_synapses : int, optional, the default is 3.
        The minimum number of synapses when making the minimized weight map.
    exclude : list-like, optional
        The neurons to be excluded from being counted in the plot.

    Returns
    -------
    syn_map : np.array
        Connectivity map of number of synapses.
    weight_map : np.array
        Weight map of the percent of connections from the presynaptic neuron
        to the given postsynaptic neuron.
    min_map : np.array
        Same as weight map but only if the neurons have greater than min_synapses
        total.
    """
    pre_count = len(pre_neurs)
    post_count = len(post_neurs)
    syn_map = np.zeros((pre_count, post_count), dtype=float)
    weight_map = np.zeros((pre_count, post_count), dtype=float)
    min_map = np.zeros((pre_count, post_count), dtype=float)
    
    for indi, i in enumerate(post_neurs):
        temp_df = df[df["post"]==i]
        total_connections = len(temp_df)
        if total_connections==0: #Prevents from dividing by zero
            print(f"{i!s} does not have connections.")
            continue
        for indj, j in enumerate(pre_neurs):
            partner_count = 0 if (i in exclude or j in exclude) else \
                    len(temp_df[temp_df["pre"]==j])
            synaptic_weight = partner_count/total_connections
            min_weight = synaptic_weight if total_connections>=min_synapses else 0
            syn_map[indj][indi] = partner_count
            weight_map[indj][indi] = synaptic_weight
            min_map[indj][indi] = min_weight
    return syn_map, weight_map, min_map


def make_type_maps(df, pre_types, post_types):
    """Makes type connectivity and weight maps given a dataframe from Flywire.

    Parameters
    ----------
    df : pandas.DataFrame
        A neuron dataframe given by flywire.get_neurons.
    pre_types : list-like
        The pre-synaptic types.
    post_types : list-like
        The post-synaptic types.

    Returns
    -------
    type : np.array
        Connectivity map of total number of synapses for the given types.
    weight_map : np.array
        Weight map of the percent of connections from the presynaptic type
        to the given postsynaptic type.
    """
    add_types(pre_types)
    add_types(post_types)
    type_map = np.zeros((len(pre_types), len(post_types)), dtype=np.uint64)
    weight_map = np.zeros((len(pre_types), len(post_types)), dtype=float)
    for indi, i in enumerate(post_types):
        temp_df = df[df["post"].isin(neur_ids[i])]
        total_connections = len(temp_df)
        for indj, j in enumerate(pre_types):
            type_connections = len(temp_df[temp_df["pre"].isin(neur_ids[j])])
            type_map[indj][indi] = type_connections
            weight_map[indj][indi] = type_connections/total_connections \
                    if total_connections!=0 else 0
    return type_map, weight_map


class ConnectionMap: 
    #Creates maps between pre types and post types in a region.
    def __init__(self, pre_types, post_types, region, min_synapses=3, \
                 exclude=[], font_size=6):
        """Creates maps between pre types and post types in a region.

        Parameters
        ----------
        pre_types : list-like
            The types with presynaptic neurons.
        post_types : list-like
            The types with postsynaptic neurons.
        region : str
            The region by which to limit the map.
        min_synapses : int, optional
            The minimum number of synapses in min_maps. The default is 3.
        exclude : list-like, optional
            The neurons to exclude from weight maps. The default is [].
        font_size : int, optional
            The font size of the generated figures. The default is 6.ß
        """
        self.pre_types = pre_types
        self.post_types = post_types
        self.region = region
        
        self.pre_neurs = ids_from_types(self.pre_types)
        self.post_neurs = ids_from_types(self.post_types)
        self.syn_df = partner_df(ids=self.post_neurs, region=self.region, \
                                 pre=False)
            
        self.syn_map, self.weight_map, self.weight_map_min_3=\
            make_syn_maps(self.syn_df, self.pre_neurs, self.post_neurs,\
                         min_synapses=min_synapses, exclude=exclude)
        self.type_map, self.type_weight_map = make_type_maps(self.syn_df,
                                            self.pre_types, self.post_types)
        
        self.font_size = font_size
        self.label_pad = 1
        self.tick_pad = 0
    
    def get_major_label(self, type_list, starting_point=-0.5):
        """Gets the major label for connection maps.

        Parameters
        ----------
        type_list : list-like
            List of types.
        starting_point : float, optional
            The starting location of the label. The default is -0.5.

        Returns
        -------
        label_position : np.array
            The positions of each label.
        """
        label_position = np.array([], dtype=float)
        for i in type_list:
            label_position = np.append(label_position, \
                                       (len(neur_ids[i])/2) + starting_point)
            starting_point += len(neur_ids[i])
        return label_position
    
    def get_minor_label(self, type_list, starting_point=-0.5):
        """Gets the minor label for connection maps.

        Parameters
        ----------
        type_list : list-like
            List of types.
        starting_point : float, optional
            The starting location of the label. The default is -0.5.

        Returns
        -------
        label_position : np.array
            The positions of each label.
        """
        label_position = np.array([], dtype=float)
        for indi, i in enumerate(type_list):
            starting_point += len(neur_ids[i])
            if indi < len(type_list)-1:
                label_position = np.append(label_position, starting_point)
        return label_position
    
    def get_type_label(self, type_list):
        """Gets the major label for type maps.

        Parameters
        ----------
        type_list : list-like
            List of types.
        starting_point : float, optional
            The starting location of the label. The default is -0.5.

        Returns
        -------
        label_position : np.array
            The positions of each label.
        """
        return np.arange(len(type_list))
    
    def get_minor_type_label(self, type_list, starting_point=-0.5):
        """Gets the minor label for type maps.

        Parameters
        ----------
        type_list : list-like
            List of types.
        starting_point : float, optional
            The starting location of the label. The default is -0.5.

        Returns
        -------
        label_position : np.array
            The positions of each label.
        """
        label_position = np.array([], dtype=float)
        order = []
        type_amounts = collections.defaultdict(int)
        for i in type_list:
            if i[0:2]=="ER":
                current_type = "ER"
            elif i[0:4]=="AOTU":
                current_type = "AOTU046"
            elif i[0:4] in ["MeTu", "TuBu", "TuTu"]:
                current_type = i[0:4]
            else:
                continue
            if not current_type in order:
                order.append(current_type)
            type_amounts[current_type] += 1
        for indi, i in enumerate(order):
            if indi==len(order)-1:
                break
            starting_point += type_amounts[i]
            label_position = np.append(label_position, starting_point)
        return label_position
                
    def plot_connectivity(self, conn_map, plot_name="", cmap_color="Purples",
                          fig_size=(1.5,1.5), save_figure=True):
        """Makes a connectivity matrix that can be exported.

        Parameters
        ----------
        conn_map : np.array
            The map to plot.
        plot_name : str, optional
            The plot name to save. The default is "".
        cmap_color : matplotlib.colors.ListerColormap, optional
            The color scheme of the matrix plot. The default is "Purples".
        fig_size : tuple, optional
            The size of the figure.. The default is (1.5,1.5).
        save_figure : TYPE, optional
            Whether to save the figure file. The default is True.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            The connectivity matrix.
        """
        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot(1,1,1)
        img = ax.imshow(conn_map, cmap=plt.get_cmap(cmap_color),
                        interpolation="none")
        cax = fig.add_axes([ax.get_position().x1+0.01,
                            ax.get_position().y0,0.02,
                            ax.get_position().height])
        cbar = plt.colorbar(img, cax=cax)       
    
        cbar.ax.tick_params(labelsize=self.font_size)
        ax.tick_params(axis='both', which='major', pad=self.tick_pad)
        
        ax.set_xlabel("Post-Synaptic Neurons", fontsize=self.font_size,
                      labelpad=self.label_pad)
        ax.set_xticks(self.get_major_label(self.post_types))
        ax.set_xticklabels(self.post_types, rotation=90, fontsize=self.font_size)
        
        ax.set_ylabel("Pre-Synaptic Neurons", fontsize=self.font_size,
                      labelpad=self.label_pad)
        ax.set_yticks(self.get_major_label(self.pre_types))
        ax.set_yticklabels(self.pre_types, fontsize = self.font_size)
            
        ax.set_xticks(self.get_minor_label(self.post_types), minor=True)
        ax.set_yticks(self.get_minor_label(self.pre_types), minor=True)
        ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.2)
        ax.tick_params(which="both", bottom=False, left=False)
        
        if save_figure:
            figures.save_fig(fig, plot_name=plot_name)
        return fig
    
    def plot_type_connectivity(self, type_map, plot_name="", 
                               cmap_color=new_colors["Green"],
                               fig_size=(1.5,1.5), save_figure=True):
        """Makes a neuron type connectivity matrix that can be exported.

        Parameters
        ----------
        type_map : np.array
            The map to plot.
        plot_name : str, optional
            The plot name to save. The default is "".
        cmap_color : matplotlib.colors.ListerColormap, optional
            The color scheme of the matrix plot. The default is 
            new_colors["Green"].
        fig_size : tuple, optional
            The size of the figure. The default is (1.5,1.5).
        save_figure : TYPE, optional
            Whether to save the figure file. The default is True.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            The connectivity matrix.
        """
        fig = plt.figure(figsize = fig_size)
        ax = fig.add_subplot(1,1,1)
        img = ax.imshow(type_map, cmap = plt.get_cmap(cmap_color), 
                         interpolation='none')
        cax = fig.add_axes([ax.get_position().x1+0.01,
                            ax.get_position().y0,0.02,
                            ax.get_position().height])
        cbar = plt.colorbar(img, cax=cax)       

        cbar.ax.tick_params(labelsize=self.font_size)
        
        ax.tick_params(axis='both', which='major', pad=self.tick_pad)
    
        ax.set_xlabel("Post-Synaptic Neurons", fontsize=self.font_size, 
                      labelpad=self.label_pad)
        ax.set_xticks(self.get_type_label(self.post_types))
        ax.set_xticklabels(self.post_types, rotation=90, fontsize=self.font_size)
        
        ax.set_ylabel("Pre-Synaptic Neurons", fontsize=self.font_size,
                      labelpad=self.label_pad)
        ax.set_yticks(self.get_type_label(self.pre_types))
        ax.set_yticklabels(self.pre_types, fontsize = self.font_size)
        
        ax.set_xticks(self.get_minor_type_label(self.post_types), minor=True)
        ax.set_yticks(self.get_minor_type_label(self.pre_types), minor=True)
        ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.2)
        ax.tick_params(which="both", bottom=False, left=False)
        
        if save_figure:
            figures.save_fig(fig, plot_name=plot_name)
        return fig
    
    def make_connectivity_plots(self, plot_name="", fig_size=(1.5,1.5)):
        """Plots the connectivity matrices.

        Parameters
        ----------
        plot_name : str, optional
            The plot name to save. The default is "".
        fig_size : tuple, optional
            The size of the figure. The default is (1.5,1.5).
        """
        self.plot_connectivity(self.syn_map, \
                        plot_name=f"{plot_name} Synaptic Connections",
                        cmap_color=new_colors["Orange"],
                        fig_size=fig_size, save_figure=False)
        self.plot_connectivity(self.weight_map, \
                        plot_name=f"{plot_name} Synaptic Weights",
                        fig_size=fig_size, save_figure=False)
        self.plot_connectivity(self.weight_map_min_3, \
                        plot_name=f"{plot_name} Synaptic Weights Min 3",
                        fig_size=fig_size, save_figure=True)
    
    def make_type_plots(self, plot_name="", fig_size=(1.5,1.5)): 
        """Plots the neuron type connectivity matrices.

        Parameters
        ----------
        plot_name : str, optional
            The plot name to save. The default is "".
        fig_size : tuple, optional
            The size of the figure. The default is (1.5,1.5).
        """
        self.plot_type_connectivity(self.type_map, \
                        plot_name=f"{plot_name} Type Connections",
                        fig_size=fig_size, save_figure=False)
        self.plot_type_connectivity(self.type_weight_map, \
                        plot_name=f"{plot_name} Type Weights",
                        fig_size=fig_size)
    
    def to_excel(self, file_name):
        """Exports all maps to Excel.
        
        Parameters
        ----------
        file_name : str
            The name of the file to be exported.
        """
        pre_header = syn_header(self.pre_types)
        pre_header = np.insert(pre_header, 0, np.empty((2,2), dtype="U50"),
                               axis=1)
        post_header = syn_header(self.post_types)

        pre_type_header = np.array([self.pre_types], dtype="U50")
        pre_type_header = np.insert(pre_type_header, 0, np.empty((1,1), dtype="U50"),
                                axis=1)
        post_type_header = np.array([self.post_types], dtype="U50")
        
        absolute_path = os.path.dirname(__file__)
        relative_path = os.path.join("Excel Plots", f"{file_name}.xlsx")
        file_path = os.path.join(absolute_path, relative_path)
        sheet_dict = {"Connections": self.syn_map, 
                      "Weights": self.weight_map, 
                      "Weights Minimized": self.weight_map_min_3, 
                      "Types": self.type_map, \
                       "Type Weights": self.type_weight_map}
        with pd.ExcelWriter(file_path) as writer:
            for i in sheet_dict:
                current_arr = np.array(sheet_dict[i], dtype="U50")
                if i in ["Connections", "Weights", "Weights Minimized"]:
                    current_arr = np.concatenate((post_header, current_arr))
                    current_arr = np.insert(current_arr, 0, pre_header, axis=1)
                else:
                    current_arr = np.concatenate((post_type_header, current_arr))
                    current_arr = np.insert(current_arr, 0, pre_type_header, axis=1)
                current_df = pd.DataFrame(current_arr)
                current_df.to_excel(writer, sheet_name = i)

        
def syn_header(types):
    """Creates a np.array with the types and their ids.
    
    Parameters
    ----------
    types: list-like
        A list of types in neur_coords
        
    Returns
    -------
    header: np.array 
        The second row is all ids and the first row is the id's respective type.
    """
    all_ids = np.array(ids_from_types(types), dtype="U50")
    
    types_sequential = np.array([], dtype="U50")
    for i in types:
        neur_count = len(neur_ids[i])
        for j in range(neur_count):
            types_sequential = np.append(types_sequential, i)
    header = np.array([types_sequential, all_ids], dtype="U50")
    return header



=======
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 16:24:03 2022

@author: Dustin Garner
"""


import os
import collections
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fafbseg
from fafbseg import flywire
import flycode.reduction as reduction
import flycode.utils as utils
import flycode.readfiles as readfiles
import flycode.figures as figures
import flycode.flywire_functions as fw


pd.options.mode.chained_assignment = None  # default='warn'

neur_coords, prev_ids = readfiles.import_coords()
neur_ids = {}
new_colors = readfiles.import_colors()


def add_types(types, output=False):
    """Adds types to the dictionary neur_ids
    
    Parameters
    ----------
    types : list-like
        List of types that are in neur_coords. They will be added to neur_ids 
        in the same order.
    output : bool, Optional
        Whether to show if neurons in Neuron Spreadsheet have been updated
        since the IDs were last added.
    """
    not_updated = np.array([], dtype=np.int64)
    not_in_spread = np.array([], dtype=str)
    for i in types:
        if i in neur_ids:
            continue
        if i not in neur_coords:
            print(f"{i} is not a valid neuron type.")
            continue
        neur_ids[i] = fw.locs_to_segments(neur_coords[i])
        if not output:
            continue
        if not i in prev_ids:
            not_in_spread = np.append(not_in_spread, i)
            continue
        for j, k in zip(neur_ids[i], prev_ids[i]):
            if not j==k:
                not_updated = np.append(not_updated, k)
    if output:
        print(f"No neurons in spreadsheet: {not_in_spread}\n")
        print(f"Not updated: {not_updated}\n")


def ids_from_types(types):
    """Returns a list of ids from the given types.
    
    Parameters
    ----------
    types: list-like 
        List of types that are in neur_coords.
    
    Returns
    -------
    ids: np.array
        Array of neur ids from types.
    """
    if isinstance(types, str):
        types = [types]
    add_types(types, output=False)
    return np.concatenate(([neur_ids.get(x, np.array([])) for x in types]))


def copy_types(neur_types, output=True):
    """Copies neurons of certain types to the clipboard.

    Parameters
    ----------
    neur_types : str, list-like
        Neuron types you want copied.
    output : bool, optional
        Whether to print the number of neurons copied. The default is True.
    """
    neur_ids = ids_from_types(neur_types)
    utils.copy_ids(neur_ids, output=output)


def find_neur_type(neur_id, search=[], unidentified_type=""):
    """Finds what neuron type the neuron ID is, assuming it is in neur_ids.

    Parameters
    ----------
    neur_id : int
        The neuron ID.
    search : list-like, optional
        The types to search through. Default is [], which means it will search
        through all keys in mapping.neur_ids.
    unidentified_type : str, optional
        The type that is returned if the ID is not in neur_ids. The default is "".

    Returns
    -------
    str
        The type or unidentified type if ID is not in neur_ids.

    """
    add_types(search)
    search = neur_ids.keys() if search==[] else search
    for i in search:
        if int(neur_id) in neur_ids[i]:
            return i
    return unidentified_type


def check_existence_in_region(neur_types, regions, min_syns=10):
    """Checks whether any neurons from the given types do not have at least 
        min_syns in the given regions.
    
    Parameters
    ----------
    neur_types : list-like
        Neuron types in Neuron Spreadsheet.
    regions : list-like
        The regions to test.
    min_syns : int, Optional
        Below this number the neuron is displayed in output.
    """
    ids = ids_from_types(neur_types)
    syn_df = fw.fetch_synapses(ids)
    for i in ids:
        temp_df = syn_df[(syn_df["pre"]==i) | (syn_df["post"]==i)]
        for j in regions:
            temp_df2 = reduction.lim_region(temp_df, j)
            number_syns = len(temp_df2)
            if number_syns<min_syns:
                print(f"{i!s} has {number_syns!s} connections in {j}")


def find_partners(neurons, region, pre_or_post, min_neurs=5):
    """Finds the pre_or_post partners of all of the given neurons. Returns a 
        dictionary of the partner with number of connection and an array of 
        just the partners.
    
    Parameters
    ----------
    neurons : list-like
        Neuron IDs.
    region : str
        The region to limit the search.
    pre_or_post : str
        Whether pre or postsynaptic connections are wanted
    min_neurs : int, Optional
        The lowest number of connections to be shown in the output.
        
    Returns
    -------
    neur_dict : dict
        A dictionary. Keys are IDs and values are the number of connections from
        the input neurons.
    """
    pre = False if pre_or_post=="pre" else True
    post = not pre
    syn_df = fw.fetch_synapses(neurons, pre=pre, post=post)
    syn_df = reduction.lim_region(syn_df, region)
    all_partners = np.unique(np.array(syn_df[pre_or_post]))
    for i in all_partners:
        temp_df = syn_df[syn_df[pre_or_post]==i]
        if len(temp_df)<min_neurs:
            syn_df = syn_df[syn_df[pre_or_post]!=i]
    new_partners = np.array(syn_df[pre_or_post])
    neur_dict = utils.count_instances(new_partners, min_neurs=min_neurs)
    return neur_dict
    

def partner_df(ids, region="Connectome", pre=True, post=True, 
           autapses=False, min_synapses=1):
    """Retrieves the partner spreadsheet from Flywire.
    
    Parameters
    ----------
    ids : int or list of int
        Neurons that you want to find the partners of.
    region : str
        The region to limit the map to.
    pre : bool
        Whether you want the neurs to be pre.
    post : bool
        Whether you want the neurs to be post.
    autapses : bool
        If False, removes autapses.
    min_synapses : int
        Limits the df to this number of neurons.

    Returns
    -------
    syn_df: pd.DataFrame
        The resulting synapse dataframe.
    """
    syn_df = fw.fetch_synapses(ids, pre=pre, post=post)
    syn_df = reduction.lim_region(syn_df, region)
    if not autapses:
        syn_df = reduction.remove_autapses(syn_df)
    syn_df = reduction.remove_min_synapses(syn_df, min_syns=min_synapses)
    return syn_df


def make_syn_maps(df, pre_neurs, post_neurs, min_synapses=3, exclude=[]):
    """Makes connectivity and weight maps given a dataframe from Flywire.

    Parameters
    ----------
    df : pandas.DataFrame
        A neuron dataframe given by flywire.get_neurons.
    pre_neurs : list-like
        The pre-synaptic neurons.
    post_neurs : list-like
        The post-synaptic neurons.
    min_synapses : int, optional, the default is 3.
        The minimum number of synapses when making the minimized weight map.
    exclude : list-like, optional
        The neurons to be excluded from being counted in the plot.

    Returns
    -------
    syn_map : np.array
        Connectivity map of number of synapses.
    weight_map : np.array
        Weight map of the percent of connections from the presynaptic neuron
        to the given postsynaptic neuron.
    min_map : np.array
        Same as weight map but only if the neurons have greater than min_synapses
        total.
    """
    pre_count = len(pre_neurs)
    post_count = len(post_neurs)
    syn_map = np.zeros((pre_count, post_count), dtype=float)
    weight_map = np.zeros((pre_count, post_count), dtype=float)
    min_map = np.zeros((pre_count, post_count), dtype=float)
    
    for indi, i in enumerate(post_neurs):
        temp_df = df[df["post"]==i]
        total_connections = len(temp_df)
        if total_connections==0: #Prevents from dividing by zero
            print(f"{i!s} does not have connections.")
            continue
        for indj, j in enumerate(pre_neurs):
            partner_count = 0 if (i in exclude or j in exclude) else \
                    len(temp_df[temp_df["pre"]==j])
            synaptic_weight = partner_count/total_connections
            min_weight = synaptic_weight if total_connections>=min_synapses else 0
            syn_map[indj][indi] = partner_count
            weight_map[indj][indi] = synaptic_weight
            min_map[indj][indi] = min_weight
    return syn_map, weight_map, min_map


def make_type_maps(df, pre_types, post_types):
    """Makes type connectivity and weight maps given a dataframe from Flywire.

    Parameters
    ----------
    df : pandas.DataFrame
        A neuron dataframe given by flywire.get_neurons.
    pre_types : list-like
        The pre-synaptic types.
    post_types : list-like
        The post-synaptic types.

    Returns
    -------
    type : np.array
        Connectivity map of total number of synapses for the given types.
    weight_map : np.array
        Weight map of the percent of connections from the presynaptic type
        to the given postsynaptic type.
    """
    add_types(pre_types)
    add_types(post_types)
    type_map = np.zeros((len(pre_types), len(post_types)), dtype=np.uint64)
    weight_map = np.zeros((len(pre_types), len(post_types)), dtype=float)
    for indi, i in enumerate(post_types):
        temp_df = df[df["post"].isin(neur_ids[i])]
        total_connections = len(temp_df)
        for indj, j in enumerate(pre_types):
            type_connections = len(temp_df[temp_df["pre"].isin(neur_ids[j])])
            type_map[indj][indi] = type_connections
            weight_map[indj][indi] = type_connections/total_connections \
                    if total_connections!=0 else 0
    return type_map, weight_map


class ConnectionMap: 
    #Creates maps between pre types and post types in a region.
    def __init__(self, pre_types, post_types, region, min_synapses=5, \
                 exclude=[], font_size=6):
        """Creates maps between pre types and post types in a region.

        Parameters
        ----------
        pre_types : list-like
            The types with presynaptic neurons.
        post_types : list-like
            The types with postsynaptic neurons.
        region : str
            The region by which to limit the map.
        min_synapses : int, optional
            The minimum number of synapses in min_maps. The default is 3.
        exclude : list-like, optional
            The neurons to exclude from weight maps. The default is [].
        font_size : int, optional
            The font size of the generated figures. The default is 6.ß
        """
        self.pre_types = pre_types
        self.post_types = post_types
        self.region = region
        
        self.pre_neurs = ids_from_types(self.pre_types)
        self.post_neurs = ids_from_types(self.post_types)
        self.syn_df = partner_df(ids=self.post_neurs, region=self.region, \
                                 pre=False)
            
        self.syn_map, self.weight_map, self.weight_map_min_3=\
            make_syn_maps(self.syn_df, self.pre_neurs, self.post_neurs,\
                         min_synapses=min_synapses, exclude=exclude)
        self.type_map, self.type_weight_map = make_type_maps(self.syn_df,
                                            self.pre_types, self.post_types)
        
        self.font_size = font_size
        self.label_pad = 1
        self.tick_pad = 0
    
    def get_major_label(self, type_list, starting_point=-0.5):
        """Gets the major label for connection maps.

        Parameters
        ----------
        type_list : list-like
            List of types.
        starting_point : float, optional
            The starting location of the label. The default is -0.5.

        Returns
        -------
        label_position : np.array
            The positions of each label.
        """
        label_position = np.array([], dtype=float)
        for i in type_list:
            label_position = np.append(label_position, \
                                       (len(neur_ids[i])/2) + starting_point)
            starting_point += len(neur_ids[i])
        return label_position
    
    def get_minor_label(self, type_list, starting_point=-0.5):
        """Gets the minor label for connection maps.

        Parameters
        ----------
        type_list : list-like
            List of types.
        starting_point : float, optional
            The starting location of the label. The default is -0.5.

        Returns
        -------
        label_position : np.array
            The positions of each label.
        """
        label_position = np.array([], dtype=float)
        for indi, i in enumerate(type_list):
            starting_point += len(neur_ids[i])
            if indi < len(type_list)-1:
                label_position = np.append(label_position, starting_point)
        return label_position
    
    def get_type_label(self, type_list):
        """Gets the major label for type maps.

        Parameters
        ----------
        type_list : list-like
            List of types.
        starting_point : float, optional
            The starting location of the label. The default is -0.5.

        Returns
        -------
        label_position : np.array
            The positions of each label.
        """
        return np.arange(len(type_list))
    
    def get_minor_type_label(self, type_list, starting_point=-0.5):
        """Gets the minor label for type maps.

        Parameters
        ----------
        type_list : list-like
            List of types.
        starting_point : float, optional
            The starting location of the label. The default is -0.5.

        Returns
        -------
        label_position : np.array
            The positions of each label.
        """
        label_position = np.array([], dtype=float)
        order = []
        type_amounts = collections.defaultdict(int)
        for i in type_list:
            if i[0:2]=="ER":
                current_type = "ER"
            elif i[0:4]=="AOTU":
                current_type = "AOTU046"
            elif i[0:4] in ["MeTu", "TuBu", "TuTu"]:
                current_type = i[0:4]
            else:
                continue
            if not current_type in order:
                order.append(current_type)
            type_amounts[current_type] += 1
        for indi, i in enumerate(order):
            if indi==len(order)-1:
                break
            starting_point += type_amounts[i]
            label_position = np.append(label_position, starting_point)
        return label_position
                
    def plot_connectivity(self, conn_map, plot_name="", cmap_color="Purples",
                          fig_size=(1.5,1.5), save_figure=True):
        """Makes a connectivity matrix that can be exported.

        Parameters
        ----------
        conn_map : np.array
            The map to plot.
        plot_name : str, optional
            The plot name to save. The default is "".
        cmap_color : matplotlib.colors.ListerColormap, optional
            The color scheme of the matrix plot. The default is "Purples".
        fig_size : tuple, optional
            The size of the figure.. The default is (1.5,1.5).
        save_figure : TYPE, optional
            Whether to save the figure file. The default is True.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            The connectivity matrix.
        """
        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot(1,1,1)
        img = ax.imshow(conn_map, cmap=plt.get_cmap(cmap_color),
                        interpolation="none")
        cax = fig.add_axes([ax.get_position().x1+0.01,
                            ax.get_position().y0,0.02,
                            ax.get_position().height])
        cbar = plt.colorbar(img, cax=cax)       
    
        cbar.ax.tick_params(labelsize=self.font_size)
        ax.tick_params(axis='both', which='major', pad=self.tick_pad)
        
        ax.set_xlabel("Post-Synaptic Neurons", fontsize=self.font_size,
                      labelpad=self.label_pad)
        ax.set_xticks(self.get_major_label(self.post_types))
        ax.set_xticklabels(self.post_types, rotation=90, fontsize=self.font_size)
        
        ax.set_ylabel("Pre-Synaptic Neurons", fontsize=self.font_size,
                      labelpad=self.label_pad)
        ax.set_yticks(self.get_major_label(self.pre_types))
        ax.set_yticklabels(self.pre_types, fontsize = self.font_size)
            
        ax.set_xticks(self.get_minor_label(self.post_types), minor=True)
        ax.set_yticks(self.get_minor_label(self.pre_types), minor=True)
        ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.2)
        ax.tick_params(which="both", bottom=False, left=False)
        
        if save_figure:
            figures.save_fig(fig, plot_name=plot_name)
        return fig
    
    def plot_type_connectivity(self, type_map, plot_name="", 
                               cmap_color=new_colors["Green"],
                               fig_size=(1.5,1.5), save_figure=True):
        """Makes a neuron type connectivity matrix that can be exported.

        Parameters
        ----------
        type_map : np.array
            The map to plot.
        plot_name : str, optional
            The plot name to save. The default is "".
        cmap_color : matplotlib.colors.ListerColormap, optional
            The color scheme of the matrix plot. The default is 
            new_colors["Green"].
        fig_size : tuple, optional
            The size of the figure. The default is (1.5,1.5).
        save_figure : TYPE, optional
            Whether to save the figure file. The default is True.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            The connectivity matrix.
        """
        fig = plt.figure(figsize = fig_size)
        ax = fig.add_subplot(1,1,1)
        img = ax.imshow(type_map, cmap = plt.get_cmap(cmap_color), 
                         interpolation='none')
        cax = fig.add_axes([ax.get_position().x1+0.01,
                            ax.get_position().y0,0.02,
                            ax.get_position().height])
        cbar = plt.colorbar(img, cax=cax)       

        cbar.ax.tick_params(labelsize=self.font_size)
        
        ax.tick_params(axis='both', which='major', pad=self.tick_pad)
    
        ax.set_xlabel("Post-Synaptic Neurons", fontsize=self.font_size, 
                      labelpad=self.label_pad)
        ax.set_xticks(self.get_type_label(self.post_types))
        ax.set_xticklabels(self.post_types, rotation=90, fontsize=self.font_size)
        
        ax.set_ylabel("Pre-Synaptic Neurons", fontsize=self.font_size,
                      labelpad=self.label_pad)
        ax.set_yticks(self.get_type_label(self.pre_types))
        ax.set_yticklabels(self.pre_types, fontsize = self.font_size)
        
        ax.set_xticks(self.get_minor_type_label(self.post_types), minor=True)
        ax.set_yticks(self.get_minor_type_label(self.pre_types), minor=True)
        ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.2)
        ax.tick_params(which="both", bottom=False, left=False)
        
        if save_figure:
            figures.save_fig(fig, plot_name=plot_name)
        return fig
    
    def make_connectivity_plots(self, plot_name="", fig_size=(1.5,1.5)):
        """Plots the connectivity matrices.

        Parameters
        ----------
        plot_name : str, optional
            The plot name to save. The default is "".
        fig_size : tuple, optional
            The size of the figure. The default is (1.5,1.5).
        """
        self.plot_connectivity(self.syn_map, \
                        plot_name=f"{plot_name} Synaptic Connections",
                        cmap_color=new_colors["Orange"],
                        fig_size=fig_size, save_figure=False)
        self.plot_connectivity(self.weight_map, \
                        plot_name=f"{plot_name} Synaptic Weights",
                        fig_size=fig_size, save_figure=False)
        self.plot_connectivity(self.weight_map_min_3, \
                        plot_name=f"{plot_name} Synaptic Weights Min 3",
                        fig_size=fig_size, save_figure=True)
    
    def make_type_plots(self, plot_name="", fig_size=(1.5,1.5)): 
        """Plots the neuron type connectivity matrices.

        Parameters
        ----------
        plot_name : str, optional
            The plot name to save. The default is "".
        fig_size : tuple, optional
            The size of the figure. The default is (1.5,1.5).
        """
        self.plot_type_connectivity(self.type_map, \
                        plot_name=f"{plot_name} Type Connections",
                        fig_size=fig_size, save_figure=False)
        self.plot_type_connectivity(self.type_weight_map, \
                        plot_name=f"{plot_name} Type Weights",
                        fig_size=fig_size)
    
    def to_excel(self, file_name):
        """Exports all maps to Excel.
        
        Parameters
        ----------
        file_name : str
            The name of the file to be exported.
        """
        pre_header = syn_header(self.pre_types)
        pre_header = np.insert(pre_header, 0, np.empty((2,2), dtype="U50"),
                               axis=1)
        post_header = syn_header(self.post_types)

        pre_type_header = np.array([self.pre_types], dtype="U50")
        pre_type_header = np.insert(pre_type_header, 0, np.empty((1,1), dtype="U50"),
                                axis=1)
        post_type_header = np.array([self.post_types], dtype="U50")
        
        absolute_path = os.path.dirname(__file__)
        relative_path = os.path.join("Excel Plots", f"{file_name}.xlsx")
        file_path = os.path.join(absolute_path, relative_path)
        sheet_dict = {"Connections": self.syn_map, 
                      "Weights": self.weight_map, 
                      "Weights Minimized": self.weight_map_min_3, 
                      "Types": self.type_map, \
                       "Type Weights": self.type_weight_map}
        with pd.ExcelWriter(file_path) as writer:
            for i in sheet_dict:
                current_arr = np.array(sheet_dict[i], dtype="U50")
                if i in ["Connections", "Weights", "Weights Minimized"]:
                    current_arr = np.concatenate((post_header, current_arr))
                    current_arr = np.insert(current_arr, 0, pre_header, axis=1)
                else:
                    current_arr = np.concatenate((post_type_header, current_arr))
                    current_arr = np.insert(current_arr, 0, pre_type_header, axis=1)
                current_df = pd.DataFrame(current_arr)
                current_df.to_excel(writer, sheet_name = i)

        
def syn_header(types):
    """Creates a np.array with the types and their ids.
    
    Parameters
    ----------
    types: list-like
        A list of types in neur_coords
        
    Returns
    -------
    header: np.array 
        The second row is all ids and the first row is the id's respective type.
    """
    all_ids = np.array(ids_from_types(types), dtype="U50")
    
    types_sequential = np.array([], dtype="U50")
    for i in types:
        neur_count = len(neur_ids[i])
        for j in range(neur_count):
            types_sequential = np.append(types_sequential, i)
    header = np.array([types_sequential, all_ids], dtype="U50")
    return header



>>>>>>> Stashed changes
