# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 10:31:15 2023

@author: dusti
"""


import math
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from fafbseg import flywire
import flycode.reduction as reduction
import flycode.readfiles as readfiles
import flycode.mapping as mapping
import flycode.utils as utils
import flycode.flywire_functions as fw


regions = readfiles.import_file("Regions", sheet_name = "Boxes")

type_colors = {
    "TuBu01_R": "TBlue",
    "TuBu02_R": "TPink",
    "TuBu03_R": "TBlue",
    "TuBu04_R": "TRed",
    "TuBu05_R": "TGreen",
    "TuBu06_R": "TPink",
    "TuBu07_R": "TOrange2",
    "TuBu08_R": "TRed",
    "TuBu09_R": "TGreen",
    "TuBu10_R": "TPink",
    "MeTu1_R": "TRed",
    "MeTu2a_R": "TBlue",
    "MeTu2b_R": "TPink",
    "MeTu3a_R": "TBlue",
    "MeTu3b_R": "TOrange2",
    "MeTu3c_R": "TGreen",
    "MeTu4a_R": "TBlue",
    "MeTu4b_R": "TPink",
    "MeTu4c_R": "TGreen",
    "MeTu4d_R": "TOrange2"
    }

type_max = {}
for i in ["MeTu1_R", "TuBu08_R"]:
    type_max[i] = 0.138351
for i in ["MeTu2a_R", "MeTu2b_R","TuBu01_R", "TuBu06_R"]:
    type_max[i] = 0.18535
for i in ["MeTu3a_R", "MeTu3b_R", "MeTu3c_R", "TuBu07_R", "TuBu09_R", "TuBu10_R"]:
    type_max[i] = 0.12935
for i in ["MeTu4a_R", "MeTu4b_R", "MeTu4c_R", "MeTu4d_R", "TuBu02_R", "TuBu03_R", 
          "TuBu04_R", "TuBu05_R"]:
    type_max[i] = 0.10949



def rotate_coord(coord: np.array, axis, angle, scale=1):
    """
    Parameters
    ----------
    coord : np.array
        [x,y,z] Coordinate you want rotated around the z-axis.
    axis : np.array
        [x1,y1,z1] The point around which you want the coord rotated.
    angle : float
        The angle of the rotation (in degrees).
    scale : float
        The new scale of the point in reference to the axis.
        
    Returns
    -------
    coord : np.array
        The rotated coord.
    """
    angle = math.radians(angle)
    coord = np.matmul(np.array([[1,0,0,-axis[0]],
                               [0,1,0,-axis[1]],
                               [0,0,1,-axis[2]]]), np.append(coord, 1))
    coord = np.matmul(np.array([[scale, 0, 0],
                                [0, scale, 0],
                                [0, 0, 1]]), coord)
    coord = np.matmul(np.array([[math.cos(angle), -math.sin(angle), 0],
                               [math.sin(angle), math.cos(angle), 0],
                               [0, 0, 1]]), coord)
    coord = np.matmul(np.array([[1,0,0,axis[0]],
                               [0,1,0,axis[1]],
                               [0,0,1,axis[2]]]), np.append(coord, 1))
    return coord
    
    
def make_box(xmin, xmax, ymin, ymax, zmin, zmax):
    """Give the coordinates of a box given min and max values.
    
    Parameters
    ----------
    xmin : float
        The minimum x value of the box.
    xmax : float
        The maximum x value of the box.
    ymin : float
        The minimum y value of the box.
    ymax : float
        The maximum y value of the box.
    zmin : float
        The minimum z value of the box.
    zmax : float
        The maximum z value of the box.
    
    Returns
    -------
    coords : np.array
        The 8 corner coordinates of the box, stored as a 2-D array.
    """
    coords = np.array([[]])
    for i,j,k in itertools.product([xmin,xmax], [ymin,ymax], [zmin,zmax]):
        temp_coord = np.array([[i,j,k]])
        if coords.size==0:
            coords = temp_coord
            continue
        coords = np.append(coords, temp_coord, axis=0)
    return coords
                

def transform_box(box, region="AOTU_R", angle=30, scale=1.3):
    """Performs transformations on the box and gives the new coordinates.

    Parameters
    ----------
    box : np.array
        A box, for instance given by make_box().
    region : str, optional
        The region, used for finding the center (if more regions are used, more
        centers should be added). The default is "AOTU_R".
    angle : float, optional
        The angle by which to rotate the box. The default is 30.
    scale : float, optional
        The scale by which to resize the box. The default is 1.3.

    Returns
    -------
    box : np.array
        The new box coordinates.
    """
    if angle == 0:
        return box
    center = np.array([0, 0])
    if region == "AOTU_R":
        center = np.array([167400, 43250])
    elif region == "Full_AOTU_R":
        center = np.array([163972, 39883])
    for i in range(len(box)):
        coord = box[i]
        box[i] = rotate_coord(coord, 
                              np.concatenate((center, np.array([coord[2]]))), 
                              angle=angle, scale=scale)
    return box


def get_steps(region="AOTU_R", spacing=10):
    """Creates a set of equidistant steps along each axis depending on the 
        region.

    Parameters
    ----------
    region : str, optional
        The region by which to find the steps. The default is "AOTU_R".
    spacing : int, optional
        The distance each step should be. The default is 10.

    Returns
    -------
    x_steps : np.array
        Equally spaced steps along the x-axis.
    y_steps : np.array
        Equally spaced steps along the y-axis.
    z_steps : np.array
        Equally spaced steps along the z-axis.
    """
    temp_regions = regions[regions["Region"]==region]
    minmax = {}
    for i in temp_regions.columns:
        if i == "Region":
            continue
        minmax[i] = np.array(temp_regions[i], dtype=float)[0] + 0.5
    
    x_steps = np.arange(minmax["xmin"], minmax["xmax"], spacing)
    y_steps = np.arange(minmax["ymin"], minmax["ymax"], spacing)
    z_steps = np.arange(minmax["zmin"], minmax["zmax"], spacing//10)
    return x_steps, y_steps, z_steps


def make_map(df, x_steps, y_steps, z_steps, spacing, plane, vert=False,
             reverse={True: False, False: False}, region="AOTU_R",
             angle=30, scale=1.3):
    """Recursively makes a synapse density map given a synapse dataframe and a 
        region.

    Parameters
    ----------
    df : pd.DataFrame
        Synapse dataframe given by flywire.get_synapses().
    x_steps : np.array
        Given by get_steps().
    y_steps : np.array
        Given by get_steps().
    z_steps : np.array
        Given by get_steps().
    spacing : int, optional
        The distance each step should be. The default is 10.
    plane : str
        "xy", "zy", or "xz" depending on which plane you would like to view from.
    vert : bool, optional
        Used for recursion, should be false upon running the function. 
        The default is False.
    reverse : dict, optional
        Should be {True: False, False: True} if plane = "xz" due to the 
        mirroring of the dataset. 
        The default is {True: False, False: False}.
    region : str, optional
        The region by which the search is performed. The default is "AOTU_R".
    angle : float, optional
        The angle by which to rotate the box. The default is 30.
    scale : float, optional
        The scale by which to resize the box. The default is 1.3.

    Returns
    -------
    np.array
        A synapse density matrix given the step size and view plane.

    """
    directions = {"x": x_steps, "y": y_steps, "z": z_steps}
    
    temp_box = make_box(x_steps[0], x_steps[-1] + spacing, 
                        y_steps[0], y_steps[-1] + spacing,
                        z_steps[0], z_steps[-1] + spacing//10)
    temp_box = transform_box(temp_box, region=region, angle=angle, scale=scale)
    temp_df = reduction.lim_vol(df, include=temp_box)
    if len(temp_df) == 0:
        return np.zeros((len(directions[plane[1]]), len(directions[plane[0]])))
    
    if vert and len(directions[plane[1]])==1:
        return np.array([[len(temp_df)]])

    axis = int(not vert)
    direct = plane[int(vert)]
    if not vert and len(directions[plane[0]])==1:
        mid = len(directions[plane[1]])//2
        direct = plane[1]
        axis = 0
    else:
        mid = len(directions[direct])//2
    
    temp_direct1, temp_direct2 = {}, {}
    for i in "xyz":
        if i==direct and not reverse[bool(axis)]:
            temp_direct1[i] = directions[i][:mid]
            temp_direct2[i] = directions[i][mid:]
        elif i == direct:
            temp_direct1[i] = directions[i][mid:]
            temp_direct2[i] = directions[i][:mid]
        else:
            temp_direct1[i] = directions[i].copy()
            temp_direct2[i] = directions[i].copy()
    
    return np.concatenate((\
            make_map(df=temp_df, x_steps=temp_direct1["x"], 
                     y_steps=temp_direct1["y"], z_steps=temp_direct1["z"], 
                     spacing=spacing, plane=plane, vert=not bool(axis), 
                     reverse=reverse, region=region, angle=angle,
                     scale=scale),
            make_map(df=temp_df, x_steps=temp_direct2["x"], 
                     y_steps=temp_direct2["y"], z_steps=temp_direct2["z"], 
                     spacing=spacing, plane=plane, vert=not bool(axis), 
                     reverse=reverse, region=region, angle=angle,
                     scale=scale)
            ), axis=axis)


@utils.time_elapsed
def xy_map(df, steps, spacing, region, angle, scale):
    """
    Uses make_map() from the xy plane point of view.
    """
    return make_map(df, steps[0], steps[1], steps[2], spacing, plane="xy",
                    region=region, angle=angle, scale=scale)

@utils.time_elapsed
def zy_map(df, steps, spacing, region, angle, scale):
    """
    Uses make_map() from the zy plane point of view.
    """
    return make_map(df, steps[0], steps[1], steps[2], spacing, plane="zy",
                    region=region, angle=angle, scale=scale)

@utils.time_elapsed
def xz_map(df, steps, spacing, region, angle, scale):
    """
    Uses make_map() from the xz plane point of view.
    """
    return make_map(df, steps[0], steps[1], steps[2], spacing, plane="xz",
            reverse={True: False, False: True}, 
            region=region, angle=angle, scale=scale)


def all_maps(df, spacing=10, region="AOTU_R", angle=30, scale=1.3):
    """Creates maps by make_map along each axis.
    
    Parameters
    ----------
    df : pd.DataFrame
        Synapse df from flywire.get_synapses()
    spacing : int, optional
        The distance each step should be. The default is 10.
    region : str, optional
        The region by which to find the steps. The default is "AOTU_R".
    angle : float, optional
        The angle by which to rotate the box. The default is 30.
    scale : float, optional
        The scale by which to resize the box. The default is 1.3.
    
    Returns
    -------
    tuple
        (xy, zy, xz) maps from their respective functions.
    """
    steps = get_steps(region = region, spacing = spacing)
    xy,zy,xz = [x(df, steps, spacing, region=region, angle=angle, scale=scale) 
                  for x in (xy_map,zy_map,xz_map)]
    print(" ")
    return (xy, zy, xz)


def maps_from_types(neur_types, spacing=10, region="AOTU_R", angle=30, scale=1.3):
    """Makes maps from certain neuron types.

    Parameters
    ----------
    neur_types : list
        The neuron types by which to find the maps.
    spacing : int, optional
        The distance each step should be. The default is 10.
    region : str, optional
        The region by which to find the steps. The default is "AOTU_R".
    angle : float, optional
        The angle by which to rotate the box. The default is 30.
    scale : float, optional
        The scale by which to resize the box. The default is 1.3.

    Returns
    -------
    tuple
        all_maps() run on the neuron types.

    """
    mapping.add_types(neur_types)
    all_ids = mapping.ids_from_types(neur_types)
    df = fw.fetch_synapses(all_ids)
    df = reduction.remove_autapses(df)
    return all_maps(df, spacing = spacing, region = region, angle = angle, 
                    scale = scale)


def maps_from_mt_tb(metu, tubu):
    """Gets maps that limit the search to just synapses between MeTu and TuBu
        neurons.

    Parameters
    ----------
    metu : list
        The MeTu types to limit the search.
    tubu : list
        The TuBu types to limit the search.

    Returns
    -------
    maps : tuple
        all_maps() run on each neuron type, but limited to the synapses between
        MeTu and TuBu neurons.
    """
    mapping.add_types(metu+tubu, output=False)
    maps = {}
    all_ids = mapping.ids_from_types(metu+tubu)
    metu_ids = mapping.ids_from_types(metu)
    tubu_ids = mapping.ids_from_types(tubu)
    df = fw.fetch_synapses(all_ids, region="AOTU_R")
    for i in metu:
        temp_df = df[df["pre"].isin(mapping.neur_ids[i])]
        temp_df = temp_df[temp_df["post"].isin(tubu_ids)]
        maps[i] = all_maps(temp_df)
    for i in tubu:
        temp_df = df[df["post"].isin(mapping.neur_ids[i])]
        temp_df = temp_df[temp_df["pre"].isin(metu_ids)]
        maps[i] = all_maps(temp_df)
    return maps
    

def plot_maps(maps, color="TPurple", blur=4, max_value=0,
              plot_name="", save_figure=True):
    """Plots the synapse density maps.
    
    Parameters
    ----------
    maps : tuple
        Maps given by all_maps()
    color : str
        A color within the Colors spreadsheet.
    blur : float
        The sigma value of the Gaussian blur on the plot.
    max_value : float
        The maximum value to be put in the corner of the plot for normalization.
    plot_name : str, optional
        The name of the plot. The default is ""
    save_figure : bool, optional
        Whether to save the figure. The default is True.
    """
    colors = readfiles.import_colors()
    color_map = colors[color]
    for i, j in zip(maps, ["xy","zy","xz"]):
        i = gaussian_filter(i, sigma=blur)
        i[0][0] = max_value #For normalizing
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        img = ax.imshow(i, cmap = color_map,
                        interpolation = "none")
        cax = fig.add_axes([ax.get_position().x1+0.01, \
                            ax.get_position().y0,0.02,ax.get_position().height])
        cbar = plt.colorbar(img, cax=cax)
        cbar.ax.tick_params(labelsize = 6)
        temp_plot_name = f"Synapse Density/{plot_name} {j}"
        if save_figure:
            utils.save_fig(fig, plot_name=temp_plot_name, transparent=True)












