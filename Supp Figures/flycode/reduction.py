<<<<<<< Updated upstream
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 16:09:31 2022

@author: dusti
"""


import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import flycode.readfiles as readfiles


region_dict = readfiles.import_regions()


def lim_region(df, region):
    """Limits a synapse dataframe to a given region (via lim_vol()).
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses(ID).
    region : str
        A region in region_dict.
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe limited to the given region.
    """
    if region in region_dict:
        return lim_vol(df, include=region_dict[region]["include"], \
                       exclude=region_dict[region]["exclude"])
    else:
        print("Region not valid")
        return df
    

def lim_vol(df, include=np.array([]), exclude=np.array([]), 
            x="post_x", y="post_y", z="post_z"):
    """Takes in a synapse dataframe and returns one with synapses limited to a 
       specified volume.
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses(ID)
    include : np.array
        The volume you would like to limit the synapses to (nm). A two-dimensional 
        array that contains coordinates that are on the outer bounds of the volume. 
        For instance if you would like to limit the synapses to be contained within
        a specific tetrahedron, this parameters would be 
        np.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3],[x4,y4,z4]])
    exclude : np.array
        This is a volume you would like the dataframe to exclude (nm). 
    x,y,z : str
        These are the respective columns of the dataframe that the coordinate 
        components are in.
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe limited to synapses within include, exempting those 
        within exclude.
    """
    oob_simplex = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0]])
    
    if include.size != 0:
        include = Delaunay(points = include)
    else:
        include = Delaunay(point = oob_simplex)
    if exclude.size != 0:
        exclude = Delaunay(points = exclude)
    else:
        exclude = Delaunay(points = oob_simplex)
    
    df_limiter = pd.Series(np.zeros((len(df)), dtype = bool), index = df.index)

    for idx, x_c, y_c, z_c in zip(df_limiter.index, df[x], df[y], df[z]):
        temp_coord = np.array([int(x_c/4), int(y_c/4), int(z_c/40)])
                
        if include.find_simplex(temp_coord)!=-1 and \
                exclude.find_simplex(temp_coord)==-1:
            df_limiter[idx] = True

    df_limed = df[df_limiter]

    return df_limed


def remove_autapses(df):
    """Takes in a synapse df and removes all autapses.
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses().
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe without autapses.
    """
    return df[df["pre"]!=df["post"]]
    

def remove_min_synapses(df, minimum_syns=1, pre_or_post="post"):
    """Takes in a synapse df and removes neurons that have fewer than 
        [minimum_syns] connections to the [pre_or_post]synaptic neuron.
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses(ID).
    minimum_syns : int
        The minimum number of synapses that a connected neuron should have to 
        be part of the df.
    pre_or_post : str
        "pre" or "post". Whether ID is presynaptic or postsynaptic.
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe with only partners that have greater than 
        minimum_syns connections.
    """
    opposite = "pre" if pre_or_post=="post" else "post"
    partners = np.unique(np.array(df[opposite]))
    lim_df = df.copy()
    for i in partners:
        if len(df[df[opposite]==i])<minimum_syns:
            lim_df = lim_df[lim_df[opposite]!=i]
    return lim_df


def coords_to_region(coords, region):
    """Limits an array of coords to a given region.
    
    Parameters
    ----------
    coords: list-like
        Coordinates to limit.
    region: str
        The region to limit the box to.
        
    Returns
    -------
    lim_coords: np.array
        The same coordinates if they fall within the region.
    """
    oob_simplex = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0]])
    include = region_dict[region]["include"]
    exclude = region_dict[region]["exclude"]
    
    if include.size != 0:
        include = Delaunay(points = include)
    else:
        include = Delaunay(point = oob_simplex)
    if exclude.size != 0:
        exclude = Delaunay(points = exclude)
    else:
        exclude = Delaunay(points = oob_simplex)
    
    lim_coords = np.array([], dtype="U50")
    for i in coords:        
        temp_coord = [int(x.strip()) for x in i.split(",")]
        if include.find_simplex(temp_coord)!=-1 and \
                exclude.find_simplex(temp_coord)==-1:
            lim_coords = np.append(lim_coords, i)
    return lim_coords
            
        
def lines_to_region(coords1, coords2, region):
    """Takes in two coord sets and returns the ones that are both within the
    given region.

    Parameters
    ----------
    coords1 : list-like
        One set of neurons.
    coords2 : list-like
        The other set of neurons
    region : str
        Region to search through.

    Returns
    -------
    lim_coords1 : np.array
        The first set of coords limited.
    lim_coords2 : np.array
        The second set of coords limited.
    """
    lim_coords1 = np.array([], dtype="U50")
    lim_coords2 = np.array([], dtype="U50")
    for i, j in zip(coords1, coords2):
        temp_coord1 = coords_to_region(np.array([i], dtype="U50"), region)
        temp_coord2 = coords_to_region(np.array([j], dtype="U50"), region)
        if temp_coord1.size>0 and temp_coord2.size>0:
            lim_coords1 = np.append(lim_coords1, temp_coord1[0])
            lim_coords2 = np.append(lim_coords2, temp_coord2[0])
    return lim_coords1, lim_coords2
            
        
    
    




=======
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 16:09:31 2022

@author: Dustin Garner
"""


import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import flycode.readfiles as readfiles


region_dict = readfiles.import_regions()


def lim_region(df, region, coord_names=[f"post_{x}" for x in "xyz"]):
    """Limits a synapse dataframe to a given region (via lim_vol()).
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses(ID).
    region : str
        A region in region_dict.
    coord_names : list-like
        The names of the x, y, and z coordinate columns in the dataframe to be
        passed into lim_vol().
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe limited to the given region.
    """
    if not region in region_dict:
        print("Region not valid")
        return df
    return lim_vol(df, 
                   include=region_dict[region]["include"],
                   exclude=region_dict[region]["exclude"],
                   coord_names=coord_names)
        
    
def remove_region(df, region, coord_names=[f"post_{x}" for x in "xyz"]):
    """Removes a certain region from a synapse dataframe.
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses(ID).
    region : str
        A region in region_dict.
    coord_names : list-like
        The names of the x, y, and z coordinate columns in the dataframe to be
        passed into lim_vol().
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe with the given region removed.
    """
    if not region in region_dict:
        print("Region not valid")
        return df
    return lim_vol(df, 
                   include=region_dict["Connectome"]["include"],
                   exclude=region_dict[region]["include"],
                   coord_names=coord_names)
    

def lim_vol(df, include=np.array([]), exclude=np.array([]), 
            coord_names = [f"post_{x}" for x in "xyz"]):
    """Takes in a synapse dataframe and returns one with synapses limited to a 
       specified volume.
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses(ID)
    include : np.array
        The volume you would like to limit the synapses to (nm). A two-dimensional 
        array that contains coordinates that are on the outer bounds of the volume. 
        For instance if you would like to limit the synapses to be contained within
        a specific tetrahedron, this parameters would be 
        np.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3],[x4,y4,z4]])
    exclude : np.array
        This is a volume you would like the dataframe to exclude (nm). 
    coord_names : list-like
        The names of the x, y, and z coordinate columns in the dataframe.
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe limited to synapses within include, exempting those 
        within exclude.
    """
    if len(coord_names)!=3:
        raise ValueError("coord_names must be of length 3.")
        
    out_of_bounds_simplex = Delaunay(\
                        points=np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0]]))
    
    include = Delaunay(points=include) if include.size!=0 \
            else out_of_bounds_simplex
    exclude = Delaunay(points=exclude) if exclude.size!=0 \
            else out_of_bounds_simplex
    df_limiter = pd.Series(np.zeros((len(df)), dtype=bool), index=df.index)
    
    x, y, z = coord_names
    for idx, x_c, y_c, z_c in zip(df_limiter.index, df[x], df[y], df[z]):
        temp_coord = np.array([x_c//4, y_c//4, z_c//40])
        if include.find_simplex(temp_coord)!=-1 and \
                exclude.find_simplex(temp_coord)==-1:
            df_limiter[idx] = True

    df_limed = df[df_limiter]
    return df_limed


def remove_autapses(df):
    """Takes in a synapse df and removes all autapses.
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses().
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe without autapses.
    """
    return df[df["pre"]!=df["post"]]
    

def remove_min_synapses(df, min_syns=1, pre_or_post="post"):
    """Takes in a synapse df and removes neurons that have fewer than 
        [min_syns] connections to the partner neurons.
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses(ID).
    min_syns : int
        The minimum number of synapses that a connected neuron should have to 
        be part of the df.
    pre_or_post : str
        "pre" or "post". Whether connected IDs are presynaptic or postsynaptic.
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe with only partners that have greater than 
        minimum_syns connections.
    """
    opposite = "pre" if pre_or_post=="post" else "post"
    partners = np.unique(np.array(df[opposite]))
    relevant_partners = []
    for i in partners:
        temp_df = df[df[opposite]==i]
        if len(temp_df)<min_syns:
            continue
        relevant_partners.append(i)
    lim_df = df[df[opposite].isin(relevant_partners)]
    return lim_df


def remove_min_weight(df, min_weight=0.05, pre_or_post="post"):
    """Takes in a synapse df and removes neurons that have fewer than 
        [min_weight] synaptic weight to the partner neurons.
    
    Parameters
    ----------
    df : pd.DataFrame
        A synapse dataframe given by flywire.get_synapses(ID).
    min_weight : float
        The minimum weight of synapses that a connected neuron should have to 
        be part of the df.
    pre_or_post : str
        "pre" or "post". Whether ID is presynaptic or postsynaptic.
    
    Returns
    -------
    df : pd.DataFrame
        The same dataframe with only partners that have greater than 
        minimum_syns connections.
    """
    opposite = "pre" if pre_or_post=="post" else "post"
    partners = np.unique(np.array(df[opposite]))
    relevant_partners = []
    for i in partners:
        temp_df = df[df[opposite]==i]
        if len(temp_df)/len(df)<min_weight:
            continue
        relevant_partners.append(i)
    lim_df = df[df[opposite].isin(relevant_partners)]
    return lim_df


def coords_to_region(coords, region):
    """Limits an array of coords to a given region.
    
    Parameters
    ----------
    coords: list-like
        Coordinates to limit.
    region: str
        The region to limit the box to.
        
    Returns
    -------
    lim_coords: np.array
        The same coordinates if they fall within the region.
    """
    oob_simplex = np.array([[0,0,0], [0,0,1], [0,1,0], [1,0,0]])
    include = region_dict[region]["include"]
    exclude = region_dict[region]["exclude"]
    
    if include.size != 0:
        include = Delaunay(points = include)
    else:
        include = Delaunay(point = oob_simplex)
    if exclude.size != 0:
        exclude = Delaunay(points = exclude)
    else:
        exclude = Delaunay(points = oob_simplex)
    
    lim_coords = np.array([], dtype="U50")
    for i in coords:        
        temp_coord = [int(x.strip()) for x in i.split(",")]
        if include.find_simplex(temp_coord)!=-1 and \
                exclude.find_simplex(temp_coord)==-1:
            lim_coords = np.append(lim_coords, i)
    return lim_coords
            
        
def lines_to_region(coords1, coords2, region):
    """Takes in two coord sets and returns the ones that are both within the
    given region.

    Parameters
    ----------
    coords1 : list-like
        One set of neurons.
    coords2 : list-like
        The other set of neurons
    region : str
        Region to search through.

    Returns
    -------
    lim_coords1 : np.array
        The first set of coords limited.
    lim_coords2 : np.array
        The second set of coords limited.
    """
    lim_coords1 = np.array([], dtype="U50")
    lim_coords2 = np.array([], dtype="U50")
    for i, j in zip(coords1, coords2):
        temp_coord1 = coords_to_region(np.array([i], dtype="U50"), region)
        temp_coord2 = coords_to_region(np.array([j], dtype="U50"), region)
        if temp_coord1.size>0 and temp_coord2.size>0:
            lim_coords1 = np.append(lim_coords1, temp_coord1[0])
            lim_coords2 = np.append(lim_coords2, temp_coord2[0])
    return lim_coords1, lim_coords2
            
        
    
    




>>>>>>> Stashed changes
