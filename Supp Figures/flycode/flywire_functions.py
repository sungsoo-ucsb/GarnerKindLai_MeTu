#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 13:23:17 2024

@author: Dustin Garner
"""


import datetime
import numpy as np
import pandas as pd
from fafbseg import flywire
import flycode.utils as utils
import flycode.reduction as reduction
import flycode.readfiles as readfiles


client = flywire.get_cave_client(dataset="flat_630") #"flat_630" actually is version 783.
version = client.materialize.version


def locs_to_segments(coords):
    """
    flywire.locs_to_segments() but with the correct materialization version.
    """
    return flywire.locs_to_segments(coords, timestamp=f"mat_{version!s}")


@utils.attempt_func
def fetch_synapses(*args, **kwargs):
    """
    flywire.get_synapses() but with the correct materialization version.
    Pass in a region parameter to limit the search to a region.
    """
    region = kwargs.pop("region", "")
    syn_df = flywire.get_synapses(materialization = int(version), 
                                *args, **kwargs)
    if region!="":
        syn_df = reduction.lim_region(syn_df, region)
    return syn_df


def fetch_cave_synapses(neur_id):
    """
    Gets synapses from Caveclient of a single neuron.
    """
    synapse_table = client.info.get_datastack_info()['synapse_table']
    pre = client.materialize.live_query(synapse_table,
                            datetime.datetime.now(datetime.timezone.utc),
                            filter_equal_dict = {'pre_pt_root_id': neur_id})
    post = client.materialize.live_query(synapse_table,
                            datetime.datetime.now(datetime.timezone.utc),
                            filter_equal_dict = {'post_pt_root_id': neur_id})
    for i in [pre,post]:
        i.attrs = {}
    df = pd.concat([pre,post])
    df = df[df["pre_pt_root_id"]!=df["post_pt_root_id"]]
    df = df[(df["pre_pt_root_id"] != 0) & (df["post_pt_root_id"] != 0)]
    df = df[df["cleft_score"]>=50]
    
    return df


def fetch_boxed_synapses(region, further_limit=True):
    """
    Retrieves all the synapses within a region bounding box.
    """
    regions = readfiles.import_file("Regions", sheet_name = "Boxes")
    regions = regions[regions.Region==region]
    bounding_box = []
    for i in ["min","max"]:
        temp_bound = []
        for j,k in zip("xyz", [4,4,40]):
            temp_bound.append(np.asarray(regions[f"{j}{i}"])[0] * k)
        bounding_box.append(temp_bound)
    columns = [
        "pre_pt_root_id",
        "post_pt_root_id",
        "cleft_score",
        "pre_pt_position",
        "post_pt_position",
        "id",
        "pre_pt_supervoxel_id",
        "post_pt_supervoxel_id"
    ]
    df = client.materialize.query_view(
                    view_name="valid_synapses_nt_v2_view",
                    materialization_version=version,
                    split_positions=True,
                    select_columns=columns,
                    filter_spatial_dict={"post_pt_position": bounding_box})
    df.rename(
        {
            "post_pt_root_id": "post",
            "pre_pt_root_id": "pre",
            "post_pt_position_x": "post_x",
            "post_pt_position_y": "post_y",
            "post_pt_position_z": "post_z",
            "pre_pt_position_x": "pre_x",
            "pre_pt_position_y": "pre_y",
            "pre_pt_position_z": "pre_z",
            "pre_pt_supervoxel_id": "pre_supervoxel",
            "post_pt_supervoxel_id": "post_supervoxel",
            "idx": "id",  # this may exists if we made a join query
            "id_x": "id",  # this may exists if we made a join query
        },
        axis=1,
        inplace=True,
    )
    if further_limit:
        df = reduction.lim_region(df, region)
    return df
    

@utils.time_elapsed
def get_region_neurons(region, min_syns=1):
    """Retrieves all of the neurons within a region.

    Parameters
    ----------
    region : str
        The region to search through.
    min_syns : int, optional
        The minimum number of synapses that a neuron needs to have to be in
        the returned list.

    Returns
    -------
    np.array
        A unique list of all the neurons within the region.
    """
    df = fetch_boxed_synapses(region)
    all_neurs = np.unique(np.concatenate([np.asarray(x) for x in [df.pre,df.post]]))
    limited_neurs = np.array([], dtype=np.int64)
    for i in all_neurs:
        temp_df = df[(df.pre==i) | (df.post==i)]
        if len(temp_df)>=min_syns:
            limited_neurs = np.append(limited_neurs, i)
    return limited_neurs









