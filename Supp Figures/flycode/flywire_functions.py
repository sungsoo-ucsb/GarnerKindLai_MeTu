#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 13:23:17 2024

@author: dustin
"""


import datetime
import pandas as pd
from fafbseg import flywire
import flycode.utils as utils
import flycode.reduction as reduction


client = flywire.get_cave_client()
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
    if region=="":
        return syn_df
    return reduction.lim_region(syn_df, region)


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




