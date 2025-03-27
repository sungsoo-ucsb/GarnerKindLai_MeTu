# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 11:42:24 2023

@author: Dustin Garner
"""


import os
import navis
from cloudvolume import CloudVolume
import flycode.mapping as mapping


link = "graphene://https://prodv1.flywire-daf.com/segmentation/1.0/fly_v31"
cv = CloudVolume(link, progress=True, use_https=True)
patched = False


def get_file_path(file_name):
    """Gets a .ply file path from a file_name.
    
    Parameters
    ----------
    file_name : str
        The name of the file.
        
    Returns
    -------
    file_path : str
        The .ply path in the Meshes directory.
    """
    absolute_path = os.path.dirname(__file__)
    relative_path = os.path.join("Meshes", f"{file_name}.ply")
    file_path = os.path.join(absolute_path, relative_path)
    return file_path


def save_neurs(neur_ids, file_name="", downsample=False):
    """Saves a single neuron or multiple neuron IDs as a single mesh.

    Parameters
    ----------
    neur_ids : int or list-like
        A neuron ID or multiple IDs.
    file_name : str, optional
        The name to store the mesh as. The default is the first neuron ID.
    downsample : bool, optional
        Whether to downsample the neurons. The default is False.
    """
    global patched
    if not patched:
        navis.patch_cloudvolume()
        patched = True
    if isinstance(neur_ids, (str, int)):
        neur_ids = [neur_ids]
    mesh = cv.mesh.get(neur_ids, as_navis=True)
    file_name = str(neur_ids[0]) if file_name=="" else file_name
    if downsample:
        navis.downsample_neuron(mesh,
                               downsampling_factor=50,
                               inplace=True)
    file_path = get_file_path(file_name)
    print(file_path)
    navis.write_mesh(mesh, filepath=file_path)


def save_individually(neur_types, downsample=False):
    """Saves all neurons of each type individually.

    Parameters
    ----------
    neur_type : list-like
        The neuron types in mapping.neur_ids.
    downsample : bool, optional
        Whether to downsample the neurons. The default is False.
    """
    for i in mapping.id_from_types(neur_types):
        save_neurs(i, downsample=downsample)


def save_types(neur_types, downsample=False):
    """Saves each neuron type as a different mesh.

    Parameters
    ----------
    neur_type : list-like
        The neuron types in mapping.neur_ids.
    downsample : bool, optional
        Whether to downsample the neurons. The default is False.
    """
    for i in neur_types:
        save_neurs(mapping.ids_from_types(i),
                   file_name=i,
                   downsample=downsample)


def save_types_together(neur_types, file_name="Full Mesh", downsample=False,):
    """Saves all neuron types together as one mesh.

    Parameters
    ----------
    neur_type : list-like
        The neuron types in mapping.neur_ids.
    file_name : str, optional
        
    downsample : bool, optional
        Whether to downsample the neurons. The default is False.
    """
    save_neurs(mapping.ids_from_types(neur_types),
               file_name=file_name,
               downsample=downsample)



