from cloudvolume import CloudVolume
import flycode.mapping as mapping


link = "graphene://https://prodv1.flywire-daf.com/segmentation/1.0/fly_v31"
cv = CloudVolume(link, progress=True, use_https=True)

    
def save_neur(neur_id):
    """Saves a single neuron ID as a mesh.

    Parameters
    ----------
    neur_id : int
        A neuron ID.
    """
    cv.mesh.save(neur_id)


def save_individually(neur_types):
    """Saves all neurons of each type individually.

    Parameters
    ----------
    neur_type : list-like
        The neuron types in mapping.neur_ids.
    """
    for i in mapping.id_from_types(neur_types):
        save_neur(i)


def save_types(neur_types):
    """Saves each neuron type as a different mesh.

    Parameters
    ----------
    neur_type : list-like
        The neuron types in mapping.neur_ids.
    """
    mapping.add_types(neur_types)
    for i in neur_types:
        cv.mesh.save(mapping.neur_ids[i])


def save_types_together(neur_types):
    """Saves all neuron types together as one mesh.

    Parameters
    ----------
    neur_type : list-like
        The neuron types in mapping.neur_ids.
    """
    neur_ids = mapping.ids_from_types(neur_types)
    cv.mesh.save(neur_ids)




