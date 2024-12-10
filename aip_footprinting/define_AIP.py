"""
A set of functions that define AIP variables and thus defines the AIP class object.
Versions exist for search of AIP with extreme MEP value, or ones dependent on location
within a pi system as split into KMedoid clusters; here AIP can be defined by the extreme
value of MEP or the value of the geometric centre's MEP point.
@author: Katarzyna Joanna Zator (kz265) and Maria Chiara Storer (mcs92)
"""

import logging
import numpy as np
from scipy.spatial import ConvexHull
from scipy import spatial
from scipy.linalg import expm, norm
from sklearn_extra.cluster import KMedoids
from aip_footprinting.AIP_class import AIP as AIPclass
from aip_footprinting.linear_fit_aip import alpha_linear_002, alpha_linear_104, \
    beta_linear_002, beta_linear_03

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


def rotation_matrix(axis, theta):
    """Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians."""
    return expm(np.cross(np.eye(3), axis/norm(axis)*theta))


def excise_lone_pairs(aips, lp_size, MEPS_df_sub):
    """Excise area taken up by lone pairs in the given MEPS so that pi system can be added (at right angles)"""
    first_AIP_xyz = aips[-1].xyz
    MEPS_np_xyz = MEPS_df_sub[["x", "y", "z"]]
    # inter-surface matrix of distances for 1st AIP
    ismd1 = spatial.distance.cdist(first_AIP_xyz, MEPS_np_xyz)
    outer_xyz = MEPS_np_xyz.iloc[
        np.argmin(ismd1)].to_numpy().reshape(1, 3)
    # outer-surface matrix of distances for 1st outer
    osmd1 = spatial.distance.cdist(outer_xyz, MEPS_np_xyz)[0]
    # inter-surface far points; points which could be future AIPs
    MEPS_df_sub = MEPS_df_sub[osmd1 > lp_size]

    if len(aips) == 2:
        second_AIP_xyz = aips[-2].xyz
        MEPS_np_xyz = MEPS_df_sub[["x", "y", "z"]]
        # inter-surface matrix of distances for 2nd AIP
        ismd2 = spatial.distance.cdist(second_AIP_xyz, MEPS_np_xyz)
        outer2_xyz = MEPS_np_xyz.iloc[
            np.argmin(ismd2)].to_numpy().reshape(1, 3)
        # outer-surface matrix of distances for 1st outer
        osmd2 = spatial.distance.cdist(outer2_xyz, MEPS_np_xyz)[0]
        MEPS_df_sub = MEPS_df_sub[osmd2 > lp_size]

    return MEPS_df_sub


def get_AIP_value(atom, AIP_mepsvalue, polar=False, sigma=False, dual=False):
    """Return AIP value as fitted by linear_fit_aip matching aipAtomType and sign"""
    if sigma == True:
        if atom.atom_type in alpha_linear_104.keys():
            c0, c1 = alpha_linear_104[atom.atom_type]
            AIP_value = c0 + c1 * AIP_mepsvalue
            if AIP_value < 0:
                AIP_value = 0
        else:
            c0, c1 = alpha_linear_104["H.soft"]
            AIP_value = c0 + c1 * AIP_mepsvalue
            if AIP_value < 0:
                AIP_value = 0
    elif atom.atom_type[0] is "H":
        if atom.atom_type in alpha_linear_104.keys():
            c0, c1 = alpha_linear_104[atom.atom_type]
            AIP_value = c0 + c1 * AIP_mepsvalue
        else:
            # if we have not fitted the parameters the value defaults to
            c0, c1 = alpha_linear_104["H.soft"]
            AIP_value = c0 + c1 * AIP_mepsvalue
        if AIP_value < 0:
            AIP_value = 0
    elif (polar is False) and (AIP_mepsvalue > 0):
        if atom.atom_type in alpha_linear_002.keys():
            c0, c1 = alpha_linear_002[atom.atom_type]
            AIP_value = c0 + c1 * AIP_mepsvalue
        elif dual:
            c0, c1 = alpha_linear_002["dual"]
            AIP_value = c0 + c1 * AIP_mepsvalue
        else:
            # if we have not fitted the parameters the value defaults to "H.soft" alpha parameter
            c0, c1 = alpha_linear_002["H.soft"]
            AIP_value = c0 + c1 * AIP_mepsvalue
    elif (polar is False) and (AIP_mepsvalue < 0):
        if atom.atom_type in beta_linear_002.keys():
            c0, c1 = beta_linear_002[atom.atom_type]
        else:
            # if we have not fitted the parameters the value defaults to "C.ar" beta parameter
            c0, c1 = beta_linear_002["C.ar"]
        AIP_value = - (c0 + c1 * AIP_mepsvalue)
    elif (polar is True) and (AIP_mepsvalue < 0):
        if atom.atom_type in beta_linear_03.keys():
            c0, c1 = beta_linear_03[atom.atom_type]
        elif atom.atom_type.startswith("N.3"):
            c0, c1 = beta_linear_03["N.3.primary"]
        elif atom.atom_type.startswith("O.3"):
            c0, c1 = beta_linear_03["O.3.any"]
        else:
            c0, c1 = beta_linear_03["O.2.other"]
        AIP_value = - (c0 + c1 * AIP_mepsvalue)
        if AIP_value > 0:
            AIP_value = 0
    elif polar == "outer-polar":
        if atom.atom_type in beta_linear_002.keys():
            c0, c1 = beta_linear_002[atom.atom_type]
        else:
            c0, c1 = beta_linear_002["C.ar"]
        AIP_value = - (c0 + c1 * AIP_mepsvalue)
    else:
        # LOGGER.warning("A polar AIP with a positive MEP: {}, {}. Extrapolating AIP value".format(
        #     atom.atom_type, AIP_mepsvalue))
        if atom.atom_type in beta_linear_03.keys():
            c0, c1 = beta_linear_03[atom.atom_type]
        else:
            c0, c1 = beta_linear_03["O.2.other"]
        AIP_value = - (c0 + c1 * AIP_mepsvalue)
        if AIP_value > 0:
            AIP_value = 0
    return round(AIP_value, 2)


def define_extreme_AIP(df, atom, minimum=True, polar=False, index=None):
    """Creates an AIP class instance for the interaction point defined by extreme value. It uses the dataframe,
       df, given as input to find the maximum or minimum value (as specified) on polar or non-polar surface (as
       specified) for specific atom (here an object of the Atom class), including hydrogen which uses a middle,
       0.0104 e.Bohr-3 surface."""

    if index is not None:
        AIP_index = index
        AIP_mepsvalue = df.loc[index]["charge"]
    elif minimum:
        AIP_mepsvalue = min(df["charge"])
    else:
        AIP_mepsvalue = max(df["charge"])
    AIP_index = df[df["charge"] == AIP_mepsvalue].index[0]
    AIP_loc = df.loc[AIP_index][["x", "y", "z"]].to_numpy()
    AIP_loc = np.reshape(AIP_loc, (1, AIP_loc.size))
    if atom.atom_type[0] is "H":
        AIP_value = get_AIP_value(atom, AIP_mepsvalue)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "hydrogen",
                         atom.index, atom.atom_type, atom.atom_name])
    elif polar is False:
        AIP_value = get_AIP_value(atom, AIP_mepsvalue)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "non-polar",
                         atom.index, atom.atom_type, atom.atom_name])
    elif polar is True:
        AIP_value = get_AIP_value(atom, AIP_mepsvalue, polar=True)
        if AIP_value is None:
            return None
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "polar",
                         atom.index, atom.atom_type, atom.atom_name])


def define_average_AIP(df, atom, minimum=True, polar=False):
    """Creates an AIP class instance for the interaction point defined by extreme value. It uses the dataframe,
       df, given as input to find the maximum or minimum value (as specified) on polar or non-polar surface (as
       specified) for specific atom (here an object of the Atom class), including hydrogen which uses a middle,
       0.0104 e.Bohr-3 surface."""

    AIP_mepsvalue = df["charge"].mean()
    if minimum:
        value_for_index = min(df["charge"])
    else:
        value_for_index = max(df["charge"])
    AIP_index = df[df["charge"] == value_for_index].index[0]

    AIP_loc = df.loc[AIP_index][["x", "y", "z"]].to_numpy()
    AIP_loc = np.reshape(AIP_loc, (1, AIP_loc.size))

    if atom.atom_type[0] is "H":
        AIP_value = get_AIP_value(atom, AIP_mepsvalue)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "hydrogen",
                         atom.index, atom.atom_type, atom.atom_name])
    elif polar is False:
        AIP_value = get_AIP_value(atom, AIP_mepsvalue)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "non-polar",
                         atom.index, atom.atom_type, atom.atom_name])
    elif polar is True:
        AIP_value = get_AIP_value(atom, AIP_mepsvalue, polar=True)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "polar",
                        atom.index, atom.atom_type, atom.atom_name])


def define_extreme_geometric_AIP(df, first_AIP, atom, atom2, no_AIPs, polar=True):
    """Adds AIPs along the atom2-atom axis a regular distance away, they all equal value of the first AIP.
        no_AIPs is the total including the first_AIP"""
    O_xyz = np.array(atom.xyz)
    PS_xyz = np.array(atom2.xyz)
    PS_O_lp_angle = np.arccos((np.dot((O_xyz - PS_xyz), (O_xyz - first_AIP.xyz).reshape(3, 1)) /
                              (np.linalg.norm(O_xyz - PS_xyz) * np.linalg.norm(O_xyz - first_AIP.xyz)))[0])
    lp_circle_r = np.linalg.norm(first_AIP.xyz - O_xyz) \
        * np.sin(np.pi - PS_O_lp_angle)
    angle = 360/no_AIPs
    AIP_object_list = []
    for a in range(1, no_AIPs):
        rot = rotation_matrix((O_xyz - PS_xyz), np.radians(angle*a))
        new_lp_loc = (np.dot(rot, (first_AIP.xyz - O_xyz).T).T + O_xyz)
        closest_MEP = df.iloc[np.argmin(spatial.distance.cdist(
            new_lp_loc, df[["x", "y", "z"]]))].name
        cMEP_xyz = df.loc[closest_MEP][[
            "x", "y", "z"]].to_numpy().reshape(1, 3)
        cMEP_mepsvalue = df.loc[closest_MEP][["charge"]].to_numpy()[0]
        dist_to_aip = np.linalg.norm(new_lp_loc - cMEP_xyz)
        if dist_to_aip < (np.pi*lp_circle_r*0.25):
            if polar == True:
                cMEP_value = get_AIP_value(atom, cMEP_mepsvalue, polar=True)
                AIP_object_list.append(AIPclass(
                    [cMEP_value, cMEP_mepsvalue, cMEP_xyz, closest_MEP, "polar",
                     atom.index, atom.atom_type, atom.atom_name]))
            else:
                cMEP_value = get_AIP_value(atom, cMEP_mepsvalue, polar=False)
                AIP_object_list.append(AIPclass(
                    [cMEP_value, cMEP_mepsvalue, cMEP_xyz, closest_MEP, "non-polar",
                     atom.index, atom.atom_type, atom.atom_name]))
    return AIP_object_list


def define_single_cluster(df, atom, sigma=False):
    df_xyz = df[["x", "y", "z"]].to_numpy()
    kmeans = KMedoids(n_clusters=1, random_state=0,
                      init='k-medoids++').fit(df_xyz)
    AIP_loc = kmeans.cluster_centers_
    AIP_index = df[(df[["x", "y", "z"]] == AIP_loc).all(1)].index[0]
    AIP_mepsvalue = df["charge"].max()
    if sigma == False:
        AIP_value = get_AIP_value(atom, AIP_mepsvalue)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "outer-sigma",
                         atom.index, atom.atom_type, atom.atom_name])
    else:
        AIP_value = get_AIP_value(atom, AIP_mepsvalue, sigma=True)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "sigma",
                         atom.index, atom.atom_type, atom.atom_name])


def define_hydrogen(df, atom, atom2, Htype=True):
    H = np.array(atom.xyz)
    X = np.array(atom2.xyz)
    df_xyz = df[["x", "y", "z"]].to_numpy()

    dist = [norm(np.cross(H-p, X-H))/norm(X-H) for p in df_xyz]
    # get index of the minimum distance
    AIP_index = df.index[np.argmin(dist)]
    AIP_loc = df.loc[AIP_index][["x", "y", "z"]].to_numpy().reshape(1, 3)
    # AIP_index = df[(df[["x", "y", "z"]] == AIP_loc).all(1)].index[0]
    if Htype == True:
        AIP_mepsvalue = df["charge"].max()
        AIP_value = get_AIP_value(atom, AIP_mepsvalue, polar=False)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "hydrogen",
                        atom.index, atom.atom_type, atom.atom_name])
    else:
        AIP_mepsvalue = df["charge"].min()
        AIP_value = get_AIP_value(atom, AIP_mepsvalue, polar=True)
        return AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "polar",
                        atom.index, atom.atom_type, atom.atom_name])


def define_geom_sigma_holes(df, atom, atom2, atom3):
    S = np.array(atom.xyz)
    df_xyz = df[["x", "y", "z"]].to_numpy()
    AIP_object_list = []
    for a in [atom2, atom3]:
        X = np.array(a.xyz)
        dist = [norm(np.cross(S-p, X-S))/norm(X-S) for p in df_xyz]
        # get index of the minimum distance
        AIP_index = df.index[np.argmin(dist)]
        AIP_loc = df.loc[AIP_index][["x", "y", "z"]].to_numpy().reshape(1, 3)
        # AIP_index = df[(df[["x", "y", "z"]] == AIP_loc).all(1)].index[0]
        AIP_mepsvalue = df.charge.loc[AIP_index]
        AIP_value = get_AIP_value(atom, AIP_mepsvalue, sigma=True)
        AIP_object_list.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index, "sigma",
                                         atom.index, atom.atom_type, atom.atom_name]))
    return AIP_object_list


def define_kmedoid_AIP_value_dual(df, atom, sigma):
    """Works on pi system MEPS which is usually split into two clusters on each face of the system.
       The function uses KMedoid to find geometric centres of those clusters.
       The dualAIP value is determined by for the two extreme MEP values (if such exist)    ``for each cluster.
       Because it determines two AIPs, they are returned as a list."""
    df_xyz = df[["x", "y", "z"]].to_numpy()
    kmeans = KMedoids(n_clusters=2, random_state=0,
                      init='k-medoids++').fit(df_xyz)
    locs = kmeans.cluster_centers_

    p0 = df_xyz[kmeans.labels_ == 0]
    p1 = df_xyz[kmeans.labels_ == 1]
    p0h = ConvexHull(p0)
    p1h = ConvexHull(p1)

    # double check another initialisation won't give better separated sites
    if np.sqrt(np.sum(np.power((locs[0] - locs[1]), 2))) < (atom.vdW_radius*np.sqrt(2)) \
            or p0h.area > 4*p1h.area or p1h.area > 4*p0h.area:
        kmeans = KMedoids(n_clusters=2, random_state=0,
                          init='heuristic').fit(df_xyz)
        locs = kmeans.cluster_centers_

        p0 = df_xyz[kmeans.labels_ == 0]
        p1 = df_xyz[kmeans.labels_ == 1]
        p0h = ConvexHull(p0)
        p1h = ConvexHull(p1)
        # if the clusters still have distinctly different sizes or only one distinct patch
        if np.sqrt(np.sum(np.power((locs[0] - locs[1]), 2))) < (atom.vdW_radius*np.sqrt(2)) \
                or p0h.area > 4*p1h.area or p1h.area > 4*p0h.area:
            kmeans = KMedoids(n_clusters=1, random_state=0,
                              init='k-medoids++').fit(df_xyz)
            locs = kmeans.cluster_centers_
            LOGGER.info("Clusters of atom ", atom.index,
                        " highly uneven; using single AIP instead")

    AIP_object_list = []
    df_updated = df.copy()
    df_updated["label"] = list(kmeans.labels_)

    pos_mask = [sum(df_updated[df_updated["label"] == i]["charge"] > 0)
                for i in range(kmeans.n_clusters)]
    neg_mask = [sum(df_updated[df_updated["label"] == i]["charge"] < 0)
                for i in range(kmeans.n_clusters)]

    dual_critetion = 0.1
    for i in range(len(pos_mask)):
        ## DECIDE what the cutoff for considering atom as dual is: whether some % of surface is positive 
        ## now it's 10%
        if neg_mask[i] > dual_critetion*(pos_mask[i]+neg_mask[i]):
            meps = [df_updated[df_updated["label"] == i]["charge"].min()
                        for i in range(kmeans.n_clusters)]
            indices = [df_updated[df_updated["label"] == i]["charge"].idxmin()
                        for i in range(kmeans.n_clusters)]

            dualAIPs = [pos_mask[i] > dual_critetion*(pos_mask[i]+neg_mask[i]) for i in range(len(pos_mask))]
        else:
            meps = [df_updated[df_updated["label"] == i]["charge"].max()
                        for i in range(kmeans.n_clusters)]
            indices = [df_updated[df_updated["label"] == i]["charge"].idxmax()
                        for i in range(kmeans.n_clusters)]
            dualAIPs = [False for i in range(len(pos_mask))]
    
    meps_max = [df_updated[df_updated["label"] == i]["charge"].max()
                        for i in range(kmeans.n_clusters)]
    indices_max = [df_updated[df_updated["label"] == i]["charge"].idxmax()
                        for i in range(kmeans.n_clusters)]

    if sigma == False:
        values = [get_AIP_value(atom, i) for i in meps]
        values_max = [get_AIP_value(atom, i, dual=True) for i in meps_max]
        for s, AIP_value, AIP_mepsvalue, AIP_index, dualAIP, AIP_value_max, AIP_mepsvalue_max, AIP_index_max \
            in zip(locs, values, meps, indices, dualAIPs, values_max, meps_max, indices_max):
            AIP_loc = s.reshape(1, 3)
            AIP_object_list.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index,
                                    "non-polar", atom.index, atom.atom_type, atom.atom_name]))
            if dualAIP:
                AIP_object_list[-1].set_dual(AIP_value_max, AIP_mepsvalue_max, AIP_index_max)

    else:
        values = [get_AIP_value(atom, i, sigma=True) for i in meps]
        for s, AIP_value, AIP_mepsvalue, AIP_index in zip(locs, values, meps, indices):
            AIP_loc = s.reshape(1, 3)
            AIP_object_list.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index,
                                    "sigma", atom.index, atom.atom_type, atom.atom_name]))
    return AIP_object_list


def define_kmedoid_AIP_value_critical_point(df, atom, sigma, minimum=True, sp1=False):
    """Works on pi system MEPS which is usually split into two clusters on each face of the system.
       The function uses KMedoid to find geometric centres of those clusters.
       The AIP value is, however, determined by the extreme MEP value for each cluster.
       Because it determines two AIPs, they are returned as a list."""
    df_xyz = df[["x", "y", "z"]].to_numpy()
    kmeans = KMedoids(n_clusters=2, random_state=0,
                      init='k-medoids++').fit(df_xyz)
    locs = kmeans.cluster_centers_

    p0 = df_xyz[kmeans.labels_ == 0]
    p1 = df_xyz[kmeans.labels_ == 1]
    p0h = ConvexHull(p0)
    try:
        p1h = ConvexHull(p1)
    except:
        p1h = p0h
    # double check another initialisation won't give better separated sites
    if np.sqrt(np.sum(np.power((locs[0] - locs[1]), 2))) < (atom.vdW_radius*np.sqrt(2)) \
            or p0h.area > 4*p1h.area or p1h.area > 4*p0h.area:
        kmeans = KMedoids(n_clusters=2, random_state=0,
                          init='heuristic').fit(df_xyz)
        locs = kmeans.cluster_centers_

        p0 = df_xyz[kmeans.labels_ == 0]
        p1 = df_xyz[kmeans.labels_ == 1]
        p0h = ConvexHull(p0)
        p1h = ConvexHull(p1)
        if np.sqrt(np.sum(np.power((locs[0] - locs[1]), 2))) < (atom.vdW_radius**np.sqrt(2)) \
                or p0h.area > 4*p1h.area or p1h.area > 4*p0h.area:
            kmeans = KMedoids(n_clusters=2, random_state=0,
                              init='random').fit(df_xyz)
            locs = kmeans.cluster_centers_

            p0 = df_xyz[kmeans.labels_ == 0]
            p1 = df_xyz[kmeans.labels_ == 1]
            p0h = ConvexHull(p0)
            p1h = ConvexHull(p1)
        # if the clusters still have distinctly different sizes or only one distinct patch
        if np.sqrt(np.sum(np.power((locs[0] - locs[1]), 2))) < (atom.vdW_radius*np.sqrt(2)) \
                or p0h.area > 4*p1h.area or p1h.area > 4*p0h.area:
            kmeans = KMedoids(n_clusters=1, random_state=0,
                              init='k-medoids++').fit(df_xyz)
            locs = kmeans.cluster_centers_
            LOGGER.info("Clusters of atom ", atom.index,
                        " highly uneven; using single AIP instead")

    if sp1 == True:
        kmeans = KMedoids(n_clusters=4, random_state=0,
                          init='heuristic').fit(df_xyz)
        locs = kmeans.cluster_centers_

    AIP_object_list = []
    df_updated = df.copy()
    df_updated["label"] = list(kmeans.labels_)

    if minimum:
        meps = [df_updated[df_updated["label"] == i]["charge"].min()
                for i in range(kmeans.n_clusters)]
        indices = [df_updated[df_updated["label"] == i]["charge"].idxmin()
                for i in range(kmeans.n_clusters)]
    else:
        meps = [df_updated[df_updated["label"] == i]["charge"].max()
                for i in range(kmeans.n_clusters)]
        indices = [df_updated[df_updated["label"] == i]["charge"].idxmax()
                for i in range(kmeans.n_clusters)]    

    if sigma == False:
        values = [get_AIP_value(atom, i) for i in meps]
        for s, AIP_value, AIP_mepsvalue, AIP_index in zip(locs, values, meps, indices):
            AIP_loc = s.reshape(1, 3)
            AIP_object_list.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index,
                                             "non-polar", atom.index, atom.atom_type, atom.atom_name]))
    else:
        values = [get_AIP_value(atom, i, sigma=True) for i in meps]
        for s, AIP_value, AIP_mepsvalue, AIP_index in zip(locs, values, meps, indices):
            AIP_loc = s.reshape(1, 3)
            AIP_object_list.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index,
                                             "sigma", atom.index, atom.atom_type, atom.atom_name]))

    return AIP_object_list


def define_kmedoid_AIP_value_polar(df, atom, sigma):
    """Works on pi system MEPS which is usually split into two clusters on each face of the system.
       The function uses KMedoid to find geometric centres of those clusters.
       The AIP value is, however, determined by the extreme MEP value for each cluster.
       Because it determines two AIPs, they are returned as a list"""
    df_xyz = df[["x", "y", "z"]].to_numpy()
    kmeans = KMedoids(n_clusters=2, random_state=0,
                      init='k-medoids++').fit(df_xyz)
    locs = kmeans.cluster_centers_

    p0 = df_xyz[kmeans.labels_ == 0]
    p1 = df_xyz[kmeans.labels_ == 1]
    p0h = ConvexHull(p0)
    p1h = ConvexHull(p1)
    # double check another initialisation won't give better separated sites
    if atom.atom_type == "S.3" or atom.atom_type == "S.2.phene":
        dist_limit = (atom.vdW_radius*np.sqrt(2)/2)
    else:
        dist_limit = (atom.vdW_radius*np.sqrt(2))

    if np.sqrt(np.sum(np.power((locs[0] - locs[1]), 2))) < dist_limit \
            or p0h.area > 4*p1h.area or p1h.area > 4*p0h.area:
        kmeans = KMedoids(n_clusters=2, random_state=0,
                          init='heuristic').fit(df_xyz)
        locs = kmeans.cluster_centers_

        p0 = df_xyz[kmeans.labels_ == 0]
        p1 = df_xyz[kmeans.labels_ == 1]
        p0h = ConvexHull(p0)
        p1h = ConvexHull(p1)
        if np.sqrt(np.sum(np.power((locs[0] - locs[1]), 2))) < dist_limit \
                or p0h.area > 4*p1h.area or p1h.area > 4*p0h.area:
            kmeans = KMedoids(n_clusters=2, random_state=0,
                              init='random').fit(df_xyz)
            locs = kmeans.cluster_centers_

            p0 = df_xyz[kmeans.labels_ == 0]
            p1 = df_xyz[kmeans.labels_ == 1]
            p0h = ConvexHull(p0)
            p1h = ConvexHull(p1)
            # if the clusters still have distinctly different sizes or only one distinct patch
            if np.sqrt(np.sum(np.power((locs[0] - locs[1]), 2))) < dist_limit \
                    or p0h.area > 4*p1h.area or p1h.area > 4*p0h.area:
                kmeans = KMedoids(n_clusters=1, random_state=0,
                                  init='k-medoids++').fit(df_xyz)
                locs = kmeans.cluster_centers_
                LOGGER.info("Clusters of atom ", atom.index,
                            " highly uneven; using single AIP instead")

    AIP_object_list = []
    df_updated = df.copy()
    df_updated["label"] = list(kmeans.labels_)

    meps = [df_updated[df_updated["label"] == i]["charge"].max()
            for i in range(kmeans.n_clusters)]
    indices = [df_updated[df_updated["label"] == i]["charge"].idxmax()
               for i in range(kmeans.n_clusters)]
    if sigma == False:
        values = [get_AIP_value(atom, i) for i in meps]
        for s, AIP_value, AIP_mepsvalue, AIP_index in zip(locs, values, meps, indices):
            AIP_loc = s.reshape(1, 3)
            AIP_object_list.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index,
                                             "non-polar", atom.index, atom.atom_type, atom.atom_name]))
    else:
        values = [get_AIP_value(atom, i, sigma=True) for i in meps]
        for s, AIP_value, AIP_mepsvalue, AIP_index in zip(locs, values, meps, indices):
            AIP_loc = s.reshape(1, 3)
            AIP_object_list.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, AIP_index,
                                             "sigma", atom.index, atom.atom_type, atom.atom_name]))

    return AIP_object_list
