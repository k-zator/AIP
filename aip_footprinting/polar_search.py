import logging
import numpy as np
import pandas as pd
from scipy import spatial
from aip_footprinting.AIP_class import AIP as AIPclass
from aip_footprinting.define_AIP import define_extreme_AIP, define_extreme_geometric_AIP, get_AIP_value

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


def find_connecting_atom(MEPS, Atom, atom):
    bonds = MEPS._cml.list_bonds
    conn_atom_names = [bond.atom2 if (bond.atom1 == atom.atom_name) else
                       bond.atom1 if (bond.atom2 == atom.atom_name) else
                       np.nan for bond in bonds]
    conn_atom_names = [x for x in conn_atom_names if pd.notnull(x)]
    conn_atom = np.where(Atom.atom_name ==
                         conn_atom_names[0])[0][0]
    atom2 = Atom.Atom[conn_atom]
    return atom2


def find_local_minima(ox, atom, outer_meps=False, radius=1):
    minima = []
    AIPs = []
    ox_xyz = ox[["x", "y", "z"]].to_numpy()
    for index, point in ox.iterrows():
        xyz = point[["x", "y", "z"]].to_numpy().reshape(1, 3)
        charge = point["charge"]
        mat_dis_to_point = spatial.distance.cdist(ox_xyz, xyz).flatten()
        close_to_point = ox[mat_dis_to_point < radius]
        if (close_to_point["charge"].min() >= charge) and charge < 0:
            minima.append(index)

    for index in minima:
        AIP_mepsvalue = ox.loc[index]["charge"]
        AIP_loc = ox.loc[index][["x", "y", "z"]].to_numpy()
        AIP_loc = np.reshape(AIP_loc, (1, AIP_loc.size))
        if outer_meps == False:
            AIP_value = get_AIP_value(atom, AIP_mepsvalue, polar=True)
            AIPs.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, index, "polar",
                                  atom.index, atom.atom_type, atom.atom_name]))
        else:
            AIP_value = get_AIP_value(atom, AIP_mepsvalue, polar="outer-polar")
            AIPs.append(AIPclass([AIP_value, AIP_mepsvalue, AIP_loc, index, "outer-polar",
                                  atom.index, atom.atom_type, atom.atom_name]))

    values = [a.value for a in AIPs]
    if None in values:
        # positive valued-AIPs, so retry with the outer surface
        return None
    elif len(AIPs) == 0:
        # if it didn't find any local minima
        return None
    else:
        return AIPs


def polar_AIP_search(MEPS_p, AIP, Atom, atom, MEPS_df_sub, lp_excl_r=1.74, outer_meps=False):
    """Polar-AIP search at the 0.03 e.Bohr-3 isosurface for a selection of N, O, S, and Se atom types.
       It identifies MEP points with the minimum value. Find 1 or 2 AIPs as specified by atom type."""
    if outer_meps == True:
        lp_size = 0.95
        if atom.atom_type[0] == "S":
            lp_size = 1.29
    else:
        lp_size = lp_excl_r
    MEPS_df_p = MEPS_p.MEPS_df[MEPS_p.MEPS_owner == atom.index]

    # for non-two-well-defined systems or else
    if atom.atom_type in ["O.2.sulfone", "O.2.noxide", "O.2.po", "S.2.ps"]:
        LOGGER.info(f"The {atom.index} oxygen/sulfur ({atom.atom_type}) cannot be properly described by 2 lp, \
                    hence three geometrically placed lone pairs are used instead")
        PS = find_connecting_atom(MEPS_p, Atom, atom)
        if outer_meps == False:
            AIP.append(define_extreme_AIP(MEPS_df_p, atom, polar=True))
            AIP += define_extreme_geometric_AIP(MEPS_df_p,
                                                AIP[-1], atom, PS, 3)
        else:
            AIP.append(define_extreme_AIP(MEPS_df_p, atom))
            AIP += define_extreme_geometric_AIP(MEPS_df_p,
                                                AIP[-1], atom, PS, 3, polar=False)

        if "ps" in atom.atom_type:
            pass  # as it still needs to find the sigma hole
        else:
            MEPS_df_sub = []

    else:
        aips = find_local_minima(
            MEPS_df_p, atom, radius=1, outer_meps=outer_meps)
        if aips == None:
            return MEPS_df_sub
        else:
            AIP += aips
            # and now excise the lone pairs in the outer MEPS
            first_AIP_xyz = AIP[-1].xyz
            MEPS_np_xyz = MEPS_df_sub[["x", "y", "z"]]
            # inter-surface matrix of distances for 1st AIP
            ismd1 = spatial.distance.cdist(first_AIP_xyz, MEPS_np_xyz)
            outer_xyz = MEPS_np_xyz.iloc[
                np.argmin(ismd1)].to_numpy().reshape(1, 3)
            # outer-surface matrix of distances for 1st outer
            osmd1 = spatial.distance.cdist(outer_xyz, MEPS_np_xyz)[0]
            # inter-surface far points; points which could be future AIPs
            MEPS_df_sub = MEPS_df_sub[osmd1 > lp_size]
            if len(aips) == 2 and len(MEPS_df_sub) > 0:
                second_AIP_xyz = AIP[-2].xyz
                MEPS_np_xyz = MEPS_df_sub[["x", "y", "z"]]
                # inter-surface matrix of distances for 2nd AIP
                ismd2 = spatial.distance.cdist(second_AIP_xyz, MEPS_np_xyz)
                outer2_xyz = MEPS_np_xyz.iloc[
                    np.argmin(ismd2)].to_numpy().reshape(1, 3)
                # outer-surface matrix of distances for 1st outer
                osmd2 = spatial.distance.cdist(outer2_xyz, MEPS_np_xyz)[0]
                MEPS_df_sub = MEPS_df_sub[osmd2 > lp_size]

    return MEPS_df_sub
