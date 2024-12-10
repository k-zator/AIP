import numpy as np
import pandas as pd
from aip_footprinting.polar_search import polar_AIP_search as sulfur_AIP_search
from aip_footprinting.define_AIP import define_extreme_AIP, define_geom_sigma_holes, \
    define_kmedoid_AIP_value_critical_point, define_extreme_geometric_AIP, \
    define_kmedoid_AIP_value_dual, define_hydrogen, define_single_cluster
    
def find_connecting_atom(MEPS, Atom, atom):
    bonds = MEPS._cml.list_bonds
    conn_atom_names = [bond.atom2 if (bond.atom1 == atom.atom_name) else
                       bond.atom1 if (bond.atom2 == atom.atom_name) else
                       np.nan for bond in bonds]
    conn_atom_names = [x for x in conn_atom_names if pd.notnull(x)]
    if len(conn_atom_names) > 1:
        conn_atoms = [np.where(Atom.atom_name == c) for c in conn_atom_names]
        atom2 = [Atom.Atom[c[0][0]] for c in conn_atoms]
    else:
        conn_atom = np.where(Atom.atom_name ==
                             conn_atom_names[0])[0][0]
        atom2 = Atom.Atom[conn_atom]
    return atom2


def define_for_hydrogen(AIP, _Atom, MEPS, MEPS_df, atom):
    """For hydrogen which requires only one AIP and use of another MEPS"""
    atom2 = find_connecting_atom(MEPS, _Atom, atom)
    AIP.append(define_hydrogen(MEPS_df, atom, atom2))


def define_for_fluorine(AIP, _Atom, MEPS_p, MEPS_df, atom):
    """For fluorine, find three half-integer AIPs as a custom description, necessarily on 0.0300.
       If surface is positive, assigns value of 0"""
    atom2 = find_connecting_atom(MEPS_p, _Atom, atom)
    AIP.append(define_hydrogen(MEPS_df, atom, atom2, Htype=False))


def define_for_sp1(AIP, _Atom, MEPS, MEPS_df, atom, n1=False):
    """For C.1 carbon with its cylindrical pi system, look for AIPs of the sign of the majority of surface"""
    # place geometrically
    atom2 = find_connecting_atom(MEPS, _Atom, atom)
    if type(atom2) == list:
        atom2 = atom2[0]
    if n1:
        AIP.append(define_extreme_AIP(
            MEPS_df, atom, minimum=False))
    else:
        AIP.append(define_extreme_AIP(
            MEPS_df, atom))
    first_AIP = AIP[-1]
    AIP += define_extreme_geometric_AIP(
        MEPS_df, first_AIP, atom, atom2, 4, polar=False)


def define_for_nonpolar_sulfur(AIP, Atom, MEPS, MEPS_sigma_df, atom):
    """Determine sulfur AIPs depending on atom type"""

    if atom.atom_type == "S.O":
        # exception as highly polarised, thus AIP placed at cluster centre
        A1 = define_single_cluster(MEPS_sigma_df, atom, sigma=True)
        if A1.value is None:
            MEPS_df = MEPS.MEPS_df[MEPS.MEPS_owner == atom.index]
            A1 = define_single_cluster(MEPS_df, atom, sigma=False)
        AIP.append(A1)

    else:  # atom_type == "S.3 or "S.2.phene""
        MEPS_sigma_df = sulfur_AIP_search(
            MEPS, AIP, Atom, atom, MEPS_sigma_df, lp_excl_r=1.74, outer_meps=True)
        # and now add the sigma holes
        [atom2, atom3] = find_connecting_atom(MEPS, Atom, atom)
        aips = define_geom_sigma_holes(MEPS_sigma_df, atom, atom2, atom3)
        Avals = [a.value for a in aips]
        if None in Avals:
            no_polar_AIPs = [s.atom_name for s in AIP].count(atom.atom_name)
            [AIP.pop(-1) for _ in range(no_polar_AIPs)]
            MEPS_df = MEPS.MEPS_df[MEPS.MEPS_owner == atom.index]
            MEPS_df = sulfur_AIP_search(
                MEPS, AIP, Atom, atom, MEPS_df, outer_meps=True)
            aips = define_geom_sigma_holes(MEPS_df, atom, atom2, atom3)
        AIP += aips


def define_for_halides(AIP, Atom, MEPS, MEPS_sigma_df, atom):
    """For remaining halides, look for whole-integer AIP to represent the sigma hole, then ones
       for the electronegative belt of lone pairs represented by half-integer AIPs"""
    # identify the sigma hole - now on middle surface
    A1 = define_single_cluster(MEPS_sigma_df, atom, sigma=True)
    if A1.value is None:
        MEPS_df = MEPS.MEPS_df[MEPS.MEPS_owner == atom.index]
        A1 = define_single_cluster(MEPS_df, atom, sigma=False)
    AIP.append(A1)

    MEPS_df = MEPS.MEPS_df[MEPS.MEPS_owner == atom.index]
    atom2 = find_connecting_atom(MEPS, Atom, atom)
    AIP.append(define_extreme_AIP(MEPS_df, atom))
    AIP += define_extreme_geometric_AIP(MEPS_df,
                                        AIP[-1], atom, atom2, 4, polar=False)


def define_for_pi_system(AIP, MEPS_df, atom, dualAIP, sigma=False):
    """For pi system-containing atoms, place AIP according to the KMedoid at the geometric centre of pi cluster"""
    if dualAIP:
        AIP = define_kmedoid_AIP_value_dual(MEPS_df, atom, sigma)
    else:
        AIP = define_kmedoid_AIP_value_critical_point(MEPS_df, atom, sigma)
    return AIP


def define_for_pi_system_polar(AIP, MEPS_df, atom, dualAIP, sigma=False):
    """For pi system-containing polar atoms, look for AIPs with most positive MEP value placed
       at its pi cloud cluster's geometrical centre"""
    if dualAIP:
        AIP = define_kmedoid_AIP_value_dual(MEPS_df, atom, sigma)
    else:
        AIP = define_kmedoid_AIP_value_critical_point(MEPS_df, atom, sigma)
    return AIP


def non_polar_AIP_search(AIP, Atom, MEPS, MEPS_p, MEPS_df, atom, dualAIP):
    """Surveys the 0.002 e.Bohr-3 isosurface of each atom, but 0.0104 for H, and identified extreme values as AIPs.
        Each distinctive atom type has a separate description, including pi systems of two cluster which are
        identied by KMedoid and where AIP is placed in geometric centre of each cluster instead."""
    # as they have already defined and added for this atom
    no_polar_AIPs = [s.atom_name for s in AIP].count(atom.atom_name)

    if atom.atom_type[0] == "H":
        define_for_hydrogen(AIP, Atom, MEPS, MEPS_df, atom)

    elif atom.atom_type == "F":
        define_for_fluorine(AIP, Atom, MEPS, MEPS_df, atom)

    elif atom.atom_type in ["C.1", "N.1"]:
        define_for_sp1(AIP, Atom, MEPS, MEPS_df, atom)

    elif atom.atom_type == "S.2.ps":
        # remaining sigma hole but did not have the polar site excised
        # hence defined as centre to ensure correct positioning
        A1 = define_single_cluster(MEPS_df, atom, sigma=True)
        if A1.value is None:
            MEPS_df = MEPS.MEPS_df[MEPS.MEPS_owner == atom.index]
            A1 = define_single_cluster(MEPS_df, atom, sigma=False)
        AIP.append(A1)

    elif no_polar_AIPs > 0 \
            and ("2" in atom.atom_type or ".ar" in atom.atom_type) \
            and ("one_lp" not in atom.atom_type):
        AIP += define_for_pi_system_polar(AIP, MEPS_df, atom, dualAIP)
        # IS THIS JUSTIFIABLE? OR DO WE WANT ONE_LP TO ALSO BE POLARISED?
    elif no_polar_AIPs > 0 and "one_lp" not in atom.atom_type:
        pass

    elif atom.atom_type in ["S.2.phene", "S.3", "S.O"]:
        define_for_nonpolar_sulfur(AIP, Atom, MEPS, MEPS_df, atom)

    elif atom.atom_type in ["Cl", "Br", "I"]:
        define_for_halides(AIP, Atom, MEPS, MEPS_df, atom)

    else:
        AIP += define_for_pi_system(AIP, MEPS_df, atom, dualAIP)
