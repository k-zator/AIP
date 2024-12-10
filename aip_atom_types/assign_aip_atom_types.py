import networkx as nx
import networkx.algorithms.isomorphism as iso
from aip_atom_types.fg_graphs import water, amide, vinylogous_amide, carbonyl, aldehyde,\
    ammonia, primary_amine, secondary_amine, alcohol, any_O3, epoxide, noxide, one_lp_1, one_lp_2,\
    tertiary_amine, planar_aniline, primary_pl3_nitrogen, secondary_pl3_nitrogen, tertiary_pl3_nitrogen,  quaternary_nitrogen, quaternary_sulfur, n_sulfoxide, \
    phos, silicon, not_pyridine, nitro, amide_n, sulfone, sulfoxide, selenone, selenoxide, po,  oh, nh, softh, ps, phene, special_car, special_car2, \
    special_c2, special_c1, special_nam, incorrectly_assigned_pl3
import types


def assign_aip_atom_types(network, aromatic=False):
    dict_atoms = {}
    if aromatic:
        match = match_atoms_by_fg_description(network, special_car, dict_atoms)
        for k, v in match.items():
            dict_atoms[k] = v
        match = match_atoms_by_fg_description(
            network, special_car2, dict_atoms)
        for k, v in match.items():
            dict_atoms[k] = v
        match = match_atoms_by_fg_description(network, special_c2, dict_atoms)
        for k, v in match.items():
            dict_atoms[k] = v
        match = match_atoms_by_fg_description(network, special_c1, dict_atoms)
        for k, v in match.items():
            dict_atoms[k] = v
        match = match_atoms_by_fg_description(network, special_nam, dict_atoms)
        for k, v in match.items():
            dict_atoms[k] = v

    match = match_atoms_by_fg_description(network, water, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, ammonia, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, one_lp_1, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, one_lp_2, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, amide_n, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, amide, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(
        network, vinylogous_amide, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, aldehyde, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, carbonyl, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, nitro, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, noxide, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, po, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, alcohol, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, any_O3, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, epoxide, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, sulfone, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, sulfoxide, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(
        network, quaternary_sulfur, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(
        network, n_sulfoxide, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, phos, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, silicon, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, selenone, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, selenoxide, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(
        network, incorrectly_assigned_pl3, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, not_pyridine, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, primary_amine, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, secondary_amine, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, tertiary_amine, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, planar_aniline, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(
        network, primary_pl3_nitrogen, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(
        network, secondary_pl3_nitrogen, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(
        network, tertiary_pl3_nitrogen, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(
        network, quaternary_nitrogen, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, ps, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, phene, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, nh, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, oh, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    match = match_atoms_by_fg_description(network, softh, dict_atoms)
    for k, v in match.items():
        dict_atoms[k] = v
    dict_atoms = update_dict_atoms_with_sybyl(network, dict_atoms)
    return dict_atoms


def match_atoms_by_fg_description(network1, network2, dict_atoms={}):
    nm = node_match_by_category_list(['sybyl', 'elementType'], "ignore_key")
    em = edge_match_by_category('bondOrder', "ignore_key")
    GM = iso.GraphMatcher(network1, network2, node_match=nm, edge_match=em)
    ssip_atom_type_dict = nx.get_node_attributes(network2, 'aipAtomType')
    return update_dict_atoms(GM, dict_atoms, ssip_atom_type_dict)


def update_dict_atoms(GM, dict_atoms, ssip_atom_type_dict):
    for subgraph in GM.subgraph_isomorphisms_iter():
        for k, v in subgraph.items():
            if ssip_atom_type_dict[v] != None and k not in dict_atoms.keys():
                dict_atoms[k] = ssip_atom_type_dict[v]
    return dict_atoms


def update_dict_atoms_with_sybyl(network, dict_atoms):
    sybyl_atom_type_dict = nx.get_node_attributes(network, 'sybyl')
    for atom in network.nodes:
        if atom not in dict_atoms.keys():
            dict_atoms[atom] = sybyl_atom_type_dict[atom]
    return dict_atoms


def node_match_by_category_list(list_attr, default, ignore_key="ignore_key"):
    def match(data1, data2):
        list_rules = []
        for attr in list_attr:
            attribute1 = data1.get(attr, default)
            attribute2 = data2.get(attr, default)
            list_rules.append(matching_rules(
                attribute1, attribute2, ignore_key))
        return all(list_rules)
    return match


def node_match_by_category(attr, default, ignore_key="ignore_key"):
    def match(data1, data2):
        attribute1 = data1.get(attr, default)
        attribute2 = data2.get(attr, default)
        return matching_rules(attribute1, attribute2, ignore_key)
    return match


def matching_rules(attribute1, attribute2, ignore_key):
    if attribute2 == ignore_key:
        return True
    elif type(attribute2) == list or type(attribute2) == tuple:
        list_no = []
        list_yes = []
        for a in attribute2:
            if "!" in a:
                list_no.append(a.replace("!", ""))
            elif "!" not in a:
                list_yes.append(a)
        if len(list_yes) == 0:
            return attribute1 not in list_no
        else:
            return attribute1 in list_yes
    elif "!" in attribute2:
        attribute2_new = attribute2.replace("!", "")
        return attribute1 != attribute2
    else:
        if attribute1 == attribute2:
            return True
        else:
            return False


def edge_match_by_category(attr, default, ignore_key="ignore_key"):
    def match(data1, data2):
        values1 = data1.get(attr, default)
        values2 = data2.get(attr, default)
        if values2 == ignore_key:
            return True
        if str(values1) == str(values2):
            return True
        else:
            return False
    return match
