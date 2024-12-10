from networkx.algorithms import shortest_path_length
import numpy as np
import sympy as sym
import logging
import time


logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


def compute_cross(c1, c2, c3):
    v1 = c2 - c1
    v2 = c3 - c1
    vcross = np.cross(v1, v2)
    mcross = np.linalg.norm(vcross)
    return mcross


def compute_transformation_matrix(r12, r13, rcross):
    m = np.column_stack((r12, r13, rcross))
    return np.linalg.inv(m)


def compute_weights(ssip_coord, r1, r12, r13, rcross):
    a = ssip_coord[0] - r1
    transf_matrix = compute_transformation_matrix(r12, r13, rcross)
    return np.dot(transf_matrix, a)


def get_anchors(neigh, network, dict_atom):
    anchor_list = []
    h_num = 0
    j = 0
    if neigh.element == 'H':
        h_num += 1
    #start = time.time()
    while len(anchor_list) < 2 and j < 3:
        j += 1
        for i in list(network.nodes):
            if shortest_path_length(network, source=neigh.aname, target=i, weight=None, method='dijkstra') == j:
                if dict_atom[i].element == 'H':
                    h_num += 1
                if h_num < 2 and dict_atom[i].element == 'H':
                    anchor_list.append(i)
                elif dict_atom[i].element != 'H':
                    anchor_list.append(i)
    if len(anchor_list) < 2:
        j = 0
        while len(anchor_list) < 2 and j < 3:
            j += 1
            for i in list(network.nodes):
                if shortest_path_length(network, source=neigh.aname, target=i, weight=None, method='dijkstra') == j and i not in anchor_list:
                    anchor_list.append(i)
    if len(anchor_list) >= 2 and 0.0 < compute_cross(neigh.xyz, dict_atom[anchor_list[0]].xyz, dict_atom[anchor_list[1]].xyz) < 0.01:
        for i in list(network.nodes):
            short_path = shortest_path_length(
                network, source=neigh.aname, target=i, weight=None, method='dijkstra')
            if short_path > 1 and i not in anchor_list:
                anchor_list.append(i)
        anchor_list.remove(anchor_list[1])
    elif len(anchor_list) < 2:
        LOGGER.error("Could not find all the anchors")
        pass
    anchor1 = neigh
    anchor2 = dict_atom[anchor_list[0]]
    anchor3 = dict_atom[anchor_list[1]]
    return anchor1, anchor2, anchor3


def get_weights(ssip, anchor1, anchor2, anchor3):
    r1 = anchor1.xyz
    r12 = anchor2.xyz-anchor1.xyz
    r13 = anchor3.xyz-anchor1.xyz
    rcross = np.cross(r12, r13)
    weights = compute_weights(ssip.xyz, r1, r12, r13, rcross)
    return weights


def get_average_anchors(dict_atom):
    anchor1 = dict_atom['a1']
    anchor2 = dict_atom['a2']
    anchor3 = dict_atom['a3']
    return anchor1, anchor2, anchor3


def get_average_weights(ssip, anchor1, anchor2, anchor3):
    A, B, C = sym.symbols('A,B,C')
    eq1 = sym.Eq(A*anchor1.xyz[0] + B*anchor2.xyz[0] +
                 C*anchor3.xyz[0], ssip.xyz[0][0])
    eq2 = sym.Eq(A*anchor1.xyz[1] + B*anchor2.xyz[1] +
                 C*anchor3.xyz[1], ssip.xyz[0][1])
    eq3 = sym.Eq(A*anchor1.xyz[2] + B*anchor2.xyz[2] +
                 C*anchor3.xyz[2], ssip.xyz[0][2])
    result = sym.solve([eq1, eq2, eq3], (A, B, C))
    weights = np.array([v for v in result.values()])

    return weights
