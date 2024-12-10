#
from pep562 import pep562


"""
The linear fits are discussed in the PhD thesis by M. C. Storer
Parameters for calculation of AIP values, the c0 and c1 coefficients in the linear equations for specific atom types
The fits for the pi systems and alkyl hydrogens depend on the water model that is selected for the n-hexadecane to water phase transfer predictions.
Currently, the descriptions selected assume that solvent water has an alpha_s = 2.8 and a beta_s = 3.9, however, as discussed in the thesis, the description of solvent water may need further optimisation.
"""


use_dnn = False


water_model_A = {"soft_CH_alpha": [-0.6359573747831553, 20.480692125613984], "soft_CH_alpha_dnn": [
    -0.019354775078593345, 36.876772801702586], "dual_pi": [0, 37.24], "positive_pi": [0, 65], "negative_pi": [0, -65]}
water_model_C = {"soft_CH_alpha": [-1.0699570085894077, 23.888092751655694], "soft_CH_alpha_dnn": [
    -0.349159683161984, 42.890707095854104], "positive_pi": [0, 50], "negative_pi": [0, -55]}

selected_water = water_model_A

dft_constants = {
    "alpha_linear_002": {
        "H.soft":     selected_water["positive_pi"],
        "H.O":        [-1.12762720990721, 48.2496376196178],
        "H.O.water":  [2.8,                0],
        "H.N":        [-0.364029,         39.386340],
        "H.N.group2": [0.234542,          36.550231],
        "dual":     selected_water["dual_pi"],},

    "alpha_linear_104": {
        # sigma holes default to H.soft
        "H.soft":     selected_water["soft_CH_alpha"],
        "H.O":        [-2.801581219562319, 34.63350687],
        "H.O.water":  [2.8,                0],
        "H.N":        [-2.801581219562319, 34.63350687],
        "H.N.group2": [-2.801581219562319, 34.63350687],
        "Cl":         [-2.801581219562319, 34.63350687],
        "Br":         [-2.801581219562319, 34.63350687],
        "I":          [-2.801581219562319, 34.63350687],
        "S.3":        [-2.801581219562319, 34.63350687],
        "S.2.phene":  [-2.801581219562319, 34.63350687],
        "S.2.ps":     [-2.801581219562319, 34.63350687],
        "S.O":        [-2.801581219562319, 34.63350687]},

    "beta_linear_all_002": {
        # S.2.phene, S.2.ps, Cl, Br, I default to C.ar
        "C.ar":           selected_water["negative_pi"],
        "N.1":            [-4.131606006121907,  -130.99818527620562],
        "N.2":            [-1.6029221625513292, -106.78643006265129], ##1
        "N.3.aniline":    [-1.6029221625513292, -106.78643006265129], 
        "N.3.primary":    [-5.527027045359226,  -172.62029457370295],
        "N.3.secondary":  [-1.6819870070452145, -140.839394571341],
        "N.3.tertiary":   [-2.4919584688506173, -160.70610443222984],
        "N.ar":           [-1.6029221625513292, -106.78643006265129], ##1
        "O.2.carbonyl":   [-0.42,   -88.21], ##2 
        "O.2.po":         [2.03,	-88.21], ##4
        "O.2.sulfone":    [0.44,	-88.21], ##4
        "O.2.selenone":   [0.44,	-88.21], ##4
        "O.2.am":         [0.44, 	-88.21], ##3
        "O.2.other":      [-0.42,	-88.21], ##5
        "O.2.aldehyde":   [-0.42,	-88.21], ##5
        "O.2.sulfoxide":  [0.44,	-88.21], ##5
        "O.2.selenoxide": [0.44,	-88.21], ##5
        "O.2.nitro":      [-1.33,	-88.21], ##5
        "O.2.noxide":     [1.51,	-88.21], ##6
        "O.3.alcohol":    [-4.668547114931645,  -156.90103580145487],
        "O.3.any":        [-2.190290028287709,  -123.86546946514677],
        "S.2":            [-0.424585,           -102.781660],
        "S.2.phene":      [-0.424585,           -102.781660],
        "S.2.ps":         [-5.183053,           -177.653591]},

    "beta_linear_002": {
        # Cl, Br, I, and all other default to C.ar
        # for the sake of simplicity, S.3 defaults also rather than "S.3":  [0.792843,  -57.11]
        "C.ar": selected_water["negative_pi"]},

    "beta_linear_03": {
        "N.1":            [-0.64,	-88.21],
        "N.2":            [-1.08,	-88.21],
        "N.3.ammonia":    [6.8,       0],
        "N.3.aniline":    [-2.57,	-88.21],
        "N.3.primary":    [-2.94, 	-88.21],
        "N.3.secondary":  [-1.63,	-88.21],
        "N.3.tertiary":   [-0.98, 	-88.21],
        "N.ar":           [-1.08,	-88.21],
        "O.2.carbonyl":   [-0.42,   -88.21],
        "O.2.po":         [2.03,	-88.21],
        "O.2.sulfone":    [0.44,	-88.21],
        "O.2.selenone":   [0.44,	-88.21],
        "O.2.am":         [0.44, 	-88.21],
        "O.2.other":      [-0.42,	-88.21],
        "O.2.aldehyde":   [-0.42,	-88.21],
        "O.2.sulfoxide":  [0.44,	-88.21],
        "O.2.selenoxide": [0.44,	-88.21],
        "O.2.nitro":      [-1.33,	-88.21],
        "O.2.noxide":     [1.51,	-88.21],
        "O.3.water":      [4.50,     0],
        "O.3.alcohol":    [-1.57,   -88.21],
        "O.3.any":        [-0.66,   -88.21],
        "S.2":            [2.94,    -88.21],
        "S.2.phene":      [2.94,    -88.21],
        "S.2.ps":         [4.52,    -88.21],
        "F":              [0.31,    -88.21]}
}

# parameters for calculation of AIP values, the c0 and c1 coefficients in the linear equations for specific atom types
dnn_constants = {"alpha_linear_002": {"H.soft":    selected_water["positive_pi"],
                                      "H.O":       [-1.67, 58.3],
                                      "H.O.water": [2.8,    0],
                                      "H.N":       [-1.06, 58.3],
                                      "Br":        [-1.20, 58.3],
                                      "I":         [-1.20, 58.3]},

                 "alpha_linear_104": {"H.soft":    selected_water["soft_CH_alpha_dnn"],
                                      "H.O":       [-2.9308567162353967, 72.11751794003942],
                                      "H.O.water": [2.8, 0],
                                      "H.N":       [0.316370907988456, 39.113487146716665]},
                 "beta_linear_03": {
    "C.ar":           selected_water["negative_pi"],
    "N.3.ammonia":    [6.8,       0],
    "N.1":            [-3.04,  -117.2],
    "N.ar":           [-0.159, -117.2],
    "N.2":            [-0.159, -117.2],
    "N.3.aniline":    [-1.70,  -117.2],
    "N.3.primary":    [-1.23,  -117.2],
    "N.3.secondary":  [-0.208, -117.2],
    "N.3.tertiary":   [0.844,  -117.2],
    "O.3.water":      [4.5,       0],
    "O.3.alcohol":    [-2.38,  -117.2],
    "O.3.any":        [-2.230, -117.2],
    "S.2":            [-1.746, -117.2],
    "S.2.ps":         [-0.868, -117.2],
    "S.2.phene":      selected_water["negative_pi"],
    "Cl":             selected_water["negative_pi"],
    "Br":             selected_water["negative_pi"],
    "I":              selected_water["negative_pi"],
    "S.3":            [-0.755, -117.2],
    "Se.3":           [-0.755, -117.2],
    "F":              [-2.788, -117.2],
    "O.2.po":         [-1.162, -117.2],
    "O.2.sulfone":    [-2.724, -117.2],
    "O.2.selenone":   [-2.724, -117.2],
    "O.2.am":         [-1.331, -117.2],
    "O.2.other":      [-2.174, -117.2],
    "O.2.aldehyde":   [-2.174, -117.2],
    "O.2.sulfoxide":  [-2.724, -117.2],
    "O.2.selenoxide": [-2.724, -117.2],
    "O.2.nitro":      [-3.831, -117.2],
    "O.2.noxide":     [0.599,  -117.2],
    "O.2.carbonyl":   [-2.724, -117.2]},

    "beta_linear_all_002": {
    "C.ar":           selected_water["negative_pi"],
    "N.3.ammonia":    [6.8,       0],
    "N.1":            [-3.04,  -117.2],
    "N.ar":           [-0.159, -117.2],
    "N.2":            [-0.159, -117.2],
    "N.3.aniline":    [-1.70,  -117.2],
    "N.3.primary":    [-1.23,  -117.2],
    "N.3.secondary":  [-0.208, -117.2],
    "N.3.tertiary":   [0.844,  -117.2],
    "O.3.water":      [4.5,       0],
    "O.3.alcohol":    [-2.38,  -117.2],
    "O.3.any":        [-2.230, -117.2],
    "S.2":            [-1.746, -117.2],
    "S.2.ps":         [-0.868, -117.2],
    "S.2.phene":      selected_water["negative_pi"],
    "Cl":             selected_water["negative_pi"],
    "Br":             selected_water["negative_pi"],
    "I":              selected_water["negative_pi"],
    "S.3":            [-0.755, -117.2],
    "Se.3":           [-0.755, -117.2],
    "F":              [-2.788, -117.2],
    "O.2.po":         [-1.162, -117.2],
    "O.2.sulfone":    [-2.724, -117.2],
    "O.2.selenone":   [-2.724, -117.2],
    "O.2.am":         [-1.331, -117.2],
    "O.2.other":      [-2.174, -117.2],
    "O.2.aldehyde":   [-2.174, -117.2],
    "O.2.sulfoxide":  [-2.724, -117.2],
    "O.2.selenoxide": [-2.724, -117.2],
    "O.2.nitro":      [-3.831, -117.2],
    "O.2.noxide":     [0.599,  -117.2],
    "O.2.carbonyl":   [-2.724, -117.2]},

    "beta_linear_002": {
    "C.ar":           selected_water["negative_pi"],
    "N.3.ammonia":    [6.8,       0],
    "N.1":            [-3.04,  -117.2],
    "N.ar":           [-0.159, -117.2],
    "N.2":            [-0.159, -117.2],
    "N.3.aniline":    [-1.70,  -117.2],
    "N.3.primary":    [-1.23,  -117.2],
    "N.3.secondary":  [-0.208, -117.2],
    "N.3.tertiary":   [0.844,  -117.2],
    "O.3.water":      [4.5,       0],
    "O.3.alcohol":    [-2.38,  -117.2],
    "O.3.any":        [-2.230, -117.2],
    "S.2":            [-1.746, -117.2],
    "S.2.ps":         [-0.868, -117.2],
    "S.2.phene":      selected_water["negative_pi"],
    "Cl":             selected_water["negative_pi"],
    "Br":             selected_water["negative_pi"],
    "I":              selected_water["negative_pi"],
    "S.3":            [-0.755, -117.2],
    "Se.3":           [-0.755, -117.2],
    "F":              [-2.788, -117.2],
    "O.2.po":         [-1.162, -117.2],
    "O.2.sulfone":    [-2.724, -117.2],
    "O.2.selenone":   [-2.724, -117.2],
    "O.2.am":         [-1.331, -117.2],
    "O.2.other":      [-2.174, -117.2],
    "O.2.aldehyde":   [-2.174, -117.2],
    "O.2.sulfoxide":  [-2.724, -117.2],
    "O.2.selenoxide": [-2.724, -117.2],
    "O.2.nitro":      [-3.831, -117.2],
    "O.2.noxide":     [0.599,  -117.2],
    "O.2.carbonyl":   [-2.724, -117.2]},

    "sulfur_lp_002": [-1.746,   -117.2]
}


def __getattr__(name):
    if name in dnn_constants and name in dft_constants:
        if use_dnn:
            return dnn_constants[name]
        else:
            return dft_constants[name]
    else:
        raise AttributeError


pep562(__name__)
