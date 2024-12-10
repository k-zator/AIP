# parameters for calculation of AIP values, the c0 and c1 coefficients in the linear equations for specific atom types
alpha_linear_002 = {"H.soft": [0., 36.497381],
                    "H.O": [-2.9308567162353967,72.11751794003942],
                    "H.O.water": [2.8, 0],
                    "H.N": [0.316370907988456,39.113487146716665],
                    "H.N.group2": [0.316370907988456,39.113487146716665],
                    # for description of sigma holes
                    "Cl": [-0.7449993596273781,51.80929114708699],
                    "Br": [-0.7449993596273781,51.80929114708699],
                    "I": [-0.7449993596273781,51.80929114708699]}

alpha_linear_104 = {"H.soft": [0., 36.497381],
                    "H.O": [-2.9308567162353967,72.11751794003942],
                    "H.O.water": [2.8, 0],
                    "H.N": [0.316370907988456,39.113487146716665],
                    "H.N.group2": [0.316370907988456,39.113487146716665],
                    # for description of sigma holes
                    "Cl": [-0.7449993596273781,51.80929114708699],
                    "Br": [-0.7449993596273781,51.80929114708699],
                    "I": [-0.7449993596273781,51.80929114708699]}



beta_linear_all_002 = {
    #'C.ar': [0.,             -73.9],
    #polarised water
    'C.ar': [0.,             -64.7],
    'N.3.ammonia': [6.8,                    0],
    'N.1': [-1.6037415351236737,-98.21389950294099],
    'N.ar': [-2.249192, -148.7796],
    'N.3.aniline': [2.909980850554726,-30.57694435270371],
    'N.3.aniline': [2.909980850554726,-30.57694435270371],
    'N.3.primary': [-2.4841467424507595,-133.9918349636577],
    'N.3.secondary': [4.4237806954237575,-51.589591782067444],
    'N.3.tertiary': [6.294446743013814,-17.66046193348477],
    'O.3.water': [4.5,                    0],
    'O.3.alcohol': [-3.55301979108733,-136.81545867373666],
    'O.3.any': [0.016417461095928232,-73.8203084790694],
    'S.2':[-2.1419460628892493,-135.82417864244235],
    'S.2.ps': [3.495461917175738,-38.94304823910551],
    'S.2.phene':[0.,             -52.8],
    'Cl': [-0.15093923358001793,-82.67959155130272],
    'Br': [0.6971005133783847,-61.17322260248237],
    'I': [1.1340786126415336,-48.40306984420406],
    'S.3': [2.0963343381150894,-39.469990784264674],
    'Se.3': [2.0963343381150894,-39.469990784264674],
    'F': [-1.973147578679959,-100.92972433747371],
    'O.2.po': [-3.898359649571889,-147.1844726832033],
    'O.2.sulfone': [-3.898359649571889,-147.1844726832033],
    'O.2.selenone': [-3.898359649571889,-147.1844726832033],
    'O.2.am': [-2.1150506478248046,-127.62805272616036],
    'O.2.other': [-2.3217175997458446,-118.50730624231623],
    'O.2.aldehyde': [-2.3217175997458446,-118.50730624231623],
    'O.2.sulfoxide': [-3.898359649571889,-147.1844726832033],
    'O.2.selenoxide': [-3.898359649571889,-147.1844726832033],
    'O.2.nitro':[-2.3217175997458446,-118.50730624231623],
    'O.2.noxide': [0.06359096764516181,-118.50730624231623],
    'O.2.carbonyl': [-2.3217175997458446,-118.50730624231623],

}
beta_linear_all_03 = {
    #'C.ar': [0.,             -73.9],
    #polarised water
    'C.ar': [0.,             -64.7],
    'N.3.ammonia': [6.8,                    0],
    'N.1': [-1.6037415351236737,-98.21389950294099],
    'N.ar': [-2.249192, -148.7796],
    'N.3.aniline': [2.909980850554726,-30.57694435270371],
    'N.3.aniline': [2.909980850554726,-30.57694435270371],
    'N.3.primary': [-2.4841467424507595,-133.9918349636577],
    'N.3.secondary': [4.4237806954237575,-51.589591782067444],
    'N.3.tertiary': [6.294446743013814,-17.66046193348477],
    'O.3.water': [4.5,                    0],
    'O.3.alcohol': [-3.55301979108733,-136.81545867373666],
    'O.3.any': [0.016417461095928232,-73.8203084790694],
    'S.2':[-2.1419460628892493,-135.82417864244235],
    'S.2.ps': [3.495461917175738,-38.94304823910551],
    'S.2.phene':[0.,             -52.8],
    'Cl': [-0.15093923358001793,-82.67959155130272],
    'Br': [0.6971005133783847,-61.17322260248237],
    'I': [1.1340786126415336,-48.40306984420406],
    'S.3': [2.0963343381150894,-39.469990784264674],
    'Se.3': [2.0963343381150894,-39.469990784264674],
    'F': [-1.973147578679959,-100.92972433747371],
    'O.2.po': [-3.898359649571889,-147.1844726832033],
    'O.2.sulfone': [-3.898359649571889,-147.1844726832033],
    'O.2.selenone': [-3.898359649571889,-147.1844726832033],
    'O.2.am': [-2.1150506478248046,-127.62805272616036],
    'O.2.other': [-2.3217175997458446,-118.50730624231623],
    'O.2.aldehyde': [-2.3217175997458446,-118.50730624231623],
    'O.2.sulfoxide': [-3.898359649571889,-147.1844726832033],
    'O.2.selenoxide': [-3.898359649571889,-147.1844726832033],
    'O.2.nitro':[-2.3217175997458446,-118.50730624231623],
    'O.2.noxide': [0.06359096764516181,-118.50730624231623],
    'O.2.carbonyl': [-2.3217175997458446,-118.50730624231623],

}


beta_linear_002 = {
    #'C.ar': [0.,             -73.9],
    #polarised water
    'C.ar': [0.,             -64.7],
    'N.3.ammonia': [6.8,                    0],
    'N.1': [-1.6037415351236737,-98.21389950294099],
    'N.ar': [-2.249192, -148.7796],
    'N.3.aniline': [2.909980850554726,-30.57694435270371],
    'N.3.aniline': [2.909980850554726,-30.57694435270371],
    'N.3.primary': [-2.4841467424507595,-133.9918349636577],
    'N.3.secondary': [4.4237806954237575,-51.589591782067444],
    'N.3.tertiary': [6.294446743013814,-17.66046193348477],
    'O.3.water': [4.5,                    0],
    'O.3.alcohol': [-3.55301979108733,-136.81545867373666],
    'O.3.any': [0.016417461095928232,-73.8203084790694],
    'S.2':[-2.1419460628892493,-135.82417864244235],
    'S.2.ps': [3.495461917175738,-38.94304823910551],
    'S.2.phene':[0.,             -64.7],
    'Cl': [-0.15093923358001793,-82.67959155130272],
    'Br': [0.6971005133783847,-61.17322260248237],
    'I': [1.1340786126415336,-48.40306984420406],
    'S.3': [2.0963343381150894,-39.469990784264674],
    'Se.3': [2.0963343381150894,-39.469990784264674],
    'F': [-1.973147578679959,-100.92972433747371],
    'O.2.po': [-3.898359649571889,-147.1844726832033],
    'O.2.sulfone': [-3.898359649571889,-147.1844726832033],
    'O.2.selenone': [-3.898359649571889,-147.1844726832033],
    'O.2.am': [-2.1150506478248046,-127.62805272616036],
    'O.2.other': [-2.3217175997458446,-118.50730624231623],
    'O.2.aldehyde': [-2.3217175997458446,-118.50730624231623],
    'O.2.sulfoxide': [-3.898359649571889,-147.1844726832033],
    'O.2.selenoxide': [-3.898359649571889,-147.1844726832033],
    'O.2.nitro':[-2.3217175997458446,-118.50730624231623],
    'O.2.noxide': [0.06359096764516181,-118.50730624231623],
    'O.2.carbonyl': [-2.3217175997458446,-118.50730624231623],

}


sulfur_lp_002 = [0.792843,   -57.11]

beta_linear_03 = beta_linear_all_03