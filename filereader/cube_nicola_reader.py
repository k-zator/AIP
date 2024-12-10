import numpy as np
import pandas as pd

class CubeNicolaReader:
    def __init__(self, cube_file, scale = 1):
        self.df = self.read_cube_file(cube_file)
        self.scale = scale
        self.list_points = self.get_list_points(self.df, self.scale)
        self.xyz = np.array([i.xyz for i in self.list_points])
        self.colors = np.array([i.color for i in self.list_points])
        self.values = np.array([i.value for i in self.list_points]) 
   
    @staticmethod
    def get_list_points(df, scale):
        list_points = []
        for n, i in df.iterrows():
            p = MEPPoint()
            p.set_xyz(i[0], i[1], i[2], scale=scale)
            p.set_value(i[3])
            p.set_color(i[3])
            list_points.append(p)
        return list_points
   
    @staticmethod
    def read_cube_file(cube_file):
        df =pd.read_csv(cube_file, delimiter=" ", header = None, skiprows=3)
        return df
    

class MEPPoint:
    def __init__(self):
        self.xyz = None
        self.value = None
        self.color = "white"
    def set_xyz(self, x, y, z, scale):
        self.xyz = np.array([x,y,z])*scale
    def set_value(self, value):
        self.value = value
    def set_color(self, value):
        #threshold = 0.04571*2
        threshold = 0.04571*4

        pos_thr = threshold/2
        neg_thr = -threshold
        if value > pos_thr :
            self.color = "rgb(0, 0, 255)"
        elif value < neg_thr:
            self.color = "rgb(255, 0, 0)"
        elif value > 0:
            k = int(round(value*255/(pos_thr)))
            self.color = 'rgb({},{},{})'.format(255-k, 255-k, 255)
        elif value <0:
            k = int(round(value*255/(neg_thr)))
            self.color = 'rgb({},{},{})'.format(255, 255-k, 255-k)
