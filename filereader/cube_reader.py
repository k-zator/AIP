import numpy as np
import pandas as pd
from filereader.cml_reader import CmlReader
import mendeleev


class CubeReader:
    def __init__(self, cube_file, scale = 0.529177):
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
            try:
                p = MEPPoint()
                x = float(i[0])
                y = float(i[1])
                z = float(i[2])
                value = float(i[3])
                if np.isnan(x)!=True and np.isnan(y)!=True and np.isnan(z)!=True and np.isnan(value)!=True:
                    p.set_xyz(x, y, z, scale=scale)
                    p.set_value(value)
                    p.set_color(value)
                    list_points.append(p)
            except:
                pass
        return list_points
   
    @staticmethod
    def read_cube_file(cube_file):
        df1 =pd.read_csv(cube_file, delimiter=" ", header = None, skiprows=2, nrows=1)
        df1.dropna(axis=1, how='all', inplace=True)
        df1=df1.T.reset_index().T
        skiplines = int(df1.loc[0][0]+6)
        df =pd.read_fwf(cube_file, header = None, skiprows=skiplines, widths=[13, 13, 13, 13])
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
        threshold = 0.05
        #pos_thr = threshold/2
        pos_thr = threshold
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
        
class TMesh:
    def __init__(self, tmesh_file, scale = 1/0.529177):
        self.df = self.read_tmesh(tmesh_file)
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
            p.set_value(i[7])
            p.set_color(i[7])
            list_points.append(p)
        return list_points
   
    @staticmethod
    def read_tmesh(tmesh_file):
        df1 =pd.read_csv(tmesh_file, delimiter=' ', header = None, nrows=1)
        number_rows = df1[0][0]
        df = pd.read_csv(tmesh_file, delimiter=' ', header = None, nrows=number_rows, skiprows=1)
        df = df.dropna()
        df.astype('float64').dtypes
        return df
   
    @staticmethod
    def get_atoms(cml_file, scale):
        cml = CmlReader(cml_file)
        string_total = ""
        for n, i in enumerate(cml.list_atoms):
            atomic_number = str(mendeleev.element(i.element).atomic_number)
            x = "{:.6f}".format(i.xyz[0]*scale)
            y = "{:.6f}".format(i.xyz[1]*scale)
            z = "{:.6f}".format(i.xyz[2]*scale)
            if len(x) == 8:
                x = " "+x
            if len(y) == 8:
                y=" "+y
            if len(z) == 8:
                z = " " + z
            if len(atomic_number) == 1:
                atomic_number = " "+atomic_number
            string_total += f"   {atomic_number}    0.000000   {x}   {y}   {z}"
            string_total += "\n"
        return string_total, len(cml.list_atoms)
    
    def write_cube(self, out_file, cml_file=None):
        if cml_file != None:
            string_atoms, len_atoms = self.get_atoms(cml_file, self.scale)
            len_atoms = str(len_atoms)
        else:
            string_atoms = None
            len_atoms = 0
        with open(out_file, "w") as writefile:
            writefile.write(f""" NWChem Enclosed Volume: 0 atomic_units
     Gaussian Cube file
    {len_atoms}    0.000000    0.000000    0.000000
    0    0.000000    0.000000    0.000000
    0    0.000000    0.000000    0.000000
    0    0.000000    0.000000    0.000000
""")
            if string_atoms != None:
                writefile.write(string_atoms)

            for i in self.list_points:
                for j in i.xyz:
                    len_1 = 13
                    number_1 = round(j, 7)
                    string_1 = "{:.6f}".format(number_1)
                    len_string_1 = len(string_1)
                    len_spaces_1 = len_1-len_string_1
                    out_1 = " "*len_spaces_1+str(string_1)
                    writefile.write(out_1)
                len_1 = 13
                number_1 = round(i.value, 7)
                string_1 = "{:.6f}".format(number_1)
                len_string_1 = len(string_1)
                len_spaces_1 = len_1-len_string_1
                out_1 = " "*len_spaces_1+str(string_1)
                writefile.write(out_1)
                writefile.write("\n")
