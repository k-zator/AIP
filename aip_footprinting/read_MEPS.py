"""
Creates instance of MEPS class for each cube file. Identifies division of the surface based on atom (type).
@author: Katarzyna Joanna Zator (kz265) and Maria Chiara Storer (mcs92)
"""

import logging
import numpy as np
import pandas as pd
from mendeleev import element
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

from aip_footprinting.atom_type_reader import ATReader
from aip_footprinting.constants import bohr_to_Angstrom

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


class MEPS():
    """Class describing MEPS on a specified isosurface, splitting it on atom ownership and determining areas,
       those will be handed over to the Footprinting class for AIP calculation
       It uses input surface from a cube(esque) file - it only contains one isosurface and in Cartesians."""

    def __init__(self, cube_file, cml_file, own_dist_scaled=True):
        """
        Parameters
        ----------
        cube_file : string
            path to the .cube file of NWChem origin for the MEPS coordinates
        cml_file : string
            path to the .cml file for atom type information of the ligand
        own_dist_scaled : boolean
            default: True; if True the closest atom neighbour is decided based a relative distance, 
            scaled by the vdW surface as defined by vdW radii, and otherwise uses absolute distances
        """
        self._cube_file = cube_file
        self._cml_file = cml_file
        self._get_MEPS_ownership(own_dist_scaled)
        self.get_areas()

    def _get_Atoms_df_and_MEPS_df(self):
        """
        Reads the cube file and returns two dataframes

        Atoms_df: pandas dataframe with info of the atoms
                  atomic mass, charge, x, y, z, atom type, atom name 

        MEPS_df : pandas dataframe with info of the MEPS
                  x, y, z, MEPS charge
        """
        open_file = open(self._cube_file, "r")
        read_file = open_file.readlines()
        if len(read_file[0].split(" ")) > 4:
            if type(read_file[0].split(" ")[4]) == float:
                self.vdw_volume = float(read_file[0].split(" ")[4]) * np.power(0.529177, 3)
            else:
                self.vdw_volume = 0
        else:
            self.vdw_volume = 0
        data = pd.DataFrame([[float(x) for x in e.split()] for e in read_file[6:]],
                            columns=['amass', 'charge', 'x', 'y', 'z'])
        no_atoms = int(read_file[2].split()[0])
        open_file.close()
        Atoms_df = data[:no_atoms].reset_index(drop=True)
        Atoms_df.loc[:, ['x', 'y', 'z']] = Atoms_df[[
            'x', 'y', 'z']] * bohr_to_Angstrom

        # add atom type information
        self._cml = ATReader(self._cml_file)
        Atoms_df["atype"] = self._cml.aipAtomType
        Atoms_df["aname"] = self._cml.atoms

        MEPS_df = data[no_atoms:].reset_index(drop=True)
        MEPS_df.columns = ['x', 'y', 'z', 'charge', 'e']
        MEPS_df = MEPS_df.drop('e', axis=1)
        MEPS_df.loc[:, ['x', 'y', 'z']] = MEPS_df[[
            'x', 'y', 'z']] * bohr_to_Angstrom

        if len(MEPS_df) == 0:
            LOGGER.error(
                "\n The cube file does not cointain any MEPS points. Exiting")
            exit(1)
        return Atoms_df, MEPS_df

    @staticmethod
    def _find_MEPS_owners(MEPS_df, Atoms_df, own_dist_scaled=True):
        '''Finds the closest atom (by index) for all the MEPS points'''
        atom_xyz = Atoms_df[['x', 'y', 'z']].to_numpy()
        meps_xyz = MEPS_df[['x', 'y', 'z']].to_numpy()
        radii = [element(int(a)).vdw_radius/100 for a in Atoms_df.amass]
        #radii = [1.3 if round(r,2)==1.10 else r for r in radii] #extended H by 0.1 A
        if own_dist_scaled == True:
            dists = cdist(meps_xyz, atom_xyz) - radii
        else:
            dists = cdist(meps_xyz, atom_xyz)
        indexes = np.array([np.argmin(i) for i in np.array(dists)])
        return indexes

    def _get_MEPS_ownership(self, own_dist_scaled):
        """Reads in .cube and .cml file into MEPS_df and Atoms_df and defines ownership of 
           each MEPS point based on its proximity to the van der Waals surface of the atom.
           It moreover ensures that no negative MEPS points are assigned to hydrogen atoms"""

        # define dataframes
        self.Atoms_df, self.MEPS_df = self._get_Atoms_df_and_MEPS_df()

        self.no_atoms = len(self.Atoms_df)
        self.no_MEPS = len(self.MEPS_df)

        # assign ownership of MEPS fragments
        indexes = self._find_MEPS_owners(self.MEPS_df, self.Atoms_df, own_dist_scaled)
        self.MEPS_owner = indexes


    def get_areas(self):
        """Calculation of the surface area MEPS using ConvexHull: looks for the convex figure 
           with edges joining the coordinate points"""
        hull = ConvexHull(np.array(self.MEPS_df[['x', 'y', 'z']]))
        self.total_area = hull.area

    def visualise(self, method):
        """Visualisation of the MEPS according to the method: "charge", "atom_owner", "sign" """
        from mayavi import mlab

        if method == "charge":
            limit = np.max((np.abs(self.MEPS_df["charge"].min()),
                            np.abs(self.MEPS_df["charge"].max())))
            nodes = mlab.points3d(self.MEPS_df["x"].T, self.MEPS_df["y"].T,
                                  self.MEPS_df["z"].T, self.MEPS_df["charge"].T, colormap="RdBu",
                                  scale_mode="none", vmin=-limit, vmax=limit)
            mlab.colorbar(title="MEPS value", nb_labels=3)
            mlab.show()

        elif method == "atom_owner":
            hyd = []
            car = []
            ox = []
            f = []
            nit = []
            cl = []
            sul = []
            p = []
            for h, i in enumerate(self.MEPS_owner):
                point = self.MEPS_df[['x', 'y', 'z']].iloc[h]
                if self.Atoms_df.iloc[i]["amass"] == 1:
                    hyd.append(point)
                elif self.Atoms_df.iloc[i]["amass"] == 6:
                    car.append(point)
                elif self.Atoms_df.iloc[i]["amass"] == 7:
                    nit.append(point)
                elif self.Atoms_df.iloc[i]["amass"] == 8:
                    ox.append(point)
                elif self.Atoms_df.iloc[i]["amass"] == 9:
                    f.append(point)
                elif self.Atoms_df.iloc[i]["amass"] == 17:
                    cl.append(point)
                elif self.Atoms_df.iloc[i]["amass"] == 16:
                    sul.append(point)
                else:
                    pass

            if len(hyd) > 0:
                hyd = np.array(hyd)
                mlab.points3d(hyd.T[0], hyd.T[1], hyd.T[2], color=(1, 1, 1))
            if len(car) > 0:
                car = np.array(car)
                mlab.points3d(car.T[0], car.T[1], car.T[2],
                              color=(0.5, 0.5, 0.5))
            if len(nit) > 0:
                nit = np.array(nit)
                mlab.points3d(nit.T[0], nit.T[1], nit.T[2],
                              color=(0.1, 0.4, 1))
            if len(ox) > 0:
                ox = np.array(ox)
                mlab.points3d(ox.T[0], ox.T[1], ox.T[2], color=(1, 0.1, 0.1))
            if len(cl) > 0:
                cl = np.array(cl)
                mlab.points3d(cl.T[0], cl.T[1], cl.T[2], color=(0.4, 1, 0.4))
            if len(sul) > 0:
                sul = np.array(sul)
                mlab.points3d(sul.T[0], sul.T[1], sul.T[2], color=(1, 1, 0.2))
            if len(f) > 0:
                f = np.array(f)
                mlab.points3d(f.T[0],  f.T[1],  f.T[2], color=(0.1, 0.9, 0.7))
            mlab.show()

        elif method == "sign":
            positive_points = []
            negative_points = []
            for m in self.MEPS_df.iloc:
                if m['charge'] > 0:
                    positive_points.append(m[['x', 'y', 'z']])
                else:
                    negative_points.append(m[['x', 'y', 'z']])
            positive_points = np.array(positive_points)
            negative_points = np.array(negative_points)
            pts = mlab.points3d(
                positive_points.T[0], positive_points.T[1], positive_points.T[2], color=(0, 0.5, 1))
            pts = mlab.points3d(
                negative_points.T[0], negative_points.T[1], negative_points.T[2], color=(1, 0.1, 0.1))
            mlab.show()
        else:
            print("Please enter valid method")
