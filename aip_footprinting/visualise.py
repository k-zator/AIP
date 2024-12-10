import numpy as np
import pandas as pd
from mayavi import mlab
from scipy import spatial
from aip_footprinting.constants import subset_r_nn


def visualise(self, method="charge", opacity=1, isosurface=0.002):
    """Visualisation of the MEPS and AIPs according to the method.
    Choices include "charge", "percentile", "both", "points" """

    if method == "charge":
        if isosurface == "0.0300":
            MEPS_df = self.MEPS_p.MEPS_df
        elif isosurface == "0.0104":
            MEPS_df = self.MEPS_m.MEPS_df
        else:
            MEPS_df = self.MEPS.MEPS_df
        limit = np.max((np.abs(MEPS_df["charge"].min()), np.abs(
            MEPS_df["charge"].max())))
        nodes = mlab.points3d(MEPS_df["x"].T, MEPS_df["y"].T, MEPS_df["z"].T,
                              MEPS_df["charge"].T, colormap="RdBu", scale_mode="none",
                              vmin=-limit, vmax=limit, opacity=opacity)
        mlab.colorbar(title="MEPS value", nb_labels=3)
        for s in self.AIP:
            pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                color=(0.5, 0.2, 0.7), scale_factor=0.7, vmin=-1, vmax=1)
        mlab.show()

    elif method == "percentile":
        for i in range(self.MEPS.no_atoms):
            subset = self.MEPS.MEPS_df[[
                "x", "y", "z"]][self.MEPS.MEPS_owner == i]
            subset_mat_dis = spatial.distance.cdist(subset, subset)
            LDA = sum(subset_mat_dis < subset_r_nn)/(subset_r_nn**2*np.pi)
            colours = (LDA - LDA.min())/(LDA.max() - LDA.min())
            limit = np.max((0, np.abs(LDA.max())))
            nodes = mlab.points3d(subset["x"].T, subset["y"].T, subset["z"].T, LDA, colormap="coolwarm",
                                  scale_mode="none", vmin=0, vmax=limit)
        mlab.colorbar(title="Local MEPS density value", nb_labels=3)
        for s in self.AIP:
            pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                color=(0.5, 0.2, 0.7), scale_factor=0.7, vmin=-1, vmax=1)
        mlab.show()

    elif method == "both":
        limit = np.max((np.abs(self.MEPS.MEPS_df["charge"].min()), np.abs(
            self.MEPS.MEPS_df["charge"].max())))
        self.non_edge_df = pd.DataFrame()
        for atom in self._Atom.Atom:
            df = self.EdgeDetection(self.MEPS, atom, self.csp)
            self.non_edge_df = self.non_edge_df.append(df, ignore_index=True)

        pts = mlab.points3d(self.non_edge_df['x'], self.non_edge_df['y'], self.non_edge_df['z'],
                            self.non_edge_df["charge"].T, colormap="RdBu", scale_mode="none", vmin=-limit, vmax=limit, opacity=opacity)
        pts = mlab.points3d(
            self.MEPS.Atoms_df['x'], self.MEPS.Atoms_df['y'], self.MEPS.Atoms_df['z'], color=(0.5, 0.5, 0.5))
        for s in self.AIP:
            if s.type == "polar":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.5, 0.7, 0.2), scale_factor=0.7, vmin=-1, vmax=1)
            elif s.type == "non-polar":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.5, 0.2, 0.7), scale_factor=0.7, vmin=-1, vmax=1)
            elif s.type == "outer-polar":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.5, 0.2, 0.7), scale_factor=0.7, vmin=-1, vmax=1)
            elif s.type == "hydrogen":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.3, 0.5, 0.6), scale_factor=0.7, vmin=-1, vmax=1)
            elif s.type == "sigma":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.3, 0.5, 0.6), scale_factor=0.7, vmin=-1, vmax=1)
        mlab.colorbar(title="MEPS value", nb_labels=3)
        mlab.show()

    elif method == "points":
        pts = mlab.points3d(
            self.MEPS.Atoms_df['x'], self.MEPS.Atoms_df['y'], self.MEPS.Atoms_df['z'], color=(0.5, 0.5, 0.5))
        for s in self.AIP:
            if s.type == "polar":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.5, 0.7, 0.2), scale_factor=0.7, vmin=-1, vmax=1)
            elif s.type == "non-polar":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.5, 0.2, 0.7), scale_factor=0.7, vmin=-1, vmax=1)
            elif s.type == "outer-polar":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.5, 0.2, 0.7), scale_factor=0.7, vmin=-1, vmax=1)
            elif s.type == "hydrogen":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.3, 0.5, 0.6), scale_factor=0.7, vmin=-1, vmax=1)
            elif s.type == "sigma":
                pts = mlab.points3d(s.xyz[0][0], s.xyz[0][1], s.xyz[0][2],
                                    color=(0.3, 0.5, 0.6), scale_factor=0.7, vmin=-1, vmax=1)
        mlab.show()
    else:
        print("Please enter valid method")
