#    nwchemcmlutils can generate MEPS files with NWChem and SSIP descriptions.
#    Copyright (C) 2019  Mark D. Driver
#
#    nwchemcmlutils is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
# -*- coding: utf-8 -*-
"""
Script containing a class for representing the unformated Gaussian cube file
produced to store information on MEPS points.

@author: Mark Williamson
@author: Mark Driver
"""

import logging
import copy
import math
import numpy as np
from nwchemcmlutils.DataClasses.DistanceUnits import DistanceUnits
from nwchemcmlutils.DataClasses.Vector3D import Vector3D
from nwchemcmlutils.DataClasses.Cartesian3D import Cartesian3D
from nwchemcmlutils.DataClasses.PropertyValue import PropertyValue
from nwchemcmlutils.DataClasses.PropertyValueList import PropertyValueList
from nwchemcmlutils.DataClasses.Atom import Atom
from nwchemcmlutils.DataClasses.Molecule import Molecule
import nwchemcmlutils.nwUtils.Trig as Trig

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


class GaussianCube(object):
    """A representation of the unformatted Gaussian CUBE format"""

    def __init__(self, **kwargs):
        self.title_line_one = kwargs.get("title_line_one", "Title")
        self.title_line_two = kwargs.get(
            "title_line_two", "Description of property stored in cubefile"
        )
        self.atomic_units_to_angstroms = 0.529177249
        self.dist_units = kwargs.get("dist_units", DistanceUnits.atomic_units)
        self.number_of_atoms = kwargs.get("number_of_atoms")
        self.grid_origin = kwargs.get("grid_origin")
        self.number_of_x_grid_points = kwargs.get("number_of_x_grid_points")
        self.x_grid_vector = kwargs.get("x_grid_vector")
        self.number_of_y_grid_points = kwargs.get("number_of_y_grid_points")
        self.y_grid_vector = kwargs.get("y_grid_vector")
        self.number_of_z_grid_points = kwargs.get("number_of_z_grid_points")
        self.z_grid_vector = kwargs.get("z_grid_vector")
        # New create the atomCoordinates array
        self.atom_coordinates = kwargs.get("atom_coordinates")
        self.property_values = kwargs.get("property_values")

    def __repr__(self):
        """overload __repr__ so we can reproduce object.
        """
        repr_string = """GaussianCube(title_line_one={self.title_line_one!r},
                        title_line_two={self.title_line_two!r},
                        dist_units={self.dist_units!r},
                        number_of_atoms={self.number_of_atoms!r},
                        grid_origin={self.grid_origin!r},
                        number_of_x_grid_points={self.number_of_x_grid_points!r},
                        number_of_y_grid_points={self.number_of_y_grid_points!r},
                        number_of_z_grid_points={self.number_of_z_grid_points!r},
                        x_grid_vector={self.x_grid_vector!r},
                        y_grid_vector={self.y_grid_vector!r},
                        z_grid_vector={self.z_grid_vector!r},
                        atom_coordinates={self.atom_coordinates!r},
                        property_values={self.property_values!r})
""".format(
            self=self
        )
        return repr_string

    def __str__(self):
        """Overload string representation to produce list of values.
        """
        lines = []
        lines.append("numberOfAtoms: " + str(self.number_of_atoms))
        lines.append("Distance Units: " + str(self.dist_units))
        lines.append("gridOrigin: " + str(self.grid_origin))

        lines.append("numberOfXGridPoints: " + str(self.number_of_x_grid_points))
        lines.append("xGridVector: " + str(self.x_grid_vector))
        lines.append("numberOfYGridPoints: " + str(self.number_of_y_grid_points))
        lines.append("yGridVector: " + str(self.y_grid_vector))
        lines.append("numberOfZGridPoints: " + str(self.number_of_z_grid_points))
        lines.append("zGridVector: " + str(self.z_grid_vector))

        lines.append("atomCoordinates: " + str(self.atom_coordinates))

        lines.append("grid: " + str(self.property_values))
        return str("\n".join(lines))

    def __eq__(self, other_cube):
        """Overload __eq__ so that all components need to be equal.
        """
        return (
            self.title_line_one == other_cube.title_line_one
            and self.title_line_two == other_cube.title_line_two
            and self.dist_units == other_cube.dist_units
            and self.number_of_atoms == other_cube.number_of_atoms
            and self.grid_origin == other_cube.grid_origin
            and self.number_of_x_grid_points == other_cube.number_of_x_grid_points
            and self.x_grid_vector == other_cube.x_grid_vector
            and self.number_of_y_grid_points == other_cube.number_of_y_grid_points
            and self.y_grid_vector == other_cube.y_grid_vector
            and self.number_of_z_grid_points == other_cube.number_of_z_grid_points
            and self.z_grid_vector == other_cube.z_grid_vector
            and self.atom_coordinates == other_cube.atom_coordinates
            and self.property_values == other_cube.property_values
        )

    def _read_title_lines(self, title_lines):
        """reads in the title lines from a list of len == 2.
        """
        if len(title_lines) == 2:
            LOGGER.info("Reading in title lines")
            self.title_line_one = title_lines[0].strip()
            self.title_line_two = title_lines[1].strip()
        else:
            LOGGER.info("incorrect dimensions.")
            raise ValueError("Incorrect number of title lines.")

    def _read_atoms_and_origin(self, atom_origin_line):
        """reads in the number of atoms and grid origin.
        """
        value_list = atom_origin_line.split()
        if len(value_list) == 4:
            if value_list[0][0] == "-":
                LOGGER.debug("We have a cube file in angstroms.")
                self.dist_units = DistanceUnits.angstroms
            self.number_of_atoms = abs(int(value_list[0]))
            grid_origin_vals = [float(num) for num in value_list[1:]]
            self.grid_origin = Cartesian3D(grid_origin_vals, unit=self.dist_units)
        else:
            LOGGER.info("Incorrect number of values in #atoms, grid origin line.")
            raise ValueError("Incorrect number of values in line.")

    def _read_x_grid_points_vector(self, x_grid_line):
        """This reads in the x vector and the number of x grid points
        """
        value_list = x_grid_line.split()
        if len(value_list) == 4:
            self.number_of_x_grid_points = int(value_list[0])
            x_vector_vals = [float(num) for num in value_list[1:]]
            self.x_grid_vector = Cartesian3D(x_vector_vals, unit=self.dist_units)
        else:
            LOGGER.info("Incorrect number of values in x grid line")
            raise ValueError("Incorrect number of values in line.")

    def _read_y_grid_points_vector(self, y_grid_line):
        """This reads in the x vector and the number of x grid points
        """
        value_list = y_grid_line.split()
        if len(value_list) == 4:
            self.number_of_y_grid_points = int(value_list[0])
            y_vector_vals = [float(num) for num in value_list[1:]]
            self.y_grid_vector = Cartesian3D(y_vector_vals, unit=self.dist_units)
        else:
            LOGGER.info("Incorrect number of values in y grid line")
            raise ValueError("Incorrect number of values in line.")

    def _read_z_grid_points_vector(self, z_grid_line):
        """This reads in the x vector and the number of x grid points
        """
        value_list = z_grid_line.split()
        if len(value_list) == 4:
            self.number_of_z_grid_points = int(value_list[0])
            z_vector_vals = [float(num) for num in value_list[1:]]
            self.z_grid_vector = Cartesian3D(z_vector_vals, unit=self.dist_units)
        else:
            LOGGER.info("Incorrect number of values in z grid line")
            raise ValueError("Incorrect number of values in line.")

    def _read_atom_line(self, atom_line, id_num=0):
        """This reads in an atom line and returns an array entry.
        """
        value_list = atom_line.split()
        if len(value_list) == 5:
            LOGGER.info("Reading in atom")
            proton_number = int(value_list[0])
            LOGGER.info("Ignoring values_list[1] as this is always 0")
            coord_list = [float(num) for num in value_list[2:]]
            atom = Atom(
                proton_number,
                Cartesian3D([copy.deepcopy(coord_list)], unit=self.dist_units),
                id_num=id_num,
            )
            return atom
        else:
            LOGGER.debug("Incorrect number of values.")
            raise ValueError("Incorrect number of values.")

    def _read_atom_lines(self, atom_line_list):
        """This reads in the atom lines and returns an array.
        """
        if self.number_of_atoms != None:
            if len(atom_line_list) == self.number_of_atoms:
                id_num = 0
                atom_list = []
                LOGGER.debug("Reading Atom lines")
                for atom_line in atom_line_list:
                    atom = self._read_atom_line(atom_line, id_num=id_num)
                    atom_list.append(atom)
                    id_num += 1
                self.atom_coordinates = Molecule(atom_list, unit=self.dist_units)
            else:
                LOGGER.debug("Incorrect number of atom_lines.")
                raise ValueError("Incorrect number of atom_lines.")
        else:
            LOGGER.info("No number of atoms are specified.")
            raise AttributeError("No number of atoms")

    def _read_property_value_line(self, value_line):
        """This reads in a value line.
        """
        value_line_list = value_line.split()
        if len(value_line_list) == 4:
            LOGGER.info("Line has correct dimensions")
            value_list = [float(num) for num in value_line_list]
            position = Cartesian3D(value_list[:-1], unit=self.dist_units)
            property_value = PropertyValue(value_list[-1], position)
            return property_value
        else:
            LOGGER.info("Line has incorrect number of entries.")
            raise ValueError("Incorrect number of values.")

    def _read_property_values_lines(self, property_lines):
        """This reads in all the property lines to the property_values array.
        """
        prop_val_list = []
        for property_line in property_lines:
            property_value = self._read_property_value_line(property_line)
            prop_val_list.append(property_value)
        self.property_values = PropertyValueList(prop_val_list, unit=self.dist_units)

    def read(self, file_name):
        """Generate a GaussianCube object from a .cube file"""
        # TODO- deal with units- au or Ang.
        with open(file_name, "r") as file_in:

            file_contents = file_in.readlines()
            # Read Title
            self._read_title_lines(file_contents[0:2])

            # Number of Atoms, Grid Origin x,y,z
            # 12   -9.933301   -8.667271   -7.381550
            self._read_atoms_and_origin(file_contents[2])

            # number of X grid points, X grid vector
            # 106    0.166257    0.000000    0.000000
            self._read_x_grid_points_vector(file_contents[3])

            # number of Y grid points, Y grid vector
            # 108    0.000000    0.166313    0.000000
            self._read_y_grid_points_vector(file_contents[4])

            # number of Z grid points, Z grid vector
            #  93    0.000000    0.000000    0.165508
            self._read_z_grid_points_vector(file_contents[5])

            # Read atoms in
            self._read_atom_lines(file_contents[6 : 6 + self.number_of_atoms])

            self._read_property_values_lines(file_contents[6 + self.number_of_atoms :])

    def _write_title_lines_to_list(self):
        """function outputs the title lines in a list, ready for writing out
        the file.
        """
        title_lines_out = [" " + self.title_line_one, " " + self.title_line_two]
        return title_lines_out

    def _write_atoms_and_origin_line(self):
        """This writes the line containing the number of atoms and grid origin.
        It returns a string.
        """
        vector_fortran_format = "{:5d}{:12.6f}{:12.6f}{:12.6f}"
        number_of_atoms = self.number_of_atoms
        if self.dist_units == DistanceUnits.angstroms:
            number_of_atoms = -number_of_atoms
        atom_grid_line = vector_fortran_format.format(
            number_of_atoms,
            self.grid_origin.vector_3d[0][0],
            self.grid_origin.vector_3d[1][0],
            self.grid_origin.vector_3d[2][0],
        )
        return atom_grid_line

    def _write_x_grid_points_vector(self):
        """This writes the line containing the number of x points and x grid
        vector. It returns a string.
        """
        vector_fortran_format = "{:5d}{:12.6f}{:12.6f}{:12.6f}"
        x_grid_line = vector_fortran_format.format(
            self.number_of_x_grid_points,
            self.x_grid_vector.vector_3d[0][0],
            self.x_grid_vector.vector_3d[1][0],
            self.x_grid_vector.vector_3d[2][0],
        )
        return x_grid_line

    def _write_y_grid_points_vector(self):
        """This writes the line containing the number of y points and y grid
        vector. It returns a string.
        """
        vector_fortran_format = "{:5d}{:12.6f}{:12.6f}{:12.6f}"
        y_grid_line = vector_fortran_format.format(
            self.number_of_y_grid_points,
            self.y_grid_vector.vector_3d[0][0],
            self.y_grid_vector.vector_3d[1][0],
            self.y_grid_vector.vector_3d[2][0],
        )
        return y_grid_line

    def _write_z_grid_points_vector(self):
        """This writes the line containing the number of z points and z grid
        vector. It returns a string.
        """
        vector_fortran_format = "{:5d}{:12.6f}{:12.6f}{:12.6f}"
        z_grid_line = vector_fortran_format.format(
            self.number_of_z_grid_points,
            self.z_grid_vector.vector_3d[0][0],
            self.z_grid_vector.vector_3d[1][0],
            self.z_grid_vector.vector_3d[2][0],
        )
        return z_grid_line

    def _write_atom_lines(self):
        """This writes the lines for the atom entries.
        This returns a list of the lines.
        """
        mol_line_list = self.atom_coordinates.writeCubeFileLines()
        return mol_line_list

    def _write_property_values_lines(self):
        """This writes the lines for the value entries.
        """
        prop_val_list = self.property_values.writecubeFileLines()
        return prop_val_list

    def _write_file_lines(self):
        """This writes the lines for the file.
        """
        file_lines = self._write_title_lines_to_list()
        file_lines.append(self._write_atoms_and_origin_line())
        file_lines.append(self._write_x_grid_points_vector())
        file_lines.append(self._write_y_grid_points_vector())
        file_lines.append(self._write_z_grid_points_vector())
        atom_lines = self._write_atom_lines()
        for atom_line in atom_lines:
            file_lines.append(atom_line)
        property_lines = self._write_property_values_lines()
        for property_line in property_lines:
            file_lines.append(property_line)
        return file_lines

    def _write_file_lines_string(self):
        """This writes the file lines and then merges them to output a single
        string.
        """
        file_lines = self._write_file_lines()
        file_lines_string = "\n".join(file_lines)
        file_lines_string += "\n"
        return file_lines_string

    def write(self, file_name):
        """Write a Gaussian cube file from a GaussianCube"""
        out_file_lines_string = self._write_file_lines_string()

        with open(file_name, "w") as out_file:
            out_file.write(out_file_lines_string)

    def rotateByMatrix(self, rotation_matrix):
        """Apply a rotation matrix to the GaussianCube, return resulting
        GaussianCube.
        """
        LOGGER.info("rotate atoms about molecule centroid.")
        rot_atom_coords = self.atom_coordinates.rotateMoleculeAboutCentroid(
            rotation_matrix
        )
        # Use atom coordinates to calculate centroid
        trans_vector = self.atom_coordinates.calculateCentroidVector()
        LOGGER.debug("Rotate molecule")
        rot_prop_vals = self.property_values.rotatePropertyValueListAboutVector(
            rotation_matrix, trans_vector
        )
        LOGGER.debug("Take deep copy of original cube")
        # we only rotate the atoms and the property values. This is so it
        # doesn't mess with the cube file visualisation in VMD, JMOL, etc.
        rot_cube = copy.deepcopy(self)
        rot_cube.atom_coordinates = rot_atom_coords
        rot_cube.property_values = rot_prop_vals
        return rot_cube

    def rotateAboutAxis(self, angle_in_degrees, axis_of_rotation):
        """create the rotation matrix, then apply rotation.
        """
        LOGGER.info("Calculating rotation matrix")
        theta_in_radians = math.radians(angle_in_degrees)
        rotation_matrix = Trig.generateRotationMatrix(
            axis_of_rotation, theta_in_radians
        )
        LOGGER.debug("Roation Matrix : %s", rotation_matrix)
        rotated_cube = self.rotateByMatrix(rotation_matrix)
        return rotated_cube

    def addGausianCube(self, other_cube):
        """Function appends the property values of a second cube file to this,
        if the Molecules are the same, i.e. they have the same units and
        orientation.
        """
        LOGGER.debug("This molecules atoms:")
        LOGGER.debug(self.atom_coordinates)
        LOGGER.debug("Other molecules atoms:")
        LOGGER.debug(other_cube.atom_coordinates)
        if self.dist_units == other_cube.dist_units:
            LOGGER.info("Units match")
            if self.atom_coordinates == other_cube.atom_coordinates:
                LOGGER.info("Molecules also match")
                merged_property_vals = self.property_values.appendPropertyValueList(
                    other_cube.property_values
                )
                merged_cube = copy.deepcopy(self)
                merged_cube.property_values = merged_property_vals
                return merged_cube
            else:
                hash_atom_dict = {
                    "self": [
                        (hash(atom), atom) for atom in self.atom_coordinates.atom_list
                    ],
                    "other_cube": [
                        (hash(atom), atom)
                        for atom in other_cube.atom_coordinates.atom_list
                    ],
                }
                LOGGER.debug(hash_atom_dict)
                raise ValueError("Molecules do not match")
        else:
            raise ValueError("Units do not match")

    def calculateVoxelVolume(self):
        """This calculates the volume of a voxel using the grid vectors.
        """
        return (
            np.dot(
                np.cross(
                    self.x_grid_vector.vector_3d.T, self.y_grid_vector.vector_3d.T
                ),
                self.z_grid_vector.vector_3d,
            )[0][0],
            self.dist_units,
        )

    def returnMinValue(self):
        """Function returns the PropertyValue with the minimum value.
        """
        return self.property_values.returnMinValue()

    def returnMaxValue(self):
        """Function returns the PropertyValue with the maximum value.
        """
        return self.property_values.returnMaxValue()

    def returnPlotValues(self, units):
        """Function returns the values required for plotting the 
        """
        if units == self.dist_units:
            plot_values = self.property_values.returnPlotValues()
        elif units == DistanceUnits.angstroms:
            converted_vals = self.property_values.convertPropertyValueListToAng()
            plot_values = converted_vals.returnPlotValues()
        elif units == DistanceUnits.atomic_units:
            converted_vals = self.property_values.convertPropertyValueListToAU()
            plot_values = converted_vals.returnPlotValues()
        return plot_values

    def plot(self, units=DistanceUnits.angstroms):
        """Method generates plot of the cube file values, using a 3D plot
        and colormap.
        
        units is the distance units used in the plots.
        """
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        import matplotlib.pyplot as plt

        plt.rc("text", usetex=True)
        fig = plt.figure()
        ax = fig.gca(projection="3d")
        plot_values = self.returnPlotValues(units=units)
        x_coord = plot_values[0]
        y_coord = plot_values[1]
        z_coord = plot_values[2]
        w_coord = plot_values[3]
        ax.scatter(
            x_coord,
            y_coord,
            z_coord,
            c=w_coord,
            marker="8",
            cmap=cm.coolwarm_r,
            linewidth=0,
        )
        if units == DistanceUnits.angstroms:
            unit_string = r"\AA"
        elif units == DistanceUnits.atomic_units:
            unit_string = r"AU"
        x_label = r"x/" + unit_string
        y_label = r"y/" + unit_string
        z_label = r"z/" + unit_string
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_zlabel(z_label)
        return fig

    def savePlot(self, filename_stem, file_format="eps", **kwargs):
        """Function saves the figure to file.
        """
        figure = self.plot(**kwargs)
        filename = filename_stem + "." + file_format
        figure.savefig(filename, format=file_format)
        return 0
