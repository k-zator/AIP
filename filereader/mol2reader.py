import numpy as np

class Mol2Reader():
    def __init__(self, string, from_file):
        if from_file==True:
            self._mol2_string_ = self._get_file_as_string_(string)
        else:
            self._mol2_string_ = string
        self._list_lines_ = self._get_list_lines_(self._mol2_string_)
        self._main_lines_ = self._get_main_lines_(self._list_lines_)
        self.list_atoms = self._get_list_atoms_(self._main_lines_)

    @staticmethod
    def _get_list_lines_(mol2_string):
        return mol2_string.split("\n")
    
    
    @staticmethod
    def _get_line_index_(list_lines, search_string):
        for index, line in enumerate(list_lines):
            if search_string in line:
                return index
    
    @staticmethod
    def _get_file_as_string_(filepath):
        with open(filepath, "r") as readfile:
            return "".join(readfile.readlines())
    
    def _get_start_end_(self, list_lines):
        start = self._get_line_index_(list_lines, "@<TRIPOS>ATOM")+1
        end = self._get_line_index_(list_lines, "@<TRIPOS>BOND")
        return start, end
    
    def _get_main_lines_(self, list_lines):
        start, end = self._get_start_end_(list_lines)
        return list_lines[start:end]
    
    
    @staticmethod
    def _get_list_atoms_(main_lines):
        list_atoms = []
        for line in main_lines:
            length_line = len(line.split())
            if length_line == 9:
                atom = Mol2Atom(line)
                list_atoms.append(atom)
            else:
                return list_atoms
        return list_atoms


class Mol2Atom():
    def __init__(self, line):
        split = self._get_split_line_(line)
        self.element = split[1]
        self.xyz = np.array([
            float(split[2]), 
            float(split[3]), 
            float(split[4])
            ])
        self.sybyl = split[5]


    @staticmethod
    def _get_split_line_(line):
        return line.split()
