import os

import numpy as np
from scipy import linalg
from ase import Atoms

from Species import Species, Element

class Config (object):
    """Config class stores attributes relating to the crystal structure

    In addition to storing the crystal structure the class should perform all manipulations
    on the structure such as swapping positions, converting to and from ase Atoms object.
    """

    def __init__(self):
        
        self.title = 'untitled config'
        self.levcfg = 0
        self.pbc = True
        self.natoms = 0
        self.numghost = 0
        self.symbol = []
        self.vectors = np.zeros((3,3))
        self.charge = None
        self.mass = None
        self.pos = None
        self.label = None
        self.total_energy = 0.0

    def reset_config(self):

        self.title = 'untitled config'
        self.levcfg = 0
        self.pbc = 3
        self.natoms = 0
        self.numghost = 0
        self.vectors = np.zeros((3,3))
        self.total_energy = 0.0
        del self.symbol [:]
        self.pos = None
        self.label = None
        self.charge = None
        self.mass = None

    def create_atoms_object(self):
        """create an atoms object removing all ghost atoms and returns the object"""

        nreal = self.natoms - self.numghost
        #print("create ",nreal,self.natoms,self.numghost)
        
        new_pos = np.zeros((nreal, 3), dtype=np.float64)
        symbols = []
        cell = np.zeros((3,3), dtype=np.float64)
        
        #np.copyto(old_pos, basin.positions)
        for i in range(nreal):  # remove the ghost atoms from the atoms object
            symbols.append(self.symbol[i])
            for j in range(3):
                new_pos[i,j] = self.pos[i,j]
        for i in range(3):
            for j in range(3):
                cell[i,j] = self.vectors[i,j]
        return Atoms(symbols = symbols, positions=new_pos, cell= cell, pbc=True)
    
    def update_from_atoms(self, basin: Atoms):
        """receives an atoms object and takes positions and vectors to upadte those of the
           config object. It assumes the ordering of chemical symbols has not been changed
        """

        symbols = basin.get_chemical_symbols()
        positions = basin.get_positions()
        cell = basin.get_cell()
        
        for i in range(len(basin)):
            self.symbol[i] = symbols[i]
            for j in range(3):
               self.pos[i,j] = positions[i,j]
        for i in range(3):
            for j in range(3):
                self.vectors[i,j] = cell[i,j]

        #for i in range(self.natoms):
        #    print("update", self.symbol[i], self.pos[i,:])
    def setup_configuration(self, spec: Species):
        """Takes a list of the element types and puts the mass and charge on each atom
        The function also checks that all the species defined in the basin.xyz are present in the
        element list
        """

        self.charge = np.zeros(self.natoms, dtype=np.float64)
        self.mass = np.zeros(self.natoms, dtype=np.float64)
        self.label = np.zeros(self.natoms, dtype=np.int32)

        for i in range(self.natoms):

            for k in range(spec.get_num_species()):
                element = spec.get_species(k)
                if self.symbol[i] == element.name:
                    self.label[i] = k 
                    self.charge[i] = element.charge
                    self.mass[i] = element.mass

            if "ghost" in self.symbol[i]:
                self.numghost += 1
            elif self.numghost > 1:
                print("the ghost atoms must come last!")
                exit()
                
    
    def get_number_of_atoms(self): #returns the number of atoms in a list to controling program
        return self.natoms 
        
    def get_volume(self):
        """ calculates volume using numpy """
        axb = np.cross(self.vectors[0][:], self.vectors[1][:])
        self.volume = np.dot(axb, self.vectors[2][:])

        return self.volume
        
                  
    #create a copy of the config
    def copy_config(self):
        """ creates a copy of the Config object. The function creates a new object 
            copies arrays using numpy.copyto and then returns the object
        """

        c = Config()
        c.title = self.title
        c.levcfg = self.levcfg
        c.pbc = self.pbc
        c.natoms = self.natoms
        c.numghost = self.numghost
        c.symbol = self.symbol

        c.pos = np.zeros((self.natoms,3), dtype=np.float64)
        c.charge = np.zeros(self.natoms, dtype=np.float64)
        c.mass = np.zeros(self.natoms, dtype=np.float64)
        c.label = np.zeros(self.natoms, dtype=np.int32)
        
        np.copyto(c.vectors, self.vectors)

        c.total_energy = self.total_energy
    
        np.copyto(c.pos, self.pos)
        np.copyto(c.charge, self.charge)
        np.copyto(c.mass, self.mass)

        np.copyto(c.label, self.label)

        return c

    def select_atom(self) -> int:
        """ selects an atom at random and reurns its index number"""
        
        choice = -1
            
        choice = int(self.get_number_of_atoms * np.random.random())

        
        return choice     
        
    def select_atom_of_type(self, typ: str) -> int:
        """ selects an atom of given type at random and reurns its index number"""

        atm = -1 
        choice = -1
            
        found = False

        atm_list = []
    
        for i in range(self.natoms):
            if self.symbol[i] == typ:
                found = True
                atm_list.append(i)
        
        
        if found == False:
            return atm 

        choice = int(len(atm_list) * np.random.random())

        atm = atm_list[choice]

        return atm
    
    def randomise(self, typ1, typ2, positions):
        """ creates a rondom distribution of two types of atoms (e.g. Si/Al) """
        swap_list = []
        shuffle_list = []

        for i in range(self.natoms):
            if self.symbol[i] == typ1 or self.symbol[i] == typ2:
                swap_list.append(i)
                shuffle_list.append(i)

        print("unshuffled ", shuffle_list)

        np.random.shuffle(shuffle_list)
        print("shuffled ", shuffle_list)
        for i in range(len(swap_list)):
            na = swap_list[i]
            nb = shuffle_list[i]
            self.pos[na,:] = positions[nb,:]

    def find_num_types(self, typ: str) -> int:
        """ determines the number of a given type of atom and returns it as an integer """
        num_typ = 0
        for i in range(self.natoms):
            if typ == self.symbol[i]:
                num_typ += 1

        return num_typ
    

    def swap_atom_positions(self, atm1: int, atm2: int):
        """ swaps the positions of two atoms """

        for j in range(3):
            tmp = self.pos[atm1,j]
            self.pos[atm1,j] = self.pos[atm2,j]
            self.pos[atm2,j] = tmp

    
    def swap_atom_types(self, atm1: int, atm2: int):
        """ swaps the atom types - not does not work if ghosts are present """

        tmp = self.symbol[atm1]
        self.symbol[atm1] = self.symbol[atm2]
        self.symbol[atm2] = tmp

    def mutate_atom(self, atm: int, typ: int, spec: Species):
        ele = spec.get_species(typ)
        self.symbol[atm] = ele.name
        self.mass[atm] = ele.mass
        self.charge[atm] = ele.charge
        self.label[atm] = typ         
        
    def displace_atoms(self, delta):
        for i in range(self.natoms):
            self.pos[i,0] += delta[i,0]
            self.pos[i,1] += delta[i,1]
            self.pos[i,2] += delta[i,2]
        
    def read_config(self, instream):
        """
        reads in the simplified xyz file. very simplified for at the moment
        """
        restart_iteration = 0
        restart_time = 0.0
        restart_energy = 0.0

        line = instream.readline()
        words = line.split()
        self.natoms = int(words[0])

        self.pos = np.zeros((self.natoms, 3), dtype=np.float64)
        self.label = np.zeros(self.natoms, dtype=np.int32)
        self.frozen = np.zeros(self.natoms, dtype=np.int32)

        #read in the information line containing vectors etc
        info = instream.readline()
        info = info.lower()
        words = info.split()

        idx = 0
        for w in words:
            if "lattice=" in w:
                break
            idx +=1

        tmp = words[idx].replace('lattice=','')
        tmp = tmp.replace('"','')
        self.vectors[0,0] = float(tmp)
        self.vectors[0][1] = float(words[idx+1])
        self.vectors[0][2] = float(words[idx+2])
        self.vectors[1][0] = float(words[idx+3])
        self.vectors[1][1] = float(words[idx+4])
        self.vectors[1][2] = float(words[idx+5])
        self.vectors[2][0] = float(words[idx+6])
        self.vectors[2][1] = float(words[idx+7])
        tmp = words[idx+8].replace('"','')
        self.vectors[2][2] = float(tmp)

        for w in words:
            if "iteration" in w:
                restart_iteration = int(w.replace('iteration=',''))
                break

        for w in words:
            if "time" in w:
                restart_time = float(w.replace('time=',''))
                break

        for w in words:
            if "energy" in w:
                restart_energy = float(w.replace('energy=',''))
                break
    
        #next read in the atom pos
        for i in range(self.natoms):
            
            line = instream.readline()
            
            words = line.split()
            self.symbol.append(words[0])
            
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
            
            self.pos[i][0] = x
            self.pos[i][1] = y
            self.pos[i][2] = z

       
        return restart_iteration, restart_time, restart_energy
        
    def write_config(self, outstream, total_energy=None, iteration = None, time = None):
        """
        writes out a simplified form of the extended xyz format used by ASE
        """

        outstream.write(f"{self.natoms} \n")
        
        #prepare the info line
        info = 'Lattice="'
        tmp = self.vectors.flatten()
        for i in range(9):
            info = info + str(tmp[i]) + ' '
        info = info + '"  '

        if total_energy != None:
            info = info + "energy=" + str(total_energy) + "  "

        if iteration != None:
            info = info + "iteration=" + str(iteration) + " "

        if time != None:
            info = info + "time=" + str(time) + " "

        info = info + 'pbc="T T T"  '  #not used within this prog but may be externally so included

        outstream.write(f"{info} \n")
    
        for i in range(self.natoms):
            outstream.write(f"{self.symbol[i]}    {self.pos[i][0]}  {self.pos[i][1]}  {self.pos[i][2]}\n")


        outstream.flush()
        outstream.close()
        
    def cell_properties(self):
        """ calculates the properties of a cell and returns a vector """
        cell_prop = np.zeros(6)
    
        
        # Calculate lengths of cell vectors
        cell_prop[0] = np.linalg.norm(self.vectors[0, :])  # Length of first row vector
        cell_prop[1] = np.linalg.norm(self.vectors[1, :])  # Length of second row vector
        cell_prop[2] = np.linalg.norm(self.vectors[2, :])  # Length of third row vector

        # Calculate cosines of cell angles
        cell_prop[3] = np.dot(self.vectors[0, :], self.vectors[1, :]) / (cell_prop[0] * cell_prop[1]) #alpha
        cell_prop[4] = np.dot(self.vectors[0, :], self.vectors[2, :]) / (cell_prop[0] * cell_prop[2]) #beta
        cell_prop[5] = np.dot(self.vectors[1, :], self.vectors[2, :]) / (cell_prop[1] * cell_prop[2]) #gamma
    
        return cell_prop
