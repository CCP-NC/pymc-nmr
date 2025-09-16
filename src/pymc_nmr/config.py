"""
The positions and vectors are used in this module. It is similar to ASE atoms object but much more simple, but it allows the system to have fictional atoms
which must be at the end of the file. Also it will possibly allow an interface with LAMMPS at a later date.
"""
import os

import numpy as np
from scipy import linalg
from ase import Atoms

from species import Species, Element

class Config (object):
    """Config class stores attributes relating to the crystal structure

    In addition to storing the crystal structure the class should perform all manipulations
    on the structure such as swapping positions, converting to and from ase Atoms object.
    """

    def __init__(self):
        
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
        """
        Resets all the parameters within a Config object
        """

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
        """
        create an atoms object removing all ghost atoms and returns the object. Any ghost atoms are filtered out
        as ASE does not like them.
        
        Returns
        -------

        Atoms object for use with ASE calculators
        """

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
        """
        Receives an atoms object and takes positions and vectors to upadte those of the
           config object. It assumes the ordering of chemical symbols has not been changed
           and allows for ghost atoms to be kept at the end.

        Parameters
        ----------

        basin : Atoms object
            new set of positions
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
    def setup_configuration(self, spec: Species, out_io):
        """
        Takes a list of the element types and puts the mass and charge on each atom
        The function also checks that all the species defined in the basin.xyz are present in the
        element list.

        Parameters
        ----------

        spec : Species object
            container for elemenal informations

        out_io : IO
            output stream
        """

        self.charge = np.zeros(self.natoms, dtype=np.float64)
        self.mass = np.zeros(self.natoms, dtype=np.float64)
        self.label = np.zeros(self.natoms, dtype=np.int32)

        for i in range(self.natoms):

            found = False

            for k in range(spec.get_num_species()):
                element = spec.get_species(k)
                if self.symbol[i] == element.name:
                    self.label[i] = k 
                    self.charge[i] = element.charge
                    self.mass[i] = element.mass
                    found = True

            if "ghost" in self.symbol[i]:
                self.numghost += 1
            elif self.numghost > 1:
                print("the ghost atoms must come last!")
                exit()

            if found == False:
                out_io.write(f"\n the atom type {self.symbol[i]} was not found in the potential file")
                out_io.flush()
                exit()
                
    
    def get_number_of_atoms(self): #returns the number of atoms in a list to controling program
        """ 
        returns the number of atoms 
        
        Returns
        -------
        
        natoms : int
            The number of atoms
        """
        return self.natoms 
        
    def get_volume(self):
        """ 
        calculates volume using numpy 
        
        Returns
        -------
        
        volume : float
            The cell volume
        """
        axb = np.cross(self.vectors[0][:], self.vectors[1][:])
        self.volume = np.dot(axb, self.vectors[2][:])

        return self.volume
        
                  
    #create a copy of the config
    def copy_config(self):
        """ 
        creates a copy of the Config object. The function creates a new object 
            copies arrays using numpy.copyto and then returns the object

        Returns
        -------

        c : Config object
            a new configuration
        """

        c = Config()
        
        c.natoms = self.natoms
        c.numghost = self.numghost

        c.symbol = []
        for i in range(c.natoms):
            c.symbol.append(self.symbol[i])

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
        """ 
        chhoses an atom at random 
        
        Returns
        -------
        
        choice : int
            The chosen atom
        """
        
        choice = -1
            
        choice = int(self.natoms * np.random.random())

        
        return choice     
        
    def select_atom_of_type(self, typ: str) -> int:
        """ 
        selects an atom of given type at random and reurns its index number

        Parameters
        ----------

        typ : str
            The atom type that is bing sought after
        
        Returns
        -------
        
        atm : int
            The chosen atom index
        """

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
        """ 
        creates a rondom distribution of two types of atoms (e.g. Si/Al)

        Parameters
        ----------

        typ1 : str
            The first atom type
        typ2 : str
            The second atoms type to be randomised with the first
        
        """
    
        swap_list = []
        shuffle_list = []

        for i in range(self.natoms):
            if self.symbol[i] == typ1 or self.symbol[i] == typ2:
                swap_list.append(i)
                shuffle_list.append(i)

        for i in range(self.natoms):
            if self.symbol[i] == "K" or self.symbol[i] == "ghost":
                print(i, self.pos[i,:])

        print("unshuffled ", shuffle_list)

        np.random.shuffle(shuffle_list)
        print("shuffled ", shuffle_list)
        for i in range(len(swap_list)):
            na = swap_list[i]
            nb = shuffle_list[i]
            self.pos[na,:] = positions[nb,:]

        for i in range(self.natoms):
            if self.symbol[i] == "K" or self.symbol[i] == "ghost":
                print(i, self.pos[i,:])

    def find_num_types(self, typ: str) -> int:
        """ 
        determines the number of a given type of atom and returns it as an integer

        Parameters
        ----------

        typ : str
            The atom type that is being counted
        
        Returns
        -------
        
        num_typ : int
            The number of atoms of a given type

        """
        num_typ = 0
        for i in range(self.natoms):
            if typ == self.symbol[i]:
                num_typ += 1

        return num_typ
    

    def swap_atom_positions(self, atm1: int, atm2: int):
        """ swaps the positions of two atoms 

        Parameters
        ----------

        atm1 : int
            The index of the first atom
        atm2 : int
            The index of the second atom
        
        """

        for j in range(3):
            tmp = self.pos[atm1,j]
            self.pos[atm1,j] = self.pos[atm2,j]
            self.pos[atm2,j] = tmp

    
    def swap_atom_types(self, atm1: int, atm2: int):
        """ 
        swaps the atom types - note does not work if ghosts are present 
        
        
        Parameters
        ----------

        atm1 : int
            The index of the first atom
        atm2 : int
            The index of the second atom
            
        """

        tmp = self.symbol[atm1]
        self.symbol[atm1] = self.symbol[atm2]
        self.symbol[atm2] = tmp

    def mutate_atom(self, atm: int, typ: int, spec: Species):

        """ 
        turns an atom into a different type
        
        Parameters
        ----------
        atm1 : int
            The index of the atom to be altered

        typ : int
            The index of the atom type in the element list

        spec : Species object
            Contains all the parameters for the new species

        """
        ele = spec.get_species(typ)
        self.symbol[atm] = ele.name
        self.mass[atm] = ele.mass
        self.charge[atm] = ele.charge
        self.label[atm] = typ         
        
    def displace_atoms(self, delta):
        """
        displaces all the atoms by a specified amount
        
        Parameters
        ----------
        
        delta : vector(float) 
            displacement vector of size 3
        """
        for i in range(self.natoms):
            self.pos[i,0] += delta[i,0]
            self.pos[i,1] += delta[i,1]
            self.pos[i,2] += delta[i,2]

    def make_atom_move(self, atm: int, dist_max: float):
        """
        displaces a single atom by a random amount upto a specified distance
        
        Parameters
        ----------
        
        atm : int
            the index of the atom to be removed
            
        dist_max : float
            maximum distance for translation
            
        Returns
        -------
        
        old_pos : vector(float)
            the original position of the atom and can be used to reset the atom
            
        """

        old_pos = np.zeros(3, dtype=np.float64)
        for j in range(3):
            old_pos[j] = self.pos[atm,j]

        #r = np.zeros(3, dtype=np.float64)
        self.pos[atm,0] += (np.random.random() - 0.5) * dist_max
        self.pos[atm,1] += (np.random.random() - 0.5) * dist_max
        self.pos[atm,2] += (np.random.random() - 0.5) * dist_max

        return old_pos
    
    def reject_atom_move(self, atm: int, old_pos: np.ndarray):
        """
        Reset the atom positions after a displacement

        Parameters
        ----------
        
        atm : int
            the index of the atom to be removed
            
        old_pos : vector(float)
            the positions used to reset the atom

        """
        for j in range(3):
            self.pos[atm][j] = old_pos[j]

    def read_config(self, instream):
        """
        reads in the simplified xyz file. very simplified for at the moment

        Parameters
        ----------
        
        instream : IO
            the input stream for the file

        Returns
        -------

        restart_iteration : int
            the current simulation iteration used to continue a calculations

        restart_time : float
            the time from the internal clock 

        restart_energy : float
            the current energy of the configuration

        """

        restart_iteration = 1
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
    
    def expand_cell_cubic(self, bulks, max_vol_change):
        """
        Isotropic expansion of the unit cell in MC simulations

        Parameters
        ----------

        bulks : vector(float)
            the scaling parameter for each direction in x,y,z
        
        max_vol_change: float
            the maximum amount the unit cell vectors can be altered

        Returns
        -------
    
        volume : float
            the new volume of the unit cell
        """
        r = np.random.random()
        
        scale = 1.0 + (r - 0.5) * max_vol_change

        cell = self.vectors
        cell[0][0] *= scale
        cell[1][1] *= scale
        cell[2][2] *= scale

        volume = self.get_volume()
        #print("new volume ", scale, volume, max_vol_change)
 
        bulks[0] = scale
        bulks[1] = scale
        bulks[2] = scale

        self.scale_positions(bulks)

        return volume
   
    def expand_cell_tetragonal(self, indx, bulks, max_vol_change):
        """
        Expansion of the unit cell in MC simulations with tetragonal symmetry ie a == b != c

        Parameters
        ----------

        indx : int
            the index of the cell that is going to be placed

        bulks : vector(float)
            the scaling parameter for each direction in x,y,z
        
        max_vol_change: float
            the maximum amount the unit cell vectors can be altered

        Returns
        -------
    
        volume : float
            the new volume of the unit cell
        """
        r = np.random.random()
        #bulks = np.ones(3, dtype=np.float64)
        cell = self.vectors
        
        scale = 1.0 + (r - 0.5) * max_vol_change[indx]

        if indx == 0:
            cell[0][0] *= scale
            cell[1][1] *= scale
            bulks[0] = scale
            bulks[1] = scale
            bulks[2] = 1.0
        else:
            cell[2][2] *= scale
            bulks[0] = 1.0
            bulks[1] = 1.0
            bulks[2] = scale

        volume = self.get_volume()
        
        self.scale_positions(bulks)

        return volume
       
    def expand_cell_orthorhombic(self, indx, bulks, max_vol_change):
        """
        Expansion of the unit cell in MC simulations with orthorohmbic symmetry ie a != b != c

        Parameters
        ----------

        indx : int
            the index of the cell that is going to be placed

        bulks : vector(float)
            the scaling parameter for each direction in x,y,z
        
        max_vol_change: float
            the maximum amount the unit cell vectors can be altered

        Returns
        -------
    
        volume : float
            the new volume of the unit cell
        """
        r = np.random.random()
        cell = self.vectors

        scale = 1.0 + (r - 0.5) * max_vol_change[indx]

        if indx == 0:
            cell[0][0] *= scale
            bulks[0] = scale
            bulks[1] = 1.0
            bulks[2] = 1.0
        elif indx == 1:
            cell[1][1] *= scale
            bulks[0] = 1.0
            bulks[1] = scale
            bulks[2] = 1.0
        else:
            cell[2][2] *= scale
            bulks[0] = 1.0
            bulks[1] = 1.0
            bulks[2] = scale

        volume = self.get_volume()
        
        self.scale_positions(bulks)

        return volume
    
    def distort_cell(self, indx, bulks, max_vol_change):
        """
        Expansion of the unit cell in MC simulations with no assumed symmetry ie a != b != c

        Parameters
        ----------

        indx : int
            the index of the cell that is going to be placed

        bulks : vector(float)
            the scaling parameter for each direction in x,y,z
        
        max_vol_change: float
            the maximum amount the unit cell vectors can be altered

        Returns
        -------
    
        volume : float
            the new volume of the unit cell
        """
        r = np.random.random()
        
        cell = np.zeros((3,3), dtype=np.float64)
        np.copyto(cell, self.vectors)

        scale = (r - 0.5) * max_vol_change[indx]
        bulks[:] = 0.0
        bulks[indx] = scale

        for j in range(3):

            self.vectors[j,0] = (1.0 + bulks[0])*cell[j,0] + 0.5*bulks[5]*cell[j,1] + 0.5*bulks[4]*cell[j,2]
            self.vectors[j,1] = (1.0 + bulks[1])*cell[j,1] + 0.5*bulks[5]*cell[j,0] + 0.5*bulks[3]*cell[j,2]
            self.vectors[j,2] = (1.0 + bulks[2])*cell[j,2] + 0.5*bulks[4]*cell[j,0] + 0.5*bulks[3]*cell[j,1]

        volume = self.get_volume()
        
        self.scale_positions_strain(bulks)

        #print(cell)
        #print(self.pos)
        #print(bulks)
        #print(self.vectors)
        #print(self.pos)
        #print(volume)
        
        #exit()
        return volume

    def scale_positions(self, bulks):
        """
        Adjustment of the atoms in the unit cell of a MC simulations to coincide with a volume displacement.
        Assumes alpha = beta = gamma = 90.0

        Parameters
        ----------

        bulks : vector(float)
            the scaling parameter for each direction in x,y,z
        
        """

        for i in range(self.natoms):
           self.pos[i,:] *= bulks[:]

    def scale_positions_strain(self, bulks):
        """
        Adjustment of the atoms in the unit cell of a MC simulations to coincide with a volume displacement. There
        is no assumption on angles.

        Parameters
        ----------

        bulks : vector(float)
            the scaling parameter for each direction in x,y,z
        
        """
        cm = np.zeros(3,dtype=np.float64)

        for i in range(self.natoms):

            cm[:] = self.pos[i,:]
            
            self.pos[i,0] = (1.0 + bulks[0]) * cm[0] + 0.5 * bulks[5] * cm[1] + 0.5 * bulks[4]* cm[2]
            self.pos[i,1] = (1.0 + bulks[1]) * cm[1] + 0.5 * bulks[5] * cm[0] + 0.5 * bulks[3]* cm[2]
            self.pos[i,2] = (1.0 + bulks[2]) * cm[2] + 0.5 * bulks[4] * cm[0] + 0.5 * bulks[3]* cm[1]

    def get_positions(self):
        """
        returns a copy of the positions. NB a new matrix is created
        
        Returns
        -------
        
        pos : np.ndarray
            positions of the atoms
            
        """
        pos = np.zeros((self.natoms,3), dtype=np.float64)
        
        np.copyto(pos, self.pos)

        return pos
    
    def set_positions(self, pos):
        """
        resets the positions. 
        
        Parameters
        ----------
        
        pos : np.ndarray
            the new positions of the atoms to reset the config positions
            
        """
        np.copyto(self.pos, pos)

    def  get_atom_positions(self, atm):
        """
        returns a copy of the positions of an individual atom. NB a new matrix is created

        Parameters
        ----------

        atm : int
            the index of the atom to get the positions
        
        Returns
        -------
        
        pos : np.ndarray
            positions of the atoms
            
        """
        pos = np.zeros(3, dtype=np.float64)
        
        pos[:] = self.pos[atm,:]

        return pos
    
    def set_atom_positions(self, atm, pos):
        """
        set the positions of an individual atom. 

        Parameters
        ----------

        atm : int
            the index of the atom to get the positions

        pos : np.ndarray
            the new positions of the atom to reset the config positions
            
        """

        self.pos[atm,:] = pos[:] 

    def get_vectors(self):
        """
        returns a copy of the cell vectors. NB a new matrix is created
        
        Returns
        -------
        
        vec : np.ndarray
            new copy of the cell vectors
            
        """

        vec = np.zeros((3,3), dtype=np.float64)
        
        np.copyto(vec, self.vectors)

        return vec
    
    def set_vectors(self, vec):
        """
        resets the cell vectors. 
        
        Parameters
        ----------
        
        vec : np.ndarray
            the new cell vectors
            
        """
        np.copyto(self.vectors, vec)

    def restore_cell(self, basin, bulks, indx):
        cell = basin.get_cell()

        if indx == 0:
            scale = 1.0 / bulks[0]
            cell[0][0] *= scale
            bulks[0] = scale
            bulks[1] = 1.0
            bulks[2] = 1.0
        elif indx == 1:
            scale = 1.0 / bulks[1]
            cell[1][1] *= scale
            bulks[0] = 1.0
            bulks[1] = scale
            bulks[2] = 1.0
        else:
            scale = 1.0 / bulks[2]
            cell[2][2] *= scale
            bulks[0] = 1.0
            bulks[1] = 1.0
            bulks[2] = scale

        self.scale_positions(basin, bulks)
        
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
        """ 
        calculates the dimensions of a cell and angles
        
        Returns
        -------
         
        cell_prop : vector 
            vector containing cell lenghts and angles
            
        """
        cell_prop = np.zeros(6)
    
        
        # Calculate lengths of cell vectors
        cell_prop[0] = np.linalg.norm(self.vectors[0, :])  # Length of first row vector
        cell_prop[1] = np.linalg.norm(self.vectors[1, :])  # Length of second row vector
        cell_prop[2] = np.linalg.norm(self.vectors[2, :])  # Length of third row vector

        # Calculate cosines of cell angles
        fact = 1.0 / (3.14159265358979323846264338328 / 180.0)
        cell_prop[5] = fact * np.arccos(np.dot(self.vectors[0, :], self.vectors[1, :]) / (cell_prop[0] * cell_prop[1])) #gamma
        cell_prop[4] = fact * np.arccos(np.dot(self.vectors[0, :], self.vectors[2, :]) / (cell_prop[0] * cell_prop[2])) #beta
        cell_prop[3] = fact * np.arccos(np.dot(self.vectors[1, :], self.vectors[2, :]) / (cell_prop[1] * cell_prop[2])) #galpha

        return cell_prop
    
    def check_overlap(self, x , y, z, radius):
        """
        checks whether a given position overlaps with atom positions. Decided if the distance is less than a buffer zone radius
        
        Parameters
        ----------
        
        x, y, z : 3*float
            the cartesian positions of the point in space
            
        radius : float
            the minimum distance that is used to trigger an overlap
            
        Returns
        -------
        
        overlap : bool
            true if an overlap is found otherwise false is returned
            
        """

        overlap = False

        for na in range(self.natoms):
            
            ax, ay, az = self.pos[na,0],self.pos[na,1], self.pos[na,2]

            for nx in range(-1,2):
                for ny in range(-1,2):
                    for nz in range(-1,2):
                        xx = ax + nx * self.vectors[0,0] + ny * self.vectors[1,0] + nz * self.vectors[2,0]
                        yy = ay + nx * self.vectors[0,1] + ny * self.vectors[1,1] + nz * self.vectors[2,1]
                        zz = az + nx * self.vectors[0,2] + ny * self.vectors[1,2] + nz * self.vectors[2,2]

                        rx, ry, rz = x - xx, y - yy, z - zz
                        rsq = rx * rx + ry * ry + rz * rz
                        #print("check ",na,rx,ry,rz,np.sqrt(rsq), radius)
                        if rsq <= radius:
                            overlap = True

        return overlap
    
    def find_closest_atom(self, atm, typ):
        """
        Finds the closest atom of a given type to another atom

        Parameters
        ----------

        atm : int
            the index of the atom we want to find the closest atom to

        typ : str
            the chemical symbol of atoms being searched over

        Returns
        -------

        choice : int
            the index of the atom that is closest to the original atom
        """

        min_dist = 1.0e6
        choice = -1

        ax, ay, az = self.pos[atm,0], self.pos[atm,1], self.pos[atm,2]

        for i in range(self.natoms):
            if self.symbol[i] != typ or i == atm:
                continue

            bx, by, bz = self.pos[i,0],self.pos[i,1], self.pos[i,2]

            for nx in range(-1,2):
                for ny in range(-1,2):
                    for nz in range(-1,2):
                        xx = bx + nx * self.vectors[0,0] + ny * self.vectors[1,0] + nz * self.vectors[2,0]
                        yy = by + nx * self.vectors[0,1] + ny * self.vectors[1,1] + nz * self.vectors[2,1]
                        zz = bz + nx * self.vectors[0,2] + ny * self.vectors[1,2] + nz * self.vectors[2,2]

                        rx, ry, rz = ax - xx, ay - yy, az - zz
                        rsq = rx * rx + ry * ry + rz * rz
                        #print("check ",i,rx,ry,rz,np.sqrt(rsq), np.sqrt(min_dist))
                        if rsq <= min_dist:
                            min_dist = rsq
                            choice = i

        #print ("min dist", choice, min_dist)
        return choice



