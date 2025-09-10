import re

class Element:
    def __init__(self, name="", mass=0.0, charge=0.0, atomic_number=0):
        self.name = name
        self.mass = mass
        self.charge = charge
        self.atomic_number = atomic_number

class Species:
    def __init__(self):
        self.number_of_elements = 0
        self.ele_data = []

    def __del__(self):
        pass

    def load_species(self, in_stream, num):
        """
        reads in the species / element data from a file. It includes charge and mass at the moment to future proof
        incase someone uses a empirical potential library
        
        Parameters
        ----------
        
        in_stream : io
            the file stream for the data

        num : int
            the number of elements for input

        """
        self.number_of_elements = num

        for _ in range(self.number_of_elements):
            line = in_stream.readline().strip()

            if line[0] == '#':
                continue
            
            words = self.split(line)

            ele = Element()
            ele.name = words[0]
            ele.mass = float(words[1])
            ele.charge = float(words[2])
            if len(words) > 3:
                ele.atomic_number = int(words[3])

            self.ele_data.append(ele)

    def read_species(self, in_stream):
        """
        reads in the species / element data from a file. It includes charge and mass at the moment to future proof
        incase someone uses a empirical potential library
        
        Parameters
        ----------
        
        in_stream : io
            the file stream for the data

        num : int
            the number of elements for input
            
        """

        for _ in range(10000):
            line = in_stream.readline().strip()
            words = self.split(line)
            if words[0].lower() == "species":
                self.number_of_elements = int(words[1])
                break

        for _ in range(self.number_of_elements):
            line = in_stream.readline().strip()
            words = self.split(line)

            ele = Element()
            ele.name = words[0]
            ele.mass = float(words[1])
            ele.charge = float(words[2])
            if len(words) > 3:
                ele.atomic_number = int(words[3])

            self.ele_data.append(ele)

    def print_species(self, out_stream):
        """
        writes out to a file the element data

        Parameters
        ----------

        out_stream : io
            the file stream for the writing of data
        """
        out_stream.write("\n \n *********************************************************************\n")
        out_stream.write("                         atom data \n")
        out_stream.write("\n \n *********************************************************************\n")
        out_stream.write("\n        element         charge           mass\n")

        for ele in self.ele_data:
            out_stream.write(f"{ele.name:15} {ele.charge:15.3f} {ele.mass:15.3f}\n")

    def get_species(self, i) -> Element:
        """
        returns a requested element

        Parameters
        ----------

        i : int
            the index of the element to be returned

        Returns
        -------
        returns an Element object

        """
        return self.ele_data[i]
    
    def get_num_species(self) -> int:
        """
        returns the number of elements

        Returns
        -------
        returns an integer of the number of elements
        """
        return len(self.ele_data)

    def get_mass(self, i) -> float:
        """
        returns the mass of a requested element

        Parameters
        ----------

        i : int
            the index of the element to get the mass

        Returns
        -------
        returns a float of the mass 

        """
        return self.ele_data[i].mass

    def split(self, s):
        """
        utility function to get strings
        """
        return re.split(r'\s+', s.strip())
