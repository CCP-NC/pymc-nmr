import re

class Element:
    def __init__(self, name=""):
        self.name = name

class Species:
    def __init__(self):
        self.number_of_elements = 0
        self.ele_data = []

    def __del__(self):
        pass

    def load_species(self, in_stream, num):
        """
        reads in the species / element data from a file.
        """
        self.number_of_elements = num

        for _ in range(self.number_of_elements):
            line = in_stream.readline().strip()

            if line[0] == '#':
                continue
            
            words = self.split(line)

            ele = Element(name=words[0])
            self.ele_data.append(ele)

    def read_species(self, in_stream):
        """
        reads in the species / element data from a file.
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

            ele = Element(name=words[0])
            self.ele_data.append(ele)

    def print_species(self, out_stream):
        """
        writes out to a file the element data
        """
        out_stream.write("\n \n *********************************************************************\n")
        out_stream.write("                         atom data \n")
        out_stream.write("\n \n *********************************************************************\n")
        out_stream.write("\n        element\n")

        for ele in self.ele_data:
            out_stream.write(f"{ele.name:15}\n")

    def get_species(self, i) -> Element:
        return self.ele_data[i]
    
    def get_num_species(self) -> int:
        return len(self.ele_data)

    def split(self, s):
        """
        utility function to get strings
        """
        return re.split(r'\s+', s.strip())
