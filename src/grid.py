#grid with cavites
import numpy as np

from config import Config

class Grid:

    def __init__(self, numx, numy, numz, rcut = 2.0):
        
        self.numx: int = numx
        self.numy: int = numy
        self.numz: int = numz

        self.grid_pos = None
        
        self.grid_occ = None
        self.rcut = rcut

        ntot = int(self.numx * self.numy * self.numz)
        self.grid_pos = np.zeros((ntot, 3), dtype=np.float64)
        self.grid_occ = np.zeros(ntot, dtype=np.int32)


    def build_grid(self, cfg:Config):

        
        #cell is the lattice vectors
        cell = cfg.get_vectors()
        self.grid_occ[:] = 0
       
        #spacing along axis
        dx = 1.0 / self.numx
        dy = 1.0 / self.numy
        dz = 1.0 / self.numz
        #create the positions
        idx = 0
        for i in range(self.numx):
            for j in range(self.numy):
                for k in range(self.numz):
                    ax = i * dx
                    ay = j * dy
                    az = k * dz

                    #convert to cartesian positions
                    rx = ax * cfg.vectors[0,0] + ay * cfg.vectors[1,0] + az * cfg.vectors[2,0]
                    ry = ay * cfg.vectors[0,1] + ay * cfg.vectors[1,1] + az * cfg.vectors[2,1]
                    rz = az * cfg.vectors[0,2] + ay * cfg.vectors[1,2] + az * cfg.vectors[2,2]

                    #store positions
                    self.grid_pos[idx,0] = rx
                    self.grid_pos[idx,1] = ry
                    self.grid_pos[idx,2] = rz
                    
                    #print("cav ",i,j,k, idx,rx,ry,rz)

                    idx += 1

        #once the grid positions have been 
        rad = self.rcut * self.rcut

        for i in range(len(self.grid_occ)):
            overlap = cfg.check_overlap(self.grid_pos[i,0], self.grid_pos[i,1], self.grid_pos[i,2], rad)
            if overlap:
                self.grid_occ[i] = 1

            #print("overlap", self.grid_occ[i])

    def find_empty_grids(self, ax, ay, az, vectors, rcut):

        radius = rcut * rcut
        grd_list = []

        for i in range(len(self.grid_occ)):

            bx, by, bz = self.grid_pos[i,0],self.grid_pos[i,1], self.grid_pos[i,2]

            for nx in range(-1,2):
                for ny in range(-1,2):
                    for nz in range(-1,2):
                        xx = bx + nx * vectors[0,0] + ny * vectors[1,0] + nz * vectors[2,0]
                        yy = by + nx * vectors[0,1] + ny * vectors[1,1] + nz * vectors[2,1]
                        zz = bz + nx * vectors[0,2] + ny * vectors[1,2] + nz * vectors[2,2]

                        rx, ry, rz = ax - xx, ay - yy, az - zz
                        rsq = rx * rx + ry * ry + rz * rz
                        #print("check ",i,rx,ry,rz,np.sqrt(rsq), radius)
                        if rsq <= radius and self.grid_occ[i] == 0:
                            print("empty grid ",i,rx,ry,rz,np.sqrt(rsq), bx, by, bz)
                            grd_list.append(i)

        return grd_list