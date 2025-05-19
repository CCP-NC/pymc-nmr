
import numpy as np
import math

from typing import List
from Species import Species  
from Energy import Energy
from config import Config

"""
    
    //     scalar sum routines
    //
    //    copyright daresbury laboratory
    //    author - w.smith october 2001
    //
    //    $Author: wl $
    //    $Date: 2007/05/08 08:55:34 $
    //    $Revision: 1.3 $
    //    $State: Exp $
    //
    """
def sclsum(n, a, i):

    k = 0
    sclsum = 0.0

    for j in range(n):
        sclsum = sclsum + a[k]
        k = k + i
  
    return sclsum

def sclsum2(n, m, a, i):

    

    k = 0
    sclsum = 0.0

    for j in range(n):
        sclsum = sclsum + a[k][m]
        k = k + i
  

    return sclsum



class Statistics:
    def __init__(self):
        # Initialization of member variables
        self.m_print_unit = 1.0
        self.m_nstk = 0
        self.m_printlat = True
        self.m_sample = 0
        
        # Variables initialized to 0.0
        self.m_stptoteng = 0.0
        self.m_stppres = 0.0
        self.m_stpvir = 0.0
        self.m_stpenthalpy = 0.0
        self.m_stpcfg = 0.0
        self.m_stprecip = 0.0
        self.m_stpreal = 0.0
        self.m_stpvdw = 0.0
        self.m_stppair = 0.0
        self.m_stpthree = 0.0
        self.m_stpmany = 0.0
        
        self.m_avetoteng = 0.0
        self.m_avevir = 0.0
        self.m_aveenthalpy = 0.0
        self.m_averecip = 0.0
        self.m_avecfg = 0.0
        self.m_avereal = 0.0
        self.m_avevdw = 0.0
        self.m_avepair = 0.0
        self.m_avethree = 0.0
        self.m_avemany = 0.0
        self.m_avepres = 0.0
        
        self.m_flctoteng = 0.0
        self.m_flcvir = 0.0
        self.m_flcpres = 0.0
        self.m_flcenthalpy = 0.0
        self.m_flccfg = 0.0
        self.m_flcrecip = 0.0
        self.m_flcreal = 0.0
        self.m_flcvdw = 0.0
        self.m_flcpair = 0.0
        self.m_flcthree = 0.0
        self.m_flcmany = 0.0
        
        self.m_ravetoteng = 0.0
        self.m_raveenthalpy = 0.0
        self.m_ravecfg = 0.0
        self.m_raverecip = 0.0
        self.m_ravereal = 0.0
        self.m_ravevdw = 0.0
        self.m_ravepair = 0.0
        self.m_ravethree = 0.0
        self.m_ravemany = 0.0
        self.m_ravepres = 0.0
        self.m_ravevir = 0.0
        
        self.m_zumtoteng = 0.0
        self.m_zumenthalpy = 0.0
        self.m_zumcfg = 0.0
        self.m_zumrecip = 0.0
        self.m_zumreal = 0.0
        self.m_zumvdw = 0.0
        self.m_zumpair = 0.0
        self.m_zumthree = 0.0
        self.m_zummany = 0.0
        self.m_zumpres = 0.0
        self.m_zumvir = 0.0
        
        # Additional member variables used in the implementation
        self.m_stktoteng = None
        self.m_stkcfg = None
        self.m_stkrecip = None
        self.m_stkreal = None
        self.m_stkvdw = None
        self.m_stkpair = None
        self.m_stkthree = None
        self.m_stkmany = None
        self.m_stkenthalpy = None
        self.m_stkpres = None
        self.m_stkvir = None
        
        self.m_stp_latvec = None
        self.m_ave_latvec = None
        self.m_rave_latvec = None
        self.m_zum_latvec = None
        self.m_flc_latvec = None
        self.m_stk_latvec = None
        self.m_stp_strs = None
        self.m_ave_strs = None
        self.m_rave_strs = None
        self.m_zum_strs = None
        self.m_flc_strs = None
    
    def zero(self, num, unit, flag):
        self.m_print_unit = "INTERNALTOEV"
        self.m_nstk = num
        self.m_printlat = True
        self.m_sample = 0
        
        # Initialize arrays with zeros
        self.m_stktoteng = np.zeros(num)
        self.m_stkcfg = np.zeros(num)
        self.m_stkrecip = np.zeros(num)
        self.m_stkreal = np.zeros(num)
        self.m_stkvdw = np.zeros(num)
        self.m_stkpair = np.zeros(num)
        self.m_stkthree = np.zeros(num)
        self.m_stkmany = np.zeros(num)
        self.m_stkenthalpy = np.zeros(num)
        self.m_stkpres = np.zeros(num)
        self.m_stkvir = np.zeros(num)
        
        self.m_stp_latvec = np.zeros(9)
        self.m_ave_latvec = np.zeros(9)
        self.m_rave_latvec = np.zeros(9)
        self.m_zum_latvec = np.zeros(9)
        self.m_flc_latvec = np.zeros(9)
        self.m_stk_latvec = np.zeros(num * 9)
        self.m_stp_strs = np.zeros(9)
        self.m_ave_strs = np.zeros(9)
        self.m_rave_strs = np.zeros(9)
        self.m_zum_strs = np.zeros(9)
        self.m_flc_strs = np.zeros(9)
        
        
        # Allocate memory for arrays
        self.m_stktoteng = np.zeros(num)
        self.m_stkcfg = np.zeros(num)
        self.m_stkrecip = np.zeros(num)
        self.m_stkreal = np.zeros(num)
        self.m_stkvdw = np.zeros(num)
        self.m_stkpair = np.zeros(num)
        self.m_stkthree = np.zeros(num)
        self.m_stkmany = np.zeros(num)
        self.m_stkenthalpy = np.zeros(num)
        self.m_stkpres = np.zeros(num)
        self.m_stkvir = np.zeros(num)
        
        self.m_stp_latvec = np.zeros(9)
        self.m_ave_latvec = np.zeros(9)
        self.m_rave_latvec = np.zeros(9)
        self.m_zum_latvec = np.zeros(9)
        self.m_flc_latvec = np.zeros(9)
        self.m_stk_latvec = np.zeros(num * 9)
        self.m_stp_strs = np.zeros(9)
        self.m_ave_strs = np.zeros(9)
        self.m_rave_strs = np.zeros(9)
        self.m_zum_strs = np.zeros(9)
        self.m_flc_strs = np.zeros(9)

    def sample(self, sysequil: int, iter: int, totalEnergy: Energy, volume: np.float64, vec: np.ndarray, outStream):

        self.m_stptoteng = totalEnergy.get_total_energy()
        self.m_stpcfg = 0.0 #totalEnergy.get_total_energy() * volume;
        
        self.m_stpenthalpy = 0.0
        self.m_stppres = 0.0
        self.m_stpvir = 0.0

        for i in range(9):
            self.m_stp_strs[i] = 1.0
            self.m_stp_latvec[i] = vec[i]
            
    
        tmp = np.zeros(self.m_nstk)
        kstk=((iter - 1) % self.m_nstk)
        if iter > self.m_nstk:
    
            if kstk == 0:
        
                self.m_zumtoteng = sclsum(self.m_nstk,self.m_stktoteng,1)
                self.m_zumcfg = sclsum(self.m_nstk,self.m_stkcfg,1)
                self.m_zumrecip = sclsum(self.m_nstk,self.m_stkrecip,1)
                self.m_zumreal = sclsum(self.m_nstk,self.m_stkreal,1)
                self.m_zumvdw = sclsum(self.m_nstk,self.m_stkvdw,1)
                self.m_zumpair = sclsum(self.m_nstk,self.m_stkpair,1)
                self.m_zumthree = sclsum(self.m_nstk,self.m_stkthree,1)
                self.m_zummany = sclsum(self.m_nstk,self.m_stkmany,1)
                self.m_zumenthalpy = sclsum(self.m_nstk,self.m_stkenthalpy,1)
                self.m_zumpres = sclsum(self.m_nstk,self.m_stkpres,1)
                self.m_zumvir = sclsum(self.m_nstk,self.m_stkvir,1)

                for i in range(9):
                    for j in range(self.m_nstk):
                        tmp[j] = self.m_stk_latvec[i*self.m_nstk + j]
                
                    self.m_zum_latvec[i] = sclsum(self.m_nstk, tmp, 1)
        
            self.m_zumtoteng -= self.m_stktoteng[kstk]
            self.m_zumcfg -= self.m_stkcfg[kstk]
            self.m_zumrecip -= self.m_stkrecip[kstk]
            self.m_zumreal -= self.m_stkreal[kstk]
            self.m_zumvdw -= self.m_stkvdw[kstk]
            self.m_zumpair -= self.m_stkpair[kstk]
            self.m_zumthree -= self.m_stkthree[kstk]
            self.m_zummany -= self.m_stkmany[kstk]
            self.m_zumenthalpy -= self.m_stkenthalpy[kstk]
            self.m_zumpres -= self.m_stkpres[kstk]
            self.m_zumvir -= self.m_stkvir[kstk]
            for i in range(9):
                self.m_zum_latvec[i] -= 1.0 * self.m_stk_latvec[i*self.m_nstk + kstk]
               #m_zum_strs[i] -= 1.0 * m_stk_strs[kstk][i];

        self.m_stktoteng[kstk] = self.m_stptoteng
        self.m_stkcfg[kstk] = self.m_stpcfg
        self.m_stkrecip[kstk] = self.m_stprecip
        self.m_stkreal[kstk] = self.m_stpreal
        self.m_stkvdw[kstk] = self.m_stpvdw
        self.m_stkpair[kstk] = self.m_stppair
        self.m_stkthree[kstk] = self.m_stpthree
        self.m_stkmany[kstk] = self.m_stpmany
        self.m_stkenthalpy[kstk] = self.m_stpenthalpy
        self.m_stkpres[kstk] = self.m_stppres
        self.m_stkvir[kstk] = self.m_stpvir

        for i in range(9):
            self.m_stk_latvec[i*self.m_nstk + kstk] = self.m_stp_latvec[i];
            #m_stk_strs[kstk][i] = m_stp_strs[i];
    
        self.m_zumtoteng += self.m_stptoteng
        self.m_zumcfg += self.m_stpcfg
        self.m_zumrecip += self.m_stprecip
        self.m_zumreal += self.m_stpreal
        self.m_zumvdw += self.m_stpvdw
        self.m_zumpair += self.m_stppair
        self.m_zumthree += self.m_stpthree
        self.m_zummany += self.m_stpmany
        self.m_zumenthalpy += self.m_stpenthalpy
        self.m_zumpres += self.m_stppres
        self.m_zumvir += self.m_stpvir

        for i in range(9):
            self.m_zum_latvec[i] += self.m_stp_latvec[i]
            #m_zum_strs[i] += m_stp_strs[i];

        #calculate rolling averages

        if self.m_nstk < iter:
            self.m_zistk = self.m_nstk
        else:
            self.m_zistk = iter
    
        self.m_ravetoteng = self.m_zumtoteng / self.m_zistk
        self.m_ravecfg = self.m_zumcfg / self.m_zistk
        self.m_raverecip = self.m_zumrecip / self.m_zistk
        self.m_ravereal = self.m_zumreal / self.m_zistk
        self.m_ravevdw = self.m_zumvdw / self.m_zistk
        self.m_ravepair = self.m_zumpair / self.m_zistk
        self.m_ravethree = self.m_zumthree / self.m_zistk
        self.m_ravemany = self.m_zummany / self.m_zistk
        self.m_raveenthalpy = self.m_zumenthalpy / self.m_zistk
        self.m_ravepres = self.m_zumpres / self.m_zistk
        self.m_ravevir = self.m_zumvir / self.m_zistk

        for i in range(9):
            self.m_rave_latvec[i] = self.m_zum_latvec[i] / self.m_zistk;
            #m_rave_strs[i] = m_zum_strs[i] / m_zistk;

        # accumulate totals over steps

        self.m_sample += 1
        sclnv1 = np.float64(self.m_sample-1) / np.float64(self.m_sample)
        sclnv2 = 1.0 / np.float64(self.m_sample);
        self.m_flctoteng  = sclnv1 * (self.m_flctoteng + sclnv2 * pow((self.m_stptoteng - self.m_avetoteng ),2))
        self.m_flccfg  = sclnv1 * (self.m_flccfg + sclnv2 * pow((self.m_stpcfg - self.m_avecfg ),2))
        self.m_flcrecip  = sclnv1 * (self.m_flcrecip + sclnv2 * pow((self.m_stprecip - self.m_averecip ),2))
        self.m_flcreal  = sclnv1 * (self.m_flcreal + sclnv2 * pow((self.m_stpreal - self.m_avereal ),2))
        self.m_flcvdw  = sclnv1 * (self.m_flcvdw + sclnv2 * pow((self.m_stpvdw - self.m_avevdw ),2))
        self.m_flcpair  = sclnv1 * (self.m_flcpair + sclnv2 * pow((self.m_stppair - self.m_avepair ),2))
        self.m_flcthree  = sclnv1 * (self.m_flcthree + sclnv2 * pow((self.m_stpthree - self.m_avethree ),2))
        self.m_flcmany  = sclnv1 * (self.m_flcmany + sclnv2 * pow((self.m_stpmany - self.m_avemany ),2))
        self.m_flcenthalpy  = sclnv1 * (self.m_flcenthalpy + sclnv2 * pow((self.m_stpenthalpy - self.m_aveenthalpy ),2))
        self.m_flcpres  = sclnv1 * (self.m_flcpres + sclnv2 * pow((self.m_stppres - self.m_avepres ),2))
        self.m_flcvir  = sclnv1 * (self.m_flcvir + sclnv2 * pow((self.m_stpvir - self.m_avevir ),2))

        for i in range(9):
            self.m_flc_latvec[i]  = sclnv1 * (self.m_flc_latvec[i] + sclnv2 * pow((self.m_stp_latvec[i] - self.m_ave_latvec[i] ),2))
            #self.m_flc_strs[i]  = sclnv1 * (self.m_flc_strs[i] + sclnv2 * pow((self.m_stp_strs[i] - self.m_ave_strs[i] ),2));

        self.m_avetoteng = sclnv1 * self.m_avetoteng + sclnv2 * self.m_stptoteng
        self.m_avecfg = sclnv1 * self.m_avecfg + sclnv2 * self.m_stpcfg
        self.m_averecip = sclnv1 * self.m_averecip + sclnv2 * self.m_stprecip
        self.m_avereal = sclnv1 * self.m_avereal + sclnv2 * self.m_stpreal
        self.m_avevdw = sclnv1 * self.m_avevdw + sclnv2 * self.m_stpvdw
        self.m_avepair = sclnv1 * self.m_avepair + sclnv2 * self.m_stppair
        self.m_avethree = sclnv1 * self.m_avethree + sclnv2 * self.m_stpthree
        self.m_avemany = sclnv1 * self.m_avemany + sclnv2 * self.m_stpmany
        self.m_aveenthalpy = sclnv1 * self.m_aveenthalpy + sclnv2 * self.m_stpenthalpy
        self.m_avepres = sclnv1 * self.m_avepres + sclnv2 * self.m_stppres
        self.m_avevir = sclnv1 * self.m_avevir + sclnv2 * self.m_stpvir

        for i in range(9):
            self.m_ave_latvec[i] = sclnv1 * self.m_ave_latvec[i] + sclnv2 * self.m_stp_latvec[i]
            #self.m_ave_strs[i] = sclnv1 * self.m_ave_strs[i] + sclnv2 * self.m_stp_strs[i];


        if iter <= sysequil:
            self.m_sample = 0
            self.m_avetoteng = 0.0
            self.m_aveenthalpy = 0.0
            self.m_averecip = 0.0
            self.m_avecfg = 0.0
            self.m_avereal = 0.0
            self.m_avevdw = 0.0
            self.m_avepair = 0.0
            self.m_avethree = 0.0
            self.m_avemany = 0.0
            self.m_avepres = 0.0
            self.m_avevir = 0.0
            for i in range(9):
                self.m_ave_latvec[i] = 0.0
                #self.m_ave_strs[i] = 0.0;

            self.m_flctoteng = 0.0
            self.m_flcpres = 0.0
            self.m_flcvir = 0.0
            self.m_flcenthalpy = 0.0
            self.m_flccfg = 0.0
            self.m_flcrecip = 0.0
            self.m_flcreal = 0.0
            self.m_flcvdw = 0.0
            self.m_flcpair = 0.0
            self.m_flcthree = 0.0
            self.m_flcmany = 0.0
            for i in range(9):
                self.m_flc_latvec[i] = 0.0
                #self.m_flc_strs[i] = 0.0;
    
    

    def check_point(self, iter, equil, timeelp, outStream):
        outStream.write("\n\n *****************************************************************************************************\n")
        outStream.write("\n iteration {}\n".format(iter))
        outStream.write("\n energy                                    instantaneous          rolling average")

        outStream.write("\n\n total energy                  {:>25.15e}{:>25.15e}".format(self.m_stptoteng, self.m_ravetoteng))
        

        if self.m_printlat:
            outStream.write("\n\n {:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}".format(
                self.m_stp_latvec[0], self.m_stp_latvec[1], self.m_stp_latvec[2],
                self.m_rave_latvec[0], self.m_rave_latvec[1], self.m_rave_latvec[2]))
            outStream.write("\n {:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}".format(
                self.m_stp_latvec[3], self.m_stp_latvec[4], self.m_stp_latvec[5],
                self.m_rave_latvec[3], self.m_rave_latvec[4], self.m_rave_latvec[5]))
            outStream.write("\n {:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}".format(
                self.m_stp_latvec[6], self.m_stp_latvec[7], self.m_stp_latvec[8],
                self.m_rave_latvec[6], self.m_rave_latvec[7], self.m_rave_latvec[8]))

        if iter == equil:
            outStream.write("\n *****************************************************************************************************\n")
            outStream.write(" equilibration period ended at step {}\n".format(equil))
            outStream.write(" *****************************************************************************************************\n")

        outStream.flush()

    def last_summary(self, outStream, nstep):
        # Calculate averages and fluctuations
        outStream.write("\n *****************************************************************************************************\n")
        outStream.write("\n run closing at step {} final averages and fluctuations over {}\n".format(nstep - 1, self.m_sample))
        outStream.write(" *****************************************************************************************************\n")
        outStream.write("\n                                                averages              fluctuations")

        if self.m_flctoteng > 0.0:
            self.m_flctoteng = math.sqrt(self.m_flctoteng)
        self.m_flcvir = math.sqrt(self.m_flcvir)
        self.m_flccfg = math.sqrt(self.m_flccfg)
        self.m_flcrecip = math.sqrt(self.m_flcrecip)
        self.m_flcreal = math.sqrt(self.m_flcreal)
        self.m_flcvdw = math.sqrt(self.m_flcvdw)
        self.m_flcpair = math.sqrt(self.m_flcpair)
        self.m_flcthree = math.sqrt(self.m_flcthree)
        self.m_flcmany = math.sqrt(self.m_flcmany)

        for i in range(9):
            self.m_flc_latvec[i] = math.sqrt(self.m_flc_latvec[i])
            self.m_flc_strs[i] = math.sqrt(self.m_flc_strs[i])

        # Write out averages and fluctuations
        outStream.write("\n avg total energy                  {:>25.15e}{:>25.15e}".format(self.m_avetoteng, self.m_flctoteng))
        

        if self.m_printlat:
            outStream.write("\n\n {:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}".format(
                self.m_ave_latvec[0], self.m_ave_latvec[1], self.m_ave_latvec[2],
                self.m_flc_latvec[0], self.m_flc_latvec[1], self.m_flc_latvec[2]))
            outStream.write("\n {:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}".format(
                self.m_ave_latvec[3], self.m_ave_latvec[4], self.m_ave_latvec[5],
                self.m_flc_latvec[3], self.m_flc_latvec[4], self.m_flc_latvec[5]))
            outStream.write("\n {:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}{:>15.7e}".format(
                self.m_ave_latvec[6], self.m_ave_latvec[7], self.m_ave_latvec[8],
                self.m_flc_latvec[6], self.m_flc_latvec[7], self.m_flc_latvec[8]))

class TypeStatistics:
    def __init__(self):
        self.numSpec = 0
        self.m_nstk = 10  # Example value, replace with appropriate initialization
        self.m_sample = 1000  # Example value, replace with appropriate initialization
        self.m_zistk = 100  # Example value, replace with appropriate initialization
        self.m_stp_nspc = None
        self.m_ave_nspc = None
        self.m_rave_nspc = None
        self.m_zum_nspc = None
        self.m_flc_nspc = None
        self.m_stk_nspc = None

    def zero_types(self, num_s: int, flag: bool):
        self.numSpec = num_s

        if num_s != 0:
            self.m_stp_nspc = np.zeros(num_s)
            self.m_ave_nspc = np.zeros(num_s)
            self.m_rave_nspc = np.zeros(num_s)
            self.m_zum_nspc = np.zeros(num_s)
            self.m_flc_nspc = np.zeros(num_s)
            self.m_stk_nspc = np.zeros(self.m_nstk * num_s)

    def find_num_types(self, bas:Config, typ) -> int:
        num_typ = 0
        for i in range(bas.natoms):
            if typ == bas.symbol[i]:
                num_typ += 1

        return num_typ

    def sample_types(self, iter: int, equil: int, bas: Config, spec:Species):
        tmp = [0.0] * self.m_nstk
        sclnv1 = float(self.m_sample - 1) / float(self.m_sample)
        sclnv2 = 1.0 / float(self.m_sample)

        for i in range(self.numSpec):
            ele = spec.get_species(i)
            self.m_stp_nspc[i] = float(self.find_num_types(bas, ele.name))

        kstk = ((iter - 1) % self.m_nstk)

        if iter > self.m_nstk:
            if kstk == 0:
                for i in range(self.numSpec):
                    for j in range(self.m_nstk):
                        tmp[j] = self.m_stk_nspc[i * self.m_nstk + j]
                    self.m_zum_nspc[i] = sclsum(self.m_nstk, tmp, 1)

            for i in range(self.numSpec):
                self.m_zum_nspc[i] -= 1.0 * self.m_stk_nspc[i * self.m_nstk + kstk]

        for i in range(self.numSpec):
            self.m_stk_nspc[i * self.m_nstk + kstk] = self.m_stp_nspc[i]
            self.m_zum_nspc[i] += self.m_stp_nspc[i]

        for i in range(self.numSpec):
            self.m_rave_nspc[i] = self.m_zum_nspc[i] / self.m_zistk

        for i in range(self.numSpec):
            self.m_flc_nspc[i] = sclnv1 * (self.m_flc_nspc[i] + sclnv2 * pow((self.m_stp_nspc[i] - self.m_ave_nspc[i]), 2))

        for i in range(self.numSpec):
            self.m_ave_nspc[i] = sclnv1 * self.m_ave_nspc[i] + sclnv2 * self.m_stp_nspc[i]

        if iter <= equil:
            for i in range(self.numSpec):
                self.m_ave_nspc[i] = 0.0
                self.m_flc_nspc[i] = 0.0

    def check_point_types(self, spec: Species, out_stream):
        out_stream.write("\n\n number of atom types ")
        for i in range(self.numSpec):
            ele = spec.get_species(i)
            symbol = ele.name
            out_stream.write(f"\n {symbol} {self.m_stp_nspc[i]:.5e} {self.m_rave_nspc[i]:.5e}")
        out_stream.flush()

    def last_summary_types(self, spec: Species, out_stream):
        out_stream.write("\n\n number of atom types ")
        for i in range(self.numSpec):
            self.m_flc_nspc[i] = math.sqrt(self.m_flc_nspc[i])
            ele = spec.get_species(i)
            symbol = ele.name
            out_stream.write(f"\n {symbol} {self.m_ave_nspc[i]:.5e} {self.m_flc_nspc[i]:.5e}")
        out_stream.flush()

    def get_average_type(self, i: int) -> float:
        return self.m_ave_nspc[i]

