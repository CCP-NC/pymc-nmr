# pymc-nmr

pymc-nmr generates multiple configurations allowing the user to weight calculated NMR energies.

Pymc-nmr has three modes of operation that are based on swap moves. A swap move is where the positions of two atoms, chosen at random, are interchanged.The difference in energy of the old (*o*) and new (*n*) cells is calculated and employed within the Metropolis-Hastings sampling. 

This first is the standard Monte Carlo method described in the books by Allen and Tildesley  and Frenkel and Smit.
DL_MONTE. Simulations can be performed in the canonical (:math:`NVT`) ensemble and isobaric-isothermal (:math:`NPT`) (where MC moves attempting variations in 
the volume of the system are applied). In addition, the swap moves described above can be employed.


The second is the basin hopping approach is useful when the difference in energy between swapped particles is large. Each swap of particles is followed by an energy minimisation. The 
difference in energy is taken after the minimisation and the Metropolis selection rule applied. Although the calculation of the energy consumes more resources
the probability of the swap being accepted is significantly greater. 
In some cases the position of a particle is not known and we have created the possibility of using a NVT molecular dynamics calculation followed by energy minimisation
to further reduce the total energy of the cell.


The last method is an ab initio random structure search method. In this approach all the particles of two specific type are swapped and an energy minimisation performed.


