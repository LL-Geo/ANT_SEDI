
_________________________________________________________________________________________________________
Folder Retreat_rate_SSB and Retreat_rate_Basement

Each sub folder contains the simulation result of a 3 km sedimentary basin or basement only with 19 ka ice sheet advancing,
1 ka glacial maximum, and ice sheet retreat at different times.

Folder name: Retreat time _ Vertical Permeability of basement, confining unit and aquifer

For example 10k_964 indicate glacier retreat in 10 ka, with the permeability of basement, confining unit and 
aquifer as 10-19m2, 10-16m2 and 10-14m2 respectively.  

cvfem_input_100ka.txt show geomechanical simulation result with 100 year time step, which is then called to simulate 
groundwater flow
* Note althought the name of this file is similar in each folder. The result is different in each folder. The retreat 
rate is consistent with the folder name.

rift2d_964.dat shows SSB parameters including its geometry

sfc_tec_964.dat shows surface water flux


_________________________________________________________________________________________________________
Folder Permeability 

Each sub folder contains a 3 km sedimentary basin with different aquifer permability


_________________________________________________________________________________________________________
Folder SSB_thickness

Each sub folder contains sedimentary basin with different thickness


_________________________________________________________________________________________________________
vz_output

show vertical water flux in each simulation

The result can be plotted by Figure_4.ipynb