# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:32:56 2020

@author: seb18121
"""

import numpy as np
# import matplotlib.pyplot as plt
import tables
import gc
import sys

# python /mnt/d/My_python_script/feedback_aperp_1D.py XXX_aperp_NO feedbackfactor detuning 
filename1 = sys.argv[1] # retreive the base name
# filename1 = 'CEP_aperp_41'
h5name1 = filename1 + ".h5"

# Read the HDF5 file (Puffin_aperp file)
f1 = tables.open_file(h5name1, 'r')
aperps = f1.root.aperp.read()
Aperp_x = aperps[0] # x-polarised field
Aperp_y = aperps[1] # y-polarised field
f1.close()
del aperps
gc.collect()

# cavity loss = 1 - (E1/E0)**2
Ff = float(sys.argv[2])

Ap0_x = Aperp_x
Ap1_x = np.sqrt(Ff)*Ap0_x

Ap0_y = Aperp_y
Ap1_y = np.sqrt(Ff)*Ap0_y

del Aperp_x, Aperp_y
gc.collect()

h5f = tables.open_file(h5name1, mode='r')

wavelength = h5f.root.runInfo._v_attrs.lambda_r
nx = h5f.root.runInfo._v_attrs.nX
ny = h5f.root.runInfo._v_attrs.nY
nz = h5f.root.runInfo._v_attrs.nZ2
Lc = h5f.root.runInfo._v_attrs.Lc
Lg = h5f.root.runInfo._v_attrs.Lg
meshsizeX = h5f.root.runInfo._v_attrs.sLengthOfElmX
meshsizeY = h5f.root.runInfo._v_attrs.sLengthOfElmY
meshsizeZ2 = h5f.root.runInfo._v_attrs.sLengthOfElmZ2
meshsizeXSI = meshsizeX*np.sqrt(Lc*Lg)
meshsizeYSI = meshsizeY*np.sqrt(Lc*Lg)
meshsizeZSI = meshsizeZ2*Lc
zsep = meshsizeZSI/wavelength
h5f.close()

# detuning in no. of wavelength 
detuning = float(sys.argv[3])
d = int(np.round(detuning*wavelength/meshsizeZSI))

Ax = np.pad(Ap1_x, (0,d), 'constant')
Ax = Ax[d:]
Ay = np.pad(Ap1_y, (0,d), 'constant')
Ay = Ay[d:]

detuned_Ap = np.concatenate((Ax,Ay))
detuned_Ap = np.reshape(detuned_Ap, (2,nz))
del Ax, Ay
gc.collect()

print ("Saving to h5 file ...\n")
#save h5 file
hf = tables.open_file("entrance.h5",'w')
saperp = hf.create_array('/','aperp', detuned_Ap)
runInfo = hf.create_group('/','runInfo','')
saperp.attrs['iCsteps'] = 0
runInfo._v_attrs.lambda_r = wavelength
runInfo._v_attrs.nX = nx
runInfo._v_attrs.nY = ny
runInfo._v_attrs.nZ2 = nz
runInfo._v_attrs.Lc = Lc
runInfo._v_attrs.Lg = Lg
runInfo._v_attrs.sLengthOfElmX = meshsizeX
runInfo._v_attrs.sLengthOfElmY = meshsizeY
runInfo._v_attrs.sLengthOfElmZ2 = meshsizeZ2

hf.close()
del detuned_Ap
gc.collect

print ("Done\n")
sys.exit()

# plt.xlim(0,1000)
# plt.plot(Ax)