# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:32:56 2020

@author: seb18121
"""

import numpy as np
# import matplotlib.pyplot as plt
import h5py
import tables
import sys

# python /mnt/d/My_python_script/feedback_aperp_1D.py XXX_aperp_NO feedbackfactor detuning 
filename1 = sys.argv[1] # retreive the base name
# filename1 = 'CEP_aperp_41'
h5name1 = filename1 + ".h5"

# Read the HDF5 file (Puffin_aperp file)
f1 = h5py.File(h5name1, 'r')
# List all groups
a_group_key1 = list(f1.keys())[0]
# Get the aperp data
data1 = list(f1[a_group_key1])
Aperp_x = data1[0] # x-polarised field
Aperp_y = data1[1] # y-polarised field
f1.close()

# cavity loss = 1 - (E1/E0)**2
Ff = float(sys.argv[2])

Ap0_x = Aperp_x
Ap1_x = np.sqrt(Ff)*Ap0_x

Ap0_y = Aperp_y
Ap1_y = np.sqrt(Ff)*Ap0_y


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

print ("Saving to h5 file ...\n")
#save h5 file
hf = h5py.File("entrance.h5",'w')
saperp = hf.create_dataset('aperp', data=detuned_Ap)
runInfo = hf.create_group('runInfo')
saperp.attrs['iCsteps'] = 0
runInfo.attrs['lambda_r'] = wavelength
runInfo.attrs['nX'] = nx
runInfo.attrs['nY'] = ny
runInfo.attrs['nZ2'] = nz
runInfo.attrs['Lc'] = Lc
runInfo.attrs['Lg'] = Lg
runInfo.attrs['sLengthOfElmX'] = meshsizeX
runInfo.attrs['sLengthOfElmY'] = meshsizeY
runInfo.attrs['sLengthOfElmZ2'] = meshsizeZ2

hf.close()
print ("Done\n")
sys.exit()

# plt.xlim(0,1000)
# plt.plot(Ax)