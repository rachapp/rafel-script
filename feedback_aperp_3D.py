# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:32:56 2020

@author: seb18121
"""

import numpy as np
# import matplotlib.pyplot as plt
import tables, sys, gc

# python /mnt/d/My_python_script/feedback_aperp_1D.py XXX_aperp_NO feedbackfactor detuning 
basename = sys.argv[1] # retreive the base name
# filename1 = 'CEP_aperp_41'
filename = basename + ".h5"
outfilename = "entrance.h5"

print ("Reading the field file ...\n")

h5 = tables.open_file(filename, 'r')

wavelength = h5.root.runInfo._v_attrs.lambda_r
nx = h5.root.runInfo._v_attrs.nX
ny = h5.root.runInfo._v_attrs.nY
nz = h5.root.runInfo._v_attrs.nZ2
Lc = h5.root.runInfo._v_attrs.Lc
Lg = h5.root.runInfo._v_attrs.Lg
meshsizeX = h5.root.runInfo._v_attrs.sLengthOfElmX
meshsizeY = h5.root.runInfo._v_attrs.sLengthOfElmY
meshsizeZ2 = h5.root.runInfo._v_attrs.sLengthOfElmZ2
meshsizeXSI = meshsizeX*np.sqrt(Lc*Lg)
meshsizeYSI = meshsizeY*np.sqrt(Lc*Lg)
meshsizeZSI = meshsizeZ2*Lc

fieldin = h5.root.aperp.read()

# cavity loss = 1 - (E1/E0)**2
Ff = float(sys.argv[2]) # feedback factor => 1 - total cavity loss
# detuning in no. of wavelength 
detuning = float(sys.argv[3])

d = int(np.round(detuning*wavelength/meshsizeZSI))
print ("Reading file ... DONE \n")
print ("Appying the feedback factor and detune length ...")
print ("Total cavity loss = " + str((1.0-Ff)*100) + "%, " + "Cavity detuning = " + str(detuning) + " wavelengths\n")

fieldout = np.pad(fieldin*np.sqrt(Ff), ((0,0),(0,d),(0,0),(0,0)), 'constant')[:, d:, :, :]


del(fieldin)
h5.close()
gc.collect()

print ("Saving to new field file ...\n")
#save h5 file
hf = tables.open_file(outfilename,'w')
saperp = hf.create_array('/','aperp', fieldout)
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
print ("Return entrance.h5 ... Done\n")
sys.exit()

# plt.xlim(0,1000)
# plt.plot(Ax)