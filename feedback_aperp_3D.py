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

tables.copy_file(filename, outfilename, overwrite=1)
h5 = tables.open_file(outfilename, 'r+')

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

fieldout = np.pad(fieldin*np.sqrt(Ff), ((0,0),(0,d),(0,0),(0,0)), 'constant')[:, d:, :, :]

del(fieldin)
gc.collect()

h5.create_array('/','aperpN',fieldout)
h5.copy_node_attrs('/aperp','/aperpN')
h5.remove_node('/aperp')
h5.rename_node('/aperpN','aperp')

print ("Done\n")
sys.exit()

# plt.xlim(0,1000)
# plt.plot(Ax)