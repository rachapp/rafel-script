# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:03:45 2019

@author: seb18121
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:38:37 2019

@author: seb18121
"""
# noted only x polarization of the Aperp field will be converted to OPC format
import numpy as np
import itertools
import tables, gc
from scipy.signal import hilbert
import sys

filename = sys.argv[1] # retreive the base name
h5name = filename + ".h5"
binname_x = filename + "_x.dfl"
paramname_x = filename + "_x.param"

print ("Reading aperp file ..." + h5name + "\n")
h5f = tables.open_file(h5name, mode='r')

# Read the HDF5 file (Puffin_aperp file)
aperps = h5f.root.aperp.read()
Aperp_x = np.array(aperps[0]) # x-polarised field
# Aperp_y = np.array(aperps[1]) # y-polarised field

print ("Getting file attributes ... \n")
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

print("Getting the complex envelope from x-field ..." + "\n")
Aperp_x_complex = hilbert(Aperp_x,None,0)
print("Performing Hilbert transform ... DONE \n")
envelope_x = np.abs(Aperp_x_complex) # get the envelope
print("Envelope DONE \n")
phase_x = np.unwrap(np.angle(Aperp_x_complex)) # get the phase information
print("Phase DONE \n")
del(Aperp_x)
gc.collect()

print("Re-ordering the envelope into the OPC format" + "\n")
# reorder the data to binary file for X-polarised field
bin_x = envelope_x*np.exp(1j*-(phase_x))
# bin_x = Aperp_x_complex
bin_x = np.reshape(bin_x,nx*ny*nz)

dfl_x = list(itertools.chain(*zip(np.real(bin_x), np.imag(bin_x)))) # interleave real and imaginary part

del(bin_x)
gc.collect()

print("Saving x-field to binary file ..." + "\n")
# save binary file
np.array(dfl_x).tofile(binname_x)

# write the parameter file for physical interpretation
param_x = open(paramname_x, 'w')
param_x.write(" $optics\n")
param_x.write(" nslices = " + str(nz) +"\n")
param_x.write(" zsep = " + str(zsep) +"\n")
# for 1D data or 3D
if nx-1 == 0:
    param_x.write(" mesh_x = " + str(1) +"\n")
    param_x.write(" mesh_y = " + str(1) +"\n")
else:
    param_x.write(" mesh_x = " + str(meshsizeXSI) +"\n")
    param_x.write(" mesh_y = " + str(meshsizeYSI) +"\n")

param_x.write(" npoints_x = " + str(nx) +"\n")
param_x.write(" npoints_y = " + str(ny) +"\n")
param_x.write(" Mx = " + str(1) +"\n")
param_x.write(" My = " + str(1) +"\n")
param_x.write(" lambda = " + str(wavelength) +"\n")
param_x.write(" field_next = 'none'" + "\n")
param_x.write(" /\n")
param_x.write(" $puffin\n")
param_x.write(" Lc = " + str(Lc) +"\n")
param_x.write(" Lg = " + str(Lg) +"\n")
param_x.write(" sLengthOfElmX = " + str(meshsizeX) +"\n")
param_x.write(" sLengthOfElmY = " + str(meshsizeY) +"\n")
param_x.write(" sLengthOfElmZ2 = " + str(meshsizeZ2) +"\n")
param_x.write(" /\r")
param_x.close()

del(dfl_x)
gc.collect()

print("DONE\n")