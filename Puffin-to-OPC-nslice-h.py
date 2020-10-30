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
import tables
from array import array
from scipy.signal import hilbert
import sys
# import matplotlib.pyplot as plt

# python /mnt/d/My_python_script/Puffin-to-OPC.py Puffin_aperp_file model_length_in_mm 
filename = sys.argv[1] # retreive the base name
h5name = filename + ".h5"
binname_x = filename + "_x.dfl"
paramname_x = filename + "_x.param"
# paramname_y = filename + "_y.param"
# binname_y = filename + "_y.dfl"

h5f = tables.open_file(h5name, mode='r')

# Read the HDF5 file (Puffin_aperp file)
aperps = h5f.root.aperp.read()
Aperp_x = aperps[0] # x-polarised field
Aperp_y = aperps[1] # y-polarised field

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

# define the grid size
lengthX = (nx-1)*meshsizeXSI
lengthY = (nx-1)*meshsizeXSI

x = np.linspace(-lengthX/2.0, lengthX/2.0, nx)
y = np.linspace(-lengthY/2.0, lengthY/2.0, ny)
z = np.linspace(0.0,1.0,nz)
y, x, z = np.meshgrid(y, x, z)

print("Getting the complex envelope from x-field ..." + "\n")
Aperp_x_complex = hilbert(Aperp_x,None,0)
envelope_x = np.abs(Aperp_x_complex) # get the envelope
phase_x = np.angle(Aperp_x_complex) # get the phase information

"""
print("Getting the complex envelope from y-field ..." + "\n")
# Get the envelope for y-polarised field
Aperp_y_complex = hilbert(Aperp_y,None,0)
envelope_y = np.abs(Aperp_y_complex)
phase_y = np.angle(Aperp_y_complex)
"""

print("Saving x-field to binary file ..." + "\n")
# reorder the data to binary file for X-polarised field
bin_x_n = envelope_x*np.exp(1j*-(phase_x - np.pi))
real_x = np.real(bin_x_n) # take real part to variable
real_x = np.reshape(real_x,nx*ny*nz) # reshape array to 1D
imag_x = np.imag(bin_x_n) # take imaginary part to variable
imag_x = np.reshape(imag_x,nx*ny*nz) # reshape array to 1D
dfl_x = list(itertools.chain(*zip(real_x, imag_x))) # interleave real and imaginary part
dfl_x = np.array(dfl_x)
# save binary file
output_file_x = open(binname_x, 'wb')
float_array_x = array('d', dfl_x)
float_array_x.tofile(output_file_x)
output_file_x.close()
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


"""
print("Saving y-field to binary file ..." + "\n")
# reorder the data to binary file for Y-polarised field
bin_y_n = envelope_y*np.exp(1j*-(phase_y - np.pi))
real_y = np.real(bin_y_n)
real_y = np.reshape(real_y,nx*ny*nz)
imag_y = np.imag(bin_y_n)
imag_y = np.reshape(imag_y,nx*ny*nz)
dfl_y = list(itertools.chain(*zip(real_y, imag_y)))
dfl_y = np.array(dfl_y)
# save binary file
output_file_y = open(binname_y, 'wb')
float_array_y = array('d', dfl_y)
float_array_y.tofile(output_file_y)
output_file_y.close()
# write the parameter file for physical interpretation
param_y = open(paramname_y, 'w')
param_y.write(" $optics\n")
param_y.write(" nslices = " + str(nz) +"\n")
param_y.write(" zsep = str(zsep)\n")
param_y.write(" mesh_x = " + str(lengthX/1000.0/(nx-1)) +"\n")
param_y.write(" mesh_y = " + str(lengthY/1000.0/(ny-1)) +"\n")
param_y.write(" npoints_x = " + str(nx) +"\n")
param_y.write(" npoints_y = " + str(ny) +"\n")
param_y.write(" lambda = " + str(wavelength) +"\n")
param_y.write(" field_next = 'none'" + "\n")
param_y.write(" $end\r")
param_y.close()
"""

f.close()

print("DONE\n")