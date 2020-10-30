# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:52:01 2019

@author: seb18121
"""

import struct
import numpy as np
import tables
import sys
import f90nml

paramname = sys.argv[1] + ".param"
binname = sys.argv[1] + ".dfl"
h5name = sys.argv[1] + ".h5"

#get parameters from .param file 
with open(paramname, "r") as paramfile:
    nml = f90nml.read(paramfile)
wlambda = nml['optics']['lambda']
npoints_x = nml['optics']['npoints_x']
npoints_y = nml['optics']['npoints_y']
nslices = nml['optics']['nslices']
mesh_x = nml['optics']['mesh_x']
mesh_y = nml['optics']['mesh_x']
zsep = nml['optics']['zsep']

sLengthOfElmX = nml['puffin']['sLengthOfElmX']
sLengthOfElmY = nml['puffin']['sLengthOfElmY']
sLengthOfElmZ2 = nml['puffin']['sLengthOfElmZ2']
Lc = nml['puffin']['Lc']
Lg = nml['puffin']['Lg']
#######get parameters from .param file done
print ("Reading binary file ...\n")

f = open(binname,'rb')    # open the binary file
f = f.read()        # read entire file
count = len(f)/8    # check the byte size, 8-byte floating number

i = 0
dt = np.dtype('f8')
raw_data = []   # define array for data

while (i < count):  # insert each byte into array
    dt = struct.unpack('d',f[i*8:i*8+8])[0] # unpack data every 8-byte
    raw_data.append(dt)
    i = i + 1

u = np.array(raw_data)
real = u[0:][::2]   # even index
real = np.array(real)
imag = u[1:][::2]   # odd index
imag = np.array(imag)

if npoints_x == 1:
    Aperp_x = real
else:
    Aperp_x = -np.reshape(real, (nslices,npoints_y, npoints_x))

lengthX = mesh_x*(npoints_x - 1.0)
lengthY = mesh_y*(npoints_y - 1.0)
nz = nslices
x = np.linspace(-lengthX/2.0, lengthX/2.0, npoints_x)
y = np.linspace(-lengthY/2.0, lengthY/2.0, npoints_y)
z = np.linspace(0.0,1.0,nz)
y, x, z = np.meshgrid(y, x, z)


Aperp_y = Aperp_x*0.0

aperp = np.concatenate((Aperp_x,Aperp_y))

if npoints_x == 1:
    aperp = np.reshape(aperp, (2,nz))
else:
    aperp = np.reshape(aperp, (2,nz,npoints_y,npoints_x))

print ("Saving to h5 file ...\n")
#save h5 file
hf = tables.open_file(h5name,'w')
saperp = hf.create_array('/','aperp', aperp)
runInfo = hf.create_group('/','runInfo','')
saperp.attrs['iCsteps'] = 0
runInfo._v_attrs.lambda_r = wlambda
runInfo._v_attrs.nX = np.int32(npoints_x)
runInfo._v_attrs.nY = np.int32(npoints_y)
runInfo._v_attrs.nZ2 = np.int32(nslices)
runInfo._v_attrs.Lc = Lc
runInfo._v_attrs.Lg = Lg
runInfo._v_attrs.sLengthOfElmX = sLengthOfElmX
runInfo._v_attrs.sLengthOfElmY = sLengthOfElmY
runInfo._v_attrs.sLengthOfElmZ2 = sLengthOfElmZ2

hf.close()

print ("DONE\n")