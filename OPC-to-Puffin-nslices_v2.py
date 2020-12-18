# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 06:42:31 2020

@author: seb18121
"""

import numpy as np
import sys, tables, gc

f = sys.argv[1] + ".dfl"
p = sys.argv[1] + ".param"
h5name = sys.argv[1] + ".h5"

print ("Reading parameter from .param file ..." + p + "\n")
with open(p, 'r') as pfile:
    # param = np.array(pfile.read().replace(' ', '').splitlines())
    param = pfile.read().replace(' ', '')

def substring_after(s, delim):
    return s.partition(delim)[2]

def substring_before(s, delim):
    return s.partition(delim)[0]

def search_param(param_name, param_array):
    return substring_before(substring_after(param_array,param_name + '='),'\n')

mesh_x = np.float64(search_param('mesh_x', param))
mesh_y = np.float64(search_param('mesh_y', param))
nslices = np.int64(search_param('nslices', param))
npoints_x = np.int32(search_param('npoints_x', param))
npoints_y = np.int32(search_param('npoints_y', param))
wavelength = np.float64(search_param('lambda',param))
Lc = np.float64(search_param('Lc', param))
Lg = np.float64(search_param('Lg', param))
sLengthOfElmX = np.float64(search_param('sLengthOfElmX', param))
sLengthOfElmY = np.float64(search_param('sLengthOfElmY', param))
sLengthOfElmZ2 = np.float64(search_param('sLengthOfElmZ2', param))

print ("Reading binary file ..." + f + "\n")
data = np.fromfile(f, dtype='f8')
complex_field = np.array(data[0:][::2]) + 1j*np.array(data[1:][::2])
mag = np.abs(complex_field)
phase = np.unwrap(np.angle(complex_field) + 1.25*np.pi) # phase correction while converting Puffin/OPC
complex_field = mag*np.exp(1j*(-phase))

pfile.close()

del(data)
gc.collect()

print ("Converting to Puffin format ...\n")
Aperp_x = np.real(complex_field)
Aperp_y = Aperp_x*0.0
aperp = np.concatenate((Aperp_x,Aperp_y))
aperp = np.reshape(aperp, (2, nslices, npoints_y, npoints_x))

del(Aperp_x,Aperp_y,complex_field)
gc.collect()

print ("Saving to h5 file ...\n")
#save h5 file
a = tables.Float64Atom()
shape = (2, nslices, npoints_y, npoints_x)

hf = tables.open_file(h5name,'w')
saperp = hf.create_array('/','aperp', obj = aperp)

runInfo = hf.create_group('/','runInfo','')
saperp.attrs['iCsteps'] = 0
runInfo._v_attrs.lambda_r = wavelength
runInfo._v_attrs.nX = np.int32(npoints_x)
runInfo._v_attrs.nY = np.int32(npoints_y)
runInfo._v_attrs.nZ2 = np.int32(nslices)
runInfo._v_attrs.Lc = Lc
runInfo._v_attrs.Lg = Lg
runInfo._v_attrs.sLengthOfElmX = sLengthOfElmX
runInfo._v_attrs.sLengthOfElmY = sLengthOfElmY
runInfo._v_attrs.sLengthOfElmZ2 = sLengthOfElmZ2

hf.close()
print ("Saving done ...." + h5name +"\n" )
