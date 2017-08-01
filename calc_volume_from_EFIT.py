#! /usr/bin/python

from read_EFIT import *
#from read_iterdb_file import *
from calc_fields_from_EFIT import *
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
import sys
import numpy as np

def areaTriangle(vec1, vec2):
    return abs(0.5 * np.cross(vec1, vec2))

def dVolume(R_in, Z_in, R_out, Z_out, plot_flux_surface = False):
#    area = 0.
#    vol = 0
    area = []
    vol = []
    for i in range(len(R_in) - 1):
#        print(i)
        vec1 = [R_out[i] - R_in[i], Z_out[i] - Z_in[i]]
        vec2 = [R_out[i] - R_out[i + 1], Z_out[i] - Z_out[i + 1]]
        vec3 = [R_in[i + 1] - R_out[i + 1], Z_in[i + 1] - Z_out[i + 1]]
        vec4 = [R_in[i + 1] - R_in[i], Z_in[i + 1] - Z_in[i]]
        this_area = areaTriangle(vec1, vec4) + areaTriangle(vec2, vec3)
        area.append(this_area)
        R_center = 0.25 * (R_in[i] + R_in[i + 1] + R_out[i] + R_out[i + 1])
        this_vol = this_area * 2. * np.pi * R_center
        vol.append(this_vol)
        if 1 == 0:
            plt.scatter(vec1[0], vec1[1], label = 'vec1', s = 1)
            plt.scatter(vec2[0], vec2[1], label = 'vec2', s = 1)
            plt.scatter(vec3[0], vec3[1], label = 'vec3', s = 1)
            plt.scatter(vec4[0], vec4[1], label = 'vec4', s = 1)
            plt.legend()
            plt.show()
    if 1 == 0:
        plt.plot(area)
        plt.show()
        plt.plot(vol)
        plt.show()
      
    if plot_flux_surface:
        plt.plot(R_in, Z_in)
        plt.plot(R_out, Z_out)
        plt.axis('equal')
        plt.xlabel('R (m)')
        plt.ylabel('Z (m)')
#        plt.title(EFIT_file_name+', psip ='+str(fs_psipn))
        plt.show()
    return sum(vol)

def totalV(EFIT_file_name, ntheta):
    EFITdict = read_EFIT(EFIT_file_name)
    dV = []
    psipn = []
    for i in range(len(EFITdict['psipn']) - 1):
        print(i)
        R_in, Z_in, B_pol_in, B_tor_in, B_tot_in = BfieldsFS(EFITdict, EFITdict['psipn'][i], ntheta)
        R_out, Z_out, B_pol_out, B_tor_out, B_tot_out = BfieldsFS(EFITdict, EFITdict['psipn'][i + 1], ntheta)
        dV.append(dVolume(R_in, Z_in, R_out, Z_out, False))
        psipn.append(EFITdict['psipn'][i])
    return sum(dV)
    

