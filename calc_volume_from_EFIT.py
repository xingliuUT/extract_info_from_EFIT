#! /usr/bin/python

from read_EFIT import *
from calc_fields_from_EFIT import *
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
import sys
import numpy as np

def areaTriangle(vec1, vec2):
    return abs(0.5 * np.cross(vec1, vec2))

def dVolume(R_in, Z_in, R_out, Z_out, plot_flux_surface = False):
    area = []
    vol = []
    for i in range(len(R_in) - 1):
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

def totalV(EFIT_file_name, OUT_file_name, ntheta):
    EFITdict = read_EFIT(EFIT_file_name)
    psipn_vec = EFITdict['psipn']
    f = open(OUT_file_name, 'w')
    sys.stdout = f
    dV = []
    for i in range(len(psipn_vec) - 2):
#        print(i)
        R_in, Z_in, B_pol_in, B_tor_in, B_tot_in = BfieldsFS(EFITdict, psipn_vec[i], ntheta)
        R_out, Z_out, B_pol_out, B_tor_out, B_tot_out = BfieldsFS(EFITdict, psipn_vec[i + 1], ntheta)
        this_dV = dVolume(R_in, Z_in, R_out, Z_out, False)
        dV.append(this_dV)
        print('{:.8f}    {:.8f}'.format(psipn_vec[i], this_dV))
   
    return sum(dV)
    
def surfaceArea(EFITdict, ntheta, psipn_fs = 1.):
    psipn_fs = float(psipn_fs)
    psipn_vec = EFITdict['psipn']
    psipnInd = np.argmin(abs(psipn_vec - psipn_fs))
    R_fs, Z_fs, B_pol_fs, B_tor_fs, B_tot_fs = BfieldsFS(EFITdict, psipn_fs, ntheta)
    area = []
    for i in range(len(R_fs) - 1):
        this_seg = np.sqrt( (R_fs[i] - R_fs[i + 1])**2 + (Z_fs[i] - Z_fs[i + 1])**2 )
        R_center = 0.5 * (R_fs[i] + R_fs[i + 1])
        area.append(this_seg * 2. * np.pi * R_center)
    
    return sum(area)

def totalN(EFITdict, ITERDBdict, OUT_file_name, ntheta):

    psipn_vec = EFITdict['psipn']
    rhotn_vec = EFITdict['rhotn']

    nz = interp(ITERDBdict['rhot_nz'], ITERDBdict['nz'], rhotn_vec)

    f = open(OUT_file_name, 'w')
    sys.stdout = f
    dV = []
    dn = []
    for i in range(len(psipn_vec) - 2):
        #print(i)
        R_in, Z_in, B_pol_in, B_tor_in, B_tot_in = BfieldsFS(EFITdict, psipn_vec[i], ntheta)
        R_out, Z_out, B_pol_out, B_tor_out, B_tot_out = BfieldsFS(EFITdict, psipn_vec[i + 1], ntheta)

        this_dV = dVolume(R_in, Z_in, R_out, Z_out, False)
        this_dn = this_dV * nz[i]

        dV.append(this_dV)
        dn.append(this_dn)

        print('{:.8f}    {:.8f}    {:.8f}'.format(psipn_vec[i], this_dV, this_dn))
   
    return sum(dV), sum(dn)

def totalE(EFITdict, ITERDBdict, OUT_file_name, ntheta):

    psipn_vec = EFITdict['psipn']
    rhotn_vec = EFITdict['rhotn']

    ne = interp(ITERDBdict['rhot_ne'], ITERDBdict['ne'], rhotn_vec)
    ni = interp(ITERDBdict['rhot_ni'], ITERDBdict['ni'], rhotn_vec)
    nz = interp(ITERDBdict['rhot_nz'], ITERDBdict['nz'], rhotn_vec)
    te = interp(ITERDBdict['rhot_te'], ITERDBdict['te'], rhotn_vec)
    ti = interp(ITERDBdict['rhot_ti'], ITERDBdict['ti'], rhotn_vec)

    f = open(OUT_file_name, 'w')
    sys.stdout = f
    dV = []
    de = []
    for i in range(len(psipn_vec) - 2):
        #print(i)
        R_in, Z_in, B_pol_in, B_tor_in, B_tot_in = BfieldsFS(EFITdict, psipn_vec[i], ntheta)
        R_out, Z_out, B_pol_out, B_tor_out, B_tot_out = BfieldsFS(EFITdict, psipn_vec[i + 1], ntheta)

        this_dV = dVolume(R_in, Z_in, R_out, Z_out, False)
        this_de = this_dV * ((ni[i] + nz[i]) * ti[i] + ne[i] * te[i]) * 1.6E-19

        dV.append(this_dV)
        de.append(this_de)

        print('{:.8f}    {:.8f}    {:.8f}'.format(psipn_vec[i], this_dV, this_de))
   
    return sum(dV), sum(de)

def alt_totalE(EFITdict, OUT_file_name, ntheta):

    psipn_vec = EFITdict['psipn']
    rhotn_vec = EFITdict['rhotn']
    pressure = EFITdict['Pres']

    f = open(OUT_file_name, 'w')
    sys.stdout = f
    dV = []
    de = []
    for i in range(len(psipn_vec) - 2):
        #print(i)
        R_in, Z_in, B_pol_in, B_tor_in, B_tot_in = BfieldsFS(EFITdict, psipn_vec[i], ntheta)
        R_out, Z_out, B_pol_out, B_tor_out, B_tot_out = BfieldsFS(EFITdict, psipn_vec[i + 1], ntheta)

        this_dV = dVolume(R_in, Z_in, R_out, Z_out, False)
        this_de = this_dV * pressure[i]

        dV.append(this_dV)
        de.append(this_de)

        print('{:.8f}    {:.8f}    {:.8f}'.format(psipn_vec[i], this_dV, this_de))
   
    return sum(dV), sum(de)
