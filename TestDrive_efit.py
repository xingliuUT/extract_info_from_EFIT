#! /usr/bin/python

from read_EFIT import *
from calc_fields_from_EFIT import *
from calc_volume_from_EFIT import *
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
import sys

EFIT_file_name = sys.argv[1]
EFITdict = read_EFIT(EFIT_file_name)
print(list(EFITdict.keys()))

ntheta = 512
if 1 == 1:
    psip_fs = 1.
    R_fs0, Z_fs0, B_pol, B_tor, B_tot = BfieldsFS(EFITdict, psip_fs, ntheta, False)
    psip_fs = 0.995
    R_fs1, Z_fs1, B_pol, B_tor, B_tot = BfieldsFS(EFITdict, psip_fs, ntheta, False)
    psip_fs = 0.997
    R_fs2, Z_fs2, B_pol, B_tor, B_tot = BfieldsFS(EFITdict, psip_fs, ntheta, False)
    if 1 == 1:
        plt.plot(R_fs0, Z_fs0, label = 'psi_pol = 1.')
        plt.plot(R_fs1, Z_fs1, label = 'psi_pol = 0.995')
        plt.plot(R_fs2, Z_fs2, label = 'psi_pol = 0.997')
        plt.axis('equal')
        plt.xlabel('R (m)')
        plt.ylabel('Z (m)')
#        plt.title(EFIT_file_name+', psip ='+str(fs_psipn))
        plt.legend()
        plt.show()

if 1 == 1:
    uni_rhot, shat, Ls = magneticShear(EFITdict)
    plt.plot(EFITdict['rhotn'], shat, label = 'shat')
    plt.plot(uni_rhot, shat)
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 0:
   plt.plot(EFITdict['psipn'], EFITdict['qpsi'], label = 'q')
   plt.xlabel('psip')
   plt.ylabel('q')
   plt.show()
if 1 == 1:
   plt.plot(EFITdict['rhotn'], EFITdict['jtor'], label = 'jtor')
   plt.xlabel('rhotn')
   plt.legend()
   plt.show()

if 1 == 0:
    V = totalV(EFIT_file_name, sys.argv[2], ntheta)
    print('# Total volume = {}'.format(V))
if 1 == 0:
    print(surfaceArea(EFITdict, ntheta, psipn_fs = 0.975))
