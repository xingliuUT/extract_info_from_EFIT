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
print list(EFITdict.keys())

ntheta = 128
if 1 == 0:
    psip_fs = 1.
    R_fs, Z_fs, B_pol, B_tor, B_tot = BfieldsFS(EFITdict, psip_fs, ntheta, True)
if 1 == 0:
    V = totalV(EFIT_file_name, sys.argv[2], ntheta)
    print('# Total volume = {}'.format(V))
if 1 == 0:
    print(surfaceArea(EFITdict, ntheta))
if 1 == 0:
    uni_rhot, shat, Ls = magneticShear(EFITdict)
    plt.plot(EFITdict['rhotn'], shat)
    plt.plot(uni_rhot, shat)
    plt.show()
    plt.plot(uni_rhot, Ls)
    plt.show()
