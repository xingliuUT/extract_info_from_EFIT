#! /usr/bin/python

from read_EFIT import *
#from read_iterdb_file import *
from calc_fields_from_EFIT import *
from calc_volume_from_EFIT import *
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
import sys

EFIT_file_name = sys.argv[1]
EFITdict = read_EFIT(EFIT_file_name)
print list(EFITdict.keys())
#psip_fs = 1.
ntheta = 1024
#R_fs, Z_fs, B_pol, B_tor, B_tot = BfieldsFS(EFITdict, psip_fs, ntheta, True)
#V = totalV(EFIT_file_name, sys.argv[2], ntheta)
#print('# Total volume = {}'.format(V))
print(surfaceArea(EFITdict, ntheta))
