#! /usr/bin/python

from read_EFIT import *
#from read_iterdb_file import *
from calc_fields_from_EFIT import *
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
import sys

EFIT_file_name = sys.argv[1]
EFITdict = read_EFIT(EFIT_file_name)
#print list(EFITdict.keys())
plt.plot(EFITdict['psipn'], '.')
plt.ylabel('psip_n')
plt.show()
psip_fs = EFITdict['psipn'][1]
R_fs, Z_fs, B_pol, B_tor, B_tot = BfieldsFS(EFIT_file_name, psip_fs, True)
