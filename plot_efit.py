#! /usr/bin/python

from read_EFIT import *
from calc_fields_from_EFIT import *
from calc_volume_from_EFIT import *
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
import sys

fontsize0 = 18
figsize0 = (9, 6)
linewidth0 = 3

matplotlib.rcParams.update({'font.size': fontsize0})

EFIT_file_name = sys.argv[1]
case_file_name = sys.argv[2]
output_format = sys.argv[3]
EFITdict = read_EFIT(EFIT_file_name)
print(list(EFITdict.keys()))

if 1 == 1:
    fig, ax = plt.subplots(figsize = figsize0)
    ax.plot(EFITdict['rhotn'], EFITdict['Pres'], label = r'$P (N/m^2)$', color = 'black', linewidth = linewidth0)
    ax.set_xlabel(r'$\rho_t$', fontsize = fontsize0)
    ax.set_xlim([0.92, 1.])
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, 1.01, 0.01))
    ax.set_ylim([0, 10000])
    ax.tick_params(axis='both',which='major',labelsize=fontsize0,length=5,width=2,direction='out')
    ax.legend()
    fig.tight_layout()
    if output_format == 'pdf':
        plt.savefig(case_file_name + '_Pvsrhot.pdf', format = 'pdf')
    else:
        plt.show()

if 1 == 0:
   plt.plot(EFITdict['rhotn'], EFITdict['qpsi'], label = 'q')
   plt.xlabel('rhot')
   plt.ylabel('q')
   plt.show()
if 1 == 0:
    uni_rhot, shat, Ls = magneticShear(EFITdict)
    plt.plot(EFITdict['rhotn'], shat, label = 'shat')
    #plt.plot(uni_rhot, shat)
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
