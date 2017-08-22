#! /usr/bin/python

from scipy import interpolate
from read_EFIT_file_new import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
import numpy as np
from calc_fields_from_EFIT import *
#from interp import *


def psirzSpl(efit_file_name):
    psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file_new(efit_file_name)
    psirz_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,psirz)
    return psirz_spl,psiax,psisep,zmag,F,psip_n,ffprime,pprime,Rgrid

def theta_grid_old(R,Z,q,psirz_spl,F,psip_n,flux_surface):

    fs = float(flux_surface)
    
    psi = np.empty(len(R))
    dpsidz = np.empty(len(R))
    dpsidr = np.empty(len(R))

    for i in range(len(R)):
        psi[i] = psirz_spl(Z[i],R[i])[0][0]
        dpsidz[i] = psirz_spl(Z[i],R[i],0,1,0)[0][0]
        dpsidr[i] = psirz_spl(Z[i],R[i],0,0,1)[0][0]
    gradpsi = np.sqrt(dpsidz**2+dpsidr**2)

    F_spl = interpolate.UnivariateSpline(psip_n,F)
    F_fs = F_spl(fs)
    Bt = abs(F_fs/R)

    theta_fs = np.empty(0,dtype='float')
    theta0 = 0.
    theta_fs = np.append(theta_fs,theta0)
    for i in range(len(R)-1):
        dl = np.sqrt((R[i+1]-R[i])**2+(Z[i+1]-Z[i])**2)
        dtheta = dl*(Bt[i]/q/gradpsi[i]+Bt[i+1]/q/gradpsi[i+1])/2.
        theta_fs = np.append(theta_fs,theta_fs[i]+dtheta)
        
    return theta_fs
