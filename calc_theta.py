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

def originalSurf(efit_file_name,flux_surface,ntheta):

    fs=float(flux_surface)

    psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file_new(efit_file_name)
    rdim,zdim,rctr,rmin,zmid,Bctr,curr,nh = read_EFIT_parameters(efit_file_name)
    print 'rdim = ', rdim
    print 'rctr = ', rctr
    print 'rmag = ', rmag

    psirz = (psirz - psiax)/(psisep-psiax)
    psirz_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,psirz)

    tol = 1.0E-7
    l_1 = (rmin+rdim-rmag)*fs*0.99
    l_2 = (rmin+rdim-rmag)*fs*1.01
    theta_grid = np.empty(ntheta)
    l_grid = np.empty(ntheta)
    R_fs = np.empty(ntheta)
    Z_fs = np.empty(ntheta)

    for i in range(ntheta):
        theta = 2.*np.pi*i/(ntheta-1)
        theta_grid[i] = theta
        R_1 = rmag+l_1*np.cos(theta)
        Z_1 = zmag-l_1*np.sin(theta)
        psi_1 = psirz_spl(Z_1,R_1)[0][0] - fs
        R_2 = rmag+l_2*np.cos(theta)
        Z_2 = zmag-l_2*np.sin(theta)
        psi_2 = psirz_spl(Z_2,R_2)[0][0] - fs
        istep = 0
        while abs(psi_1)>tol and abs(psi_2)>tol and istep<20:
            this_l = l_1-psi_1*(l_1-l_2)/(psi_1-psi_2)
	    this_R = rmag+this_l*np.cos(theta)
	    this_Z = zmag-this_l*np.sin(theta)
	    this_psi = psirz_spl(this_Z,this_R)[0][0] - fs
	    psi_2 = psi_1
	    psi_1 = this_psi
	    l_2 = l_1
	    l_1 = this_l
	    istep = istep +1
        l_grid[i] = l_1
        R_fs[i] = rmag+l_grid[i]*np.cos(theta)
        Z_fs[i] = zmag-l_grid[i]*np.sin(theta)

    dpsidz = np.empty(len(R_fs))
    dpsidr = np.empty(len(R_fs))

    for i in range(len(R_fs)):
        dpsidz[i] = psirz_spl(Z_fs[i],R_fs[i],0,1,0)[0][0]
        dpsidr[i] = psirz_spl(Z_fs[i],R_fs[i],0,0,1)[0][0]

    F_spl = interpolate.UnivariateSpline(psip_n,F)
    F_fs = F_spl(fs)
    Bt = F_fs/R_fs

    MU = 4.E-7

    jtor = - Rgrid*pprime - ffprime/Rgrid
    jtor = jtor*np.pi*MU

    Bp = np.sqrt(dpsidz**2+dpsidr**2)/R_fs

    return zmag,R_fs,Z_fs,dpsidz,dpsidr,Bp,Bt,qpsi,jtor,psip_n
