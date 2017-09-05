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
def theta_grid(zmag,R,Z,Bp,Bt,gradpsi,flux_surface):
    
    fs = float(flux_surface)

    theta_fs = np.empty(len(R),dtype='float')
    #print 'len(theta_fs)=', len(theta_fs)
    theta_fs[0] = 0.
    theta_fs[len(R)-1] = 0.
    theta_lower = np.empty(0,dtype='float')
    theta_upper = np.empty(0,dtype='float')
    R_lower = np.empty(0,dtype='float')
    Z_lower = np.empty(0,dtype='float')
    R_upper = np.empty(0,dtype='float')
    Z_upper = np.empty(0,dtype='float')
    Bt_lower = np.empty(0,dtype='float')
    Bt_upper = np.empty(0,dtype='float')
    Bp_lower = np.empty(0,dtype='float')
    Bp_upper = np.empty(0,dtype='float')
    theta_lower = np.append(theta_lower,theta_fs[0])
    R_lower = np.append(R_lower,R[0])
    Z_lower = np.append(Z_lower,Z[0])
    Bt_lower = np.append(Bt_lower,Bt[0])
    Bp_lower = np.append(Bp_lower,Bp[0])
    
    if 1 == 0:
        plt.plot(R,label='R')
        plt.plot(Z,'.',label='Z')
        plt.axhline(y=zmag,label='zmag')
        plt.legend()
        plt.show()
    # be careful which (R,Z) each theta correspond to
    for i in range(len(R)-1):
        if (Z[i]+Z[i+1])/2. <= zmag:
            dl = np.sqrt((R[i+1]-R[i])**2+(Z[i+1]-Z[i])**2)
            dtheta = dl*(Bt[i]/gradpsi[i]+Bt[i+1]/gradpsi[i+1])/2.
            theta_fs[i+1] = theta_fs[i]-np.abs(dtheta)
            theta_lower = np.append(theta_lower,theta_fs[i+1])
            R_lower = np.append(R_lower,R[i+1])
            Z_lower = np.append(Z_lower,Z[i+1])
            Bt_lower = np.append(Bt_lower,Bt[i+1])
            Bp_lower = np.append(Bp_lower,Bp[i+1])
    for i in range(len(R)-1,0,-1):
        if (Z[i]+Z[i-1])/2. > zmag:
            dl = np.sqrt((R[i-1]-R[i])**2+(Z[i-1]-Z[i])**2)
            dtheta = dl*(Bt[i]/gradpsi[i]+Bt[i-1]/gradpsi[i-1])/2.
            theta_fs[i-1] = theta_fs[i]+ np.abs(dtheta)
            theta_upper = np.append(theta_upper,theta_fs[i-1])
            R_upper = np.append(R_upper,R[i-1])
            Z_upper = np.append(Z_upper,Z[i-1])
            Bt_upper = np.append(Bt_upper,Bt[i-1])
            Bp_upper = np.append(Bp_upper,Bp[i-1])

    qtheta_new = np.concatenate((np.flipud(theta_lower),theta_upper))
    R_new = np.concatenate((np.flipud(R_lower),R_upper))
    Z_new = np.concatenate((np.flipud(Z_lower),Z_upper))
    Bt_new = np.concatenate((np.flipud(Bt_lower),Bt_upper))
    Bp_new = np.concatenate((np.flipud(Bp_lower),Bp_upper))

    q_new = (np.max(qtheta_new)-np.min(qtheta_new))/2./np.pi
    #q_new is the q that makes theta goes 2pi
    #which could be a little different than q from EFIT file

    return qtheta_new, q_new, R_new, Z_new, Bp_new, Bt_new
