#! /usr/bin/python

from scipy import interpolate
#from read_EFIT_file_new import *
from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
import numpy as np
#from interp import *
#from read_rbsProfs import *
from read_EFIT import *
from finite_differences_x import *

#def BfieldsFS(EFIT_file_name, fs_psipn, plot_flux_surface = False):
def BfieldsFS(EFITdict, fs_psipn, ntheta, plot_flux_surface = False):

    fs_psipn = float(fs_psipn)
    ntheta = int(ntheta)

    #EFITdict = read_EFIT(EFIT_file_name)

    psirz = EFITdict['psirz']
    Zgrid = EFITdict['Zgrid']
    Rgrid = EFITdict['Rgrid']
    rleft = EFITdict['rleft']
    rdim = EFITdict['rdim']
    rmaxis = EFITdict['rmaxis']
    zmaxis = EFITdict['zmaxis']
    simag = EFITdict['simag']
    sibry = EFITdict['sibry']
    psipn = EFITdict['psipn']
    Fpol = EFITdict['Fpol']

    Z0_ind = np.argmin(abs(Zgrid-zmaxis))
    R0_ind = np.argmin(abs(Rgrid-rmaxis))
    nR = len(Rgrid)
    nZ = len(Zgrid)

    psirz_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,psirz)

    Bp_Z_grid = np.empty(np.shape(psirz))
    for i in range(nZ):
        Bp_Z_grid_this = first_derivative(psirz[i,:].flatten(),Rgrid)/Rgrid
        Bp_Z_grid[i,:] = Bp_Z_grid_this.copy()

    Bp_R_grid = np.empty(np.shape(psirz))
    for i in range(nR):
        Bp_R_grid_this = - first_derivative(psirz[:,i].flatten(),Zgrid)/Rgrid[i]
        Bp_R_grid[:,i] = Bp_R_grid_this.copy()

    Bp_R_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,Bp_R_grid)
    Bp_Z_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,Bp_Z_grid)

    Bp_tot_grid = np.sqrt(Bp_R_grid**2+Bp_Z_grid**2)

    tol = 1.0E-7
    l_1 = (rleft + rdim - rmaxis)*fs_psipn*0.99
    l_2 = (rleft + rdim - rmaxis)*fs_psipn*1.01

    theta_grid = np.empty(ntheta, dtype = float128)
    l_grid = np.empty(ntheta, dtype = float128)
    R_fs = np.empty(ntheta, dtype = float128)
    Z_fs = np.empty(ntheta, dtype = float128)
    Bp_R = np.empty(ntheta, dtype = float128)
    Bp_Z = np.empty(ntheta, dtype = float128)
    Bp = np.empty(ntheta, dtype = float128)
    kZ = np.empty(ntheta, dtype = float128)
    kR = np.empty(ntheta, dtype = float128)

    for i in range(ntheta):
        theta = 2.*np.pi*i/(ntheta-1)
        theta_grid[i] = theta
        R_1 = rmaxis + l_1 * np.cos(theta)
        Z_1 = zmaxis + l_1 * np.sin(theta)
        psi_1 = (psirz_spl(Z_1,R_1)[0][0]-simag)/(sibry-simag) - fs_psipn
        R_2 = rmaxis + l_2 * np.cos(theta)
        Z_2 = zmaxis + l_2 * np.sin(theta)
        psi_2 = (psirz_spl(Z_2,R_2)[0][0]-simag)/(sibry-simag) - fs_psipn
        istep = 0
        while abs(psi_1)>tol and abs(psi_2)>tol and istep<20:
            this_l = l_1-psi_1*(l_1-l_2)/(psi_1-psi_2)
            this_R = rmaxis + this_l*np.cos(theta)
            this_Z = zmaxis + this_l*np.sin(theta)
            this_psi = (psirz_spl(this_Z,this_R)[0][0]-simag)/(sibry-simag) - fs_psipn
            psi_2 = psi_1
            psi_1 = this_psi
            l_2 = l_1
            l_1 = this_l
            istep = istep +1
        l_grid[i] = l_1
        R_fs[i] = rmaxis + l_grid[i] * np.cos(theta)
        Z_fs[i] = zmaxis + l_grid[i] * np.sin(theta)
        Bp_R[i] = Bp_R_spl(Z_fs[i],R_fs[i])[0][0]
        Bp_Z[i] = Bp_Z_spl(Z_fs[i],R_fs[i])[0][0]
        Bp[i] = np.sqrt(Bp_R[i]**2+Bp_Z[i]**2)
        kZ[i] = -np.sqrt(Bp_Z[i]**2/(Bp_R[i]**2+Bp_Z[i]**2))
        kR[i] = np.sqrt(Bp_R[i]**2/(Bp_Z[i]**2+Bp_R[i]**2))

    F_spl = interpolate.UnivariateSpline(psipn,Fpol)
    F_fs = F_spl(fs_psipn)
    Bt = F_fs/R_fs
    B_tot = np.sqrt(Bp**2+Bt**2)
    
    if plot_flux_surface:
        plt.plot(R_fs,Z_fs)
        plt.axis('equal')
        plt.xlabel('R (m)')
        plt.ylabel('Z (m)')
        #plt.title(EFIT_file_name+', psip ='+str(fs_psipn))
        plt.title('psip ='+str(fs_psipn))
        plt.show()

    return R_fs,Z_fs,Bp,abs(Bt),B_tot

def Bfields(EFIT_file_name, fs_psipn):

    fs_psipn = float(fs_psipn)

    EFITdict = read_EFIT(EFIT_file_name)

    psirz = EFITdict['psirz']
    Zgrid = EFITdict['Zgrid']
    Rgrid = EFITdict['Rgrid']
    rleft = EFITdict['rleft']
    rdim = EFITdict['rdim']
    rmaxis = EFITdict['rmaxis']
    zmaxis = EFITdict['zmaxis']
    simag = EFITdict['simag']
    sibry = EFITdict['sibry']
    psipn = EFITdict['psipn']
    Fpol = EFITdict['Fpol']

    Z0_ind = np.argmin(abs(Zgrid-zmaxis))
    R0_ind = np.argmin(abs(Rgrid-rmaxis))
    nR = len(Rgrid)
    nZ = len(Zgrid)

    R_obmp = Rgrid[R0_ind:]
    psirz_obmp = psirz[Z0_ind,R0_ind:]
    psipn_obmp = (psirz_obmp - simag) / (sibry - simag)

    sepInd = np.argmin(abs(psipn_obmp - 1.))
    psipn_obmp = psipn_obmp[:sepInd + 1]
    R_obmp = R_obmp[:sepInd + 1]

    psirz_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,psirz)

    Bp_Z_grid = np.empty(np.shape(psirz))
    for i in range(nZ):
        Bp_Z_grid_this = first_derivative(psirz[i,:].flatten(),Rgrid)/Rgrid
        Bp_Z_grid[i,:] = Bp_Z_grid_this.copy()

    Bp_R_grid = np.empty(np.shape(psirz))
    for i in range(nR):
        Bp_R_grid_this = - first_derivative(psirz[:,i].flatten(),Zgrid)/Rgrid[i]
        Bp_R_grid[:,i] = Bp_R_grid_this.copy()

    Bp_R_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,Bp_R_grid)
    Bp_Z_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,Bp_Z_grid)

    Bp_tot_grid = np.sqrt(Bp_R_grid**2+Bp_Z_grid**2)
    Bp_obmp = Bp_tot_grid[Z0_ind,R0_ind:R0_ind + sepInd + 1]

    tol = 1.0E-7
    l_1 = (rleft + rdim - rmaxis)*fs_psipn*0.99
    l_2 = (rleft + rdim - rmaxis)*fs_psipn*1.01
    ntheta = 1000
    theta_grid = np.empty(ntheta, dtype = float128)
    l_grid = np.empty(ntheta, dtype = float128)
    R_fs = np.empty(ntheta, dtype = float128)
    Z_fs = np.empty(ntheta, dtype = float128)
    Bp_R = np.empty(ntheta, dtype = float128)
    Bp_Z = np.empty(ntheta, dtype = float128)
    Bp = np.empty(ntheta, dtype = float128)
    kZ = np.empty(ntheta, dtype = float128)
    kR = np.empty(ntheta, dtype = float128)

    for i in range(ntheta):
        theta = 2.*np.pi*i/(ntheta-1)
        theta_grid[i] = theta
        R_1 = rmaxis + l_1 * np.cos(theta)
        Z_1 = zmaxis + l_1 * np.sin(theta)
        psi_1 = (psirz_spl(Z_1,R_1)[0][0]-simag)/(sibry-simag) - fs_psipn
        R_2 = rmaxis + l_2 * np.cos(theta)
        Z_2 = zmaxis + l_2 * np.sin(theta)
        psi_2 = (psirz_spl(Z_2,R_2)[0][0]-simag)/(sibry-simag) - fs_psipn
        istep = 0
        while abs(psi_1)>tol and abs(psi_2)>tol and istep<20:
            this_l = l_1-psi_1*(l_1-l_2)/(psi_1-psi_2)
            this_R = rmaxis + this_l*np.cos(theta)
            this_Z = zmaxis + this_l*np.sin(theta)
            this_psi = (psirz_spl(this_Z,this_R)[0][0]-simag)/(sibry-simag) - fs_psipn
            psi_2 = psi_1
            psi_1 = this_psi
            l_2 = l_1
            l_1 = this_l
            istep = istep +1
        l_grid[i] = l_1
        R_fs[i] = rmaxis + l_grid[i] * np.cos(theta)
        Z_fs[i] = zmaxis + l_grid[i] * np.sin(theta)
        Bp_R[i] = Bp_R_spl(Z_fs[i],R_fs[i])[0][0]
        Bp_Z[i] = Bp_Z_spl(Z_fs[i],R_fs[i])[0][0]
        Bp[i] = np.sqrt(Bp_R[i]**2+Bp_Z[i]**2)
        kZ[i] = -np.sqrt(Bp_Z[i]**2/(Bp_R[i]**2+Bp_Z[i]**2))
        kR[i] = np.sqrt(Bp_R[i]**2/(Bp_Z[i]**2+Bp_R[i]**2))

    F_spl = interpolate.UnivariateSpline(psipn,Fpol)
    F_fs = F_spl(fs_psipn)
    Bt = F_fs/R_fs
    B_tot = np.sqrt(Bp**2+Bt**2)
    Bt_obmp = F_spl(psipn_obmp)/R_obmp
    
    if 1 == 0:
        plt.plot(R_fs,Z_fs)
        plt.axis('equal')
        plt.xlabel('R (m)')
        plt.ylabel('Z (m)')
        plt.title(EFIT_file_name+', psip ='+str(fs_psipn))
        plt.show()

    return R_fs,Z_fs,Bp,abs(Bt),B_tot,psipn_obmp,R_obmp,Bp_obmp,abs(Bt_obmp)

def flux_surface_B_fields(efit_file_name,flux_surface):
#efit_file_name = 'efit_base'
#flux_surface = '0.97'

#parser = op.OptionParser()
#options,args = parser.parse_args()
#efit_file_name = args[0]
#flux_surface = args[1]
    fs=float(flux_surface)

    psip_n, Rgrid, Zgrid, F, p, ffprime, pprime, psirz, qpsi, rmag, zmag, nw, psiax, psisep = read_EFIT_file_new(efit_file_name)

    Z0_ind = np.argmin(abs(Zgrid-zmag))
    nR = len(Rgrid)
    nZ = len(Zgrid)

#print 'psiax = ', psiax
#print 'psisep = ', psisep
#print 'len(Rgrid)', 'len(Zgrid)'
#print nR, nZ
#print 'rmag','zmag'
#print rmag,zmag
#print Zgrid[Z0_ind]
#print 'shape of psirz'
#print np.shape(psirz)

    psirz_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,psirz)

    Bp_Z_grid = np.empty(np.shape(psirz))
    for i in range(nZ):
        Bp_Z_grid_this=first_derivative(psirz[i,:].flatten(),Rgrid)/Rgrid
        Bp_Z_grid[i,:]=Bp_Z_grid_this.copy()

    Bp_R_grid = np.empty(np.shape(psirz))
    for i in range(nR):
        Bp_R_grid_this=-first_derivative(psirz[:,i].flatten(),Zgrid)/Rgrid[i]
        Bp_R_grid[:,i]=Bp_R_grid_this.copy()

    Bp_R_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,Bp_R_grid)
    Bp_Z_spl = interpolate.RectBivariateSpline(Zgrid,Rgrid,Bp_Z_grid)

    Bp_tot_grid = np.sqrt(Bp_R_grid**2+Bp_Z_grid**2)
    Bp_obmp = Bp_tot_grid[Z0_ind,:]
    #f = open('Bp_obmp.dat','w')
    #f.write('# R Bp_t\n # \n')
    #np.savetxt(f,np.column_stack((Rgrid,Bp_obmp)))
    #f.close()

    #plt.plot(Bp_obmp)
    #plt.show()

    tol = 1.0E-7         #1.0E-06
    l_1 = 0.5        #0.5#(rmin+rdim-rmag)*fs*0.99
    l_2 = 1.5        #0.1#(rmin+rdim-rmag)*fs*1.01
    ntheta = 10000
    theta_grid = np.empty(ntheta)
    l_grid = np.empty(ntheta)
    R_fs = np.empty(ntheta)
    Z_fs = np.empty(ntheta)
    Bp_R = np.empty(ntheta)
    Bp_Z = np.empty(ntheta)
    Bp = np.empty(ntheta)
    kZ = np.empty(ntheta)
    kR = np.empty(ntheta)

    for i in range(ntheta):
        theta = 2.*np.pi*i/ntheta
        theta_grid[i] = theta
        R_1 = rmag+l_1*np.cos(theta)
        Z_1 = zmag-l_1*np.sin(theta)
        psi_1 = (psirz_spl(Z_1,R_1)[0][0]-psiax)/(psisep-psiax) - fs
        R_2 = rmag+l_2*np.cos(theta)
        Z_2 = zmag-l_2*np.sin(theta)
        psi_2 = (psirz_spl(Z_2,R_2)[0][0]-psiax)/(psisep-psiax) - fs
        istep = 0
        while abs(psi_1)>tol and abs(psi_2)>tol and istep<20:
            this_l = l_1-psi_1*(l_1-l_2)/(psi_1-psi_2)
            this_R = rmag+this_l*np.cos(theta)
            this_Z = zmag-this_l*np.sin(theta)
            this_psi = (psirz_spl(this_Z,this_R)[0][0]-psiax)/(psisep-psiax) - fs
            psi_2 = psi_1
            psi_1 = this_psi
            l_2 = l_1
            l_1 = this_l
            istep = istep +1
        l_grid[i] = l_1
        R_fs[i] = rmag+l_grid[i]*np.cos(theta)
        Z_fs[i] = zmag-l_grid[i]*np.sin(theta)
        Bp_R[i] = Bp_R_spl(Z_fs[i],R_fs[i])[0][0]
        Bp_Z[i] = Bp_Z_spl(Z_fs[i],R_fs[i])[0][0]
        Bp[i] = np.sqrt(Bp_R[i]**2+Bp_Z[i]**2)
        kZ[i] = -np.sqrt(Bp_Z[i]**2/(Bp_R[i]**2+Bp_Z[i]**2))
        kR[i] = np.sqrt(Bp_R[i]**2/(Bp_Z[i]**2+Bp_R[i]**2))

#plt.plot(R_fs,Bp_R/Bp,label='theta_Z')
#plt.plot(R_fs,Bp_Z/Bp,label='theta_R')
#plt.legend()
#plt.show()
#plt.plot(R_fs,kZ,'.',label='r_Z')
#plt.plot(R_fs,kR,'.',label='r_R')
#plt.legend()
#plt.show()

#tz = kZ*Bp_R+kR*Bp_Z
#print max(tz)

    F_spl = interpolate.UnivariateSpline(psip_n,F)
    F_fs = F_spl(fs)
    #print fs,F_fs
    Bt = F_fs/R_fs
    B_tot = np.sqrt(Bp**2+Bt**2)
#plt.plot(data[:,0],data[:,2])
#plt.plot(R_fs,Bt)
#plt.show()
    plt.plot(R_fs,Z_fs)
    plt.axis('equal')
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.title(efit_file_name+', psip ='+flux_surface)
    plt.show()

#f = open('fs.dat','w')
#f.write('# theta kZ Bp_R_n\n # \n')
#np.savetxt(f,np.column_stack((theta_grid,kZ,Bp_R/Bp)))
#f.close()
#rbsProfs_file_name = 'rbsProfs'
#R_obmp, tprime_obmp, fprime_obmp, gamE_obmp, B_pol_obmp, B_tot_obmp = read_rbsProfs(rbsProfs_file_name,flux_surface,efit_file_name)
#tprime_fs = R_fs*Bp*tprime_obmp/R_obmp/B_pol_obmp
#fprime_fs = R_fs*Bp*fprime_obmp/R_obmp/B_pol_obmp
#gammaE_fs = (R_fs*Bp)**2/B_tot*gamE_obmp*B_tot_obmp/(R_obmp*B_pol_obmp)**2
#f=open('tpfp_iterbase.dat','w')
#for i in np.arange(0,len(R_fs)):
#    f.write('%12.4f%12.4f%12.4f%12.4f%12.4f\n' %(tprime_fs[i],fprime_fs[i],theta_grid[i],kZ[i],Bp_R[i]/Bp[i])) 
#f.close()
    #q_spl = interpolate.UnivariateSpline(psip_n,qpsi)
    #shat_fs = q_spl(fs,nu=1)*fs/q_spl(fs)
    #q_fs = q_spl(fs)
    #print fs,q_fs
    return R_fs,Z_fs,Bp,Bt,B_tot
