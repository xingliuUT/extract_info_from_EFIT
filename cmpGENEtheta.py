from calc_fields_from_EFIT import *
from scipy import interpolate
from calc_theta import *

efit_file_name = 'g098889.04530'
psi0 = 0.9854333
pnts = np.genfromtxt('pnts_DIIID_psi09854.txt')
q_gene = 4.4762795

R = pnts[:,0]
Z = pnts[:,1]

psirz_spl0, psiax, psisep, zmag, F, psip_n, ffprime, pprime, R0 = \
psirzSpl(efit_file_name)
theta_fs = theta_grid_old(R,Z,q_gene,psirz_spl0,F,psip_n,psi0)
theta_fs = theta_fs - np.pi

if 1 == 1:
    plt.plot(theta_fs/np.pi,Z,'.',label='gene')
    plt.xlabel('theta/pi')
    plt.ylabel('Z (m)')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(R,Z,'.',label='gene')
    plt.axis('equal')
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.title(efit_file_name)
    plt.legend()
    plt.show()
