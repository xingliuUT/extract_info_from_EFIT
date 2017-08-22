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

ntheta = 401
zmag,R_fs00,Z_fs00,dPsidZ0,dPsidR0,Bp_fs00,Bt_fs00,q_fs00,jtor0,psip_n0 = originalSurf(efit_file_name,psi0,ntheta)
jtor_spl = interpolate.UnivariateSpline(psip_n0,jtor0)
gradPsi20 = dPsidZ0**2+dPsidR0**2
qtheta_fs00,q_new0,R_new0,Z_new0,Bp_new0,Bt_new0 = theta_grid(zmag,R_fs00,Z_fs00,Bp_fs00,Bt_fs00,np.sqrt(gradPsi20),psi0)

if 1 == 1:
    plt.plot(theta_fs/np.pi,Z,'.',label='gene')
    plt.plot(qtheta_fs00/q_new0/np.pi,Z_new0,'.',label='x')
    plt.xlabel('theta/pi')
    plt.ylabel('Z (m)')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(R,Z,'.',label='gene')
    plt.plot(R_new0,Z_new0,'.',label='x')
    plt.axis('equal')
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.title(efit_file_name)
    plt.legend()
    plt.show()
