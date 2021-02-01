import numpy as np
import matplotlib.pyplot as pl
import sys
def dtrap_poy(fname_in,nx,ny):
    poy_data = np.loadtxt(fname_in)
    x_export2D = poy_data[:,2].reshape((nx,ny))
    y_export2D = poy_data[:,3].reshape((nx,ny))
    poy_y_grid_re  = poy_data[:,10].reshape((nx,ny))
    temp=np.trapz(poy_y_grid_re,x_export2D[:,0],axis=0)
    return np.trapz(temp,y_export2D[0,:]) #[x_export2D,y_export2D,poy_y_grid_re] #

myDir = sys.argv[1]

intpoyz_tot = abs(dtrap_poy(myDir+'/Poy_tot_gd.pos',50,50))
intpoyz_ref = abs(dtrap_poy(myDir+'/Poy_ref_gd.pos',50,50))
intpoyz_inc = abs(dtrap_poy(myDir+'/Poy_inc_gd.pos',50,50))

Rnm = np.loadtxt(myDir+'/eff_r.txt',ndmin=2)[:,1] + 1j*np.loadtxt(myDir+'/eff_r.txt',ndmin=2)[:,2]
Tnm = np.loadtxt(myDir+'/eff_t.txt',ndmin=2)[:,1] + 1j*np.loadtxt(myDir+'/eff_t.txt',ndmin=2)[:,2]
Q = [np.loadtxt(myDir+'/temp-Q_L_%g.txt'%k,ndmin=2)[:,1] for k in range(2,7)]
Q.append(np.loadtxt(myDir+'/temp-Q_scat.txt',ndmin=2)[:,1])
Q=np.array(Q)
TOT = Rnm.real.sum()+Tnm.real.sum()+Q.sum()

if myDir[6:]=='solarcell':
    print('cf pdf')
    Nmax=2
    tab_lambdas=np.loadtxt(myDir+'/temp_lambda_step.txt',ndmin=2)[:,8]
    nb_lambdas=len(tab_lambdas)
    Rnm = Rnm.reshape((nb_lambdas,2*Nmax+1,2*Nmax+1))
    Tnm = Tnm.reshape((nb_lambdas,2*Nmax+1,2*Nmax+1))
    Rtot = [Rnm[i].real.sum() for i in range(nb_lambdas)]
    Ttot = [Tnm[i].real.sum() for i in range(nb_lambdas)]
    Abs_rods = Q[-1]
    Abs_ITO  = Q[0]
    Abs_subs  = Q[2]+Q[3]+Q[4]+Ttot
    TOT = Rtot+Abs_rods+Abs_ITO+Abs_subs
    pl.figure()
    pl.plot(tab_lambdas,Abs_ITO,label='absorption ITO electrode')
    pl.plot(tab_lambdas,Abs_rods,label='absorption in Si rods')
    pl.plot(tab_lambdas,Abs_subs,label='absorption in Si subs')
    pl.plot(tab_lambdas,Rtot,label='reflection')
    pl.plot(tab_lambdas,TOT,label='total')
    pl.legend()
    pl.xlabel(r'$\lambda$ [nm]')
    pl.ylabel('fraction of incident energy')
    pl.savefig('fig_solar_balance.pdf')
elif myDir[6:]=='conv':
    print('cf pdf')
    pl.figure()
    data_abs = np.loadtxt(myDir+'/temp-Q_scat.txt',ndmin=2)[:,1]
    pl.plot(data_abs[0:int(len(data_abs)/2)],label='linear elements')
    pl.plot(data_abs[int(len(data_abs)/2):] ,label='curved elements')
    pl.legend()
    pl.savefig('fig_convergence_absorption.pdf')
else:
    print('===> Computed from diffraction efficiencies')
    print('Rtot ',Rnm.real.sum())
    print('Ttot ',Tnm.real.sum())
    print('Atot ',Q.sum())
    print('TOT  ',TOT)
    print('===> Computed from trapz on Poynting vector')
    print('(expected to be less precise)')
    print('Rtot2',intpoyz_ref/intpoyz_inc)
    print('Ttot2',intpoyz_tot/intpoyz_inc)
    print('Atot ',Q.sum())
    print('TOT2 ',Q.sum()+(intpoyz_tot+intpoyz_ref)/intpoyz_inc)
