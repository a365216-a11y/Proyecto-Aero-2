import numpy as np
from teoria import *
D='/mnt/user-data/uploads/tobera/_analisis/'
def load(c):
    z=np.load(D+f'axis_{c}.npz'); return {k:z[k] for k in z.files}
def win(z,t0,t1):
    m=(z['t']>=t0-1e-9)&(z['t']<=t1+1e-9); 
    P=z['P'][m].mean(0);T=z['T'][m].mean(0);Ux=z['Ux'][m].mean(0)
    return P,T,Ux,m.sum()
def Mach(P,T,Ux): return Ux/np.sqrt(g*R*T)
def mdots(c,t0,t1):
    o=[]
    for q in ('gastoEntrada','gastoSalida'):
        a=np.loadtxt(D+f'{q}_caso_{c}.csv',delimiter=','); m=(a[:,0]>=t0)&(a[:,0]<=t1); o.append(abs(a[m,1].mean()))
    return o
def th_S():
    pB=290e3/p01; 
    # Me subsonico de pe/p0
    Me=np.sqrt(2/(g-1)*(pB**(-(g-1)/g)-1)); a_e=A_M(Me); Aast_At=AR/a_e  # Ae/A*=a_e -> A*/At=AR/a_e
    return Me,a_e/AR  # At/A*
if __name__=="__main__":
    for c,tw in (('S_r2',[(0.02,0.03),(0.03,0.04),(0.04,0.05),(0.05,0.06)]),('A_r2',[(0.03,0.04),(0.04,0.05),(0.05,0.06)]),('A_r1',[(0.03,0.04),(0.04,0.05),(0.05,0.06)])):
        z=load('caso_'+c); x=z['x']
        print('==',c)
        for t0,t1 in tw:
            P,T,Ux,n=win(z,t0,t1); M=Mach(P,T,Ux); mi,mo=mdots(c,t0,t1)
            print(f' ventana {t0}-{t1} (n={n}): M_exit {M[-1]:.4f} p_exit {P[-1]/1e3:.2f} kPa T_exit {T[-1]:.1f} M_garg(x~0) {M[np.argmin(abs(x))]:.4f} mdot_in {mi:.5f} mdot_out {mo:.5f}')
