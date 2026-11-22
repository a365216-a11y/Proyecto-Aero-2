import numpy as np
from ana import *
from ana2 import theory_C, shock_pos, relerr
def secwin(z,t0,t1):
    m=(z['t']>=t0-1e-9)&(z['t']<=t1+1e-9)
    return z['xs'],z['PS'][m].mean(0),z['TS'][m].mean(0),z['US'][m].mean(0),z['MS'][m].mean(0)
if __name__=="__main__":
    thC=caso_C(200e3/p01)
    for c,(t0,t1) in (('A_r1',(0.032,0.06)),('A_r2',(0.032,0.06)),('C_r1',(0.032,0.06)),('C_r1',(0.012,0.02)),('C_r4',(0.012,0.02))):
        z=load('caso_'+c); x,P,T,U,mf=secwin(z,t0,t1); M=U/np.sqrt(g*R*T)
        mdot_x=mf*h(x)*0.01
        print(f'== {c} sec-avg ventana {t0}-{t1}: gasto por seccion (min/max) {mdot_x.min():.5f}/{mdot_x.max():.5f}')
        if c.startswith('A'):
            Mt=perfil_A(x); 
        else:
            xs,l=shock_pos(x,M); Mt,pt=theory_C(x,thC); print(f'  xs(sec) {xs:.4f} xs/Ld {xs/Ld:.4f} err {abs(xs-thC["xs"])/Ld*100:.2f}%Ld')
        pt=pr(Mt) if c.startswith('A') else pt; Tt=Tr(Mt) if c.startswith('A') else None
        sel=np.ones_like(x,bool) if c.startswith('A') else ((x<xs-0.03)|(x>xs+0.03))
        for nm,cv,tv in (('M',M,Mt),('p/p01',P/p01,pt)):
            e=relerr(cv,tv); print(f'  {nm}: error medio {e[sel].mean():.2f}% RMS {np.sqrt(np.mean(((cv-tv)/tv)[sel]**2))*100:.2f}% max {e[sel].max():.2f}% (x={x[sel][e[sel].argmax()]:.3f}); fuera del choque')
        if c.startswith('C'):
            i1=np.argmin(abs(x-(xs-0.012))); i2=np.argmin(abs(x+0-(xs+0.012)))
            print(f'  M1 {M[i1]:.3f} (teo {thC["Ms1"]:.3f} en xs teorico) M2 {M[i2]:.3f}; p2/p1 {P[i2]/P[i1]:.3f}; salida M {M[-1]:.3f} (teo {thC["Me"]:.3f}), p {P[-1]/1e3:.1f}, T {T[-1]:.1f} (teo 284.2)')
