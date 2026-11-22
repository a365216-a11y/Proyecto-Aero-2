import numpy as np
from ana import *
def relerr(a,b): return np.abs(a-b)/np.abs(b)*100
def perfilA(x): return perfil_A(x)
def theory_C(x,th):
    M=np.empty_like(x); pp=np.empty_like(x); s=th['p02p01']
    for k,xi in enumerate(x):
        if xi<0: M[k]=M_from_area(h(xi)/ht,False); pp[k]=pr(M[k])
        elif xi<=th['xs']: M[k]=M_from_area(h(xi)/ht,True); pp[k]=pr(M[k])
        else: M[k]=M_from_area(h(xi)/(ht/s),False); pp[k]=s*pr(M[k])
    return M,pp
def shock_pos(x,M):
    sup=np.where((M>1)&(x>0))[0]; l=sup[-1]
    return x[l]+(1-M[l])*(x[l+1]-x[l])/(M[l+1]-M[l]), l
if __name__=="__main__":
    ea=estado_A(); print('Teoria A:',ea,'mdot',gasto_semi())
    for c,(t0,t1) in (('A_r1',(0.032,0.06)),('A_r2',(0.032,0.06))):
        z=load('caso_'+c); x=z['x']; P,T,Ux,n=win(z,t0,t1); M=Mach(P,T,Ux)
        Mt=perfil_A(x); 
        print(f'== {c} ventana {t0}-{t1} n={n}')
        for nm,cv,tv in (('M',M,Mt),('p/p01',P/p01,pr(Mt)),('T/T01',T/T01,Tr(Mt))):
            e=relerr(cv,tv); print(f'  {nm}: error medio {e.mean():.2f}%  RMS {np.sqrt(np.mean(((cv-tv)/tv)**2))*100:.2f}%  max {e.max():.2f}% en x={x[e.argmax()]:.3f}')
        i=-1; print(f'  celda salida x={x[i]:.4f}: M {M[i]:.4f} (teo {Mt[i]:.4f}), p {P[i]/1e3:.2f} kPa (teo {p01*pr(Mt[i])/1e3:.2f}), T {T[i]:.1f} (teo {T01*Tr(Mt[i]):.1f}), V {Ux[i]:.1f}')
        ig=np.argmin(abs(x)); print(f'  garganta x={x[ig]:.4f}: M {M[ig]:.4f}')
        mi,mo=mdots(c,t0,t1); print(f'  gasto in {mi:.5f} out {mo:.5f} teo {gasto_semi():.5f}: err {abs(mo-gasto_semi())/gasto_semi()*100:.3f}%, balance {abs(mi-mo)/mi*100:.3f}%')
    thC=caso_C(200e3/p01); print('Teoria C:',{k:round(float(v),5) for k,v in thC.items()})
    for c,(t0,t1) in (('C_r1',(0.012,0.02)),('C_r1',(0.032,0.06)),('C_r4',(0.012,0.02)),('C_r4',(0.008,0.014))):
        z=load('caso_'+c); x=z['x']; P,T,Ux,n=win(z,t0,t1); M=Mach(P,T,Ux)
        xs,l=shock_pos(x,M); Mt,pt=theory_C(x,thC)
        print(f'== {c} ventana {t0}-{t1} n={n}: xs CFD {xs:.4f} (teo {thC["xs"]:.4f}) xs/Ld {xs/Ld:.4f} (teo {thC["xs_Ld"]:.4f}) error pos {abs(xs-thC["xs"])/Ld*100:.2f}% Ld')
        dx=np.gradient(x)[l]; 
        up=x<xs-0.03; dn=(x>xs+0.03)
        eM=relerr(M[up],Mt[up]); ep=relerr(P[up]/p01,pt[up]); print(f'  aguas arriba (x<xs-0.03): error M medio {eM.mean():.2f}% max {eM.max():.2f}%, p {ep.mean():.2f}% max {ep.max():.2f}%')
        # estados 1 y 2 cerca del choque: 3 celdas antes / despues del cruce sonico? usar xs-0.01, xs+0.015
        i1=np.argmin(abs(x-(xs-0.012))); i2=np.argmin(abs(x-(xs+0.012)))
        print(f'  estados: M1 {M[i1]:.3f} M2 {M[i2]:.3f} p2/p1 {P[i2]/P[i1]:.3f} (celdas x={x[i1]:.3f},{x[i2]:.3f})')
        eM2=relerr(M[dn],Mt[dn]); ep2=relerr(P[dn]/p01,pt[dn]); print(f'  aguas abajo (x>xs+0.03): error M medio {eM2.mean():.2f}% max {eM2.max():.2f}%, p medio {ep2.mean():.2f}% max {ep2.max():.2f}%')
        print(f'  salida: M {M[-1]:.3f} (teo {Mt[-1]:.3f}) p {P[-1]/1e3:.1f} kPa (teo {pt[-1]*p01/1e3:.1f}) T {T[-1]:.1f}')
        mi,mo=mdots(c,t0,t1); print(f'  gasto in {mi:.5f} out {mo:.5f} teo {gasto_semi():.5f}; err {abs(mo-gasto_semi())/gasto_semi()*100:.2f}% balance {abs(mi-mo)/mi*100:.2f}%')
