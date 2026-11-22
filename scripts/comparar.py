import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from cfd_axis import read_field
from teoria import *
def axis(case,t='0.02'):
    b=f'/mnt/user-data/uploads/tobera/{case}/{t}/'
    C=read_field(b+'C',True); p=read_field(b+'p'); T=read_field(b+'T'); U=read_field(b+'U',True)
    x,y=C[:,0],C[:,1]; sel=np.where(y/h(x)<0.04)[0]; sel=sel[np.argsort(x[sel])]
    Um=np.linalg.norm(U[sel],axis=1); M=Um/np.sqrt(g*R*T[sel])
    return x[sel],p[sel],T[sel],Um,M
def relerr(a,b): return np.abs(a-b)/np.abs(b)*100
out={}
# ---------- CASO A ----------
x,p,T,U,M=axis('caso_A_r1'); Mt=perfil_A(x)
pt=pr(Mt); Tt=Tr(Mt)
print("== A_r1 (snapshot t=0.02 s, fila de celdas junto al eje)")
for nm,c,t in (("M",M,Mt),("p/p01",p/p01,pt),("T/T01",T/T01,Tt)):
    e=relerr(c,t); print(f"{nm}: error medio {e.mean():.2f}% max {e.max():.2f}% ; RMS relativo {np.sqrt(np.mean(((c-t)/t)**2))*100:.2f}%")
ea=estado_A(); print("estado salida teoria",ea)
print("CFD ultima celda: M=%.3f p/p01=%.4f T/T01=%.4f V=%.1f"%(M[-1],p[-1]/p01,T[-1]/T01,U[-1]))
fig,ax=plt.subplots(1,3,figsize=(13,3.6))
for a,(nm,c,t) in zip(ax,(("M",M,Mt),("p/p01",p/p01,pt),("T/T01",T/T01,Tt))):
    a.plot(x,t,'k-',label='Teoría 1D (A(x) real)'); a.plot(x,c,'r.',ms=4,label='CFD A r1'); a.set_xlabel('x [m]'); a.set_ylabel(nm); a.grid(alpha=.3)
ax[0].legend(); plt.tight_layout(); plt.savefig('A_r1_perfiles.png',dpi=150)
# ---------- CASO C ----------
x,p,T,U,M=axis('caso_C_r1'); th=caso_C(200e3/p01)
print("\n== C_r1; teoria:",{k:round(float(v),4) for k,v in th.items()})
dp=np.gradient(p,x); i=np.argmax(np.abs(dp)*(x>0)); 
# posicion choque: x donde M cruza por debajo de 1 despues de supersonico (interp) 
sup=np.where((M>1)&(x>0))[0]; last=sup[-1]; 
xs_cfd=x[last]+(1-M[last])*(x[last+1]-x[last])/(M[last+1]-M[last])
print("x_s CFD (M=1 post-choque) %.4f m -> xs/Ld=%.3f (teoria %.3f) ; error posicion normalizado %.2f%% de Ld"%(xs_cfd,xs_cfd/Ld,th['xs_Ld'],abs(xs_cfd-th['xs'])/Ld*100))
print("M alrededor del choque:",np.round(M[last-2:last+4],3),"x:",np.round(x[last-2:last+4],4))
Ms1=M[last-1]; Ms2=M[last+2]; p1=p[last-1]; p2=p[last+2]
print("Ms1~%.3f (teo %.3f) Ms2~%.3f (teo %.3f) p2/p1~%.3f (teo %.3f)"%(Ms1,th['Ms1'],Ms2,th['Ms2'],p2/p1,th['p2p1']))
def p0(pp,MM): return pp*(1+(g-1)/2*MM**2)**(g/(g-1))
print("p02/p01 CFD~%.3f (teo %.3f)"%(p0(p2,Ms2)/p0(p1,Ms1)*p0(p1,Ms1)/p01, th['p02p01']))
print("exit CFD: M=%.3f p=%.0f T=%.1f ; teoria Me=%.3f p=%.0f T=%.1f"%(M[-1],p[-1],T[-1],th['Me'],th['pe_p01']*p01,th['Te_T01']*T01))
# perfil teorico por tramos con choque en xs
Mt=np.empty_like(x)
for k,xi in enumerate(x):
    a=h(xi)/ht
    if xi<0: Mt[k]=M_from_area(a,False)
    elif xi<=th['xs']: Mt[k]=M_from_area(a,True)
    else: Mt[k]=M_from_area(a*(1/th['p02p01'])*1.0 if False else a/ (1/th['p02p01']) ,False) if False else M_from_area(h(xi)/(ht/th['p02p01']),False)
pt=np.where(x<=th['xs'],[pr(m) for m in Mt],[th['p02p01']*pr(m) for m in Mt])
fig,ax=plt.subplots(1,2,figsize=(9,3.6))
ax[0].plot(x,Mt,'k-',label='Teoría 1D'); ax[0].plot(x,M,'r.',ms=4,label='CFD C r1'); ax[0].set_ylabel('M'); ax[0].legend()
ax[1].plot(x,pt,'k-'); ax[1].plot(x,p/p01,'r.',ms=4); ax[1].set_ylabel('p/p01')
for a in ax: a.set_xlabel('x [m]'); a.grid(alpha=.3); a.axvline(th['xs'],color='b',ls=':')
plt.tight_layout(); plt.savefig('C_r1_perfiles.png',dpi=150)
e=relerr(M[x<th['xs']-0.02],Mt[x<th['xs']-0.02]); print("C: error M medio aguas arriba del choque %.2f%%"%e.mean())
