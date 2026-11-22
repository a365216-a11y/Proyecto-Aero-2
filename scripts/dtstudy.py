import numpy as np, json
from ana import load, mdots
from ana3 import secwin
from ana2 import theory_C, shock_pos, relerr
from teoria import *
thC=caso_C(200e3/p01)
pB_S=276e3; Me_S=np.sqrt(2/(g-1)*((pB_S/p01)**(-(g-1)/g)-1)); AtA_S=A_M(Me_S)/AR; mdot_S=gasto_semi()/AtA_S
def thM(kind,x):
    if kind=='A': return perfil_A(x)
    if kind=='S': return np.array([M_from_area(h(xi)/ht*AtA_S,False) for xi in x])
    return theory_C(x,thC)[0]
def metr(kind,name,win):
    z=load('caso_'+name); x,P,T,U,mf=secwin(z,*win); M=U/np.sqrt(g*R*T); Mt=thM(kind,x)
    sel=np.ones_like(x,bool); out={}
    if kind=='C':
        xs,_=shock_pos(x,M); out['xs_err']=abs(xs-thC['xs'])/Ld*100; sel=(x<xs-0.03)|(x>xs+0.03); out['xs']=xs/Ld
    out['M_err']=relerr(M,Mt)[sel].mean()
    mo=mdots(name,*win)[1]; mt=mdot_S if kind=='S' else gasto_semi(); out['mdot_err']=100*(mo-mt)/mt
    return out
runs={'A':[(0.05,'dt_A_r1_co005',(0.03,0.06)),(0.1,'A_r1_co01',(0.03,0.06)),(0.2,'dt_A_r1_co02',(0.03,0.06)),(0.3,'A_r1',(0.07,0.1))],
 'C':[(0.05,'dt_C_r1_co005',(0.03,0.06)),(0.1,'C_r1_co01',(0.03,0.06)),(0.2,'dt_C_r1_co02',(0.03,0.06)),(0.5,'C_r1',(0.07,0.1))],
 'S':[(0.05,'dt_S_r1_co005',(0.07,0.1)),(0.1,'S_r1_p276_co01',(0.07,0.1)),(0.2,'dt_S_r1_co02',(0.07,0.1)),(0.3,'S_r1_p276_base',(0.04,0.06))]}
res={}
for k,L in runs.items():
    res[k]=[]
    for co,nm,w in L:
        m=metr(k,nm,w); m['co']=co; res[k].append(m); print(k,co,{a:round(float(b),4) for a,b in m.items()})
u=metr('S','dt_S_r1_unif',(0.07,0.1)); n=metr('S','S_r1_p276_co01',(0.07,0.1)); print('unif',u,'\nnonunif',n)
res['unif']=u; res['nonunif']=n
# pendientes log-log del error vs CFL (puntos con CFL<=0.2)
for k,key in (('A','M_err'),('C','xs_err'),('S','M_err'),('S','mdot_err')):
    L=[(m['co'],abs(m[key])) for m in res[k] if m['co']<=0.2]
    c=np.log([a for a,b in L]); e=np.log([b for a,b in L]); s=np.polyfit(c,e,1)[0]; print(k,key,'pendiente',round(s,2),L)
json.dump(res,open('dtstudy.json','w'),default=float)
