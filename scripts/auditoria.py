"""Auditoria: ventanas consecutivas equivalentes (5 ms), definicion de Mach, S a 276 kPa."""
import numpy as np, ana
from ana import D
from ana2 import shock_pos
from teoria import g,R,p01,T01
from resultados import FINAL, pB_S
def win(z,a,b):
    m=(z['t']>a-1e-9)&(z['t']<=b+1e-9)
    return z['PS'][m].mean(0),z['TS'][m].mean(0),z['US'][m].mean(0),z['MS'][m].mean(0)
def g_(n,q,a,b):
    d=np.loadtxt(D+f'{q}_caso_{n}.csv',delimiter=','); m=(d[:,0]>a)&(d[:,0]<=b); return abs(d[m,1].mean())
rows=[]
for k,(n,w,m_,_) in FINAL.items():
    z=np.load(D+f'axis_caso_{n}.npz'); tend=z['t'].max(); x=z['xs']
    out={}
    for L in (0.005,0.010):
        if tend-2*L<0: continue
        w1=(tend-2*L,tend-L); w2=(tend-L,tend)
        P1,T1,U1,M1=win(z,*w1); P2,T2,U2,M2=win(z,*w2)
        Ma=U1/np.sqrt(g*R*T1); Mb=U2/np.sqrt(g*R*T2)
        sel=np.ones_like(x,bool)
        if k.startswith('C'):
            xs1,_=shock_pos(x,Ma); xs2,_=shock_pos(x,Mb); sel=(x<min(xs1,xs2)-0.03)|(x>max(xs1,xs2)+0.03)
        dM=np.abs(Mb-Ma)/np.abs(Ma)*100
        o1=g_(n,'gastoSalida',*w1);o2=g_(n,'gastoSalida',*w2);i2=g_(n,'gastoEntrada',*w2)
        r=dict(L=L*1e3,w=(round(w1[0],4),round(w2[1],4)),dM_mean=dM[sel].mean(),dM_max=dM[sel].max(),dmdot=abs(o2-o1)/o1*100,bal=abs(i2-o2)/i2*100)
        if k.startswith('C'): r['dxs_Ld']=abs(xs2-xs1)/0.2*100
        out[L]=r
    # definicion de Mach
    P,T,U,MS=win(z,*w); Md=U/np.sqrt(g*R*T)
    dd=np.abs(Md-MS)/MS*100
    pout=P[-1]; 
    print(k,n); [print('  ',v) for v in out.values()]
    print(f'   Mach def: |M(U_avg,T_avg)-avg(M_local)|/M media {dd.mean():.3f}% max {dd.max():.3f}%   p_salida={pout:.0f} Pa')

# --- variacion del gasto en los ultimos 100 pasos (el gasto se escribe cada 10 pasos: ultimas 11 muestras) ---
print('\nultimos 100 pasos: variacion (max-min)/media del gasto')
for k,(n,w,m_,_) in FINAL.items():
    r=[]
    for q in ('gastoEntrada','gastoSalida'):
        a=np.loadtxt(D+f'{q}_caso_{n}.csv',delimiter=','); v=np.abs(a[-11:,1]); r.append((v.max()-v.min())/v.mean()*100)
    print(f'{k}: entrada {r[0]:.3f}%  salida {r[1]:.3f}%  pasos ~{len(a)*10}')
