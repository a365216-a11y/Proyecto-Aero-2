import numpy as np
from cfd_axis import read_field
from teoria import *
TIMES=['0.012','0.014','0.016','0.018','0.02']
def axis_avg(case,times=TIMES):
    base=f'/mnt/user-data/uploads/tobera/{case}/'
    C=read_field(base+'0.02/C',True); x,y=C[:,0],C[:,1]
    sel=np.where(y/h(x)<0.04)[0]; sel=sel[np.argsort(x[sel])]
    P=[];Tt=[];Uu=[]
    for t in times:
        P.append(read_field(base+f'{t}/p')[sel]); Tt.append(read_field(base+f'{t}/T')[sel]); Uu.append(np.linalg.norm(read_field(base+f'{t}/U',True)[sel],axis=1))
    P,Tt,Uu=map(np.array,(P,Tt,Uu))
    return x[sel],P,Tt,Uu     # (ntimes,ncells)
if __name__=="__main__":
    for case in ('caso_A_r1','caso_C_r1'):
        x,P,T,U=axis_avg(case); M=U/np.sqrt(g*R*T)
        print(case,"variacion temporal (desv.est./media) en celdas de salida: p %.1f%%  M %.1f%%"%(100*P[:,-5:].std(0).mean()/P[:,-5:].mean(),100*M[:,-5:].std(0).mean()/M[:,-5:].mean()))
        print(" p_exit por instante:",np.round(P[:,-1]/1e3,1),"kPa ; M_exit:",np.round(M[:,-1],3))
