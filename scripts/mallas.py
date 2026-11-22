import re, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
U='/mnt/user-data/uploads/tobera/'
def load(case):
    d=U+case+'/constant/polyMesh/'
    t=open(d+'points').read(); i=t.index('('); 
    pts=np.array(re.findall(r'\(([^()\s]+) ([^()\s]+) ([^()\s]+)\)',t[i:]),float)
    f=open(d+'faces').read()
    F=re.findall(r'4\((\d+) (\d+) (\d+) (\d+)\)',f)
    F=np.array(F,int)
    z0=pts[:,2].min()
    m=np.all(np.isclose(pts[F][:,:,2],z0),axis=1)
    F=F[m]
    segs=[]
    for q in F:
        P=pts[q][:,:2]
        for k in range(4): segs.append([P[k],P[(k+1)%4]])
    return pts,np.array(segs)
casos=[('caso_A_r1_co01','r1: 2 600 celdas'),('caso_A_r2_co01','r2: 10 400 celdas'),('caso_C_r4_co01','r4: 41 600 celdas')]
D={c:load(c) for c,_ in casos}
fig=plt.figure(figsize=(11,8))
gs=fig.add_gridspec(2,3,height_ratios=[1,1.1])
ax=fig.add_subplot(gs[0,:])
pts,segs=D['caso_A_r1_co01']
for sg in (1,-1):
    s=segs.copy(); s[:,:,1]*=sg
    ax.add_collection(LineCollection(s,colors='0.25',linewidths=0.4))
ax.set_xlim(-0.26,0.21); ax.set_ylim(-0.085,0.085); ax.set_aspect('equal')
ax.set_title('Geometría y malla r1 (semidominio reflejado respecto al plano de simetría)',fontsize=10)
ax.set_xlabel('$x$ [m]'); ax.set_ylabel('$y$ [m]')
for k,(c,lab) in enumerate(casos):
    a=fig.add_subplot(gs[1,k]); pts,segs=D[c]
    a.add_collection(LineCollection(segs,colors='0.25',linewidths=0.5))
    a.set_xlim(-0.03,0.03); a.set_ylim(0,0.04); a.set_aspect('equal')
    a.set_title('Detalle en la garganta, malla '+lab,fontsize=9); a.set_xlabel('$x$ [m]')
    if k==0: a.set_ylabel('$y$ [m]')
for a in fig.axes:
    for s in ('top','right'): a.spines[s].set_visible(False)
fig.tight_layout()
for e in ('png','pdf'): fig.savefig('/mnt/user-data/outputs/entrega/figuras/geometria_mallas.'+e,dpi=200)
print('ok')
