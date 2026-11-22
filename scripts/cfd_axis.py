import numpy as np, re, sys
def read_field(path, vec=False):
    s=open(path).read()
    m=re.search(r'internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(',s)
    n=int(m.group(2)); body=s[m.end():]
    if vec:
        vals=re.findall(r'\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)',body)[:n]
        return np.array(vals,float)
    out=[]
    for l in body.splitlines():
        l=l.strip()
        if l==')' : break
        if l: out.append(float(l))
    return np.array(out[:n])
def axis(case,t='0.02'):
    b=f'/mnt/user-data/uploads/tobera/{case}/{t}/'
    C=read_field(b+'C',True); p=read_field(b+'p'); T=read_field(b+'T'); U=read_field(b+'U',True)
    x,y=C[:,0],C[:,1]
    xr=np.round(x,6); xs=np.unique(xr)
    idx=[np.where(xr==xi)[0][np.argmin(y[xr==xi])] for xi in xs]   # fila de celdas junto al eje (simetria)
    idx=np.array(idx)
    Um=np.linalg.norm(U[idx],axis=1)
    M=Um/np.sqrt(1.4*287*T[idx])
    return x[idx],p[idx],T[idx],Um,M,y[idx]
if __name__=="__main__":
    x,p,T,U,M,y=axis(sys.argv[1]); print(len(x),x.min(),x.max(),y.max()); print(M[:3],M[-3:])
