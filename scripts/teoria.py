"""Teoria cuasiunidimensional con la A(x) REAL de tobera.geo (perfil coseno), 2D planar."""
import numpy as np
from scipy.optimize import brentq
g,R=1.4,287.0; p01,T01=300e3,300.0
ht=0.025; AR=1.6; he=AR*ht; hi=3*ht; Lc=0.25; Ld=0.20
def h(x):
    x=np.asarray(x,float)
    return np.where(x<0, ht+(hi-ht)*(1+np.cos(np.pi*(x+Lc)/Lc))/2,
                         ht+(he-ht)*(1-np.cos(np.pi*x/Ld))/2)
def A_M(M): return (1/M)*((2/(g+1))*(1+(g-1)/2*M**2))**((g+1)/(2*(g-1)))
def pr(M): return (1+(g-1)/2*M**2)**(-g/(g-1))
def Tr(M): return 1/(1+(g-1)/2*M**2)
def M_from_area(a,sup):
    if abs(a-1)<1e-12: return 1.0
    return brentq(lambda M:A_M(M)-a, 1+1e-9, 8) if sup else brentq(lambda M:A_M(M)-a, 1e-6, 1-1e-9)
def sigma(Ms):   # p02/p01 a traves de choque normal
    t1=((g+1)*Ms**2/((g-1)*Ms**2+2))**(g/(g-1)); t2=((g+1)/(2*g*Ms**2-(g-1)))**(1/(g-1)); return t1*t2
def M2(Ms): return np.sqrt(((g-1)*Ms**2+2)/(2*g*Ms**2-(g-1)))
def p2p1(Ms): return 1+2*g/(g+1)*(Ms**2-1)

def perfil_A(x):
    a=h(x)/ht; sup=np.asarray(x)>0
    M=np.array([M_from_area(ai,s) for ai,s in zip(np.atleast_1d(a),np.atleast_1d(sup))])
    return M

def caso_C(pB_p01):
    """Devuelve dict con Ms1,Ms2,xs,As/At,p2/p1,p02/p01 para choque en divergente."""
    def f(Ms):
        s=sigma(Ms); a_e=AR*s            # Ae/A2*
        Me=M_from_area(a_e,False)
        return s*pr(Me)-pB_p01
    Ms=brentq(f,1.0001,M_from_area(AR,True))
    As_At=A_M(Ms)                        # A_s/At = A_s/A*  (aguas arriba A*=At)
    xs=Ld/np.pi*np.arccos(1-2*(As_At-1)/(AR-1))   # h(xs)/ht=As_At
    s=sigma(Ms)
    Me=M_from_area(AR*s,False)
    return dict(Ms1=Ms,Ms2=M2(Ms),As_At=As_At,xs=xs,xs_Ld=xs/Ld,p2p1=p2p1(Ms),p02p01=s,
                Me=Me,pe_p01=s*pr(Me),Te_T01=Tr(Me))   # T0 constante

def estado_A():
    Me=M_from_area(AR,True)
    return dict(Me=Me,pe_p01=pr(Me),Te_T01=Tr(Me),Ve=Me*np.sqrt(g*R*T01*Tr(Me)))
def gasto_semi(dz=0.01):
    return ht*dz*p01*np.sqrt(g/(R*T01))*(2/(g+1))**((g+1)/(2*(g-1)))
if __name__=="__main__":
    print("A:",estado_A()); print("C pB/p01=0.6667:",caso_C(200e3/p01)); print("m semi",gasto_semi())
