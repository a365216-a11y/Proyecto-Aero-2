# -*- coding: utf-8 -*-
"""
teoria_tobera.py
Modelo cuasiunidimensional (gas perfecto, adiabatico, isentropico salvo choque normal)
de la tobera convergente-divergente 2D planar del proyecto (Ae/At = 1.6).
Todas las funciones se usan en la memoria de calculo. Unidades SI.
"""
import numpy as np

# ---------------- Gas y alimentacion (valores provisionales justificados en el reporte) -------------
g, R = 1.4, 287.0            # aire: gamma, constante especifica [J/(kg K)]
p01, T01 = 300e3, 300.0      # presion y temperatura totales de alimentacion [Pa], [K]
k = (g - 1) / 2              # constante auxiliar (gamma-1)/2

# ---------------- Geometria (identica a tobera.geo) ----------------
ht = 0.025                   # semialtura de garganta [m]
AR = 1.6                     # Ae/At asignada
he = AR * ht                 # semialtura de salida [m]
hi = 3.0 * ht                # semialtura de entrada [m]
Lc, Ld = 0.25, 0.20          # longitudes convergente y divergente [m]
dz = 0.01                    # espesor en z del dominio numerico [m]

def h(x):
    """Semialtura de la pared h(x) [m]; garganta en x = 0."""
    x = np.asarray(x, float)
    conv = ht + (hi - ht) * (1 + np.cos(np.pi * (x + Lc) / Lc)) / 2
    div = ht + (he - ht) * (1 - np.cos(np.pi * x / Ld)) / 2
    return np.where(x < 0, conv, div)

def x_de_area(a):
    """Inversa en el divergente: posicion x tal que A(x)/At = a (1 <= a <= AR)."""
    return Ld / np.pi * np.arccos(1 - 2 * (a - 1) / (AR - 1))

# ---------------- Relaciones isentropicas ----------------
def A_M(M):                  # A/A* en funcion de M
    return (1 / M) * ((1 + k * M**2) / (1 + k)) ** ((g + 1) / (2 * (g - 1)))
def dA_M(M):                 # derivada d(A/A*)/dM
    return A_M(M) * (M**2 - 1) / (M * (1 + k * M**2))
def pr(M):  return (1 + k * M**2) ** (-g / (g - 1))      # p/p0
def Tr(M):  return 1 / (1 + k * M**2)                     # T/T0
def rr(M):  return (1 + k * M**2) ** (-1 / (g - 1))      # rho/rho0
def M_de_p(pratio):          # inversa de p/p0 (subsonico o supersonico)
    return np.sqrt((pratio ** (-(g - 1) / g) - 1) / k)

def newton_M(a, M0, tol=1e-13, maxit=60):
    """Resuelve A_M(M) = a con Newton-Raphson. Devuelve (M, historial)."""
    M, hist = M0, []
    for it in range(maxit):
        f = A_M(M) - a
        dM = -f / dA_M(M)
        hist.append((it, M, f, dM))
        M += dM
        if abs(dM) < tol: break
    return M, hist

# ---------------- Choque normal ----------------
def r_choque(M1):  return (g + 1) * M1**2 / ((g - 1) * M1**2 + 2)         # rho2/rho1
def Pi_choque(M1): return 1 + 2 * g / (g + 1) * (M1**2 - 1)                # p2/p1
def M2_choque(M1): return np.sqrt(((g - 1) * M1**2 + 2) / (2 * g * M1**2 - (g - 1)))
def sigma(M1):               # p02/p01 = r^(g/(g-1)) * Pi^(-1/(g-1))
    return r_choque(M1) ** (g / (g - 1)) * Pi_choque(M1) ** (-1 / (g - 1))

# ---------------- Gasto ----------------
def mdot_garganta(p0=p01, T0=T01, At=None):
    """Gasto bloqueado: m = At p0 sqrt(g/(R T0)) (2/(g+1))^((g+1)/(2(g-1)))."""
    At = ht * dz if At is None else At      # semidominio: ht*dz
    return At * p0 * np.sqrt(g / (R * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))

# ---------------- Casos ----------------
def estado_A():
    Me, _ = newton_M(AR, 2.0)
    return dict(Me=Me, pe=pr(Me) * p01, Te=Tr(Me) * T01,
                Ve=Me * np.sqrt(g * R * Tr(Me) * T01))

def caso_S(pB):
    """Subsonico no estrangulado: p_e = pB, M_e < 1, garganta no sonica."""
    Me = M_de_p(pB / p01)
    Ae_Astar = A_M(Me)
    At_Astar = Ae_Astar / AR
    Mt, _ = newton_M(At_Astar, 0.5)
    return dict(Me=Me, Ae_Astar=Ae_Astar, At_Astar=At_Astar, Mt=Mt,
                mdot=mdot_garganta() / At_Astar, Te=Tr(Me) * T01,
                Ve=Me * np.sqrt(g * R * Tr(Me) * T01))

def F_choque(Ms, pB_p01):
    s = sigma(Ms)
    Me2, _ = newton_M(AR * s, 0.5)
    return s * pr(Me2) - pB_p01

def caso_C(pB, tol=1e-12):
    """Choque normal en el divergente (biseccion sobre Ms1)."""
    lo, hi_ = 1.0 + 1e-6, newton_M(AR, 2.0)[0]
    hist = []
    for it in range(80):
        mid = 0.5 * (lo + hi_)
        Fm = F_choque(mid, pB / p01)
        hist.append((it, lo, hi_, mid, Fm))
        if Fm > 0: lo = mid           # F decrece con Ms: raiz a la derecha
        else: hi_ = mid
        if hi_ - lo < tol: break
    Ms = 0.5 * (lo + hi_)
    s = sigma(Ms)
    Me2, _ = newton_M(AR * s, 0.5)
    As_At = A_M(Ms)
    xs = float(x_de_area(As_At))
    return dict(Ms1=Ms, Ms2=M2_choque(Ms), As_At=As_At, xs=xs, xs_Ld=xs / Ld,
                p2p1=Pi_choque(Ms), r=r_choque(Ms), p02p01=s, Me=Me2, pe=s * pr(Me2) * p01,
                Te=Tr(Me2) * T01, Ve=Me2 * np.sqrt(g * R * Tr(Me2) * T01), hist=hist)

def perfil(case, x, pB=None):
    """Perfiles teoricos M, p/p01, T/T01 a lo largo del eje (misma A(x) real)."""
    x = np.atleast_1d(x); M = np.empty_like(x); P = np.empty_like(x); T = np.empty_like(x)
    if case == 'C': C = caso_C(pB); s = C['p02p01']
    for i, xi in enumerate(x):
        a = float(h(xi)) / ht
        if case == 'A':
            Mi = newton_M(a, 2.0 if xi > 0 else 0.4)[0] if abs(a - 1) > 1e-12 else 1.0
            M[i], P[i], T[i] = Mi, pr(Mi), Tr(Mi)
        elif case == 'C':
            if xi <= C['xs']:
                Mi = newton_M(a, 2.0 if xi > 0 else 0.4)[0] if abs(a - 1) > 1e-12 else 1.0
                M[i], P[i], T[i] = Mi, pr(Mi), Tr(Mi)
            else:
                Mi = newton_M(a * s, 0.4)[0]     # A/A2* = (A/At) * (p02/p01)
                M[i], P[i], T[i] = Mi, s * pr(Mi), Tr(Mi)
    return M, P, T

if __name__ == "__main__":
    print(estado_A()); print(caso_S(276e3)); C = caso_C(200e3); C.pop('hist'); print(C)
    print("mdot semidominio [kg/s] =", mdot_garganta())
