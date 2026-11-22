# -*- coding: utf-8 -*-
"""Genera datos.tex (macros y tablas) y figuras teoricas para la memoria de calculo."""
import numpy as np, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from teoria_tobera import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .3, 'figure.dpi': 150})
OUT = os.path.join(os.path.dirname(__file__), '..')
M_ = []
def mac(n, v, f="{:.4f}"): M_.append(r"\newcommand{\%s}{%s}" % (n, f.format(v) if not isinstance(v, str) else v))
def tex_table(path, header, rows, fmt, colspec):
    L = [r"\begin{tabular}{%s}" % colspec, r"\toprule", " & ".join(header) + r" \\", r"\midrule"]
    for r in rows: L.append(" & ".join(f.format(v) if not isinstance(v, str) else v for f, v in zip(fmt, r)) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]; open(os.path.join(OUT, 'tablas', path), 'w').write("\n".join(L))
def sci(v, n=2):
    m, e = f"{v:.{n}e}".split('e'); return r"$%s\times10^{%d}$" % (m, int(e))

# --- parametros ---
for n, v, f in [("gam", g, "{:.1f}"), ("Rgas", R, "{:.0f}"), ("pu", p01/1e3, "{:.0f}"), ("Tu", T01, "{:.0f}"),
                ("hth", ht, "{:.3f}"), ("heh", he, "{:.3f}"), ("hih", hi, "{:.3f}"), ("Lcc", Lc, "{:.2f}"), ("Ldd", Ld, "{:.2f}"),
                ("dzz", dz, "{:.2f}"), ("kk", k, "{:.1f}"), ("ARr", AR, "{:.1f}")]: mac(n, v, f)
# --- caso A: Newton rama supersonica y subsonica ---
Ms, hS = newton_M(AR, 2.0); Mb, hB = newton_M(AR, 0.4)
mac("MeSup", Ms, "{:.6f}"); mac("MeSub", Mb, "{:.6f}")
for nm, hh in (("newton_sup", hS), ("newton_sub", hB)):
    tex_table(nm + ".tex", [r"$n$", r"$M_n$", r"$f(M_n)=A/A^*(M_n)-1.6$", r"$f'(M_n)$", r"$\Delta M_n=-f/f'$"],
              [(i, M, f, dA_M(M), dM) for i, M, f, dM in hh], ["{}", "{:.10f}", "{:.3e}", "{:.6f}", "{:.3e}"], "rcccc")
pS, TS = pr(Ms), Tr(Ms); eA = estado_A()
mac("peSup", pS, "{:.6f}"); mac("TeSup", TS, "{:.6f}"); mac("peSupPa", eA['pe']/1e3, "{:.3f}"); mac("TeSupK", eA['Te'], "{:.2f}")
mac("VeSup", eA['Ve'], "{:.1f}"); mac("aeSup", np.sqrt(g*R*eA['Te']), "{:.1f}")
mac("peSubr", pr(Mb), "{:.6f}"); mac("TeSubr", Tr(Mb), "{:.6f}")
# --- gasto ---
m1 = mdot_garganta(); mac("mdotsemi", m1, "{:.5f}"); mac("mdotfull", 2*m1, "{:.5f}")
mac("rhostar", p01/(R*T01)*(2/(g+1))**(1/(g-1)), "{:.4f}"); mac("astar", np.sqrt(2*g*R*T01/(g+1)), "{:.2f}")
mac("rhocero", p01/(R*T01), "{:.4f}"); mac("coefm", np.sqrt(g/(R*T01))*(2/(g+1))**((g+1)/(2*(g-1))), "{:.6f}")
# --- rangos de contrapresion ---
sg = sigma(Ms); pbshock = sg * pr(newton_M(AR*sg, 0.5)[0])     # choque exactamente en la salida
# choque en la salida: p_e despues = pe_sup*p2/p1
pbshock2 = pS * Pi_choque(Ms)
mac("pbSubLim", pr(Mb), "{:.4f}"); mac("pbShockExit", pbshock2, "{:.4f}"); mac("pbAdapt", pS, "{:.4f}")
mac("pbSubLimkPa", pr(Mb)*p01/1e3, "{:.1f}"); mac("pbShockExitkPa", pbshock2*p01/1e3, "{:.1f}"); mac("pbAdaptkPa", pS*p01/1e3, "{:.1f}")
mac("pIIs", Pi_choque(Ms), "{:.4f}"); mac("sgExit", sg, "{:.4f}")
# --- caso S ---
pB_S = 276e3; S = caso_S(pB_S)
mac("pBS", pB_S/1e3, "{:.0f}"); mac("pBSr", pB_S/p01, "{:.4f}"); mac("MeS", S['Me'], "{:.5f}"); mac("AeAsS", S['Ae_Astar'], "{:.5f}")
mac("AtAsS", S['At_Astar'], "{:.5f}"); mac("MtS", S['Mt'], "{:.5f}"); mac("mdotS", S['mdot'], "{:.5f}"); mac("TeS", S['Te'], "{:.2f}"); mac("VeS", S['Ve'], "{:.2f}")
mac("AsAtS", 1/S['At_Astar'], "{:.5f}")
_, hSS = newton_M(S['At_Astar'], 0.5)
tex_table("newton_S.tex", [r"$n$", r"$M_n$", r"$f(M_n)$", r"$\Delta M_n$"], [(i, M, f, dM) for i, M, f, dM in hSS],
          ["{}", "{:.10f}", "{:.3e}", "{:.3e}"], "rccc")
# --- caso C ---
pB_C = 200e3; C = caso_C(pB_C)
for n, key, f in [("MsI", 'Ms1', "{:.6f}"), ("MsII", 'Ms2', "{:.6f}"), ("AsAt", 'As_At', "{:.6f}"), ("xsC", 'xs', "{:.5f}"), ("xsLd", 'xs_Ld', "{:.5f}"),
                  ("pIIpI", 'p2p1', "{:.5f}"), ("rrC", 'r', "{:.5f}"), ("pdospu", 'p02p01', "{:.6f}"), ("MeC", 'Me', "{:.6f}"), ("TeC", 'Te', "{:.2f}"), ("VeC", 'Ve', "{:.2f}")]:
    mac(n, C[key], f)
mac("pBC", pB_C/1e3, "{:.0f}"); mac("pBCr", pB_C/p01, "{:.5f}"); mac("pdos", C['p02p01']*p01/1e3, "{:.3f}")
mac("AstarC", 1/C['p02p01'], "{:.5f}"); mac("AeAstarC", AR*C['p02p01'], "{:.5f}")
MsC = C['Ms1']; MII = C['Ms2']
mac("TIIT", (C['p2p1']/C['r']), "{:.5f}"); mac("hxs", float(h(C['xs'])), "{:.6f}")
mac("pIkPa", pr(MsC)*p01/1e3, "{:.3f}"); mac("pIIkPa", pr(MsC)*C['p2p1']*p01/1e3, "{:.3f}"); mac("pIrel", pr(MsC), "{:.5f}")
mac("TIK", Tr(MsC)*T01, "{:.2f}"); mac("TIIK", Tr(MsC)*T01*C['p2p1']/C['r'], "{:.2f}")
mac("MsDen", (g-1)*MsC**2+2, "{:.5f}"); 
rows = [(i, lo, hi_, mid, Fm) for i, lo, hi_, mid, Fm in C['hist']]
tex_table("biseccion_C.tex", [r"$n$", r"$a_n$", r"$b_n$", r"$M_{s1}^{(n)}=(a_n+b_n)/2$", r"$F(M_{s1}^{(n)})$"], rows,
          ["{}", "{:.8f}", "{:.8f}", "{:.8f}", "{:+.3e}"], "rcccc")
mac("nbis", len(C['hist']), "{:d}")
# --- estaciones representativas ---
def est(xv, M0):
    a = float(h(xv)) / ht; M, hh = newton_M(a, M0); return a, M, hh
xA = 0.10; aA, MA, hA = est(xA, 1.5); xB = -0.10; aB, MB, hB2 = est(xB, 0.3)
mac("xAs", xA, "{:.2f}"); mac("hAs", float(h(xA)), "{:.5f}"); mac("aAs", aA, "{:.5f}"); mac("MAs", MA, "{:.6f}")
mac("pAs", pr(MA), "{:.6f}"); mac("TAs", Tr(MA), "{:.6f}"); mac("pAsPa", pr(MA)*p01/1e3, "{:.3f}"); mac("TAsK", Tr(MA)*T01, "{:.2f}")
mac("xBs", xB, "{:.2f}"); mac("hBs", float(h(xB)), "{:.6f}"); mac("aBs", aB, "{:.6f}"); mac("MBs", MB, "{:.6f}"); mac("pBs", pr(MB), "{:.6f}"); mac("TBs", Tr(MB), "{:.6f}")
mac("cosA", np.cos(np.pi*xA/Ld), "{:.1f}")
tex_table("newton_estacion.tex", [r"$n$", r"$M_n$", r"$f(M_n)$", r"$\Delta M_n$"], [(i, M, f, dM) for i, M, f, dM in hA],
          ["{}", "{:.10f}", "{:.3e}", "{:.3e}"], "rccc")
# --- tabla de geometria ---
xs_ = [-0.25, -0.20, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20]
tex_table("geometria.tex", [r"$x$ [m]", r"$h(x)$ [m]", r"$A/A_t=h/h_t$"], [(x, float(h(x)), float(h(x))/ht) for x in xs_], ["{:+.2f}", "{:.5f}", "{:.4f}"], "ccc")
# --- tabla de limites ---
tex_table("limites.tex", [r"R\'egimen", r"Condici\'on", r"$p_B/p_{01}$", r"$p_B$ [kPa]"],
          [("Subs\\'onico (S)", r"$p_B/p_{01} > %.4f$" % pr(Mb), "$>%.4f$" % pr(Mb), "$>%.1f$" % (pr(Mb)*p01/1e3)),
           ("Choque interno (C)", r"$%.4f < p_B/p_{01} < %.4f$" % (pbshock2, pr(Mb)), "$%.4f$ a $%.4f$" % (pbshock2, pr(Mb)), "%.1f a %.1f" % (pbshock2*p01/1e3, pr(Mb)*p01/1e3)),
           ("Choque en la salida", "---", "%.4f" % pbshock2, "%.1f" % (pbshock2*p01/1e3)),
           ("Adaptado (A)", "$p_B=p_{e,is}$", "%.4f" % pS, "%.1f" % (pS*p01/1e3)),
           ("Sobreexpandido (O)", r"$p_{e,is}<p_B<p_{\rm choque\ salida}$", "%.4f a %.4f" % (pS, pbshock2), "%.1f a %.1f" % (pS*p01/1e3, pbshock2*p01/1e3)),
           ("Subexpandido (U)", r"$p_B<p_{e,is}$", "$<%.4f$" % pS, "$<%.1f$" % (pS*p01/1e3))],
          ["{}", "{}", "{}", "{}"], "llcc")

# --- evaluacion manual de F en un valor de prueba Ms1 = 1.5 ---
Mt_ = 1.5; sg_ = sigma(Mt_); Me2_, hh_ = newton_M(AR*sg_, 0.5)
mac("MtrialA", Mt_, "{:.1f}"); mac("rTr", r_choque(Mt_), "{:.6f}"); mac("PiTr", Pi_choque(Mt_), "{:.6f}"); mac("MIITr", M2_choque(Mt_), "{:.6f}")
mac("sgTr", sg_, "{:.6f}"); mac("AeAdosTr", AR*sg_, "{:.6f}"); mac("MeTr", Me2_, "{:.6f}"); mac("prTr", pr(Me2_), "{:.6f}"); mac("FTr", F_choque(Mt_, pB_C/p01), "{:+.6f}")
mac("pBeTr", sg_*pr(Me2_), "{:.6f}")
mac("rpow", r_choque(Mt_)**(g/(g-1)), "{:.6f}"); mac("Pipow", Pi_choque(Mt_)**(-1/(g-1)), "{:.6f}")
# F en los extremos del intervalo
mac("FaExt", F_choque(1.0+1e-6, pB_C/p01), "{:+.5f}"); mac("FbExt", F_choque(Ms, pB_C/p01), "{:+.5f}")
# evaluacion numerica del choque a Ms1 de C
mac("rpowC", C['r']**(g/(g-1)), "{:.6f}"); mac("PipowC", C['p2p1']**(-1/(g-1)), "{:.6f}")
mac("PiCnum", 2*g/(g+1)*(MsC**2-1), "{:.6f}"); mac("MsCsq", MsC**2, "{:.6f}")
mac("numMdos", (g-1)*MsC**2+2, "{:.6f}"); mac("denMdos", 2*g*MsC**2-(g-1), "{:.6f}")
mac("MeShock", Ms, "{:.6f}")
open(os.path.join(OUT, 'datos.tex'), 'w').write("\n".join(M_) + "\n")
# --- figuras teoricas ---
x = np.linspace(-Lc, Ld, 600)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].plot(x, h(x), 'k'); ax[0].plot(x, -h(x), 'k'); ax[0].axhline(0, color='gray', lw=.6); ax[0].set_xlabel('x [m]'); ax[0].set_ylabel('±h(x) [m]'); ax[0].set_aspect('equal', 'datalim')
ax[0].set_title('Contorno de la tobera')
ax[1].plot(x, h(x)/ht, 'k'); ax[1].axhline(AR, color='r', ls='--', lw=.8, label=r'$A_e/A_t=1.6$'); ax[1].axvline(0, color='b', ls=':', lw=.8)
ax[1].set_xlabel('x [m]'); ax[1].set_ylabel(r'$A(x)/A_t$'); ax[1].legend(); ax[1].set_title(r'Distribución real $A(x)/A_t$')
plt.tight_layout(); plt.savefig(OUT + '/figuras/geometria_Ax.pdf'); plt.close()
xx = np.linspace(-Lc, Ld, 400); xx = xx[np.abs(xx) > 1e-9]
MA_, PA_, TA_ = perfil('A', xx)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.3))
for a, y, l in zip(ax, (MA_, PA_, TA_), ('M', r'$p/p_{01}$', r'$T/T_{01}$')): a.plot(xx, y, 'k'); a.set_xlabel('x [m]'); a.set_ylabel(l)
plt.suptitle('Perfiles teóricos casos A (rama supersónica en el divergente)'); plt.tight_layout(); plt.savefig(OUT + '/figuras/teoria_perfiles_A.pdf'); plt.close()
MC_, PC_, TC_ = perfil('C', xx, pB_C)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.3))
for a, y, l in zip(ax, (MC_, PC_, TC_), ('M', r'$p/p_{01}$', r'$T/T_{01}$')): a.plot(xx, y, 'k'); a.axvline(C['xs'], color='r', ls=':'); a.set_xlabel('x [m]'); a.set_ylabel(l)
plt.suptitle(r'Perfiles teóricos caso C ($p_B/p_{01}=%.4f$); línea roja: choque en $x_s=%.4f$ m' % (pB_C/p01, C['xs'])); plt.tight_layout(); plt.savefig(OUT + '/figuras/teoria_perfiles_C.pdf'); plt.close()
Ms_grid = np.linspace(1.0005, Ms - 1e-4, 300); Fv = [F_choque(m, pB_C/p01) + pB_C/p01 for m in Ms_grid]
fig, ax = plt.subplots(figsize=(5, 3.4)); ax.plot(Ms_grid, Fv, 'k'); ax.axhline(pB_C/p01, color='r', ls='--', label=r'$p_B/p_{01}=%.4f$' % (pB_C/p01))
ax.plot(C['Ms1'], pB_C/p01, 'ro'); ax.set_xlabel(r'$M_{s1}$'); ax.set_ylabel(r'$p_B/p_{01}$ resultante'); ax.legend(); plt.tight_layout(); plt.savefig(OUT + '/figuras/F_choque.pdf'); plt.close()
print("OK", len(M_), "macros")
