"""Resultados finales: figuras y tablas (CFD vs teoria cuasi-1D con la A(x) real).
Uso: python3 resultados.py   -> escribe en ./figuras y ./tablas
"""
import os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from ana import load, mdots, D
from ana3 import secwin
from ana2 import theory_C, shock_pos, relerr
from teoria import *

OUT_F = 'figuras'; OUT_T = 'tablas'
os.makedirs(OUT_F, exist_ok=True); os.makedirs(OUT_T, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': .25, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.dpi': 120, 'savefig.dpi': 200, 'legend.frameon': False})
K, B, O, G, SK = '#222222', '#0072B2', '#D55E00', '#009E73', '#56B4E9'
MESH_COL = {'r1': SK, 'r2': B, 'r4': O}

pB_S = 276e3
Me_S = np.sqrt(2/(g-1)*((pB_S/p01)**(-(g-1)/g)-1)); AtA_S = A_M(Me_S)/AR; mdot_S = gasto_semi()/AtA_S
pB_C = 200e3
thC = caso_C(pB_C/p01)

# caso final -> (nombre en disco, ventana de promedio, malla)
FINAL = {
    'A_r1': ('A_r1_co01', (0.03, 0.06), 'r1', 1),
    'A_r2': ('A_r2_co01', (0.02, 0.04), 'r2', 2),
    'S_r2': ('S_r2_p276_co01', (0.03, 0.05), 'r2', 2),
    'S_r1': ('S_r1_p276_co01', (0.04, 0.10), 'r1', 1),
    'C_r1': ('C_r1_co01', (0.03, 0.06), 'r1', 1),
    'C_r2': ('C_r2_co01', (0.02, 0.04), 'r2', 2),
    'C_r4': ('C_r4_co01', (0.012, 0.02), 'r4', 4),
}
def exists(name): return os.path.exists(D+f'axis_caso_{name}.npz')

def theory_profile(kind, x):
    if kind == 'A':
        M = perfil_A(x); p = pr(M); T = Tr(M)
    elif kind == 'S':
        M = np.array([M_from_area(h(xi)/ht*AtA_S, False) for xi in x]); p = pr(M); T = Tr(M)
    else:
        M, p = theory_C(x, thC)
        s = thC['p02p01']
        T = np.where(x <= thC['xs'], Tr(M), Tr(M))   # T0 constante en el choque
    return M, p, T

def cfd_profile(name, win):
    z = load('caso_'+name); x, P, T, U, mf = secwin(z, *win)
    M = U/np.sqrt(g*R*T)
    return x, M, P/p01, T/T01, mf

def metrics(kind, name, win):
    x, M, p, T, mf = cfd_profile(name, win)
    Mt, pt, Tt = theory_profile(kind, x)
    sel = np.ones_like(x, bool)
    out = {}
    if kind == 'C':
        xs, l = shock_pos(x, M); out['xs'] = xs
        sel = (x < xs-0.03) | (x > xs+0.03)
    for nm, c, t in (('M', M, Mt), ('p', p, pt), ('T', T, Tt)):
        e = relerr(c, t)[sel]
        out[nm] = (e.mean(), np.sqrt(np.mean(((c-t)/t)[sel]**2))*100, e.max())
    mi, mo = mdots(name, *win)
    ref = gasto_semi() if kind != 'S' else mdot_S
    out.update(mi=mi, mo=mo, ref=ref, err_m=abs(mo-ref)/ref*100, bal=abs(mi-mo)/mi*100)
    i = -1
    out['exit'] = dict(x=x[i], M=M[i], Mt=Mt[i], p=p[i]*p01, pt=pt[i]*p01, T=T[i]*T01, Tt=Tt[i]*T01,
                       V=M[i]*np.sqrt(g*R*T[i]*T01), Vt=Mt[i]*np.sqrt(g*R*Tt[i]*T01))
    return out

def savefig(fig, name):
    fig.savefig(f'{OUT_F}/{name}.pdf', bbox_inches='tight'); fig.savefig(f'{OUT_F}/{name}.png', bbox_inches='tight'); plt.close(fig)

xs_fine = np.linspace(-0.25, 0.2, 600)

# ------------------------------------------------------------ perfiles A y S
def fig_perfiles(kind, name, win, title, fname):
    x, M, p, T, mf = cfd_profile(name, win)
    Mt, pt, Tt = theory_profile(kind, xs_fine)
    Mtx, ptx, Ttx = theory_profile(kind, x)
    fig, ax = plt.subplots(2, 3, figsize=(10.5, 5.2), sharex=True, gridspec_kw={'height_ratios': [3, 1.3]})
    for j, (lab, c, t, tf) in enumerate((('M', M, Mtx, Mt), ('$p/p_{01}$', p, ptx, pt), ('$T/T_{01}$', T, Ttx, Tt))):
        ax[0, j].plot(xs_fine, tf, '-', color=K, lw=1.6, label='Teoría cuasi-1D, $A(x)$ real')
        ax[0, j].plot(x, c, 'o', color=B, ms=3, mfc='none', mew=.9, label=f'CFD ({name}, prom. {win[0]}–{win[1]} s)')
        ax[0, j].set_ylabel(lab)
        ax[1, j].plot(x, (c-t)/t*100, '-', color=O, lw=1.2); ax[1, j].axhline(0, color=K, lw=.6)
        ax[1, j].set_ylabel('error rel. [%]'); ax[1, j].set_xlabel('$x$ [m]')
        ax[0, j].axvline(0, color='gray', lw=.6, ls=':')
    ax[0, 0].legend(loc='best', fontsize=8)
    fig.suptitle(title, y=.995, fontsize=10)
    fig.tight_layout(); savefig(fig, fname)

# ------------------------------------------------------------ perfiles C (tres mallas)
def fig_perfiles_C(meshes, fname):
    fig, ax = plt.subplots(2, 3, figsize=(10.5, 5.2), sharex=True, gridspec_kw={'height_ratios': [3, 1.3]})
    Mt, pt = theory_C(xs_fine, thC); Tt = Tr(Mt)
    for j, (lab, tf) in enumerate((('M', Mt), ('$p/p_{01}$', pt), ('$T/T_{01}$', Tt))):
        ax[0, j].plot(xs_fine, tf, '-', color=K, lw=1.6, label='Teoría cuasi-1D, $A(x)$ real')
        ax[0, j].axvline(thC['xs'], color=K, lw=.8, ls='--'); ax[0, j].set_ylabel(lab)
    for key in meshes:
        name, win, mesh, _ = FINAL[key]
        x, M, p, T, mf = cfd_profile(name, win); Mx, px = theory_C(x, thC); Tx = Tr(Mx)
        xs, _ = shock_pos(x, M)
        for j, (c, t) in enumerate(((M, Mx), (p, px), (T, Tx))):
            ax[0, j].plot(x, c, '.', color=MESH_COL[mesh], ms=3.5, label=f'CFD {mesh} (choque en $x/L_d$={xs/Ld:.3f})' if j == 0 else None)
            e = (c-t)/t*100; m = (x < xs-0.03) | (x > xs+0.03)
            ax[1, j].plot(x[m], e[m], '-', color=MESH_COL[mesh], lw=1)
    for j in range(3):
        ax[1, j].axhline(0, color=K, lw=.6); ax[1, j].set_ylabel('error rel. [%] (fuera del choque)'); ax[1, j].set_xlabel('$x$ [m]')
    ax[0, 0].legend(fontsize=7.5, loc='lower left')
    fig.suptitle(f'Caso C: choque interno, $p_B/p_{{01}}$={pB_C/p01:.3f}. Línea punteada: $x_s$ teórico ({thC["xs_Ld"]:.4f} $L_d$)', y=.995, fontsize=10)
    fig.tight_layout(); savefig(fig, fname)

# ------------------------------------------------------------ mapas
def fig_mapas(name, kind, fname, title):
    z = load('caso_'+name)
    fx, fy, fp, fT, fUx, fUy = (z[k] for k in ('fx', 'fy', 'fp', 'fT', 'fUx', 'fUy'))
    M = np.hypot(fUx, fUy)/np.sqrt(g*R*fT)
    X = np.r_[fx, fx]; Y = np.r_[fy, -fy]
    tri = Triangulation(X, Y)
    cx = X[tri.triangles].mean(1); cy = np.abs(Y[tri.triangles].mean(1))
    tri.set_mask(cy > h(cx)*1.0)
    fig, ax = plt.subplots(2, 1, figsize=(8, 5.6), sharex=True)
    for a, v, lab, cm in ((ax[0], np.r_[M, M], 'Número de Mach', 'viridis'), (ax[1], np.r_[fp, fp]/p01, '$p/p_{01}$', 'magma')):
        c = a.tricontourf(tri, v, 40, cmap=cm); cb = fig.colorbar(c, ax=a, pad=.01); cb.set_label(lab)
        xx = np.linspace(-0.25, 0.2, 400); a.plot(xx, h(xx), 'k', lw=1.2); a.plot(xx, -h(xx), 'k', lw=1.2)
        a.set_aspect('equal'); a.set_ylabel('$y$ [m]'); a.grid(False)
    ax[1].set_xlabel('$x$ [m]'); ax[0].set_title(f'{title} (instantáneo en $t$={float(z["full_t"]):.3f} s; reflejado respecto al plano de simetría)', fontsize=9)
    fig.tight_layout(); savefig(fig, fname)

# ------------------------------------------------------------ historias de gasto
def fig_gasto(keys, fname):
    fig, ax = plt.subplots(1, len(keys), figsize=(4.2*len(keys), 3.2), squeeze=False)
    for a, key in zip(ax[0], keys):
        name, win, mesh, _ = FINAL[key]; kind = key[0]
        ei = np.loadtxt(D+f'gastoEntrada_caso_{name}.csv', delimiter=','); eo = np.loadtxt(D+f'gastoSalida_caso_{name}.csv', delimiter=',')
        ref = gasto_semi() if kind != 'S' else mdot_S
        a.plot(ei[:, 0]*1e3, -ei[:, 1], '-', color=SK, lw=.7, label='entrada'); a.plot(eo[:, 0]*1e3, eo[:, 1], '-', color=O, lw=.7, label='salida')
        a.axhline(ref, color=K, lw=1.2, ls='--', label='teoría'); a.axvspan(win[0]*1e3, win[1]*1e3, color=G, alpha=.12, label='ventana de promedio')
        a.set_xlabel('$t$ [ms]'); a.set_ylabel('$\\dot m$ (semidominio) [kg/s]'); a.set_title(key.replace('_', ' ') , fontsize=9)
        a.set_ylim(0, ref*1.25)
    ax[0, 0].legend(fontsize=7, loc='lower right'); fig.tight_layout(); savefig(fig, fname)

# ------------------------------------------------------------ convergencia por ventanas
def fig_ventanas(fname):
    wins = [(0.01, 0.02), (0.02, 0.03), (0.03, 0.04), (0.04, 0.05), (0.05, 0.06)]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    for key, ls in (('C_r1', 'o-'), ('C_r2', 's-'), ('C_r4', '^-')):
        name, _, mesh, _ = FINAL[key]
        if not exists(name): continue
        z = load('caso_'+name); tmax = z['t'].max(); xs_l = []; tw = []
        for w in wins:
            if w[1] > tmax+1e-9: continue
            x, M, p, T, mf = cfd_profile(name, w); xs_l.append(shock_pos(x, M)[0]/Ld); tw.append(0.5*(w[0]+w[1])*1e3)
        ax[0].plot(tw, xs_l, ls, color=MESH_COL[mesh], label=f'CFD {mesh}', ms=4)
    ax[0].axhline(thC['xs_Ld'], color=K, ls='--', label='teoría'); ax[0].set_xlabel('centro de la ventana [ms] (ancho 10 ms)'); ax[0].set_ylabel('$x_s/L_d$'); ax[0].legend(fontsize=8)
    for key, ls, col in (('A_r1', 'o-', B), ('S_r2', 's-', G), ('C_r1', '^-', O)):
        name = FINAL[key][0]
        if not exists(name): continue
        z = load('caso_'+name); tmax = z['t'].max(); e = []; tw = []
        ref = gasto_semi() if key[0] != 'S' else mdot_S
        for w in wins:
            if w[1] > tmax+1e-9: continue
            mi, mo = mdots(name, *w); e.append(abs(mo-ref)/ref*100); tw.append(0.5*(w[0]+w[1])*1e3)
        ax[1].plot(tw, e, ls, color=col, ms=4, label=key.replace('_', ' '))
    ax[1].set_xlabel('centro de la ventana [ms]'); ax[1].set_ylabel('error del gasto de salida [%]'); ax[1].legend(fontsize=8); ax[1].set_ylim(bottom=0)
    fig.tight_layout(); savefig(fig, fname)

# ------------------------------------------------------------ estudio de malla en C
def fig_malla_C(fname):
    rows = []
    for key in ('C_r1', 'C_r2', 'C_r4'):
        name, win, mesh, r = FINAL[key]
        if not exists(name): continue
        m = metrics('C', name, win); rows.append((0.2/(80*r), m['xs'], m['M'][0], m['p'][0], m['err_m'], mesh))
    if len(rows) < 2: return rows
    d = np.array([r[0] for r in rows]); fig, ax = plt.subplots(1, 3, figsize=(10, 3.2))
    exs = np.array([abs(r[1]-thC['xs'])/Ld*100 for r in rows])
    ax[0].loglog(d*1e3, exs, 'o-', color=B); ax[0].set_ylabel('error de posición del choque [% de $L_d$]')
    ax[1].loglog(d*1e3, [r[2] for r in rows], 's-', color=O, label='M'); ax[1].loglog(d*1e3, [r[3] for r in rows], '^-', color=G, label='$p/p_{01}$')
    ax[1].set_ylabel('error relativo medio fuera del choque [%]'); ax[1].legend()
    ax[2].loglog(d*1e3, [r[4] for r in rows], 'd-', color=K); ax[2].set_ylabel('error del gasto de salida [%]')
    for a in ax:
        a.set_xlabel('$\\Delta x$ en el divergente [mm]')
        for dd, r in zip(d, rows): pass
    fig.tight_layout(); savefig(fig, fname); return rows

# ------------------------------------------------------------ sensibilidad al paso de tiempo
SENS = [('C_r1 (CFL 0.5)', 'C_r1', (0.03, 0.06), 'C', O), ('C_r1 (CFL 0.1)', 'C_r1_co01', (0.03, 0.06), 'C', B)]
def fig_sens_dt(fname):
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.5))
    Mt, pt = theory_C(xs_fine, thC); ax[0].plot(xs_fine, Mt, '-', color=K, lw=1.5, label='teoría')
    for lab, name, win, kind, col in SENS:
        x, M, p, T, mf = cfd_profile(name, win); xs, _ = shock_pos(x, M); ax[0].plot(x, M, '.', color=col, ms=4, label=f'{lab}: $x_s/L_d$={xs/Ld:.3f}')
    ax[0].set_xlabel('$x$ [m]'); ax[0].set_ylabel('M'); ax[0].legend(fontsize=8)
    lab = ['A r1', 'C r1', 'S r1']; v05 = []; v01 = []
    for key, (a_, b_) in {'A r1': (('A_r1', (0.04, 0.06), 'A'), ('A_r1_co01', (0.04, 0.06), 'A')),
                          'C r1': (('C_r1', (0.04, 0.06), 'C'), ('C_r1_co01', (0.04, 0.06), 'C')),
                          'S r1': (('S_r1_p276_base', (0.04, 0.06), 'S'), ('S_r1_p276_co01', (0.04, 0.06), 'S'))}.items():
        v05.append(metrics(a_[2], a_[0], a_[1])['M'][0]); v01.append(metrics(b_[2], b_[0], b_[1])['M'][0])
    xx = np.arange(3); ax[1].bar(xx-.18, v05, .36, color=O, label='CFL 0.3–0.5'); ax[1].bar(xx+.18, v01, .36, color=B, label='CFL 0.1')
    ax[1].set_xticks(xx); ax[1].set_xticklabels(lab); ax[1].set_ylabel('error medio de M [%]'); ax[1].legend(fontsize=8); ax[1].grid(axis='x', alpha=0)
    fig.tight_layout(); savefig(fig, fname); return dict(zip(lab, zip(v05, v01)))

# ------------------------------------------------------------ tablas
def tex_escape(s): return s.replace('%', '\\%').replace('_', '\\_')
def write_tables(keys):
    rows = []
    hdr = ['Caso', 'Malla', 'Ventana [s]', 'Error M medio [%]', 'RMS M [%]', 'Error p/p01 medio [%]', 'Error T/T01 medio [%]',
           'mdot CFD [kg/s]', 'mdot teoría [kg/s]', 'Error mdot [%]', 'Balance masa [%]']
    for key in keys:
        name, win, mesh, r = FINAL[key]
        if not exists(name): continue
        m = metrics(key[0], name, win)
        rows.append([key[0], mesh, f'{win[0]}–{win[1]}', f'{m["M"][0]:.2f}', f'{m["M"][1]:.2f}', f'{m["p"][0]:.2f}', f'{m["T"][0]:.2f}',
                     f'{m["mo"]:.5f}', f'{m["ref"]:.5f}', f'{m["err_m"]:.3f}', f'{m["bal"]:.3f}'])
    import csv
    with open(f'{OUT_T}/tabla_errores_perfiles.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
    with open(f'{OUT_T}/tabla_errores_perfiles.tex', 'w') as f:
        f.write('\\begin{tabular}{llcrrrrrrrr}\n\\toprule\n' + ' & '.join(tex_escape(h_) for h_ in hdr) + ' \\\\\n\\midrule\n')
        for r_ in rows: f.write(' & '.join(r_) + ' \\\\\n')
        f.write('\\bottomrule\n\\end{tabular}\n')
    # estado de salida
    ex = []
    for key in keys:
        name, win, mesh, r = FINAL[key]
        if not exists(name): continue
        e = metrics(key[0], name, win)['exit']
        ex.append([key[0], mesh, f'{e["x"]:.4f}', f'{e["M"]:.4f}', f'{e["Mt"]:.4f}', f'{e["p"]/1e3:.2f}', f'{e["pt"]/1e3:.2f}', f'{e["T"]:.2f}', f'{e["Tt"]:.2f}', f'{e["V"]:.1f}', f'{e["Vt"]:.1f}',
                   f'{abs(e["M"]-e["Mt"])/e["Mt"]*100:.2f}', f'{abs(e["p"]-e["pt"])/e["pt"]*100:.2f}', f'{abs(e["T"]-e["Tt"])/e["Tt"]*100:.2f}', f'{abs(e["V"]-e["Vt"])/e["Vt"]*100:.2f}'])
    h2 = ['Caso', 'Malla', 'x [m]', 'Me CFD', 'Me teo', 'pe CFD [kPa]', 'pe teo [kPa]', 'Te CFD [K]', 'Te teo [K]', 'Ve CFD [m/s]', 'Ve teo [m/s]',
          'err Me [%]', 'err pe [%]', 'err Te [%]', 'err Ve [%]']
    with open(f'{OUT_T}/tabla_estado_salida.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(h2); w.writerows(ex)
    # choque
    ch = []
    for key in ('C_r1', 'C_r2', 'C_r4'):
        name, win, mesh, r = FINAL[key]
        if not exists(name): continue
        x, M, p, T, mf = cfd_profile(name, win); xs, l = shock_pos(x, M)
        # estados medidos fuera del espesor del choque (3 celdas a cada lado), extrapolados hasta xs con la pendiente local
        dx = np.gradient(x)[l]; xa, xb = xs-3*dx, xs+3*dx
        i1 = np.argmin(abs(x-xa)); i2 = np.argmin(abs(x-xb))
        Mt_, pt_ = theory_C(x, thC)
        Ms1, Ms2 = M[i1], M[i2]; p2p1 = p[i2]/p[i1]
        # Ms1 y Ms2 CFD se comparan con la teoria evaluada en las MISMAS posiciones (no con el valor en xs)
        Ms1t, Ms2t = Mt_[i1] if x[i1] < thC['xs'] else np.nan, Mt_[i2] if x[i2] > thC['xs'] else np.nan
        p0 = lambda pp, MM: pp*(1+0.2*MM**2)**3.5
        sig = p0(p[i2], Ms2)/p0(p[i1], Ms1)
        Asat = h(xs)/ht
        ch.append([mesh, f'{win[0]}–{win[1]}', f'{xs:.4f}', f'{thC["xs"]:.4f}', f'{xs/Ld:.4f}', f'{thC["xs_Ld"]:.4f}', f'{abs(xs-thC["xs"])/Ld*100:.2f}',
                   f'{Asat:.4f}', f'{thC["As_At"]:.4f}', f'{Ms1:.3f}', f'{Ms1t:.3f}', f'{Ms2:.3f}', f'{Ms2t:.3f}',
                   f'{p2p1:.3f}', f'{pt_[i2]/pt_[i1]:.3f}', f'{sig:.4f}', f'{thC["p02p01"]:.4f}', f'{x[i1]:.4f}', f'{x[i2]:.4f}'])
    h3 = ['Malla', 'Ventana [s]', 'xs CFD [m]', 'xs teo [m]', 'xs/Ld CFD', 'xs/Ld teo', 'err pos [% Ld]', 'As/At CFD', 'As/At teo', 'Ms1 CFD (xs-3dx)', 'Ms1 teo (misma x)',
          'Ms2 CFD (xs+3dx)', 'Ms2 teo (misma x)', 'p2/p1 CFD', 'p2/p1 teo (misma x)', 'p02/p01 CFD', 'p02/p01 teo', 'x1 [m]', 'x2 [m]']
    with open(f'{OUT_T}/tabla_choque_C.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(h3); w.writerows(ch)
    return rows, ex, ch

if __name__ == '__main__':
    import sys
    only = sys.argv[1:]
    keysA = [k for k in ('A_r1', 'A_r2') if exists(FINAL[k][0])]
    for k in keysA: fig_perfiles('A', FINAL[k][0], FINAL[k][1], f'Caso A: supersónico adaptado, malla {FINAL[k][2]}', f'perfiles_A_{FINAL[k][2]}'); fig_mapas(FINAL[k][0], 'A', f'mapas_A_{FINAL[k][2]}', f'Caso A, malla {FINAL[k][2]}')
    for k in ('S_r2', 'S_r1'):
        if exists(FINAL[k][0]):
            fig_perfiles('S', FINAL[k][0], FINAL[k][1], f'Caso S: subsónico no estrangulado, $p_B/p_{{01}}$={pB_S/p01:.3f}, malla {FINAL[k][2]}', f'perfiles_S_{FINAL[k][2]}'); fig_mapas(FINAL[k][0], 'S', f'mapas_S_{FINAL[k][2]}', f'Caso S, malla {FINAL[k][2]}')
    mc = [k for k in ('C_r1', 'C_r2', 'C_r4') if exists(FINAL[k][0])]
    if mc:
        fig_perfiles_C(mc, 'perfiles_C_mallas')
        for k in mc: fig_mapas(FINAL[k][0], 'C', f'mapas_C_{FINAL[k][2]}', f'Caso C, malla {FINAL[k][2]}')
    fig_gasto(['A_r2', 'S_r2', 'C_r2'], 'gasto_historias')
    fig_ventanas('convergencia_ventanas')
    print('malla C:', fig_malla_C('convergencia_malla_C'))
    print('sens dt:', fig_sens_dt('sensibilidad_dt'))
    r, e, c = write_tables(['A_r1', 'A_r2', 'S_r1', 'S_r2', 'C_r1', 'C_r2', 'C_r4'])
    for x in r: print(x)
    for x in c: print(x)
