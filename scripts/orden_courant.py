"""Figura orden_courant.pdf y ajuste e = e0 + a*C^p (caso S) a partir de dtstudy.json (generado por dtstudy.py)."""
import json, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
d = json.load(open('dtstudy.json'))
B, O, K = '#0072B2', '#D55E00', '#222222'
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': .25, 'axes.spines.top': False, 'axes.spines.right': False})
f = lambda C, e0, a, p: e0 + a*C**p
fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
Cs = np.array([m['co'] for m in d['S']]); xx = np.geomspace(0.05, 0.3, 80)
for key, c, lab in (('M_err', B, 'error de $M$'), ('mdot_err', O, 'error del gasto')):
    y = np.abs([m[key] for m in d['S']])
    ax[0].plot(Cs, y, 'o', color=c, label=lab)
    try:
        q, _ = curve_fit(f, Cs, y, p0=[.02, 10, 2], maxfev=20000); ax[0].plot(xx, f(xx, *q), '-', color=c, lw=1)
        print('S', key, 'e0=%.3f a=%.2f p=%.2f' % tuple(q))
    except Exception as e: print('ajuste fallo', key, e)
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('Courant máximo'); ax[0].set_ylabel('error [%]'); ax[0].set_title('Caso S (r1)'); ax[0].legend(frameon=False)
Cc = [m['co'] for m in d['C']]
ax[1].plot(Cc, [m['xs_err'] for m in d['C']], 'o-', color=B, label='posición del choque [% $L_d$]')
ax[1].plot(Cc, [m['M_err'] for m in d['C']], 's-', color=O, label='error de $M$ [%]')
ax[1].set_xlabel('Courant máximo'); ax[1].set_title('Caso C (r1)'); ax[1].legend(frameon=False, fontsize=8)
ax[2].plot([m['co'] for m in d['A']], [m['M_err'] for m in d['A']], 'o-', color=B); ax[2].set_xlabel('Courant máximo'); ax[2].set_ylabel('error de $M$ [%]'); ax[2].set_title('Caso A (r1)')
fig.tight_layout(); fig.savefig('orden_courant_nueva.pdf', bbox_inches='tight')
