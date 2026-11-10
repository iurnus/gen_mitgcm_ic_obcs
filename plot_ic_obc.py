#!/usr/bin/env python3
"""
Quick-look figures and sanity checks for the files written by make_ic_obc.py.

Usage:   python3 plot_ic_obc.py config_file.py [--record N] [--outdir DIR]

Reads the same config file as make_ic_obc.py, plus <OUT>/meta.npz written by it, and makes
  <OUT>/check_IC.png   initial conditions: SST, SSS, Eta, T at ~1000 m, surface speed, T section
  <OUT>/check_OBC.png  open boundaries:
                         one open side  -> T, S, normal velocity for (up to) 4 records
                         several sides  -> T, S, normal velocity on every side for one record
It also prints NaN counts and value ranges for every file.

Optional config entry (ignored by make_ic_obc.py):
  PLOT_SECTION = 'auto'          zonal or meridional section with the most ocean (default)
               = ('lat', -15.0)  zonal section at this latitude
               = ('lon', 180.0)  meridional section at this longitude (any convention)
               = 'deepest'       follow the deepest point of each row (long narrow basins)
"""
import sys, os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('config')
ap.add_argument('--record', type=int, default=0, help='OB record shown when several sides are open')
ap.add_argument('--outdir', default=None, help='where to save figures (default: OUT from config)')
args = ap.parse_args()

cfg = {'OUT': 'out', 'OUT_DTYPE': '>f4', 'PLOT_SECTION': 'auto'}
ns = {}
exec(open(args.config).read(), {'np': np}, ns)
cfg.update({k: v for k, v in ns.items() if not k.startswith('_')})
OUT = cfg['OUT']; DT = cfg['OUT_DTYPE']
FIG = args.outdir or OUT
os.makedirs(FIG, exist_ok=True)
if not os.path.exists(os.path.join(OUT, 'meta.npz')):
    sys.exit(f'{OUT}/meta.npz not found -- run make_ic_obc.py {args.config} first')

M = np.load(os.path.join(OUT, 'meta.npz'))
XC, YC, XG, YG, RC, depth = M['XC'], M['YC'], M['XG'], M['YG'], M['RC'], M['depth']
sides = str(M['sides']); dates = [str(d)[:16].replace('T', ' ') for d in M['dates']]
Nr, Ny, Nx, nt = RC.size, YC.size, XC.size, len(dates)
delR = np.asarray(cfg['delR'], float); rF = np.concatenate([[0.0], np.cumsum(delR)])
wetC = depth[None] > rF[:-1, None, None]
dW = np.minimum(depth, np.roll(depth, 1, 1)); dW[:, 0] = 0
dS = np.minimum(depth, np.roll(depth, 1, 0)); dS[0, :] = 0
wetW = dW[None] > rF[:-1, None, None]; wetS = dS[None] > rF[:-1, None, None]

def rd(name, shape):
    a = np.fromfile(os.path.join(OUT, name), DT)
    if a.size != np.prod(shape): sys.exit(f'{name}: {a.size} values, expected shape {shape}')
    return a.reshape(shape).astype('f8')

# ------------------------------------------------------------------ labels
def lon_label(v):
    v = (v + 180) % 360 - 180
    if abs(v) < 1e-6 or abs(abs(v) - 180) < 1e-6: return f'{abs(v):.0f}°'
    return f'{abs(v):g}°{"E" if v > 0 else "W"}'
def lat_label(v):
    return '0°' if abs(v) < 1e-6 else f'{abs(v):g}°{"N" if v > 0 else "S"}'
def set_lon_ticks(ax, lo, hi):
    tk = matplotlib.ticker.MaxNLocator(6).tick_values(lo, hi); tk = tk[(tk >= lo) & (tk <= hi)]
    ax.set_xticks(tk); ax.set_xticklabels([lon_label(t) for t in tk])
def set_lat_ticks(ax, lo, hi, axis='y'):
    tk = matplotlib.ticker.MaxNLocator(6).tick_values(lo, hi); tk = tk[(tk >= lo) & (tk <= hi)]
    (ax.set_yticks if axis == 'y' else ax.set_xticks)(tk)
    (ax.set_yticklabels if axis == 'y' else ax.set_xticklabels)([lat_label(t) for t in tk])

def robust(a, mask, lo=0.5, hi=99.5):
    v = a[mask]; v = v[np.isfinite(v)]
    return (None, None) if v.size == 0 else tuple(np.percentile(v, [lo, hi]))

# ------------------------------------------------------------------ IC
T = rd('T_ic.bin', (Nr, Ny, Nx)); S = rd('S_ic.bin', (Nr, Ny, Nx))
U = rd('U_ic.bin', (Nr, Ny, Nx)); V = rd('V_ic.bin', (Nr, Ny, Nx)); E = rd('Eta_ic.bin', (Ny, Nx))
print(f'{OUT}: grid {Nx}x{Ny}x{Nr}, open boundaries: {sides or "none"}, records: {nt}')
for n, a, wm in (('T_ic', T, wetC), ('S_ic', S, wetC), ('U_ic', U, wetW), ('V_ic', V, wetS)):
    print(f'  {n:7s} NaN {int(np.isnan(a).sum()):6d}   ocean range {np.nanmin(a[wm]):9.3f} .. {np.nanmax(a[wm]):9.3f}')
print(f'  Eta_ic  NaN {int(np.isnan(E).sum()):6d}   ocean range {E[wetC[0]].min():9.3f} .. {E[wetC[0]].max():9.3f}')

fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
def mapplot(a, mask, title, cmap, ax_, vmin=None, vmax=None):
    if vmin is None: vmin, vmax = robust(a, mask)
    im = ax_.pcolormesh(XC, YC, np.ma.masked_where(~mask, a), cmap=cmap, vmin=vmin, vmax=vmax,
                        shading='auto')
    plt.colorbar(im, ax=ax_, shrink=.85); ax_.set_title(title)
    ax_.set_aspect(1 / np.cos(np.deg2rad(YC.mean())))
    set_lon_ticks(ax_, XG[0], XG[-1] + (XG[-1] - XG[-2])); set_lat_ticks(ax_, YG[0], YC[-1])
d0 = dates[0][:10]
mapplot(T[0], wetC[0], f'Initial SST (°C), {d0}', 'turbo', ax[0, 0])
mapplot(S[0], wetC[0], 'Initial SSS', 'viridis', ax[0, 1])
el = max(abs(np.array(robust(E, wetC[0], 1, 99))))
mapplot(E, wetC[0], 'Initial Eta (m)', 'RdBu_r', ax[0, 2], -el, el)
k = int(np.argmin(abs(RC - min(1000, depth.max() * 0.6))))
mapplot(T[k], wetC[k], f'Initial T at {RC[k]:.0f} m (°C)', 'turbo', ax[1, 0])
# speed at C points
uc = 0.5 * (U[0] + np.roll(U[0], -1, 1)); vc = 0.5 * (V[0] + np.roll(V[0], -1, 0))
sp = np.hypot(uc, vc)
mapplot(sp, wetC[0], 'Initial surface speed (m/s)', 'magma', ax[1, 1], 0, robust(sp, wetC[0], 0, 99.5)[1])

# section
sec = cfg['PLOT_SECTION']; a3 = ax[1, 2]
if sec == 'auto':
    jz = int(np.argmax(wetC.sum(axis=(0, 2)))); im_ = int(np.argmax(wetC.sum(axis=(0, 1))))
    sec = ('lat', YC[jz]) if wetC[:, jz, :].sum() >= wetC[:, :, im_].sum() else ('lon', XC[im_])
if sec == 'deepest':
    rows = np.where(depth.max(1) > 0)[0]; ii = np.argmax(depth, axis=1)
    data = np.array([T[:, j, ii[j]] for j in rows]).T; msk = np.array([wetC[:, j, ii[j]] for j in rows]).T
    xs = YC[rows]; title = 'Initial T along the deepest point of each row'
    set_lat_ticks(a3, xs[0], xs[-1], axis='x')
elif sec[0] == 'lat':
    j = int(np.argmin(abs(YC - sec[1]))); data, msk, xs = T[:, j, :], wetC[:, j, :], XC
    title = f'Initial T along {lat_label(round(YC[j], 2))}'; set_lon_ticks(a3, XC[0], XC[-1])
else:
    lonu = (sec[1] - XG[0]) % 360 + XG[0]
    i = int(np.argmin(abs(XC - lonu))); data, msk, xs = T[:, :, i], wetC[:, :, i], YC
    title = f'Initial T along {lon_label(round(XC[i], 2))}'; set_lat_ticks(a3, YC[0], YC[-1], axis='x')
im = a3.pcolormesh(xs, -RC, np.ma.masked_where(~msk, data), cmap='turbo', shading='auto')
plt.colorbar(im, ax=a3, shrink=.85); a3.set_title(title); a3.set_ylabel('depth (m)')
a3.set_ylim(-rF[-1], 0)
plt.tight_layout(); f1 = os.path.join(FIG, 'check_IC.png'); plt.savefig(f1, dpi=80); plt.close()

# ------------------------------------------------------------------ OBs
if sides:
    info = {'E': ('U', wetC[:, :, -1], wetW[:, :, -1], YC, 'lat'),
            'W': ('U', wetC[:, :, 0], wetW[:, :, 1], YC, 'lat'),
            'N': ('V', wetC[:, -1, :], wetS[:, -1, :], XC, 'lon'),
            'S': ('V', wetC[:, 0, :], wetS[:, 1, :], XC, 'lon')}
    ob = {}
    for s in sides:
        n = Ny if s in 'EW' else Nx
        ob[s] = {v: rd(f'OB{s}_{v}.bin', (nt, Nr, n)) for v in 'TSUV'}
        bad = sum(int(np.isnan(ob[s][v]).sum()) for v in 'TSUV')
        mT = info[s][1]
        print(f'  OB{s}: NaN {bad}, T {ob[s]["T"][:, mT].min():.2f}..{ob[s]["T"][:, mT].max():.2f}, '
              f'S {ob[s]["S"][:, mT].min():.2f}..{ob[s]["S"][:, mT].max():.2f}, '
              f'normal velocity {ob[s][info[s][0]].min():+.2f}..{ob[s][info[s][0]].max():+.2f} m/s')
    if len(sides) == 1:
        s = sides; recs = list(np.unique(np.linspace(0, nt - 1, min(nt, 4)).round().astype(int)))
        panels = [(s, r) for r in recs]
    else:
        r = min(max(args.record, 0), nt - 1); panels = [(s, r) for s in sides]
    # common colour limits across panels
    lims = {}
    for row, var in enumerate(('T', 'S', 'N')):
        vals = []
        for s, r in panels:
            nv, mT, mN, _, _ = info[s]
            a = ob[s][nv if var == 'N' else var][r]; m = mN if var == 'N' else mT
            vals.append(a[m])
        v = np.concatenate(vals)
        if var == 'N':
            L = np.percentile(abs(v), 99.5) if v.size else 1; lims[var] = (-L, L)
        else:
            lims[var] = tuple(np.percentile(v, [0.5, 99.5])) if v.size else (None, None)
    nc = len(panels)
    fig, ax = plt.subplots(3, nc, figsize=(4.5 * nc, 10), squeeze=False)
    for c, (s, r) in enumerate(panels):
        nv, mT, mN, x, kind = info[s]
        for row, (var, cmap, lab) in enumerate((('T', 'turbo', 'T (°C)'), ('S', 'viridis', 'S'),
                                                ('N', 'RdBu_r', f'{nv} normal (m/s)'))):
            a = ob[s][nv if var == 'N' else var][r]; m = mN if var == 'N' else mT
            im = ax[row, c].pcolormesh(x, -RC, np.ma.masked_where(~m, a), cmap=cmap,
                                       vmin=lims[var][0], vmax=lims[var][1], shading='auto')
            plt.colorbar(im, ax=ax[row, c])
            ax[row, c].set_title(f'OB{s} {lab}   {dates[r][:10]}', fontsize=10)
            wc = np.where(mT.any(0))[0]                    # zoom to the ocean part of the edge
            x0, x1 = (x[wc[0]], x[wc[-1]]) if wc.size else (x[0], x[-1])
            pad = 0.03 * (x1 - x0) + 0.5 * abs(x[1] - x[0])
            ax[row, c].set_xlim(x0 - pad, x1 + pad)
            edge = {'E': depth[:, -1], 'W': depth[:, 0], 'N': depth[-1], 'S': depth[0]}[s]
            zmax = max(min(edge.max(), rF[-1]), RC[0])     # deepest point on this edge
            ax[row, c].set_ylim(-1.03 * zmax, 0)
            if kind == 'lon': set_lon_ticks(ax[row, c], x0, x1)
            else: set_lat_ticks(ax[row, c], x0, x1, axis='x')
        ax[0, c].set_ylabel('depth (m)') if c == 0 else None
        ax[2, c].set_xlabel({'E': f'latitude, east boundary ({lon_label(round(XC[-1], 2))})',
                             'W': f'latitude, west boundary ({lon_label(round(XC[0], 2))})',
                             'N': f'longitude, north boundary ({lat_label(round(YC[-1], 2))})',
                             'S': f'longitude, south boundary ({lat_label(round(YC[0], 2))})'}[s])
    plt.tight_layout(); f2 = os.path.join(FIG, 'check_OBC.png'); plt.savefig(f2, dpi=80); plt.close()
    print('wrote', f1, 'and', f2)
else:
    print('wrote', f1, '(no open boundaries)')
