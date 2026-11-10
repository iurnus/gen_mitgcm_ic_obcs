#!/usr/bin/env python3
"""
Build MITgcm initial conditions and open-boundary conditions from GLORYS (CMEMS) NetCDF files.

Usage:   python3 make_ic_obc.py config_file.py
         (the config is a small Python file; see config_redsea.py / config_pacific.py)

Handles
  * domains that cross the date line (or any longitude convention: -180..180, 0..360,
    model xgOrigin negative or > 180). GLORYS longitudes are unwrapped onto the model's
    longitude range and sorted, so a file that comes back split at +-180 is fine.
  * GLORYS downloaded as several files per date (e.g. one west and one east of the
    date line): files with the same time stamp are stitched together along longitude.
  * any combination of open boundaries N, S, E, W along the domain edges.
  * uniform or variable delX / delY (spherical polar grid).

Method
  1. Crop GLORYS to the model domain + MARGIN degrees.
  2. Vertical linear interpolation of each GLORYS column onto model RC.
  3. NaN-aware bilinear interpolation onto model C/U/V points.
  4. Gap filling at each model level, top to bottom: lateral diffusion-fill only through
     model ocean cells at that level; anything left takes the value from the level above;
     land: T/S/Eta nearest-filled, U/V = 0.
  5. Optional: uniform correction to the boundary-normal velocity so the net volume
     transport into the domain is zero at every record.
  6. Eta: domain-mean of the initial zos removed from IC and OBs (arbitrary reference).

Output (big-endian float32 by default):
  T_ic.bin S_ic.bin U_ic.bin V_ic.bin (Nr,Ny,Nx), Eta_ic.bin (Ny,Nx)
  OB<side>_T/S/U/V.bin  (nt,Nr,Ny) for E/W,  (nt,Nr,Nx) for N/S;   OB<side>_Eta.bin
  namelist_snippets.txt with the matching data / data.obcs / data.exf lines.
"""
import sys, glob, os, re, datetime as dt
import numpy as np
from scipy import ndimage
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h5read import H5File

# ------------------------------------------------------------------ config
DEFAULTS = dict(
    OUT='out', OPEN_BOUNDARIES='E', HFACMIN=1.0, HFACMINDR=0.0, ZERO_NET_TRANSPORT=True,
    MARGIN=0.5, BATHY_DTYPE='>f4', OUT_DTYPE='>f4', Rearth=6371.0e3,
    VARS=dict(lon='longitude', lat='latitude', depth='depth', time='time',
              T='thetao', S='so', U='uo', V='vo', Eta='zos'),
)
if len(sys.argv) < 2:
    sys.exit(__doc__)
cfg = dict(DEFAULTS)
_ns = {}
exec(open(sys.argv[1]).read(), {'np': np}, _ns)
cfg.update({k: v for k, v in _ns.items() if not k.startswith('_')})
for k in ['GLORYS_FILES', 'BATHY', 'Nx', 'Ny', 'xgOrigin', 'ygOrigin', 'delX', 'delY', 'delR']:
    if k not in cfg: sys.exit(f'config is missing {k}')
VN = dict(DEFAULTS['VARS']); VN.update(cfg['VARS'])
Nx, Ny = int(cfg['Nx']), int(cfg['Ny'])
OUT = cfg['OUT']; Rearth = cfg['Rearth']
SIDES = [s for s in 'NSEW' if s in cfg['OPEN_BOUNDARIES'].upper()]

# -------------------------------------------------------------------- grid
delR = np.asarray(cfg['delR'], float); Nr = delR.size
dX = np.broadcast_to(np.asarray(cfg['delX'], float), (Nx,)).copy()
dY = np.broadcast_to(np.asarray(cfg['delY'], float), (Ny,)).copy()
XG = cfg['xgOrigin'] + np.concatenate([[0], np.cumsum(dX)[:-1]]); XC = XG + 0.5 * dX
YG = cfg['ygOrigin'] + np.concatenate([[0], np.cumsum(dY)[:-1]]); YC = YG + 0.5 * dY
LX = dX.sum()
if LX + 2 * cfg['MARGIN'] >= 360:
    sys.exit('zonally periodic/global domains are not supported by this script')
rF = np.concatenate([[0.0], np.cumsum(delR)])
RC = 0.5 * (rF[:-1] + rF[1:])
LON0 = cfg['xgOrigin'] - cfg['MARGIN'] - 1e-3     # unwrap all longitudes to [LON0, LON0+360)

def unwrap(lon):
    return (np.asarray(lon, float) - LON0) % 360.0 + LON0

files = sorted(set(sum([sorted(glob.glob(p)) for p in np.atleast_1d(cfg['GLORYS_FILES'])], [])))
if not files: sys.exit('no GLORYS files match ' + str(cfg['GLORYS_FILES']))

depth = -np.fromfile(cfg['BATHY'], cfg['BATHY_DTYPE']).astype('f8')
if depth.size != Nx * Ny: sys.exit(f'bathymetry has {depth.size} values, expected {Nx*Ny}')
depth = depth.reshape(Ny, Nx); depth[depth < 0] = 0
print(f'grid {Nx}x{Ny}x{Nr}: lon {XG[0]:.3f}..{XG[0]+LX:.3f}  lat {YG[0]:.3f}..{YG[-1]+dY[-1]:.3f}'
      f'  open: {"".join(SIDES) or "none"}   max depth {depth.max():.0f} / grid {rF[-1]:.0f} m')

def hfac(d):
    h = np.clip((d[None] - rF[:-1, None, None]) / delR[:, None, None], 0, 1)
    mn = np.maximum(cfg['HFACMIN'], np.minimum(cfg['HFACMINDR'] / delR, 1.0))[:, None, None]
    return np.where(h < mn, np.where(h < 0.5 * mn, 0.0, mn), h)

hFacC = hfac(depth)
hFacW = np.minimum(hFacC, np.roll(hFacC, 1, axis=2)); hFacW[:, :, 0] = 0
hFacS = np.minimum(hFacC, np.roll(hFacC, 1, axis=1)); hFacS[:, 0, :] = 0
depthW = np.minimum(depth, np.roll(depth, 1, axis=1)); depthW[:, 0] = 0
depthS = np.minimum(depth, np.roll(depth, 1, axis=0)); depthS[0, :] = 0
wetC = depth[None] > rF[:-1, None, None]          # permissive masks for filling
wetW = depthW[None] > rF[:-1, None, None]
wetS = depthS[None] > rF[:-1, None, None]

# ------------------------------------------------------------- read GLORYS
def parse_time(val, units):
    m = re.match(r'\s*(\w+)\s+since\s+([0-9-]+)(?:[ T]([0-9:.]+))?', units or '')
    if not m:
        print('  WARNING: no time units found, assuming hours since 1950-01-01')
        m = re.match(r'(\w+) since ([0-9-]+)', 'hours since 1950-01-01')
    unit, day = m.group(1).lower(), m.group(2)
    hms = (m.group(3) or '00:00:00').split(':') + ['0', '0']
    ref = dt.datetime.strptime(day, '%Y-%m-%d') + dt.timedelta(
        hours=float(hms[0]), minutes=float(hms[1]), seconds=float(hms[2]))
    sec = {'seconds': 1, 'second': 1, 'minutes': 60, 'minute': 60, 'hours': 3600,
           'hour': 3600, 'days': 86400, 'day': 86400}[unit]
    return [ref + dt.timedelta(seconds=float(v) * sec) for v in np.atleast_1d(val)]

_eps = 1e-6
lat_lo, lat_hi = YG[0] - cfg['MARGIN'] - _eps, YG[-1] + dY[-1] + cfg['MARGIN'] + _eps
lon_lo, lon_hi = XG[0] - cfg['MARGIN'] - _eps, XG[0] + LX + cfg['MARGIN'] + _eps

# index every file: its times, the longitude pieces it contributes
pieces = {}                      # datetime -> list of (file, time index)
grid_ref = None
for fn in files:
    f = H5File(fn)
    tv, ta = f.read(VN['time'])
    lat = f.read(VN['lat'])[0]; dep = f.read(VN['depth'])[0]
    if grid_ref is None: grid_ref = (lat, dep)
    elif lat.shape != grid_ref[0].shape or not np.allclose(lat, grid_ref[0]) \
            or not np.allclose(dep, grid_ref[1]):
        sys.exit(f'{fn}: latitude/depth grid differs from the other files')
    for it, t in enumerate(parse_time(tv, ta['units'])):
        pieces.setdefault(t, []).append((fn, it))
    f.close()
glat_all, gdep = grid_ref
iy = np.where((glat_all >= lat_lo) & (glat_all <= lat_hi))[0]
glat = glat_all[iy]
if glat.size < 2 or glat_all.min() > YC.min() or glat_all.max() < YC.max():
    sys.exit(f'GLORYS latitudes {glat_all.min():.2f}..{glat_all.max():.2f} do not cover the model')

def read_snapshot(t):
    """Stitch all pieces for time t; return dict of cropped arrays on sorted unwrapped lon."""
    lons, data = [], {v: [] for v in 'TSUVE'}
    for fn, it in pieces[t]:
        f = H5File(fn)
        lon = unwrap(f.read(VN['lon'])[0])
        ix = np.where((lon >= lon_lo) & (lon <= lon_hi))[0]
        if ix.size:
            lons.append(lon[ix])
            for v, name in zip('TSUV', ('T', 'S', 'U', 'V')):
                data[v].append(f.read(VN[name], dtype='f4')[0][it][:, iy][:, :, ix])
            data['E'].append(f.read(VN['Eta'], dtype='f4')[0][it][iy][:, ix])
        f.close()
    if not lons: sys.exit(f'no GLORYS longitudes inside the domain for {t}')
    lon = np.concatenate(lons)
    lon, keep = np.unique(lon, return_index=True)        # sorted, duplicates dropped
    out = {v: np.concatenate(data[v], axis=-1)[..., keep] for v in data}
    out['lon'] = lon; out['date'] = t
    return out

# ---------------------------------------------------- interpolation tools
kk = np.searchsorted(gdep, RC)
k1 = np.clip(kk, 0, gdep.size - 1); k0 = np.clip(kk - 1, 0, gdep.size - 1)
wz = np.where(k1 > k0, (RC - gdep[k0]) / (gdep[k1] - gdep[k0] + 1e-30), 0.0)

def vinterp(G):
    a, b = G[k0], G[k1]; w = wz[:, None, None].astype(G.dtype)
    out = (1 - w) * a + w * b
    return np.where(np.isnan(a), b, np.where(np.isnan(b), a, out))

def frac_index(p, g):
    """fractional index of points p in monotonic coordinate g (NaN outside)."""
    return np.interp(p, g, np.arange(g.size), left=np.nan, right=np.nan)

def hinterp(F, glon, px, py):
    """NaN-aware bilinear interp of F(..., nlat, nlon) onto the (py, px) mesh."""
    fx = frac_index(px, glon); fy = frac_index(py, glat)
    if np.isnan(fx).any(): sys.exit(f'GLORYS longitudes {glon[0]:.2f}..{glon[-1]:.2f} '
                                    f'(unwrapped) do not cover model {px.min():.2f}..{px.max():.2f}')
    i0 = np.clip(np.floor(fx).astype(int), 0, glon.size - 2); ax = (fx - i0)[None, :]
    j0 = np.clip(np.floor(fy).astype(int), 0, glat.size - 2); ay = (fy - j0)[:, None]
    J, I = j0[:, None], i0[None, :]
    num = 0; den = 0
    for dj, wj in ((0, 1 - ay), (1, ay)):
        for di, wi in ((0, 1 - ax), (1, ax)):
            v = F[..., J + dj, I + di]; w = wi * wj * np.isfinite(v)
            num = num + np.nan_to_num(v) * w; den = den + w
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(den > 1e-12, num / den, np.nan)

def fill_connected(a, mask):
    a = a.copy()
    while True:
        miss = np.isnan(a) & mask
        if not miss.any(): break
        v = np.isfinite(a); a0 = np.nan_to_num(a)
        s = np.zeros_like(a); n = np.zeros_like(a)
        for sh, ax in ((1, 0), (-1, 0), (1, 1), (-1, 1)):
            vv = np.roll(v, sh, ax); aa = np.roll(a0, sh, ax)
            if ax == 0: vv[0 if sh == 1 else -1, :] = False
            else:       vv[:, 0 if sh == 1 else -1] = False
            s += aa * vv; n += vv
        new = miss & (n > 0)
        if not new.any(): break
        a[new] = s[new] / n[new]
    return a

def nearest_fill(a):
    bad = np.isnan(a)
    if bad.all() or not bad.any(): return a
    idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
    return a[tuple(idx)]

def fill3d(A, wet, land_value=None):
    A = A.copy(); stats = [0, 0]
    for k in range(Nr):
        before = np.isnan(A[k]) & wet[k]
        A[k] = fill_connected(A[k], wet[k])
        rem = np.isnan(A[k]) & wet[k]
        stats[0] += (before & ~rem).sum()
        if rem.any():
            if k == 0: A[k] = nearest_fill(A[k])
            else: A[k][rem] = A[k - 1][rem]
            stats[1] += rem.sum()
        if land_value is None: A[k] = nearest_fill(A[k])
        else: A[k][~wet[k]] = land_value
    return A, stats

XCu, XGu = unwrap(XC), unwrap(XG)     # model lon on the same unwrapped axis

def to_model(d, verbose=False):
    out = {}
    for v, px, py, wet, lv in (('T', XCu, YC, wetC, None), ('S', XCu, YC, wetC, None),
                               ('U', XGu, YC, wetW, 0.0), ('V', XCu, YG, wetS, 0.0)):
        out[v], st = fill3d(hinterp(vinterp(d[v]), d['lon'], px, py), wet, lv)
        if verbose: print(f'  {v}: {st[0]} wet cells filled laterally, {st[1]} vertically')
    out['E'] = nearest_fill(fill_connected(hinterp(d['E'], d['lon'], XCu, YC), wetC[0]))
    return out

def w(fn, a):
    a.astype(cfg['OUT_DTYPE']).tofile(os.path.join(OUT, fn))

# ----------------------------------------------------------------- process
os.makedirs(OUT, exist_ok=True)
times = sorted(pieces)
print('GLORYS records:', ', '.join(f'{t:%Y-%m-%d %H:%M} ({len(pieces[t])} file(s))' for t in times))
dts = np.diff(times)
if len(dts) and any(x != dts[0] for x in dts):
    sys.exit('GLORYS records are not evenly spaced: ' + str(sorted(set(dts))))
period = dts[0].total_seconds() if len(dts) else 0.0
print('record spacing (s):', period)

models = []
for i, t in enumerate(times):
    s = read_snapshot(t)
    gap = np.diff(s['lon']).max() / np.median(np.diff(s['lon']))
    if gap > 1.5: print(f'  WARNING: {t:%Y-%m-%d} has a longitude gap ({gap:.1f}x the grid spacing)'
                        ' -- is a file piece missing?')
    print(f'interpolating {t:%Y-%m-%d %H:%M}  (GLORYS lon {s["lon"][0]:.3f}..{s["lon"][-1]:.3f}'
          f', {s["lon"].size} pts)')
    models.append(to_model(s, verbose=(i == 0)))
    del s

eta_ref = np.sum(models[0]['E'] * wetC[0]) / max(wetC[0].sum(), 1)
print(f'removing domain-mean zos {eta_ref:.4f} m from Eta (IC and OB)')

m = models[0]
w('T_ic.bin', m['T']); w('S_ic.bin', m['S']); w('U_ic.bin', m['U']); w('V_ic.bin', m['V'])
w('Eta_ic.bin', (m['E'] - eta_ref) * wetC[0])

# ------------------------------------------------------------- open boundaries
# MITgcm OBCS conventions (0-based indices here):
#   East : T at i=Nx-1, normal U at face i=Nx-1 (OB_Ie),  V at i=Nx-1
#   West : T at i=0,    normal U at face i=1    (OB_Iw+1), V at i=0
#   North: T at j=Ny-1, normal V at face j=Ny-1 (OB_Jn),  U at j=Ny-1
#   South: T at j=0,    normal V at face j=1    (OB_Js+1), U at j=0
dyG = Rearth * np.deg2rad(dY)                                 # (Ny,)  length of U faces
dxG = lambda j: Rearth * np.cos(np.deg2rad(YG[j])) * np.deg2rad(dX)   # (Nx,) V faces at row j
along = {s: np.ones(Nx if s in 'NS' else Ny, bool) for s in SIDES}
for s in SIDES:      # exclude faces lying along another open boundary (corner cells)
    if s in 'EW':
        if 'N' in SIDES: along[s][-1] = False
        if 'S' in SIDES: along[s][0] = False
    else:
        if 'E' in SIDES: along[s][-1] = False
        if 'W' in SIDES: along[s][0] = False
spec = {  # side: (slice for T, normal var, normal face slice, sign(+ = outflow), area, tangential)
    'E': (np.s_[:, :, Nx-1], 'U', np.s_[:, :, Nx-1], +1, hFacW[:, :, Nx-1] * delR[:, None] * dyG, 'V', np.s_[:, :, Nx-1]),
    'W': (np.s_[:, :, 0],    'U', np.s_[:, :, 1],    -1, hFacW[:, :, 1]    * delR[:, None] * dyG, 'V', np.s_[:, :, 0]),
    'N': (np.s_[:, Ny-1, :], 'V', np.s_[:, Ny-1, :], +1, hFacS[:, Ny-1, :] * delR[:, None] * dxG(Ny-1), 'U', np.s_[:, Ny-1, :]),
    'S': (np.s_[:, 0, :],    'V', np.s_[:, 1, :],    -1, hFacS[:, 1, :]    * delR[:, None] * dxG(1), 'U', np.s_[:, 0, :]),
}
ob = {s: {k: [] for k in 'TSUVE'} for s in SIDES}
if SIDES:
    A_tot = sum((spec[s][4] * along[s][None]).sum() for s in SIDES)
    print('net transport out of the domain through open boundaries (Sv):')
for t, m in zip(times, models):
    normals = {}; Q = 0.0
    for s in SIDES:
        _, nv, nsl, sign, area, _, _ = spec[s]
        un = m[nv][nsl] * (area > 0)
        normals[s] = un; Q += sign * (un * area * along[s][None]).sum()
    msg = f'  {t:%Y-%m-%d %H:%M}: raw {Q/1e6:+.3f}'
    if cfg['ZERO_NET_TRANSPORT'] and SIDES and A_tot > 0:
        du = Q / A_tot                         # outflow velocity to remove
        for s in SIDES:
            _, _, _, sign, area, _, _ = spec[s]
            normals[s] = np.where((area > 0) & along[s][None], normals[s] - sign * du, normals[s])
        Q2 = sum(spec[s][3] * (normals[s] * spec[s][4] * along[s][None]).sum() for s in SIDES)
        msg += f'  ->  corrected {Q2/1e6:+.4f}'
    if SIDES: print(msg)
    for s in SIDES:
        tsl, nv, _, _, _, tv, tgsl = spec[s]
        ob[s]['T'].append(m['T'][tsl]); ob[s]['S'].append(m['S'][tsl])
        ob[s][nv].append(normals[s])
        wt = (wetW if tv == 'U' else wetS)[tgsl]
        ob[s][tv].append(m[tv][tgsl] * wt)
        ob[s]['E'].append((m['E'] - eta_ref)[tsl[1:]] * wetC[0][tsl[1:]])
for s in SIDES:
    for v, nm in zip('TSUVE', ('T', 'S', 'U', 'V', 'Eta')):
        w(f'OB{s}_{nm}.bin', np.array(ob[s][v]))
    nwet = (hFacC[0][spec[s][0][1:]] > 0).sum()
    print(f'  OB{s}: {nwet} wet surface points on the boundary')

# ------------------------------------------------------------ namelist help
name = {'N': 'OB_Jnorth', 'S': 'OB_Jsouth', 'E': 'OB_Ieast', 'W': 'OB_Iwest'}
val = {'N': f'{Nx}*-1', 'S': f'{Nx}*1', 'E': f'{Ny}*-1', 'W': f'{Ny}*1'}
t0 = times[0]
lines = ["# ---- input/data, &PARM05 ----",
         " hydrogThetaFile = 'T_ic.bin',", " hydrogSaltFile  = 'S_ic.bin',",
         " uVelInitFile    = 'U_ic.bin',", " vVelInitFile    = 'V_ic.bin',",
         " pSurfInitFile   = 'Eta_ic.bin',", "",
         "# ---- input/data.obcs ----", " &OBCS_PARM01"]
lines += [f" {name[s]} = {val[s]}," for s in SIDES]
lines += [" useOBCSprescribe = .TRUE.,"]
for s in SIDES:
    lines += [f" OB{s}tFile = 'OB{s}_T.bin',", f" OB{s}sFile = 'OB{s}_S.bin',",
              f" OB{s}uFile = 'OB{s}_U.bin',", f" OB{s}vFile = 'OB{s}_V.bin',",
              f"# OB{s}etaFile = 'OB{s}_Eta.bin',   (only with NONLIN_FRSURF)"]
lines += [" &", "", "# ---- input/data.cal ----",
          f" &CAL_NML\n TheCalendar = 'gregorian',\n startDate_1 = {t0:%Y%m%d},\n startDate_2 = {t0:%H%M%S},\n &",
          "", "# ---- input/data.exf ----",
          " &EXF_NML_01\n &\n &EXF_NML_02\n &\n &EXF_NML_03\n &\n &EXF_NML_04\n &\n &EXF_NML_OBCS"]
for s in SIDES:
    lines += [f" obcs{s}startdate1 = {t0:%Y%m%d},", f" obcs{s}startdate2 = {t0:%H%M%S},",
              f" obcs{s}period     = {period:.1f},"]
lines += [" &", "", f"# records: {len(times)}, {times[0]:%Y-%m-%d %H:%M} .. {times[-1]:%Y-%m-%d %H:%M}"
          f"  -> max run length {(times[-1]-times[0]).total_seconds():.0f} s"]
open(os.path.join(OUT, 'namelist_snippets.txt'), 'w').write('\n'.join(lines) + '\n')
np.savez(os.path.join(OUT, 'meta.npz'), XC=XC, YC=YC, XG=XG, YG=YG, RC=RC, depth=depth,
         hFacC=hFacC, period=period, sides=''.join(SIDES),
         dates=np.array([t.isoformat() for t in times]))
print('wrote files to', OUT, '(see namelist_snippets.txt)')
