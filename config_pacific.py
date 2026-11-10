# Pacific, 170E - 170W, 25S - 5S, 0.08 deg, all four sides open -- an example that crosses the date line.
# Longitudes: use the same xgOrigin as in your input/data. 170 (east of Greenwich) and
# -190 describe the same place; the script handles either, and any GLORYS longitude convention.

GLORYS_FILES = 'glorys_pacific/*.nc'   # one file per period, or east/west halves -- both fine
BATHY = 'bathy_pacific.bin'            # your model bathymetry (negative = ocean), Ny x Nx
OUT = 'out_pacific'

Nx, Ny = 250, 250                       # from code/SIZE.h
xgOrigin, ygOrigin = 170.0, -25.0        # &PARM04 xgOrigin / ygOrigin
delX, delY = 0.08, 0.08                 # scalar, or a list of Nx / Ny values
delR = [5.0, 5.1, 5.3, 5.7, 6.2, 6.9, 7.7, 8.7, 9.8, 11.0, 12.5, 14.0, 15.7, 17.6, 19.6, 21.8,
        24.1, 26.5, 29.1, 31.9, 34.8, 37.9, 41.1, 44.4, 47.9, 51.6, 55.4, 59.3, 63.4, 67.7,
        72.1, 76.6, 81.3, 86.2, 91.2, 96.3, 101.6, 107.0, 112.6, 118.4, 124.2, 130.3, 136.5,
        142.8, 149.3, 155.9, 162.7, 169.6, 176.7, 183.9]

OPEN_BOUNDARIES = 'NSEW'                # any combination of N, S, E, W (full-length edges)
HFACMIN, HFACMINDR = 1.0, 0.0           # match &PARM01 in input/data
ZERO_NET_TRANSPORT = True               # uniform correction so net inflow = 0 at each record

# Optional:
# MARGIN = 0.5                          # degrees of GLORYS kept around the domain
# BATHY_DTYPE = '>f4'                   # '>f8' if your bathymetry is real*8
# OUT_DTYPE = '>f4'                     # '>f8' with readBinaryPrec = 64

# plot_ic_obc.py only: zonal section along 15S
PLOT_SECTION = ('lat', -15.0)
