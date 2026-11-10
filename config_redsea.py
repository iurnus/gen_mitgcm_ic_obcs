# Red Sea, 0.08 deg, east boundary open (OB_Ieast = 250*-1)
GLORYS_FILES = 'glorys_redsea/*.nc'   # one file per period, or east/west halves -- both fine
BATHY = 'bathy_redsea.bin'            # your model bathymetry (negative = ocean), Ny x Nx
OUT = 'out_redsea'
Nx, Ny = 250, 250
xgOrigin, ygOrigin = 30.0, 10.0
delX, delY = 0.08, 0.08
delR = [5.0, 5.1, 5.3, 5.7, 6.2, 6.9, 7.7, 8.7, 9.8, 11.0, 12.5, 14.0, 15.7, 17.6, 19.6, 21.8,
        24.1, 26.5, 29.1, 31.9, 34.8, 37.9, 41.1, 44.4, 47.9, 51.6, 55.4, 59.3, 63.4, 67.7,
        72.1, 76.6, 81.3, 86.2, 91.2, 96.3, 101.6, 107.0, 112.6, 118.4, 124.2, 130.3, 136.5,
        142.8, 149.3, 155.9, 162.7, 169.6, 176.7, 183.9]
OPEN_BOUNDARIES = 'E'
HFACMIN, HFACMINDR = 1.0, 0.0

# plot_ic_obc.py only: section that follows the deepest point of each row (long, narrow basin)
PLOT_SECTION = 'deepest'
