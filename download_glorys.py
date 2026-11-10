"""Download GLORYS12 for a box crossing the date line (170E - 170W).
pip install copernicusmarine   ;   copernicusmarine login   (once)"""
import copernicusmarine

LAT_MIN, LAT_MAX = 20.0, 36.0          # <-- your domain
MARGIN = 0.5
common = dict(
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",       # daily; "..._P1M-m" for monthly
    variables=["thetao", "so", "uo", "vo", "zos"],
    minimum_latitude=LAT_MIN - MARGIN, maximum_latitude=LAT_MAX + MARGIN,
    minimum_depth=0, maximum_depth=6000,
    start_datetime="2016-01-01T00:00:00", end_datetime="2016-01-31T00:00:00",
    output_directory="glorys_pacific",
)

try:
    # Recent toolbox versions handle the antimeridian: 170W is written as 190.
    copernicusmarine.subset(minimum_longitude=170 - MARGIN, maximum_longitude=190 + MARGIN,
                            output_filename="glorys_pacific_201601.nc", **common)
except Exception as e:
    print("single request failed (", e, ") -> downloading the two halves")
    copernicusmarine.subset(minimum_longitude=170 - MARGIN, maximum_longitude=180,
                            output_filename="glorys_pacific_201601_west.nc", **common)
    copernicusmarine.subset(minimum_longitude=-180, maximum_longitude=-170 + MARGIN,
                            output_filename="glorys_pacific_201601_east.nc", **common)
# make_ic_obc.py stitches the halves automatically (same time stamps -> one record).
