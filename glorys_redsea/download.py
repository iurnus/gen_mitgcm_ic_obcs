import copernicusmarine

copernicusmarine.subset(
  dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
  variables=["bottomT", "mlotst", "siconc", "sithick", "so", "thetao", "uo", "usi", "vo", "vsi", "zos"],
#   minimum_longitude=-180,
#   maximum_longitude=179.9166717529297,
  minimum_longitude=28,
  maximum_longitude=52,
  minimum_latitude=8,
  maximum_latitude=32,
  start_datetime="2016-01-31T00:00:00",
  end_datetime="2016-01-31T00:00:00",
  minimum_depth=0.49402499198913574,
  maximum_depth=5727.9,
  output_filename = "CMEMS_RedSea_2016-01-31.nc"
)

