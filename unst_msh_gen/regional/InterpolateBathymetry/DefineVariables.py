mesh="meshes/HawaiiTest.msh"

BathyDir="/scratch3/NCEPDEV/climate/Keston.Smith/CoastalReliefModel/"
########################################################
# processed netcdf bathymetry files with variables:
# lon, lat and z. Specify corresponding file type, gridded or scattered,
# in list "filetype" below.
########################################################
flnms=[
    "cmems_obs-sdb_glo_phy-comp_my-oa-100m-l4-s2_static.PointValues500m.nc",
    "crm_vol1_2023.nc.S250m.VB.nc",
    "crm_vol2_2023.nc.S250m.VB.nc",
    "crm_vol3_2023.nc.S250m.VB.nc",
    "crm_vol4_2023.nc.S250m.VB.nc",
    "crm_vol5_2023.nc.S250m.VB.nc",
    "crm_vol7_2024.nc.S250m.VB.nc",
    "crm_vol9_2023.nc.S250m.VB.nc",
    "crm_vol10_2023.nc.S250m.VB.nc",
    "crm_vol6_2023.nc.S250m.VB.nc",
    "crm_vol8_2023.nc.S250m.VB.nc",
    "crm_southak.CRMformat.nc",
    "hurl_bathy_60m_nwhi.CRMformat.nc.S250m.VB.nc",
    "PIX/ngdc_bathy_10m_wake.CRM.nc",
    "PIX/pibhmc_bathy_20m_jarvis.CRM.nc",
    "PIX/pibhmc_bathy_20m_johnston.CRM.nc",
    "PIX/pibhmc_bathy_20m_kingman.CRM.nc",
    "PIX/pibhmc_bathy_40m_baker.CRM.nc",
    "PIX/pibhmc_bathy_40m_howland.CRM.nc",
    "PIX/pibhmc_bathy_40m_palmyra.CRM.nc",
    "PIX/pibhmc_bathy_40m_rose.CRM.nc",
    "PIX/pibhmc_bathy_40m_swains.CRM.nc",
    "PIX/pibhmc_bathy_40m_vailuluu.CRM.nc",
    "PIX/pibhmc_bathy_5m_palmyra.CRM.nc",
    "PIX/sopac_bathy_50m_majuro_reef.CRM.nc",
    "PIX/ngdc_bathy_180m_mariana.CRM.nc",
    "ngdc_bathy_90m_amsamoa.crm.nc",
    "RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc"
    ]

########################################################
########################################################
# Set filetype corresponding to files in flnms 
# 0 : for gridded local data and 
# 1 : for gridded global data and 
# 2 : for scattered global data
########################################################
nf=len(flnms)
filetype=np.zeros(nf, dtype=int)
filetype[0]=2
filetype[nf-1]=1
########################################################

Zmin=-11000. #deepest legit ocean depth value in case of mask fail 
# nescesarry for data sets "crm_vol6_2023.nc.S250m.nc" and "crm_vol8_2023.nc.S250m.nc"
Zmax=0. # maximum value to include in interpolation - zero out land values or ignore land 
Dmin=5./1000.# 5 meter minimum distance between observation points in km, prevent singularity
lambdaLL=.025 # set deg lat, lon search width for overlapping regions
#use Approximately NxTarget**2 points per data set, 
#each linear system is approx 2*(NxTarget**2)
NxTarget=20 # N=20 ~ 10 grid points in each cardinal direction used in interpolating 
NpointsMax=1000
