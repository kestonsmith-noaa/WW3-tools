mesh="meshes/rwps.oc_500m_10km.msh"

########################################################
# Begin: define data sets used in the interpolation
BathyDir="/scratch3/NCEPDEV/climate/Keston.Smith/CoastalReliefModel/"
BathyDir="/mnt/sda/keston/CoastalReliefModel/"
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
    "PIX/ngdc_bathy_180m_mariana.CRM.nc",
    "ngdc_bathy_90m_amsamoa.crm.nc",
    "PIX/pibhmc_bathy_10m_agrihan.CRM.S200m.nc",
    "PIX/pibhmc_bathy_10m_alamagan.CRM.S200m.nc",
    "PIX/pibhmc_bathy_10m_asuncion.CRM.S200m.nc",
    "PIX/pibhmc_bathy_10m_guguan.CRM.S200m.nc",
    "PIX/pibhmc_bathy_10m_maug.CRM.S200m.nc",
    "PIX/pibhmc_bathy_10m_pagan.CRM.S200m.nc",
    "PIX/pibhmc_bathy_10m_pajaros.CRM.S200m.nc",
    "PIX/pibhmc_bathy_10m_sarigan.CRM.S200m.nc",
    "PIX/pibhmc_bathy_10m_supply.CRM.S200m.nc",
    "PIX/pibhmc_bathy_20m_jarvis.CRM.S200m.nc",
    "PIX/pibhmc_bathy_20m_johnston.CRM.S200m.nc",
    "PIX/pibhmc_bathy_20m_nebank.CRM.S200m.nc",
    "PIX/pibhmc_bathy_20m_kingman.CRM.S200m.nc",
    "PIX/pibhmc_bathy_40m_vailuluu.CRM.S200m.nc",
    "PIX/pibhmc_bathy_40m_twoperbank.CRM.S200m.nc",
    "PIX/pibhmc_bathy_40m_swains.CRM.S200m.nc",
    "PIX/pibhmc_bathy_40m_rose.CRM.S200m.nc",
    "PIX/pibhmc_bathy_40m_palmyra.CRM.S200m.nc",
    "PIX/pibhmc_bathy_40m_howland.CRM.S200m.nc",
    "PIX/pibhmc_bathy_40m_baker.CRM.S200m.nc",
    "PIX/pibhmc_bathy_5m_palmyra.CRM.S200m.nc",
    "PIX/sopac_bathy_50m_majuro_reef.CRM.S200m.nc",
    "PIX/pibhmc_bathy_60m_rota.CRM.S200m.nc",
    "PIX/pibhmc_bathy_5m_tinian.CRM.S200m.nc",
    "PIX/ngdc_bathy_10m_tutuila.CRM.S200m.nc",
    "PIX/ngdc_bathy_10m_wake.CRM.S200m.nc",
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
filetype=np.zeros(nf, dtype=int) # local gridded bathymetry files
filetype[0]=2 # scattered point data file: cmems_obs-sdb_glo_phy-comp_my-oa-100m-l4-s2_static.PointValues500m.nc
filetype[nf-1]=1 # global gridded bathymetry file: RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc
########################################################
########################################################
# Set file observational standerd error (m) for each data file 
########################################################
nf=len(flnms)
fileStdErr=np.zeros(nf)+1. # local gridded bathymetry files
fileStdErr[0]=2. # scattered point data file: cmems_obs-sdb_glo_phy-comp_my-oa-100m-l4-s2_static.PointValues500m.nc
fileStdErr[nf-1]=10. # global gridded bathymetry file: RTopo_2_0_4_GEBCO_v2023_60sec_pixel.CRMformat.nc
# Set Background stdErr
BackGroundProcStdErr=50.

# Introduce log(observed depth) dependent error 
# has problems both nin theory and practice 
# please leave LogDepthDependantObsErr=False
# Introduced to address some artifacts due to
# data set structures
LogDepthDependantObsErr=False
ZminDepthDependantObsErr=10.
########################################################
# End: define data sets used in the interpolation
########################################################

########################################################
# covariance model parameters
#covmod="Sph" # use spherical covariance
covmod="Exp" # use exponential covariance
CovExp=2 # Exponent in covariance funtion, used if covmod="Exp", 
IDExp=2 # exponent used in inverse distance estimator
CovarianceLengthscaleFactor=1.0 # covariance Length Scale = CovarianceLengthscaleFactor * (local Mesh length scale)
                               # Prefer set to 1. if covmod="Exp" and 2, if covmod="Sph" 
########################################################
TextOutput = True # output progress and stats durring interpolation

########################################################
#Behavior of data point inclusion relative to "land" 
Zmax=1. # maximum value to include in interpolation - zero out land values or ignore land 
ExcludeLand=False # if "True" data points above Zmax(==0) are not used in the interpolation 
ZeroOutLandVals=True # if "True" data points above Zmax(==0) are set to Zmax
########################################################
Zmin=-11000. #deepest legit ocean depth value in case of mask fail 
# nescesarry for data sets "crm_vol6_2023.nc.S250m.nc" and "crm_vol8_2023.nc.S250m.nc"

Dmin=5./1000. # 5 meter minimum distance between observation points in km, prevent singularity
lambdaLL=.025 # set deg lat, lon search width for overlapping regions
#use Approximately NxTarget**2 points per data set, 
#each linear system is approx 2*(NxTarget**2)
NxTarget=20 # N=20 ~ 10 grid points in each cardinal direction used in interpolating 
NpointsMax=2400 # take NpointsMax closest points if more than that are found in data search
