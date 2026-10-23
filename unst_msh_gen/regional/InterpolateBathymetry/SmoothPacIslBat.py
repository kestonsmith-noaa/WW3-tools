import numpy as np
#from astropy.convolution import convolve
import netCDF4 as nc
from scipy import signal
from scipy.signal import convolve2d

PacificIslandBatPath="/scratch3/NCEPDEV/climate/Keston.Smith/CoastalReliefModel/PIX"

flnms=[
    "pibhmc_bathy_10m_agrihan.CRM.nc",
    "pibhmc_bathy_10m_alamagan.CRM.nc",
    "pibhmc_bathy_10m_asuncion.CRM.nc",
    "pibhmc_bathy_10m_guguan.CRM.nc",
    "pibhmc_bathy_10m_maug.CRM.nc",
    "pibhmc_bathy_10m_pagan.CRM.nc",
    "pibhmc_bathy_10m_pajaros.CRM.nc",
    "pibhmc_bathy_10m_sarigan.CRM.nc",
    "pibhmc_bathy_10m_supply.CRM.nc",
    "pibhmc_bathy_20m_jarvis.CRM.nc",
    "pibhmc_bathy_20m_johnston.CRM.nc",
    "pibhmc_bathy_20m_kingman.CRM.nc",
    "pibhmc_bathy_20m_nebank.CRM.nc",
    "pibhmc_bathy_40m_baker.CRM.nc",
    "pibhmc_bathy_40m_howland.CRM.nc",
    "pibhmc_bathy_40m_palmyra.CRM.nc",
    "pibhmc_bathy_40m_rose.CRM.nc",
    "pibhmc_bathy_40m_swains.CRM.nc",
    "pibhmc_bathy_40m_twoperbank.CRM.nc",
    "pibhmc_bathy_40m_vailuluu.CRM.nc",
    "pibhmc_bathy_5m_palmyra.CRM.nc",
    "pibhmc_bathy_5m_tinian.CRM.nc",
    "pibhmc_bathy_60m_rota.CRM.nc",
    "sopac_bathy_50m_majuro_reef.CRM.nc",
    "ngdc_bathy_10m_tutuila.CRM.nc",
    "ngdc_bathy_10m_wake.CRM.nc"
    ]
#    "ngdc_bathy_180m_mariana.CRM.nc"

filval=-99999
target_dx=.200 # taget resolution km
for fl in flnms:
    flin=PacificIslandBatPath+"/"+fl
    data0 = nc.Dataset(flin,"r")
    x=np.asarray(data0["lon"][:])
    y=np.asarray(data0["lat"][:])
    
    if np.min(x) > 90. : #RWPS longitude convention
        x=x-360.
        
    dx=abs(111.*np.cos(y[0]*np.pi/180.)*(x[1]-x[0])) # dx ~ km
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print(fl)

    print("native resolution = "+str(dx)+" (km)")
    print("target resolution = "+str(target_dx)+" (km)")

    k=int(np.floor(target_dx/dx))
    m=2*k+1
    print("subsampling at stride "+str(k))
    print("smoothed resolution = "+str(dx*float(k))+" (km)")
    print("smoothing with window of width "+str(m)+" or "+str(dx*float(m))+" (km)")
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    
    filval0=data0["z"]._FillValue
    z=np.asarray(data0["z"][:])
    jfill=np.where(z==filval0)
    z[jfill]=float('nan')
    
    nx=len(x)
    ny=len(y)
    w=np.ones((m,m),dtype=float)/float(m)/float(m)
    
    z_zeroed = np.nan_to_num(z, nan=0.0)
    mask = np.astype(~np.isnan(z), float)
    numerator = convolve2d(z_zeroed, w, mode='same')
    denominator = convolve2d(mask, w, mode='same')
    zp = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator != 0)

    x=x[k+1 :: k]
    y=y[k+1 :: k]
    zpp=zp[k+1 :: k,k+1 :: k]
    
    nx=len(x)
    ny=len(y)
    flout=fl[0:-3]+'.S200m.nc'
    flout=PacificIslandBatPath+"/"+flout
    
    jb=np.where(np.isnan(z))
    z[jb] = filval

    with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
        ncout.createDimension('lon', nx)  # Unlimited dimension
        ncout.createDimension('lat', ny)
        
        lon_var=ncout.createVariable('lon', 'f8', ('lon',))
        lon_var.units         = 'degree_east'
        lon_var.long_name     = 'longitude'
        lon_var.standard_name = 'longitude'
        lon_var.axis          = 'lon'
        lon_var[:]=x[:]

        lat_var=ncout.createVariable('lat', 'f8', ('lat',))
        lat_var.units         = 'degree_north'
        lat_var.long_name     = 'latitude'
        lat_var.standard_name = 'latitude'
        lat_var.axis          = 'lat'
        lat_var[:]=y[:]

        z_var=ncout.createVariable('z', 'f4', ('lat','lon'),fill_value    = filval)
        z_var.grid_mapping  = 'crs'
        z_var.long_name     = 'z'
        z_var.units         = 'meters'
        z_var.positive      = 'up'
        z_var.standard_name = 'height'
        z_var.vert_crs_name = 'EGM2008'
        z_var. vert_crs_epsg = 'EPSG:3855'
        z_var[:,:]=zpp[:,:]
        
        ncout.close
