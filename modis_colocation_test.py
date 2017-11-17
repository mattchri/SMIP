#-----------------------------------------------------------------
# Test code for interpolating MODIS aerosol data to the SEVIRI 
# grid
#
#15/11/17, WJ: Initial upload
#17/11/17, WJ: Update to colocate three MODIS granules within time
#   period, colocate iMerg precipitation
#-----------------------------------------------------------------

import numpy as np
import scipy as sp
import scipy.interpolate as spint
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyhdf.SD as SD

# Define paths and file names to SEVIRI and MODIS data
SEV_name = r'ESACCI-L2-CLOUD-CLD-SEVIRI_CC4CL_MSG3_201503011027_fv2.0.merged.nc'
SEV_path = r'/group_workspaces/cems2/nceo_generic/satellite_data/seviri/ORAC/2015/03/01/'
MOD_name = [r'MOD04_L2.A2015060.1030.006.2015061042711.hdf', r'MOD04_L2.A2015060.1035.006.2015061040829.hdf', r'MOD04_L2.A2015060.1040.006.2015061040357.hdf']
MOD_path = r'/group_workspaces/cems2/nceo_generic/satellite_data/modis_c6/mod04_l2/2015/060/'
# Will be replaced with functional file loader for specific time/date

# Open SEVIRI dataset
SEV = nc.Dataset(SEV_path+SEV_name)

# Get SEVIRI lon/lat arrays as masked arrays
lon = SEV.variables['lon'][:]
lat = SEV.variables['lat'][:]

# Loop over MODIS files
for i, name in enumerate(MOD_name):
	MOD = SD.SD(MOD_path+name)

	# Read AOD data, remove missing values and scale
	OD = MOD.select('Optical_Depth_Land_And_Ocean')
	aod_tmp = OD.get().astype('float')
#	aod_tmp = aod.astype('float')
	aod_tmp[aod_tmp == OD.attributes()['_FillValue']] = np.nan
	aod_tmp *= OD.attributes()['scale_factor']

	# Read modis lat/lon data
	t_lat = MOD.select('Latitude')
	t_lat = t_lat.get()
	t_lon = MOD.select('Longitude')
	t_lon = t_lon.get()
	
	# Concatenate to single array for all three datasets. Flattening required for griddata input
	if i == 0:
		aod = aod_tmp.flatten()
		m_lat = t_lat.flatten()
		m_lon = t_lon.flatten()
	else:
		aod = np.concatenate((aod, aod_tmp.flatten()))
		m_lat = np.concatenate((m_lat, t_lat.flatten()))
		m_lon = np.concatenate((m_lon, t_lon.flatten()))

# Regrid to SEVIRI using the scipy.interpolate griddata routine using multivariate linear interpolation. Note, cubic fails and produces all NaNs
aod_grid = spint.griddata((m_lon.flatten(),m_lat.flatten()),aod.flatten(),(lon.data,lat.data),method='linear')

toa_swup = SEV.variables['toa_swup'][:]

# Plot ORAC SWUP and regridded AOD
plt.imshow(toa_swup[::-1,::-1], 'gray')
plt.imshow(aod_grid[::-1,::-1])
plt.show()

# Load iMerg dataset
im_path = r'/group_workspaces/cems2/nceo_generic/satellite_data/precipitation/iMERG/download/'
im_name = r'3B-HHR-E.MS.MRG.3IMERG.20150301-S100000-E102959.0600.V04A.HDF5.nc'
IM = nc.Dataset(im_path+im_name)
prec = IM.variables['IRprecipitation']
p_dat = prec[:]
p_dat[p_dat == prec._FillValue] = np.nan
p_lat = IM.variables['lat'][:]
p_lat.shape = (1,1800)
p_lat = p_lat.repeat(3600,0)
p_lon = IM.variables['lon'][:]
p_lon.shape = (3600,1)
p_lon = p_lon.repeat(1800,1)

# Regrid. Note this takes substationally longer due to the larger dataset
p_grid = sp.interpolate.griddata((p_lon.flatten(),p_lat.flatten()), p_dat.flatten(), (lon.data, lat.data), method='linear')
