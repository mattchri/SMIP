import numpy as np
import scipy as sp
import scipy.interpolate as spint
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyhdf.SD as SD

SEV_name = r'ESACCI-L2-CLOUD-CLD-SEVIRI_CC4CL_MSG3_201503011027_fv2.0.merged.nc'
SEV_path = r'/group_workspaces/cems2/nceo_generic/satellite_data/seviri/ORAC/2015/03/01/'
MOD_name = r'MOD04_L2.A2015060.1030.006.2015061042711.hdf'
MOD_path = r'/group_workspaces/cems2/nceo_generic/satellite_data/modis_c6/mod04_l2/2015/060/'

SEV = nc.Dataset(SEV_path+SEV_name)

lon = SEV.variables['lon'][:]
lat = SEV.variables['lat'][:]

MOD = SD.SD(MOD_path+MOD_name)

OD = MOD.select('Optical_Depth_Land_And_Ocean')
aod = OD.get()
aod = aod.astype('float')
aod[aod == OD.attributes()['_FillValue']] = np.nan
aod *= OD.attributes()['scale_factor']

m_lat = MOD.select('Latitude')
m_lat = m_lat.get()
m_lon = MOD.select('Longitude')
m_lon = m_lon.get()

aod_grid = spint.griddata((m_lon.flatten(),m_lat.flatten()),aod.flatten(),(lon.data,lat.data),method='linear')

toa_swup = SEV.variables['toa_swup'][:]

plt.imshow(toa_swup[::-1,::-1], 'gray')
plt.imshow(aod_grid[::-1,::-1])
plt.show()
