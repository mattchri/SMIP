"""
modis_imerg_col.py
Colocation routine to compare the MOD04 product to regions of precipitation shown in the iMerg product.

23/11/2017, WJ: Creation of routine
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import netCDF4 as nc
import pyhdf.SD as SD
import os

# Define paths to MOD04 and iMerg datasets, as well as file identifying prefixes and suffixes.
MOD04_path = r'/group_workspaces/cems2/nceo_generic/satellite_data/modis_c6/mod04_l2/2015/060/'
MOD04_pre = r'MOD04_L2.A'
MOD04_suf = r'.hdf'
"""MOD04 filename format:
prefix + '[YYYYDOY].[HHMM].006.[...]' + suffix
"""

iM_path = r'/group_workspaces/cems2/nceo_generic/satellite_data/precipitation/iMERG/download'
iM_pre = r'3B-HHR-E.MS.MRG.3IMERG.'
iM_suf = r'.V04A.HDF5.nc'
"""iMerg filename format
prefix + '[YYYYMMDD]-S[HHMMSS]-E[HHMMSS].[MNOD]' + suffix
"""

# Generate iMerg latitude/longitude grids
IM = nc.Dataset(im_path+im_name)
p_lat = IM.variables['lat'][:]
p_lat.shape = (1,1800)
p_lat = p_lat.repeat(3600,0)
p_lon = IM.variables['lon'][:]
p_lon.shape = (3600,1)
p_lon = p_lon.repeat(1800,1)

hr_min = 12
hr_max = 24

#Locate MODIS files and extract time from filename
M_dir = os.listdir(MOD04_path)
M_files = [f for f in M_dir if f.startswith(MOD04_pre)]
M_time = [(int(f[19:21]), int(f[21:23])) for f in M_files]
M_mnod = [t[0]*60 + t[1] for t in M_time]
M_ind = [t/30 for t in M_mnod]


# Loop over entire day
for i in range(48):
    # Check for correct time
    if i >= hr_min & i < hr_max:
        aod = np.empty
        tmp_cnt = 0
        for j,f in enmuerate(M_files):
            if M_ind[j] == i:
                # Open selected MOD04 file and read aerosol data
                M = SD.SD(MOD04_path+f)
                SDS = M.select('Optical_Depth_Land_And_Ocean')
                aod_tmp = SDS.get().astype('float')
                aod_tmp[aod_tmp == SDS.attributes()['_FillValue']] = np.nan
                aod_tmp *= SDS.attributes()['scale_factor']
                
                # Read modis lat/lon data
                t_lat = M.select('Latitude')
                t_lat = t_lat.get()
                t_lon = M.select('Longitude')
                t_lon = t_lon.get()
                
        
