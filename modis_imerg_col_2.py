"""
modis_imerg_col.py
Colocation routine to compare the MOD04 product to regions of precipitation shown in the iMerg product.

23/11/2017, WJ: Creation of routine
"""

import numpy as np
import scipy as sp
import scipy.interpolate as spinterpolate
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

iM_path = r'/group_workspaces/cems2/nceo_generic/satellite_data/precipitation/iMERG/download/'
iM_pre = r'3B-HHR-E.MS.MRG.3IMERG.'
iM_suf = r'.V04A.HDF5.nc'
"""iMerg filename format
prefix + '[YYYYMMDD]-S[HHMMSS]-E[HHMMSS].[MNOD]' + suffix
"""

# Locate iMERG files and find filenames and times
iM_dir = os.listdir(iM_path)
iM_files = [f for f in iM_dir if f.startswith(iM_pre)]
iM_time = [(int(f[33:35]), int(f[35:37])) for f in iM_files]
iM_mnod = [t[0]*60 + t[1] for t in iM_time]

# Generate iMerg latitude/longitude grids
IM = nc.Dataset(iM_path+iM_files[0])
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
M_files = [f for f in M_dir if f.startswith(MOD04_pre) and f[14:17] == '060']
M_time = [(int(f[18:20]), int(f[20:22])) for f in M_files]
M_mnod = [t[0]*60 + t[1] for t in M_time]
M_ind = [t/30 for t in M_mnod]

# Initialise output arrays
aod_wh = []
n_wh = []

def load_modis_aod(path_to_file, fill_value=np.nan):
    # Open MODIS hdf4 file and extract AOD dataset
    M = SD.SD(path_to_file)
    SDS = M.select('Optical_Depth_Land_And_Ocean')
    # Read Dataset values, scale and mask
    aod_tmp = SDS.get().astype('float')
    wh = (aod_tmp == SDS.attributes()['_FillValue'])
    aod_tmp *= SDS.attributes()['scale_factor']
    aod_tmp[wh] = fill_value
    # Read modis lat/lon data
    t_lat = M.select('Latitude').get()
    t_lon = M.select('Longitude').get()
    return aod_tmp.flatten(), t_lat.flatten(), t_lon.flatten()
    
#def colocate_aod_precip(aod, precip):
    

# Loop over entire day
col_aod = []
col_n = []
col_cld = []
for i in range(48):
    # Check for correct time
    if i >= hr_min*2 and i < hr_max*2:
        print i
        tmp_cnt_j = 0
        for j,f in enumerate(M_files):
            # Loop over relevant aerosol files
            if M_ind[j] == i:
                if tmp_cnt_j == 0:
                    aod, m_lat, m_lon = load_modis_aod(MOD04_path+f)
                    tmp_cnt_j = 1
                else:
                    aod_tmp, t_lat, t_lon = load_modis_aod(MOD04_path+f)
                    aod = np.concatenate((aod, aod_tmp))
                    m_lat = np.concatenate((m_lat, t_lat))
                    m_lon = np.concatenate((m_lon, t_lon))
            
        # Interpolate MODIS AOD data to iMERG grid
        aod_grid = spinterpolate.griddata((m_lon,m_lat),aod,(p_lon,p_lat),method='linear')
            
        #Loop over iMERG files for previous 12 hours
        col_tmp = []
        n_tmp = []
        cld_tmp = []
        for k,f in enumerate(iM_files):
            # Loop over previous 12 hours of precipitation (24 time steps)
            if k >= i-24 and k < i:
                # Load iMerg netcdf file
                IM = nc.Dataset(iM_path+f)
                prec = IM.variables['IRprecipitation'][:]
                # Filter low values and missing data
                prec[prec <= 2] = np.nan
                # Colocate data points and extract info
                col_val = aod_grid[np.isfinite(aod_grid) & np.isfinite(prec)]
                col_tmp.append(np.nanmean(col_val))
                n_tmp.append(col_val.size)
                cld_tmp.append(np.count_nonzero(np.isfinite(prec) & ~np.isfinite(aod_grid)))
        # Append temporary lists to full lists
        col_aod.append(col_tmp)
        col_n.append(n_tmp)
        col_cld.append(cld_tmp)
