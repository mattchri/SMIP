#---------------------------------------------------------------
#TOP LEVEL CODE
#---------------------------------------------------------------
#PROCESS_SMIP
#
#The SMIP code is a collocation routine to match various satellite 
#and model products to the SEVIRI native resolution using ORAC as 
#the primary input product footprint of 3 km resolution at nadir.
#
#Required Inputs
#year month day {in order}
#
#Output
#NetCDF file containing merged products at each timestep
#
#Example
#python2.7 -i process_smip.py 2015 3 1 fv1.0
#
# 02/11/17, MC: upload initial version of the code to the repo
#---------------------------------------------------------------
import sys
import os
from netCDF4 import Dataset
from subroutines import *


#---------------------------------------------------------------
# input variables read from command line
#---------------------------------------------------------------
year  = int(sys.argv[1])
month = int(sys.argv[2])
dom   = int(sys.argv[3])
version   = sys.argv[4]

#Extract julian day from input time field
jday = convert_time(year,month,dom)
doy = jday_to_doy(jday)
print "Input Time: ",year,month,dom,doy,jday


#---------------------------------------------------------------
#Paths (JASMIN)
#---------------------------------------------------------------
outpath = '/group_workspaces/jasmin2/aopp/mchristensen/erc/seviri/merged/'
path_seviri = '/group_workspaces/cems2/nceo_generic/satellite_data/seviri/ORAC/'+str(year).zfill(4)+'/'+str(month).zfill(2)+'/'+str(dom).zfill(2)+'/'
path_imerg = '/group_workspaces/cems2/nceo_generic/satellite_data/precipitation/iMERG/download/'
path_ecmwf = '/badc/ecmwf-era-interim/data/ga/fs/'+str(year).zfill(4)+'/'+str(month).zfill(2)+'/'+str(dom).zfill(2)+'/' 
#---------------------------------------------------------------

#Dictionary containing data products and arrays to store collocated files
DataProd=[{'prod':'orac'    ,'prefix':'merged.nc', 'path':path_seviri, 'sensor':'seviri', 'files':[]},
          {'prod':'imerg'   ,'prefix':'V04A.HDF5.nc','path':path_imerg, 'sensor':'multi_precip' , 'files':[]}]

#Find matching files for each data product
for i in range(len(DataProd)): #loop over array elements
    path = (DataProd[i])['path']
    prefix = (DataProd[i])['prefix']
    sensor = (DataProd[i])['sensor']

    tmplist = []
    for file in os.listdir(path): #fetch files for data product specified in path
        if file.endswith(prefix):
            if sensor == 'seviri':
                tmp=extract_time_seviri(file)
            if sensor == 'multi_precip':
                tmp=extract_time_imerg(file)
            if sensor == 'era_interim':
                tmp=extract_time_ecmwf(file)
            tmplist.append(tmp)
    (DataProd[i])['files'] = tmplist


# Match Files to SEVIRI bugsrad
dt = (35.*60./86400.) #nearest file within 35 minutes
MatchedFiles = collocate_files(DataProd, dt)


# Read Precipitation Data
TP = read_ecmwf_variable('TP',year,month,dom)
#import matplotlib.pyplot as plt
#plt.imshow(ARR[0,:,:])
#plt.show()


print 'Collocating Products Together'
#Merge Data - make netCDF file

#Declare netCDF output variables
fileID = [0    ,  0  ,    0     ,    0     ,      1          ]
VARS   = ['lat','lon','toa_swup','toa_lwup','IRprecipitation']
COLL   = [0    ,  0  ,    0     ,    0     ,      1          ]


varName=VARS
varLong=['latitude','longitude','top of atmosphere upwards shortwave flux','top of atmosphere upwards longwave flux','IMERG Precipitation']
varStan=VARS
varUnit=['degrees','degrees','W/m^2','W/m^2','mm/hr','W/m^2/sr','um','mm/hr']
varFill=[-999.    , -999.   , -999. , -999. , -999. , -999.    ,-999., -999.]

nVars  = len(VARS)

sz = MatchedFiles.shape
nFileTypes = sz[0]
nFiles = sz[1]

# Read Base File Geolocation Data (this is invariant)
ncfile = Dataset( MatchedFiles[0,0], mode='r')
lonSV = ncfile.variables['lon'][:]
latSV = ncfile.variables['lat'][:]

xDim = (lonSV.shape)[0]
yDim = (latSV.shape)[1]

#Loop over each SEVIRI 15 min file
for i in range(nFiles):
    #Base files to do the collocation
    BaseDict     = DataProd[0] #by default zeroth element
    BaseFileInfo = BaseDict['files']
    BT     = BaseFileInfo[i]

    varData = np.empty( [nVars, xDim, yDim ] )    
    #Loop through each variable
    for j in range(nVars):
        inFile = MatchedFiles[ fileID[j], i ]
        print(i,j,inFile)

        #No collocation required - product already on SEVIRI grid
        if COLL[j] == 0:
            ncfile = Dataset( inFile, mode='r')
            inDATA = ncfile.variables[VARS[j]][:]
        
        #IMERG to SEVIRI collocation
        if COLL[j] == 1:
            ncfileIM = Dataset( inFile, mode='r')
            inDATA = ncfileIM.variables[VARS[j]][:]

            saveFile = outpath+'collocation_seviri_ecmwf'
            #Check if File already exists
            if os.path.isfile(saveFile+'.npy'):
                imerg = np.load(saveFile+'.npy')
                LATidIM = imerg.item().get('latid')
                LONidIM = imerg.item().get('lonid')
            else:
                print('Collocating IMERG to SEVIRI grid')
                lonIM  = ncfileIM.variables['lon'][:]
                latIM  = ncfileIM.variables['lat'][:]
                #Run collocation subroutine
                IMERG_ID = collocate_imerg_seviri(lonIM,latIM,lonSV,latSV)
                np.save(saveFile, IMERG_ID)
                LONidIM = IMERG_ID['lonid']
                LATidIM = IMERG_ID['latid']

            #Extract data from IMERG at SEVIRI location
            inDATA = (LONidIM.astype(int) > -999.)*inDATA[ LONidIM.astype(int), LATidIM.astype(int) ] + \
                     (LONidIM.astype(int) == -999.)*(-999.)

        #ECMWF Collocation Routine
        if COLL[j] == 2:
            #Fetch Precipitation Data Over Day
            stop

        varData[j,:,:] = inDATA

    #----------------------------------------------------
    #make netCDF file
    #----------------------------------------------------
    ncname=outpath+'merged_'+str(int(BT['YEAR'])).zfill(4)+str(int(BT['MONTH'])).zfill(2) + \
                   str(int(BT['DOM'])).zfill(2)+str(int(BT['HOUR'])).zfill(2) + \
                   str(int(BT['MINUTE'])).zfill(2)+'_SEVIRI_disk.'+version+'.nc'

    make_netcdf_file(ncname,varName,varData,varLong,varStan,varUnit,varFill)
    #----------------------------------------------------
