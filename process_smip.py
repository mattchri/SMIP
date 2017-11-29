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
#year month day version and method {in order}
#
#Output
#NetCDF file containing derived/merged products at each timestep
#
#METHOD 1: output products at their local time
#METHOD 2: average all products to match ECMWF time
#
#Example
#python2.7 -i process_smip.py 2015 3 1 fv1.1 1
#python2.7 -i process_smip.py 2015 3 1 fv1.2 2
#
# 02/11/17, MC: upload initial version of the code to the repo
# 29/11/17, MC: added ECMWF precip product and two methods to
#               1) output at each time step and 2) average products
#               at the temporal timestep of ECMWF.
#---------------------------------------------------------------
import sys
import os
from netCDF4 import Dataset
from subroutines import *
import numpy.ma as ma
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

#---------------------------------------------------------------
# input variables read from command line
#---------------------------------------------------------------
year  = int(sys.argv[1])
month = int(sys.argv[2])
dom   = int(sys.argv[3])
version = sys.argv[4]
method  = int(sys.argv[5])

#Extract julian day from input time field
jday = convert_time(year,month,dom)
doy = jday_to_doy(jday)
print "Input Time: ",year,month,dom,doy,jday
print "Version :",version
print "Averaging Method :",method

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
# Change this to all relevant ECMWF products
TP = read_ecmwf_variable('TP',year,month,dom) #for accumuluated data you must read in 2-12 hour periods
#plt.imshow(ARR[0,:,:])
#plt.show()

print 'Collocating Products Together'
#Merge Data - make netCDF file

#Declare netCDF output variables
fileID   = [-1   ,      1          ]
INVARS   = ['tp' ,'IRprecipitation']
COLL     = [2    ,      1          ]
OUTVARS  = ['prate_ecmwf','prate_imerg']
scaleF   = [3600.*1000. , 1.0 ]
maskVal  = [    0.      , 0.0 ]

varName=OUTVARS
varLong=['ECMWF_Precipitation','IMERG Precipitation']
varStan=OUTVARS
varUnit=['mm/hr','mm/hr']
varFill=[-999., -999.]

nVars  = len(INVARS)

sz = MatchedFiles.shape
nFileTypes = sz[0]
nFiles = sz[1]

# Read Base File Geolocation Data (this is invariant)
ncfile = Dataset( MatchedFiles[0,0], mode='r')
lonSV = ncfile.variables['lon'][:]
latSV = ncfile.variables['lat'][:]

xDim = (lonSV.shape)[0]
yDim = (latSV.shape)[1]


# Fetch time indices of each ECMWF timestep that match the SEVIRI timestep
ecmwf_tid = np.zeros(nFiles)
for i in range(nFiles):
    #Base files to do the collocation
    BaseDict     = DataProd[0] #by default zeroth element
    BaseFileInfo = BaseDict['files']
    BT     = BaseFileInfo[i]
    dTMP = (BT['JDAY'])[1] - ((TP['JDAY'])[1])
    ecmwf_tid[i] = int(np.asarray( np.where(abs(dTMP) == min(abs(dTMP))) ))




#Loop over each SEVIRI 15 min file
for i in range(nFiles):
    #Base files to do the collocation
    BaseDict     = DataProd[0] #by default zeroth element
    BaseFileInfo = BaseDict['files']
    BT     = BaseFileInfo[i]

    varData = np.empty( [nVars, xDim, yDim ] )    
    #Loop through each variable
    for j in range(nVars):
        print(i,j,INVARS[j])

        #IMERG to SEVIRI collocation
        if COLL[j] == 1:

            saveFile = outpath+'collocation_seviri_imerg'
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
 
            if method == 1: #instantaneous output
                inFile = MatchedFiles[ fileID[j], i ]
                ncfileIM = Dataset( inFile, mode='r')
                inDATA = np.asarray( ncfileIM.variables[INVARS[j]][:] )
                
                #Extract data from IMERG at SEVIRI location
                inDATA = inDATA[ LONidIM.astype(int), LATidIM.astype(int) ]
                inDATA = (LONidIM.astype(int) > -999.)*inDATA + (LONidIM.astype(int) == -999.)*(-999.)
                

            if method == 2: #select all files that coincide with the ECMWF averaging period
                imerg_step_id = (np.where(ecmwf_tid == ecmwf_tid[i]))[0]
                Nsteps = len(imerg_step_id)
                ALLData = np.zeros( [ Nsteps, len(LONidIM), len(LATidIM) ] )
                for k in range(Nsteps):
                    #Read Each IMERG File
                    inFile = MatchedFiles[ fileID[j], imerg_step_id[k]]
                    ncfileIM = Dataset( inFile, mode='r')
                    inDATA = np.asarray(ncfileIM.variables[INVARS[j]][:])
                    inDATA = (LONidIM.astype(int) > -1.)*inDATA[ LONidIM.astype(int), LATidIM.astype(int) ] + (LONidIM.astype(int) == -1.)*(0.)
                    ALLData[k,:,:] = inDATA
                    #plt.imshow(ALLData[k,:,:], vmin=0.0, vmax=5.0)
                    #plt.show()
                inDATA = np.mean(ALLData, axis=0)

        #ECMWF Collocation Routine
        if COLL[j] == 2:
            #Fetch Precipitation Data Over Day
            saveFile = outpath+'collocation_seviri_ecmwf'
            #Check if File already exists
            if os.path.isfile(saveFile+'.npy'):
                ecmwf = np.load(saveFile+'.npy')
                LATidEC = ecmwf.item().get('latid')
                LONidEC = ecmwf.item().get('lonid')
            else:
                print('Collocating ECMWF to SEVIRI grid')
                #Define lat/lon grid
                lonEC  = TP['lon']+(TP['lon'][1]-TP['lon'][0])/2.
                lonEC = (lonEC <= 180.)*(lonEC) + (lonEC > 180.)*(lonEC-360.)
                latEC  = TP['lat']

                #Run collocation subroutine
                ECMWF_ID = collocate_imerg_seviri(lonEC,latEC,lonSV,latSV) #uses same code as IMERG
                np.save(saveFile, ECMWF_ID)
                LATidEC = ECMWF_ID['lonid']
                LONidEC = ECMWF_ID['latid']

            #Determine time-step in ECMWF from SEVIRI
            dTIME = (BT['JDAY'])[1]-(TP['JDAY'])[1]
            dTIMEid = (np.argsort(dTIME))[::-1]
            inDATA = (TP['data'])[dTIMEid[0],:,:] * scaleF[j]

            #could take two indices and interpolate but will use closest index in time for now
            ##data = ma.masked_less(aot,0.0)*0.001
            lonID = (LONidEC.astype(int) > -1.)*LONidEC + (LONidEC.astype(int) <= -1.)*(-1.)
            latID = (LATidEC.astype(int) > -1.)*LATidEC + (LATidEC.astype(int) <= -1.)*(-1.)

            #Extract data from ECMWF at SEVIRI location
            inDATA = (lonID.astype(int) > -1)*inDATA[ latID.astype(int), lonID.astype(int) ] + \
                     (lonID.astype(int) == -1.)*(0.)

        inDATA = (inDATA >= maskVal[j])*inDATA + (inDATA < maskVal[j])*(-999.)
        varData[j,:,:] = inDATA

    #----------------------------------------------------
    #make netCDF file
    #----------------------------------------------------
    ncname=outpath+'merged_'+str(int(BT['YEAR'])).zfill(4)+str(int(BT['MONTH'])).zfill(2) + \
                   str(int(BT['DOM'])).zfill(2)+str(int(BT['HOUR'])).zfill(2) + \
                   str(int(BT['MINUTE'])).zfill(2)+'_SEVIRI_disk.'+version+'.nc'

    make_netcdf_file(ncname,varName,varData,varLong,varStan,varUnit,varFill)
    #----------------------------------------------------
