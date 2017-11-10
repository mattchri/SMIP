#-----------------------------------------------------------------
#Module Containing all Subroutines
#required to run the SMIP code
#The first group of subroutines deal with extracting
#times and converting them to julian day or 
#collocating files together.
#The second group deals with collocation.
#
# 02/11/17, MC: upload initial version of the code to the repo
#-----------------------------------------------------------------
from jdcal import *
import sys
import numpy as np
from netCDF4 import Dataset
from datetime import datetime

def doy_to_month(doy,year):
    #fetch julian day
    jday0 = gcal2jd(year,1,1)
    jday1 = (jday0[0],jday0[1]+doy-1)
    caldat = jd2gcal(*jday1)
    month = caldat[1]
    return month

def convert_time(year,month,dom):
    #fetch julian day
    jday = gcal2jd(year,month,dom)
    return jday

def jday_to_doy(jday):
    caldat = jd2gcal(*jday)
    year = caldat[0]
    jday0 = gcal2jd(year,1,1)
    doy = jday[1]-jday0[1]+1
    return doy

def extract_time_seviri(fname):
    year  = float(fname[38:42])
    month = float(fname[42:44])
    dom   = float(fname[44:46])
    hour  = float(fname[46:48])
    minute= float(fname[48:50])
    fday = (hour+minute/60.)/24.
    jday = gcal2jd(year,month,dom)
    jdayOUT = (jday[0],jday[1] + fday)     
    OUT = { "filename" : fname, "MONTH" : month, "YEAR" : year, "DOM" : dom, "HOUR" : hour, "MINUTE" : minute, "JDAY" : jdayOUT}
    return OUT

def extract_time_imerg(fname):
    year = float(fname[23:27])
    month = float(fname[27:29])
    dom   = float(fname[29:31])
    hour  = float(fname[33:35])
    minute = float(fname[35:37]) + 15.
    fday = (hour+minute/60.)/24.
    jday = gcal2jd(year,month,dom)
    jdayOUT = (jday[0],jday[1] + fday)
    OUT = { "filename" : fname, "MONTH" : month, "YEAR" : year, "DOM" : dom, "HOUR" : hour, "MINUTE" : minute, "JDAY" : jdayOUT}
    return OUT

def extract_time_ecmwf(fname):
    year = float(fname[4:8])
    month = float(fname[8:10])
    dom   = float(fname[10:12])
    step  = float(fname[12:14])
    hour  = float(fname[14:16])+step
    minute = step
    fday = (hour+minute/60.)/24.
    jday = gcal2jd(year,month,dom)
    jdayOUT = (jday[0],jday[1] + fday)
    OUT = { "filename" : fname, "MONTH" : month, "YEAR" : year, "DOM" : dom, "HOUR" : hour, "MINUTE" : minute, "JDAY" : jdayOUT}
    return OUT


#-----------------------------------------------------------------
#COLLOCATE_FILES
#Input the base SEVIRI files (merged output in this case) and
#all other file types (e.g. imerg, ECMWF, ect) to collocate to SEVIRI
#
#Input list of files for each data product
#Output list of collocated files from each product
#-----------------------------------------------------------------
def collocate_files(infiles,dt):
    #SEVIRI Base files Information
    BaseDict     = infiles[0] #by default zeroth element
    BaseFileInfo = BaseDict['files']
    BaseFiles = [item['filename'] for item in BaseFileInfo]
    BaseJDAY = [item['JDAY'] for item in BaseFileInfo]
    BasePath = BaseDict['path']

    #Define Matched Array
    nFileTypes = len(infiles)
    nFiles = len(BaseFiles)
    print '# Data Products: ',nFileTypes
    print '# Base (collocation) Files: ',nFiles
    matchedFiles = np.empty([nFileTypes,nFiles], dtype="S512")
    for iT in range(nFiles):
        matchedFiles[0,iT] = BasePath+BaseFiles[iT]

    i=1 #initialised value required 1st file type in for loop
    #Loop over each data product
    for iT in range(nFileTypes-1):
        TMPDict     = infiles[i] #by default zeroth element
        TMPFileInfo = TMPDict['files']
        TMPFiles    = [item['filename'] for item in TMPFileInfo]
        TMPJDAY     = [item['JDAY'] for item in TMPFileInfo]
        TMPJDAY1 = [iii[1] for iii in TMPJDAY]
        TMPPath  = TMPDict['path']
        y = np.asarray(TMPJDAY1)

        #Loop over each file
        for j in range(nFiles):
            x = float((BaseJDAY[j])[1]) #julian day temp var
            fid = np.where(abs(y-x) == min(abs(y-x)))
            z = int(fid[0])
            #print i,j,BaseFiles[j],'  ',TMPFiles[z]
            matchedFiles[ i, j ] = TMPPath + TMPFiles[z]

        i=i+1 #start with first data product (+1)
    return matchedFiles

#-----------------------------------------------------------------
#COLLOCATE_IMERG_SEVIRI
#Routine collocates the IMERG lat/lon to the SEVIRI lat/lon grid
#The top-level code produces a save file so this only should be run
#one time since the lat/lon fields are invariant for both sensors
#Input
# lonIM, latIM, lonSV, latSV
#Output
#index for longitude and latitude arrays
#-----------------------------------------------------------------
def collocate_imerg_seviri(lonIM,latIM,lonSV,latSV):
    print('Collocating IMERG to SEVIRI (~1hr processing time)')
    xDim = (lonSV.shape)[0]
    yDim = (latSV.shape)[1]
    LATidIM = np.zeros( (xDim, yDim) )
    LONidIM = np.zeros( (xDim, yDim) )
    for iX in range(xDim):
        print(iX)
        for iY in range(yDim):
            tmpLat = latSV[iX,iY]
            tmpLon = lonSV[iX,iY]
            if tmpLon > -999. and tmpLat > -999.:
                LONidIM[iX,iY] = int(((np.where(abs(lonIM-tmpLon) == min(abs(lonIM-tmpLon))) )[0])[0])
                LATidIM[iX,iY] = int(((np.where(abs(latIM-tmpLat) == min(abs(latIM-tmpLat))) )[0])[0])
            else:
                LONidIM[iX,iY] = -999.
                LATidIM[iX,iY] = -999.                   
    IMERG_ID = {'latid':LATidIM, 'lonid':LONidIM, 'xdim':xDim, 'ydim':yDim}
    return IMERG_ID


#Create NetCDF file
def make_netcdf_file(ncname,vName,vData,vLong,vStan,vUnit,vFill):
    print('Creating: '+ncname)
    sz = vData.shape
    nVars= sz[0]
    xDim = sz[1]
    yDim = sz[2]

    f = Dataset(ncname,'w', format='NETCDF4_CLASSIC')
    f.createDimension('xdim', xDim)
    f.createDimension('ydim', yDim)
    
    for iV in range(len(vName)):
        tmpstr = vName[iV]
        temp = f.createVariable(tmpstr, 'f4', ('xdim','ydim'), zlib=True)
        temp[:,:] = vData[iV,:,:]
        temp.standard_name  = vStan[iV]
        temp.long_name  = vLong[iV]
        temp.units = vUnit[iV]
        temp.FillValue = vFill[iV]

    #global attributes
    today = datetime.today()
    f.description = "Satellite Model Integrated Product"
    f.history = "Created " + today.strftime("%d/%m/%y")
    f.close()


#-----------------------------------------------------------------
#Read Accumulated (gafs) ECWMF variable from /badc/ on JASMIN
#Inputs
#variable name (e.g. 'TP' total precipitation)
#year [int]
#month [int]
#dom day of the month [int]
#Output
#ARR [accumulated field]
#-----------------------------------------------------------------
def read_ecmwf_variable(vname,year,month,dom):
    path_ecmwf = '/badc/ecmwf-era-interim/data/ga/fs/'+str(year).zfill(4)+'/'+str(month).zfill(2)+'/'+str(dom).zfill(2)+'/' 
    name = 'TP'
    prefix = str(year).zfill(4)+str(month).zfill(2)+str(dom).zfill(2)
    ftimes = ['0003','0006','0009','0012','1203','1206','1209','1212']
    ARR = np.empty( [8,256,512] )
    for i in range(len(ftimes)):
        ecmwf_file = path_ecmwf+'gafs'+prefix+ftimes[i]+'.nc'
        print(ecmwf_file)
        ncfile = Dataset( ecmwf_file, mode='r')
        ARR[i,:,:] = (ncfile.variables[vname][:])[0,0,:,:]
    return ARR
