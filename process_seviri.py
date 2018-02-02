#python2.7 -i process_seviri.py

import datetime
import glob
from jdcal import *
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import trackpy as tp
import pandas as pd
from pandas import DataFrame, Series
import pims
import numpy.ma as ma
from cloudtrack import maketrack,watershedding_2D
#conda install -c conda-forge ffmpeg
from plotting import *
from html_animation import *

#---------------------------------
#Subroutines
#---------------------------------
def extract_time_seviri(fname):
    year  = int(fname[38:42])
    month = int(fname[42:44])
    dom   = int(fname[44:46])
    hour  = int(fname[46:48])
    minute= int(fname[48:50])
    jdayOUT = datetime.datetime(year,month,dom,hour,minute)
    OUT = { "filename" : fname, "MONTH" : month, "YEAR" : year, "DOM" : dom, "HOUR" : hour, "MINUTE" : minute, "JDAY" : jdayOUT}
    return OUT

def file_mkdir(dirpath):
    if not os.path.exists(os.path.dirname(dirpath)):
        try:
            os.makedirs(os.path.dirname(dirpath))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    errval=0
    return errval



#---------------------------------
# Paths
#---------------------------------
#Data path
root_path_seviri='/group_workspaces/jasmin/acpc/Data/ORAC/cloudtrack/'

#Output path
outPath = '/group_workspaces/jasmin/acpc/Data/SEVIRI/cloudtrack/'
ferr = file_mkdir(outPath)

#Output path
figPath = '/group_workspaces/cems/cloud_ecv/public/temp_transfer/tracking/v3/'
ferr = file_mkdir(figPath)


#---------------------------------
#Start Time - input
#---------------------------------
year1 = 2017
month1= 8
day1 = 15

year2 = 2017
month2 = 8
day2 = 16

jday1 = datetime.datetime(year1,month1,day1,0,0)
jday2 = datetime.datetime(year2,month2,day2,0,0)

#---------------------------------
#Select SEVIRI area of interest
#  (note could also do this by lon/lat)
#---------------------------------
dX=400
dY=400
x0=2400
y0=1600
x1=x0+dX
y1=y0+dY


#---------------------------------
# End Inputs: Starting Main code
#---------------------------------


#---------------------------------
# Fetch & sort all SEVIRI files
#---------------------------------
jdays = []
seviri_files = []
for file in os.listdir(root_path_seviri):
    if file.endswith("merged.nc"):
        tmp = extract_time_seviri(file)
        if( (jday2 - tmp['JDAY']).total_seconds() > 0. and (jday1 - tmp['JDAY']).total_seconds() < 0.):
            jdays.append( (tmp['JDAY']) )
            seviri_files.append(tmp)

#Sort by time
sortID = np.asarray(np.argsort(jdays))
seviri_files = np.asarray(seviri_files)
seviri_files = seviri_files[sortID]
fct = sortID.shape[0]
print(seviri_files)

#---------------------------------
#assign unique filename and query if save file was already created
#cloudtrack_YYYYMMDD_YYYYMMDD_X0_Y0_X1_Y1
#---------------------------------
file_prefix = 'cloudtrack_'+str(year1).zfill(4)+str(month1).zfill(2)+str(day1).zfill(2) + '_' + str(year2).zfill(4)+str(month2).zfill(2)+str(day2).zfill(2) + '_' + str(x0).zfill(4)+'_' + str(y0).zfill(4)+'_'+str(x1).zfill(4)+ '_' +str(y1).zfill(4)


#-----------------------------------
#  Read data into iris cubes
#-----------------------------------
irisFile = outPath+file_prefix+'_'+'cube.nc'
if not os.path.isfile(irisFile):
    print('read data into iris cubes')
    files = [root_path_seviri+(seviri_files[i])['filename'] for i in range(fct)]
    cubes = iris.cube.CubeList()
    for i in range(fct):
        print('loading frame: ',i)
        inputfile = files[i]
        #ncfile = Dataset( inputfile, mode='r')
        #ARR = ncfile.variables['brightness_temperature_in_channel_no_10'][x0:x1,y0:y1]
        #cube = iris.load_cube(inputfile)[x0:x1,y0:y1]
        cube = iris.load_cube(inputfile, constraint=iris.Constraint(cube_func=(lambda c: c.var_name == 'brightness_temperature_in_channel_no_10')))[x0:x1,y0:y1]
        cubes.append(cube)

        #Filter Data
        #This step depends on the type of clouds you want to track!
        #mARR = (ARR > tMin)*10000 + (ARR <= tMin)*ARR
        #frames.append(mARR)
        #data.append(ARR)

    print('merge iris cubes')
    T_cube = cubes.merge_cube()
    print('saving: ',irisFile)
    iris.save([T_cube],irisFile)

else:
    print('loading: cubes',irisFile)
    T_cube=iris.load_cube(irisFile)

#-----------------------------------
#  Run Tracking Algorithm
#-----------------------------------
hdfFile = outPath+file_prefix+'_'+'track.hdf'
if not os.path.isfile(hdfFile):
    print('run: tracking routine')
    Track=maketrack(T_cube,grid_spacing=4000,diameter=220000,target='minimum',v_max=20,min_mass=5000,min_signal=20,stubs=6,order=2,extrapolate=2)
    print('saving tracks: ',hdfFile)
    Track.to_hdf(hdfFile,'table')

else:
    print('loading tracks: ',hdfFile)
    Track=pd.read_hdf(hdfFile,'table')

#-----------------------------------
#  Run Water shedding Algorithm
#-----------------------------------
maskFile = outPath+file_prefix+'_'+'Mask.nc'
if not os.path.isfile(maskFile):
    print('run: masking routine')
    Mask=watershedding_2D(Track,T_cube,threshold=260,target='minimum',method='watershed')
    print('saving Mask')
    iris.save([Mask],maskFile)

else:
    print('loading Masks: ',maskFile)
    Mask=iris.load_cube(maskFile)


#-----------------------------------
#  Plotting Tracks
#-----------------------------------
print('plotting tracks:')
plot_tracks_mask_field_loop(Track,T_cube,Mask,axis_extent=[-5,8,15,28],text_mass_signal=False,vmin=200,vmax=300,n_levels=50,name='Tb',plot_dir=figPath,figsize=(10,10))

#-----------------------------------
#  Make HTML & gif Files
#-----------------------------------
print('animating files')
# Fetch list of image files and sort them in ascending order by time
jdays = []
imagefiles = []
for file in os.listdir(figPath):
    if file.endswith(".png"):
        print(file)
        # read time from file
        YYYY = int(file[len(file)-23:len(file)-19])
        MM = int(file[len(file)-18:len(file)-16])
        DD = int(file[len(file)-15:len(file)-13])
        HR = int(file[len(file)-12:len(file)-10])
        MN = int(file[len(file)-9:len(file)-7])
        jdays.append(datetime.datetime(YYYY,MM,DD,HR,MN))
        imagefiles.append(file)

#Sort by time
sortID = np.asarray(np.argsort(jdays))
imagefiles = np.asarray(imagefiles)
imagefiles = imagefiles[sortID]

#Generate HTML Script
html_script( imagefiles, 'Tb', figPath, make_gif='yes')            









STOP
#Animate cloud tracks
fig = plt.figure()
ims = []
for i in range(fct):
    im = plt.imshow(Mask[i].data, animated=True)
    ims.append([im])
ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,repeat_delay=1000)
ani.save('seviri_track.mp4')
plt.show()

#Animate brightness temperature
fig = plt.figure()
ims_tb = []
for i in range(fct):
    im = plt.imshow(T_cube[i].data, animated=True)
    ims_tb.append([im])
ani_tb = animation.ArtistAnimation(fig, ims_tb, interval=50, blit=True,repeat_delay=1000)
ani_tb.save('seviri_tb.mp4')
plt.show()

