#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ProjectMain.py(w)
#  
#  Copyright 2021 Farrokh A. Ghavanini <ghavanini@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

#  
#  
import netCDF4
from WW3TabFunctions import blockIdentifier
import numpy as np
import datetime as dt

strt_time = dt.datetime.now()

"""
Initializing Variables
"""
source_lons = []
source_lats = []
source_hs = []
source_wl = []
source_tr = []
source_dir = []
source_spr = []
source_fp = []
source_pdir = []
source_pspr = []
times = []
mask = []

# Setting input file path
input_file = 'tab44.ww3'

"""
Reading Source File Data
"""
print('\n---> Reading Source File Data')
# Identifying the total number of lines for block of a data between two timesteps
block_first_line, block_last_line = blockIdentifier(input_file)
block_size = block_last_line - block_first_line + 1

# Opening input file, reading its lines, counting timesteps and setting start time and date
with open (input_file, 'r') as source_file:
    source_file_lines = source_file.readlines()
    since_date_year = source_file_lines[0].split()[2].split('/')
    since_date_time = source_file_lines[0].split()[3]
    time_units = 'seconds since ' + "-".join(since_date_year) + ' ' + since_date_time
    timesteps = 0
    for line in source_file_lines:        
        if 'Time' in line:
            timesteps += 1

# Reading mask data from the corresponding file
with open ('mask.ind', 'r') as mask_file:
    mask_file_lines = mask_file.readlines()
for line in range(len(mask_file_lines)):
    mask.append(float(mask_file_lines[line].split()[5]))

# Filling variables from lines of data
first_line, last_line = block_first_line, block_last_line
for timestep in range(timesteps):
    print('  --> Processing timestep ' + str(timestep) + ' Processed. (ET: ' + str((dt.datetime.now() - strt_time).seconds) + 's)')
    if timestep != 0:
        first_line = first_line + block_size + 7 
        last_line = last_line + block_size + 7
    for line in range(first_line, last_line + 1):
        if timestep == 0:            
            source_lons.append(float(source_file_lines[line].split()[0]))
            source_lats.append(float(source_file_lines[line].split()[1]))
        source_hs.append(float(source_file_lines[line].split()[2]))
        source_wl.append(float(source_file_lines[line].split()[3]))
        source_tr.append(float(source_file_lines[line].split()[4]))
        source_dir.append(float(source_file_lines[line].split()[5]))
        source_spr.append(float(source_file_lines[line].split()[6]))
        source_fp.append(float(source_file_lines[line].split()[7]))
        source_pdir.append(float(source_file_lines[line].split()[8]))
        source_pspr.append(float(source_file_lines[line].split()[9]))
print('')

"""
Creating Destination Grid
"""
# lons = np.arange(31.000, 77.500+.125, 0.125)
# lats = np.arange(8.000, 31.00+.125, 0.125)
lons = sorted(list(set(source_lons))) # Set type is used for removing \
lats = sorted(list(set(source_lats))) # lon and lat duplicate values

"""
Creating NetCDF File
"""
# Creating NetCDF file
print('---> Creating NetCDF file')
nc_file = netCDF4.Dataset('tab44.nc','w', format='NETCDF4')

# Creating global attributes
nc_file.title='NetCDF file created from WW3 ASCII output'
nc_file.source='Wavewatch3'
nc_file.author='Farrokh A. Ghavanini'
nc_file.contact='ghavanini[at]gmail[dot]com'

# Creating dimensions
nc_file.createDimension('lon', len(lons))
nc_file.createDimension('lat', len(lats))
nc_file.createDimension('time', None)

# Creating variables
lon = nc_file.createVariable('lon', np.float64, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
lon.point_spacing = 'even'
lat = nc_file.createVariable('lat', np.float64, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lat.point_spacing = 'even'
time = nc_file.createVariable('time', np.float64, ('time',))
time.units = time_units
time.calendar = 'gregorian'
time.long_name = 'time'
hs = nc_file.createVariable('hs', np.float64,('time','lat','lon'), zlib=True)
hs.units = 'meters'
hs.long_name = 'significant_height_of_wind_and_swell_waves'
wl = nc_file.createVariable('l', np.float64,('time','lat','lon'), zlib=True)
wl.units = 'meters'
wl.long_name = 'sea_surface_wavelength'
tr = nc_file.createVariable('tr', np.float64,('time', 'lat','lon'), zlib=True)
tr.units = 'seconds'
tr.long_name = 'peak_period of waves'
dir = nc_file.createVariable('dir', np.float64,('time', 'lat','lon'), zlib=True)
dir.units = 'degrees'
dir.long_name = 'mean_wave_direction'
dir.reference = 'true_north'
spr = nc_file.createVariable('spr', np.float64,('time', 'lat','lon'), zlib=True)
spr.units = 'degrees'
spr.long_name = 'sea_surface_wave_directional_spread'
fp = nc_file.createVariable('fp', np.float64,('time', 'lat','lon'), zlib=True)
fp.units = 's-1'
fp.long_name = 'sea_surface_wave_peak_frequency'
pdir = nc_file.createVariable('pdir', np.float64,('time', 'lat','lon'), zlib=True)
pdir.units = 'degrees'
pdir.long_name = 'sea_surface_wave_partitioned_mean_direction'
pdir.reference = 'true_north'
pspr = nc_file.createVariable('pspr', np.float64,('time', 'lat','lon'), zlib=True)
pspr.units = 'degrees'
pspr.long_name = 'sea_surface_wave_partitioned_mean_directional_spread'

# Filling variables
lon[:] = lons
lat[:] = lats
time[:] = np.array(range(timesteps)) * 3600
hs[:,:,:] = np.empty((0, len(lats), len(lons)))
wl[:,:,:] = np.empty((0, len(lats), len(lons)))
tr[:,:,:] = np.empty((0, len(lats), len(lons)))
dir[:,:,:] = np.empty((0, len(lats), len(lons)))
spr[:,:,:] = np.empty((0, len(lats), len(lons)))
fp[:,:,:] = np.empty((0, len(lats), len(lons)))
pdir[:,:,:] = np.empty((0, len(lats), len(lons)))
pspr[:,:,:] = np.empty((0, len(lats), len(lons)))

mask = np.array(mask).reshape(len(lats), len(lons))
conditional_helper = np.ones((len(lats), len(lons))).astype(bool)
mask = np.logical_and(conditional_helper, mask==0)  # For inverting mask criteria of original mask file
first_line = 0
for timestep in range(timesteps):
    print('  --> Timestep ' + str(timestep) + ' Processed. (ET: ' + str((dt.datetime.now() - strt_time).seconds) + 's)')
    hs[timestep,:,:] = np.ma.masked_array(np.array(source_hs[first_line:first_line+block_size]).reshape(len(lats), len(lons)), mask, fill_value=np.NaN)
    wl[timestep,:,:] = np.ma.masked_array(np.array(source_wl[first_line:first_line+block_size]).reshape(len(lats), len(lons)), mask, fill_value=np.NaN)
    tr[timestep,:,:] = np.ma.masked_array(np.array(source_tr[first_line:first_line+block_size]).reshape(len(lats), len(lons)), mask, fill_value=np.NaN)
    dir[timestep,:,:] = np.ma.masked_array(np.array(source_dir[first_line:first_line+block_size]).reshape(len(lats), len(lons)), mask, fill_value=np.NaN)
    spr[timestep,:,:] = np.ma.masked_array(np.array(source_spr[first_line:first_line+block_size]).reshape(len(lats), len(lons)), mask, fill_value=np.NaN)
    fp[timestep,:,:] = np.ma.masked_array(np.array(source_fp[first_line:first_line+block_size]).reshape(len(lats), len(lons)), mask, fill_value=np.NaN)
    pdir[timestep,:,:] = np.ma.masked_array(np.array(source_pdir[first_line:first_line+block_size]).reshape(len(lats), len(lons)), mask, fill_value=np.NaN)
    pspr[timestep,:,:] = np.ma.masked_array(np.array(source_pspr[first_line:first_line+block_size]).reshape(len(lats), len(lons)), mask, fill_value=np.NaN)
    first_line += block_size

print('\n --> Done!\n\n     (Total Elapsed Time: ' + str((dt.datetime.now() - strt_time).seconds) + 's)\n')

# Closing NetCDF file
nc_file.close()