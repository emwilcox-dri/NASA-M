#map_cloud_distribution_v3.py
# modified to work with the nc4

import xarray as xr
import pandas as pd
import measures
import grid
import numpy as np
from datetime import datetime
import matplotlib as mpl
mpl.use('Agg')
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt

comp_root = '/Project/'
path = comp_root + 'NASA_measures/deep_v6/nc4/'

#file_root = 'modis_clouds'
file_root = 'MODIS_Aqua_L2_DC.V001'

outpath = '/Users/ewilcox/NASA_measures/test_output/'
outfile = 'cloud_distribution_2005'

latmin = -40
latmax = 40
lonmin = -180 #stick with -180 to 180 here to avoid having to change the longitude values for millions of clouds
lonmax = 180
#lonmin = 0
#lonmax = 360

nsamples_per_day = 2. # set 1 if limiting to day or night; 2 if using both day and night

dropvars = [
					"AOT_MODIS"
					,"AOT_MERRA_at_min_IR"
					,"AOT_MERRA"
					,"CAPE_at_min_PCT"
					,"CAPE_at_min_IR"
					,"CAPE"
					,"change_of_CAPE"
#					,"cloud_areas"
					,"cloud_classification"
					,"fraction_colder_than_200K"
					,"fraction_colder_than_210K"
					,"fraction_colder_than_220K"
#					,"fraction_over_land"
					,"IR_Tb_at_min_PCT"
					,"latitude_of_min_PCT"
#					,"latitude_of_min_IR"
					,"longitude_of_min_PCT"
#					,"longitude_of_min_IR"
					,"max_CAPE"
					,"max_change_of_CAPE"
					,"max_shear"
#					,"max_view_zenith"
					,"min_PCT"
					,"min_IR"
#					,"time"
					,"particle_effective_radius"
					,"PCT_at_min_IR"
					,"PCT"
					,"shear_at_min_PCT"
					,"shear_at_min_IR"
					,"shear"
					,"wind_direction_200hPa_at_min_IR"
					,"wind_direction_500hPa_at_min_IR"
					,"wind_direction_850hPa_at_min_IR"
					]

#clouds_all = xr.open_dataset(path+"modis_clouds.200501.nc", drop_variables=dropvars)
#clouds_all = xr.open_mfdataset(path+"modis_clouds.2005*.nc",
#								drop_variables=dropvars,concat_dim="number_of_clouds_1km",combine="nested")
clouds_all = xr.open_mfdataset(path+file_root+".2005*.nc",
								drop_variables=dropvars,concat_dim="number_of_clouds",combine="nested")

timemin = (datetime(2005,1,1,0,0,0) - datetime(1990,1,1,0,0,0)).days*24*60
timemax = (datetime(2005,3,31,0,0,0) - datetime(1990,1,1,0,0,0)).days*24*60
timemin1 = (datetime(2011,1,1,0,0,0) - datetime(1900,1,1,0,0,0)).days*24*60
timemax1 = (datetime(2011,3,31,0,0,0) - datetime(1990,1,1,0,0,0)).days*24*60

clouds = clouds_all.where(
													(clouds_all["latitude_of_min_IR"] >= latmin)
													& (clouds_all["latitude_of_min_IR"] <= latmax)
													& (clouds_all["longitude_of_min_IR"] >= lonmin)
													& (clouds_all["longitude_of_min_IR"] <= lonmax)
#													& ( ((clouds_all["time"] >= timemin) &
#													(clouds_all["time"] <= timemax)) |
#													((clouds_all["time"] >= timemin1) &
#													(clouds_all["time"] <= timemax1)) )
#													& (clouds_all["max_view_zenith"]<53)
													& (clouds_all["max_view_zenith"]<65)
													).dropna(dim="number_of_clouds")

clouds_df = clouds.to_dataframe()

#toffset = pd.to_datetime("19900101",format="%Y%m%d")
#utc_time = toffset+pd.to_timedelta(clouds_df['time'],unit='minutes')
#clouds_df['UTC time'] = utc_time

hour_offset = np.floor( (7.5 + clouds_df['longitude_of_min_IR'])/15. )
#local_time = utc_time + pd.to_timedelta(hour_offset, unit="hours")
local_time = clouds_df['time'] + pd.to_timedelta(hour_offset, unit="hours")
clouds_df["local time"] = local_time

print("# of clouds = {0}".format(len(clouds_df.index)))
		
dlon = 1
dlat = 1
(lons, lats, areas) = grid.get_grid(lonmin, latmin, lonmax, latmax, dlon, dlat)
nlons = lons.size
nlats = lats.size
lonbnds = np.concatenate((lons-(dlon/2.),[lons[nlons-1]+(dlon/2.)]))
latbnds = np.concatenate((lats-(dlat/2.),[lats[nlats-1]+(dlat/2.)]))

lat_b = pd.cut(clouds_df.latitude_of_min_IR, bins=latbnds)
lon_b = pd.cut(clouds_df.longitude_of_min_IR, bins=lonbnds)
grouped = clouds_df['cloud_areas'].groupby([lat_b, lon_b]).size()

nclouds = np.zeros([nlats,nlons],dtype=np.single)

#	clouds_df grouped by lons from -180 to 180. Here I transfer it to a numpy array with
# lons arranged 0 to 360
ilons_east = list(range(int(nlons/2)))
ilons_west = [x+180 for x in ilons_east]
#ilons_east = np.indices([int(nlons/2)])
#ilons_west = np.indices([int(nlons/2)])+180
for ilat in range(nlats):
	nclouds[ilat,ilons_east] = grouped[[x+(ilat*nlons)+180 for x in ilons_east]]
	nclouds[ilat,ilons_west] = grouped[[x+(ilat*nlons) for x in ilons_east]]

lons_new = np.zeros([nlons])
lons_new[ilons_east] = lons[ilons_west]
lons_new[ilons_west] = lons[ilons_east]+360

(lon2d,lat2d) = np.meshgrid(lons_new,lats)
#print(lon2d.shape,lat2d.shape)

#ndays = np.round((np.max(clouds_df['UTC time'])-np.min(clouds_df['UTC time'])).total_seconds()/(60*60*24)) #works, but will give the wrong answer if the data is discontinuous
ndays = clouds_df['time'].dt.date.unique().size
print("ndays={0}".format(ndays))
print("min(nclouds)={0} max(nclouds)={1}".format(np.min(nclouds),np.max(nclouds)))
for i in range(nlats):
	for j in range(nlons):
		nclouds[i,j] = nclouds[i,j]/(ndays*nsamples_per_day/(365./12.))/(areas[i]/10000.) # num of clouds per month per million hectares
print("min(nclouds)={0} max(nclouds)={1}".format(np.min(nclouds),np.max(nclouds)))

f1 = plt.figure()
#s1 = f1.add_subplot(121)
#ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
proj = ccrs.Mollweide(central_longitude=180)
#data_proj = ccrs.PlateCarree(central_longitude=180)
data_proj = ccrs.PlateCarree()
ax = plt.axes(projection=proj)
#ax = plt.axes(projection=ccrs.Mollweide())
#m=ax.contourf(lon2d, lat2d, nclouds, 
#m=ax.contourf(lon2d, lat2d, nclouds, transform=ccrs.PlateCarree(),
m=ax.contourf(lon2d, lat2d, nclouds, transform=data_proj,
					levels=[5,10,15,20,25], 
#					levels=[20,40,60,80,100,120], 
					cmap='YlOrRd')
ax.set_global()
ax.coastlines()
plt.colorbar(m,orientation='horizontal')

#(h,b) = np.histogram(nclouds)
#print(h)
#print(b)

plt.savefig(outpath+outfile+".png")
plt.savefig(outpath+outfile+".ps")


