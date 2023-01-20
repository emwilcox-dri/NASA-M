#map_cloud_distribution_difference_v2.py
# modified to work with the nc4

import xarray as xr
import pandas as pd
import measures
import grid
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt

path = '/project/NASA_measures/deep_v6/nc4/'
file_root = 'MODIS_Aqua_L2_DC.V001'

outpath = '/Users/ewilcox/NASA_measures/test_output/'
outfile = 'cloud_distribution_diff'

latmin = -30
latmax = 30
lonmin = -180 #stick with -180 to 180 here to avoid having to change the longitude values for millions of clouds
lonmax = 180
#lonmin = 0
#lonmax = 360
dlon = 1
dlat = 1

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
					,"fraction_over_land"
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


clouds_2005_all = xr.open_mfdataset(path+file_root+".2005*.nc",
								drop_variables=dropvars,concat_dim="number_of_clouds",combine="nested")

clouds_2005 = clouds_2005_all.where(
													(clouds_2005_all["latitude_of_min_IR"] >= latmin)
													& (clouds_2005_all["latitude_of_min_IR"] <= latmax)
													& (clouds_2005_all["longitude_of_min_IR"] >= lonmin)
													& (clouds_2005_all["longitude_of_min_IR"] <= lonmax)
#													& ( ((clouds_2005_all["time"] >= timemin) &
#													(clouds_2005_all["time"] <= timemax)) |
#													((clouds_2005_all["time"] >= timemin1) &
#													(clouds_2005_all["time"] <= timemax1)) )
#													& (clouds_2005_all["max_view_zenith"]<53)
													& (clouds_2005_all["max_view_zenith"]<65)
													).dropna(dim="number_of_clouds")

clouds_2005_df = clouds_2005.to_dataframe()
print("# of 2005 clouds = {0}".format(len(clouds_2005_df.index)))
								
clouds_2011_all = xr.open_mfdataset(path+file_root+".2011*.nc",
								drop_variables=dropvars,concat_dim="number_of_clouds",combine="nested")

clouds_2011 = clouds_2011_all.where(
													(clouds_2011_all["latitude_of_min_IR"] >= latmin)
													& (clouds_2011_all["latitude_of_min_IR"] <= latmax)
													& (clouds_2011_all["longitude_of_min_IR"] >= lonmin)
													& (clouds_2011_all["longitude_of_min_IR"] <= lonmax)
#													& ( ((clouds_2011_all["time"] >= timemin) &
#													(clouds_2011_all["time"] <= timemax)) |
#													((clouds_2011_all["time"] >= timemin1) &
#													(clouds_2011_all["time"] <= timemax1)) )
#													& (clouds_2011_all["max_view_zenith"]<53)
													& (clouds_2011_all["max_view_zenith"]<65)
													).dropna(dim="number_of_clouds")

clouds_2011_df = clouds_2011.to_dataframe()
print("# of 2011 clouds = {0}".format(len(clouds_2011_df.index)))

(lons, lats, areas) = grid.get_grid(lonmin, latmin, lonmax, latmax, dlon, dlat)
nlons = lons.size
nlats = lats.size
lonbnds = np.concatenate((lons-(dlon/2.),[lons[nlons-1]+(dlon/2.)]))
latbnds = np.concatenate((lats-(dlat/2.),[lats[nlats-1]+(dlat/2.)]))#	clouds_df grouped by lons from -180 to 180. Here I transfer it to a numpy array with
# lons arranged 0 to 360
ilons_east = list(range(int(nlons/2)))
ilons_west = [x+180 for x in ilons_east]

nclouds = np.zeros([nlats,nlons],dtype=np.single)
lat_b = pd.cut(clouds_2005_df.latitude_of_min_IR, bins=latbnds)
lon_b = pd.cut(clouds_2005_df.longitude_of_min_IR, bins=lonbnds)
grouped_2005 = clouds_2005_df['cloud_areas'].groupby([lat_b, lon_b]).size()
for ilat in range(nlats):
	nclouds[ilat,ilons_east] = grouped_2005[[x+(ilat*nlons)+180 for x in ilons_east]]
	nclouds[ilat,ilons_west] = grouped_2005[[x+(ilat*nlons) for x in ilons_east]]

nclouds1 = np.zeros([nlats,nlons],dtype=np.single)
lat_b = pd.cut(clouds_2011_df.latitude_of_min_IR, bins=latbnds)
lon_b = pd.cut(clouds_2011_df.longitude_of_min_IR, bins=lonbnds)
grouped_2011 = clouds_2011_df['cloud_areas'].groupby([lat_b, lon_b]).size()
for ilat in range(nlats):
	nclouds1[ilat,ilons_east] = grouped_2011[[x+(ilat*nlons)+180 for x in ilons_east]]
	nclouds1[ilat,ilons_west] = grouped_2011[[x+(ilat*nlons) for x in ilons_east]]

lons_new = np.zeros([nlons])
lons_new[ilons_east] = lons[ilons_west]
lons_new[ilons_west] = lons[ilons_east]+360
(lon2d,lat2d) = np.meshgrid(lons_new,lats)
#(lon2d,lat2d) = np.meshgrid(lons,lats)
#print(lon2d.shape,lat2d.shape)

#nclouds = nclouds
#nclouds = nclouds/12.

diff_nclouds = np.zeros([nlats,nlons],dtype=np.single)
ndays = clouds_2005_df['time'].dt.date.unique().size
for i in range(nlats):
	for j in range(nlons):
		diff_nclouds[i,j] = (nclouds[i,j]-nclouds1[i,j])/(ndays*2/(365./12.))/(areas[i]/10000.) # num of clouds per month per million hectares
print("min(diff_nclouds)={0} max(diff_nclouds)={1}".format(np.min(diff_nclouds),np.max(diff_nclouds)))

f1 = plt.figure()
#s1 = f1.add_subplot(121)
#ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
ax = plt.axes(projection=ccrs.Mollweide(central_longitude=180))
#ax = plt.axes(projection=ccrs.Mollweide())
m=ax.contourf(lon2d, lat2d, diff_nclouds, transform=ccrs.PlateCarree(),
#ax.contourf(lon2d, lat2d, nclouds, transform=ccrs.PlateCarree(central_longitude=180),
#ax.contourf(lon2d, lat2d, nclouds,
#ax.contourf(lon2d, lat2d, nclouds, transform=ccrs.Mollweide(central_longitude=0),
#					levels=[-8,-6,-4,-2,0,2,4,6,8], cmap='PuOr')
					levels=[-8,-6,-4,-2,0,2,4,6,8], cmap='BrBG')
ax.set_global()
ax.coastlines()
plt.colorbar(m,orientation='horizontal')

#(h,b) = np.histogram(nclouds)
#print(h)
#print(b)

plt.savefig(outpath+outfile+".png")
plt.savefig(outpath+outfile+".ps")


