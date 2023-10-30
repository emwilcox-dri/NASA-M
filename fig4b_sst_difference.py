import sst
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

sst_path = '/project/NOAA_OI_SST_v2/'
#sst_file = 'sst.mnmean.nc'
sst_file = 'sst.mnmean.180.nc'
outpath = '/Users/ewilcox/NASA_measures/test_output/'
outfile = 'sst_diff'

i_jan2005 = 277
i_dec2005 = 289
i_jan2011 = 349
i_dec2011 = 361

(lons, lats, time, sst) = sst.read_noaa_oi_v2(sst_path,sst_file)

nlons = len(lons)
nlats = len(lats)
print(np.min(lons),np.max(lons))

sst_2005 = np.zeros([nlats,nlons])
sst_2011 = np.zeros([nlats,nlons])
sst_diff = np.zeros([nlats,nlons])

for i in range(nlons):
	for j in range(nlats):
		sst_2005[j,i] = np.average(sst[i_jan2005:i_dec2005,j,i])
		sst_2011[j,i] = np.average(sst[i_jan2011:i_dec2011,j,i])
		sst_diff[j,i] = sst_2005[j,i] - sst_2011[j,i]
		
#Need to mask sst_diff with the same mask as one instance of sst to remove the land from the map

print(np.min(sst_2005),np.max(sst_2005))
print(np.min(sst_2011),np.max(sst_2011))
print(np.min(sst_diff),np.max(sst_diff))

(lon2d,lat2d) = np.meshgrid(lons,lats)

f1 = plt.figure()
ax = plt.axes(projection=ccrs.Mollweide(central_longitude=180))
m=ax.contourf(lon2d, lat2d, sst_diff, transform=ccrs.PlateCarree(),
#					levels=[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2], cmap='PuOr_r')
#					levels=[-1,-0.75,-0.5,-0.25,0,0.5,1,1.5,2], cmap='PuOr_r')
					levels=[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2], cmap='coolwarm')
#					levels=[-1,-0.75,-0.5,-0.25,0,0.5,1,1.5,2], cmap='RdBu_r')
#m=ax.contourf(lon2d, lat2d, sst[i_jan2005,:,:], transform=ccrs.PlateCarree(),
#					levels=[-5,0,5,10,15,20,25,30,35], cmap='PuOr')
ax.set_global()
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
ax.add_feature(land_50m)
ax.coastlines()
plt.colorbar(m,orientation='horizontal')

plt.savefig(outpath+outfile+".png")
plt.savefig(outpath+outfile+".ps")