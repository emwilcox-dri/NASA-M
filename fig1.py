#map_das_steps_v2.py
#
# same as the version in terraqua/tools, but updated Apr. 2021 to use cartopy
# Further updated Dec. 2021 to take netcdf file with grid and gridded variables.
#
# netcdf files needed for this code are in data_for_fig1.zip
#
#####################

import xarray as xr
import grid
import terraqua
import numpy as np
#from numpy import ma, fromfile, reshape, float32, min, max, histogram
import matplotlib as mpl
mpl.use('Agg')
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

path_root = '/Project/NASA_measures/'
#yyyyddd = '2011027'
#hhmm = '1720'
#yyyyddd = '2011010'
#hhmm = '2140'
#hhmm = '1640'
#hhmm = '0825'
#hhmm = '0645'
#hhmm = '0505'
yyyyddd = '2006199'
hhmm = '0815'
#yyyyddd = '2006189'
#hhmm = '0730'
#thresh = 'BR'
thresh = 'deep'
#thresh = 'none'

gridfile = 'cloud_labels.{0}.{1}.nc'.format(yyyyddd,hhmm)
path = "{0}das_steps_{1}.{2}/".format(path_root,yyyyddd,hhmm)
#file = 'cloud_labels.{0}.{1}.nc'.format(yyyyddd,hhmm)
#file = 'detect_step_2.{0}.{1}.nc'.format(yyyyddd,hhmm)

if (thresh=="BR"):
	path = path+"BR/"
	nstep=4
if (thresh=="deep"):
	path = path+"deep/"
	nstep=3

outpath = path
grid = xr.open_dataset(path+gridfile)

lonmin = np.min(grid.longitude.values)
lonmax = np.max(grid.longitude.values)
latmin = np.min(grid.latitude.values)
latmax = np.max(grid.latitude.values)
extents = [lonmin,lonmax,latmin,latmax]
(lon2d,lat2d) = np.meshgrid(grid.longitude.values,grid.latitude.values,indexing='ij')

Tb31 = grid.Tb_IR.values
Tb31[np.where(Tb31>305)]=-9999.
varmap = np.ma.masked_less_equal(Tb31,-9900)
print(np.min(varmap),np.max(varmap))
#	varmap = np.ma.masked_greater(Tb31,280)
#	map_range=[200,280]
map_range=[190,308]
outfile='tb31.png'
outfile2='tb31.ps'
titl='MODIS channel 31 Tb (K)'
#	colors='YlOrRd'
#	colors='gist_stern_r'
colors='gist_ncar_r'
#	colors='cubehelix_r'
#	colors='nipy_spectral'
cbarflag=1
contourflag=0

f1 = plt.figure()
#varmap = np.transpose(varmap)
ax = plt.axes(projection=ccrs.PlateCarree())

#	m=ax.pcolormesh(labels.longitude.values, labels.latitude.values,
#m=ax.pcolormesh(lon2d, lat2d,
#								varmap, transform=ccrs.PlateCarree(), 
#								vmin=map_range[0], vmax=map_range[1], cmap=colors,
#								shading='nearest'
#								)
m = ax.imshow(varmap, extent=extents, transform=ccrs.PlateCarree(), cmap=colors)
#ax.set_extent(extents, crs=ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
ax.coastlines()
plt.colorbar(m,orientation='vertical')
print("saving file to "+outpath+outfile)
plt.savefig(outpath+outfile)
plt.savefig(outpath+outfile2)

tab20 =  plt.get_cmap('tab20',20)
white = np.array([1.,1.,1.,1.])
newc = np.empty((21,4))
newc[1:21,:] = tab20(np.linspace(0,1,20))
newc[0,:] = white
colors = ListedColormap(newc)

for step in range(nstep):
	for id in ['detect','spread']:
	
		file = '{0}_step_{1}.{2}.{3}.nc'.format(id,step,yyyyddd,hhmm)

		labels = xr.open_dataset(path+file)
		
	#	varmap = np.ma.masked_less_equal(labels.cloud32.values,0)
		varmap = np.ma.masked_less_equal(labels.cloud_labels.values,0)
		varmap = np.mod(varmap,20)+1
		map_range=[1,np.max(varmap)]
		outfile='{0}_step-{1}.{2}.{3}.png'.format(id,step,yyyyddd,hhmm)
		outfile2='{0}_step-{1}.{2}.{3}.ps'.format(id,step,yyyyddd,hhmm)
		titl=''
				
		f1 = plt.figure()
#		varmap = np.transpose(varmap)
		ax = plt.axes(projection=ccrs.PlateCarree())
#		m=ax.pcolormesh(lon2d, lat2d,
#										varmap, transform=ccrs.PlateCarree(), 
#										vmin=map_range[0], vmax=map_range[1], cmap=colors,
#										shading='nearest'
#										)
		m = ax.imshow(varmap, extent=extents, transform=ccrs.PlateCarree(), cmap=colors)
#		ax.set_extent(extents, crs=ccrs.PlateCarree())
		gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
											linewidth=2, color='gray', alpha=0.5, linestyle='--')
		gl.top_labels = False
		gl.right_labels = False
		ax.coastlines()

		print("saving file to "+outpath+outfile)
		plt.savefig(outpath+outfile)
		plt.savefig(outpath+outfile2)





