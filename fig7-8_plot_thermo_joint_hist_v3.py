#plot_thermo_joint_hist.py
#version updated Sept 2020 for use with the Measures dataset

import xarray as xr
import pandas as pd
import numpy as np
import measures
import re
#import terraqua
import math
import colormaps
from datetime import datetime
from numpy import ma
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as col

#mcs_flag = 0 #0=all clouds, 1=MCS\
#aerosol_flag = 0 #0=day/night, 1=clean/polluted
#sat_model_flag = 0  # 0=satellite  1=GCE model

comp_root = '/Project/'
#comp_root = '/discover/nobackup/projects/cloudclass/measures_deep/deep_v5/'
  
#path = comp_root + 'NASA_measures/tropical_test_v4/'
path = comp_root + 'NASA_measures/deep_v6/nc4/'
file_root = 'MODIS_Aqua_L2_DC.V001'

outpath = '/Users/ewilcox/NASA_measures/test_output/'

outfile_suffix = ''

	
#global
flag=0
if (flag):
	filelist = [path+file_root+'.2005*.nc']
	latmin=-30
	latmax=30
	lonmin=-180
	lonmax=180
	outfile_suffix = '_global_30S-30N_2005'

#Indian Ocean/India Summer
flag=0
if (flag):
	yyyymm_list = [200506,200507,200508,201106,201107,201108]
	latmin=0
	latmax=35
	lonmin=65
	lonmax=95

#Indian Ocean Winter
flag=0
if (flag):
	filelist=[]
	yyyymm_list = [200501,200502,200503,200511,200512,201101,201102,201103,201111,201112]
	yrlist = [2004,2005,2006,2007,2008,2010,2011]
	month_list = [1,2,3,11,12]
	for iyr in yrlist:
		for imn in month_list:
			mnstr = str(imn)
			if (imn<10):
				mnstr = '0'+mnstr
			filelist.append("{0}{1}.{2}{3}.nc".format(path,file_root,iyr,mnstr))
	latmin=-10
	latmax=15
	lonmin=50
	lonmax=100
	outfile_suffix = '_IO_Nov-Mar_2004-2011'

#Warm Pool
flag=0
if (flag):
	latmin=-20
	latmax=20
	lonmin=110
	lonmax=180
	outfile_suffix = '_WP_2005'

#Nino 3.4
flag=0
if (flag):
	latmin=-5
	latmax=5
	lonmin=-170
	lonmax=-120

#Amazon
flag=1
if (flag):
	filelist=[]
	yyyymm_list = [200501,200502,200503,200511,200512,201101,201102,201103,201111,201112]
	yrlist = [2004,2005,2006,2007,2008,2010,2011]
	month_list = [1,2,3,11,12]
	for iyr in yrlist:
		for imn in month_list:
			mnstr = str(imn)
			if (imn<10):
				mnstr = '0'+mnstr
			filelist.append("{0}{1}.{2}{3}.nc".format(path,file_root,iyr,mnstr))
	latmin=-10
	latmax=0
	lonmin=-65
	lonmax=-50
	outfile_suffix = '_Amazon_Nov-Mar_2004-2011'
	




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
#					,"max_CAPE"
					,"max_change_of_CAPE"
#					,"max_shear"
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

	
ratio_cloudcnt_thresh = 9 #must be more than this many clouds to evaluate ratio
    
polluted_ocean = [0.3,5.]
clean_ocean = [0.,0.3]
polluted_land = [0.,0.]
clean_land = [0.,0.]

if ( len(filelist)==1 and bool(re.search('\*',filelist[0]))==False ):
	clouds_all = xr.open_dataset(filelist[0], drop_variables=dropvars)
elif (len(filelist)==1 and bool(re.search('\*',filelist[0]))==True ):
	clouds_all = xr.open_mfdataset(filelist[0],
							drop_variables=dropvars,concat_dim="number_of_clouds",combine="nested")
else:
	clouds_all = xr.open_mfdataset(filelist,
							drop_variables=dropvars,concat_dim="number_of_clouds",combine="nested")
print("# of clouds_all = {0}".format(clouds_all["cloud_areas"].values.shape[0]))

clouds = clouds_all.where(
													(clouds_all["latitude_of_min_IR"] >= latmin)
													& (clouds_all["latitude_of_min_IR"] <= latmax)
													& (clouds_all["longitude_of_min_IR"] >= lonmin)
													& (clouds_all["longitude_of_min_IR"] <= lonmax)
#													& (clouds_all["max_view_zenith"]<53)
													& (clouds_all["max_view_zenith"]<65)
													).dropna(dim="number_of_clouds")
print("# of IR clouds = {0}".format(clouds["cloud_areas"].values.shape[0]))

clouds_df = clouds.to_dataframe()

hour_offset = np.floor( (7.5 + clouds_df["longitude_of_min_IR"])/15. )
local_time = clouds_df["time"] + pd.to_timedelta(hour_offset, unit="hours")
clouds_df["local time"] = local_time

print("# of clouds = {0}".format(len(clouds_df.index)))

shear_rang = [0,80]
cape_rang = [0,5000]
#shear_rang = [0,40]
#cape_rang = [0,2400]
ncapebins = 20
nshearbins = 20
dshear = shear_rang[1]/nshearbins
dcape = cape_rang[1]/ncapebins
shear_binbnds = np.arange(nshearbins+1)*dshear
cape_binbnds = np.arange(ncapebins+1)*dcape
shear_bins = (np.arange(nshearbins)*dshear) + (dshear/2.)
cape_bins = (np.arange(ncapebins)*dcape) + (dcape/2.)
areamean_daynight = np.zeros([ncapebins,nshearbins,4])
areacnt_daynight = np.zeros([ncapebins,nshearbins,4])
areamean_clnpltd = np.zeros([ncapebins,nshearbins,4])
areacnt_clnpltd = np.zeros([ncapebins,nshearbins,4])

cape_b = pd.cut(clouds_df.max_CAPE, bins=cape_binbnds)
shear_b = pd.cut(clouds_df.max_shear, bins=shear_binbnds)
day_b = pd.cut(pd.DatetimeIndex(clouds_df["local time"]).hour, bins=[0,6,18,24],
							include_lowest=True, labels=['night','day','night'], ordered=False)
land_b = pd.cut(clouds_df.fraction_over_land, bins=[0,0.5,1], include_lowest=True,
						labels=['ocean','land'])


grouped = clouds_df['cloud_areas'].groupby([cape_b,shear_b,land_b,day_b])
areas_day = grouped.mean()
areas_day_cnt = grouped.size()

# k = 0=day/ocean; 1=night/ocean; 2=day/land; 3=night/land
# k = 0=clean/ocean; 1=polluted/ocean; 2=clean/land; 3=polluted/land
# areamean[ncapebins,nshearbins,4]
for i in range(ncapebins):
	for k in range(4):
		j = np.arange(nshearbins)*4 + (i*nshearbins*4) +k
		areamean_daynight[i,:,k] = areas_day[j]
		areacnt_daynight[i,:,k] = areas_day_cnt[j]


				
#calculate ratios of mean cloud area
ratiomeans_daynight = np.zeros([ncapebins,nshearbins,4])
ratio_k = np.array([[1,3,3,2],[0,2,1,0]])
for i in range(ncapebins):
	for j in range(nshearbins):
		for k in range(4):
		
		  if ( (areacnt_daynight[i,j,ratio_k[0,k]]>ratio_cloudcnt_thresh) and 
		        (areacnt_daynight[i,j,ratio_k[1,k]]>ratio_cloudcnt_thresh) ):
#		    print(areacnt[i,j,ratio_k[0,k]],areamean[i,j,ratio_k[0,k]],areacnt[i,j,ratio_k[1,k]],areamean[i,j,ratio_k[1,k]])
		    ratiomeans_daynight[i,j,k] = math.log10(areamean_daynight[i,j,ratio_k[0,k]]/areamean_daynight[i,j,ratio_k[1,k]])
		  else:
		    ratiomeans_daynight[i,j,k] = 9.e+20
		  


#convert areamean to log_base_10
for i in range(ncapebins):
	for j in range(nshearbins):
		for k in range(4):
		
			if (areacnt_daynight[i,j,k]>0):
				areamean_daynight[i,j,k] = math.log10(areamean_daynight[i,j,k])
			else:
				areamean_daynight[i,j,k] = 9.e+20



#create blue/green/brown stepped sequential colormap
(r1,g1,b1) = colormaps.get_stepped_sequential_colormap()
r = r1[5:20]
r = r[::-1]
g = g1[5:20]
g = g[::-1]
b = b1[5:20]
b = b[::-1]
sscmap = col.ListedColormap(np.transpose(np.array([r,g,b])),'ss_cmap')
sscmap.set_over(color='white')
sscmap.set_under(color='white')
sscmap.set_bad(color='white')

# x and y axis labeling
nxticks = 5 #must divide evenly into nshearbins
xticklocs = np.arange(nxticks,dtype=int)*int(nshearbins/nxticks)
xtickvals = shear_bins[xticklocs]
xticklabels=[]
for i in range(nxticks):
	xticklabels.append("{0:.0f}".format(shear_binbnds[xticklocs[i]]))
xticklocs = xticklocs-0.5
nyticks = 5 #must divide evenly into ncapebins
yticklocs = np.arange(nyticks,dtype=int)*int(ncapebins/nyticks)
ytickvals = cape_bins[yticklocs]
yticklabels=[]
for i in range(nyticks):
	yticklabels.append("{0:.0f}".format(cape_binbnds[yticklocs[i]]))
yticklocs = yticklocs-0.5



# colorbar labeling
flag=0
if (flag):
	amin = 1.e+3
	amax = 1.e+6
	atickl = [1.e+3, 1.e+4, 1.e+5, 1.e+6]
	#amax = 1.e+4
	#atickl = [1.e+3, 2150, 4630, 1.e+4]
	aticks = []
	for i in atickl:
		aticks[len(aticks):] = [math.log10(i)]
	atickout = []
	for i in aticks:
		atickout[len(atickout):] = ["10^{0:.0f}".format(i)]

flag=1
if (flag):
	amin = 500.
	amax = 5.e+5
	atickl = [500, 5.e+3, 5.e+4, 5.e+5]
	atickout = ['500', '5000', '5x10^4', '5x10^5']
	aticks = []
	for i in atickl:
		aticks[len(aticks):] = [math.log10(i)]



areamean = areamean_daynight
outfile = 'thermo_joint_hist'

f1 = plt.figure()

for k in range(4):

	if (k==0):
		s1 = f1.add_subplot(331)
		pos=[0.1, 0.65, 0.3, 0.3]  #[left, bottom, width, height]
		print(areamean[:,0,k])
		print(areamean[0,:,k])
		titl='ocean/day'
	if (k==1):
		s1 = f1.add_subplot(332)
		pos=[0.45, 0.65, 0.3, 0.3]  #[left, bottom, width, height]
		titl='ocean/night'
	if (k==2):
		s1 = f1.add_subplot(334)
		pos=[0.1, 0.23, 0.3, 0.3]  #[left, bottom, width, height]
		titl='land/day'
	if (k==3):
		s1 = f1.add_subplot(335)
		pos=[0.45, 0.23, 0.3, 0.3]  #[left, bottom, width, height]
		titl='land/night'
	plt.title(titl,fontsize=9)
	s1.set_position(pos)
	areaout = ma.masked_greater(areamean[:,:,k],7)
	plt.imshow(areaout,sscmap,zorder=0,interpolation='nearest',
							vmin=math.log10(amin),vmax=math.log10(amax),origin='lower')
	plt.grid()
	plt.xlabel('shear (m s^-1)',fontsize=9)
	plt.xticks(xticklocs,xticklabels,fontsize=9)
	plt.ylabel('CAPE (J kg^-1)',fontsize=9)
	plt.yticks(yticklocs,yticklabels,fontsize=9)

ax5=f1.add_axes([0.3, 0.1, 0.25, 0.03])
norm = mpl.colors.Normalize(vmin=math.log10(amin), vmax=math.log10(amax))
cbar=mpl.colorbar.ColorbarBase(ax5,cmap=sscmap,norm=norm,orientation='horizontal',
															ticks=aticks)
cbar.set_ticklabels(atickout)
cbar.set_label('cloud area (km^2)',fontsize=9)
ll,bb,ww,hh = cbar.ax.get_position().bounds
cbar.ax.set_position([0.18, 0.1, 0.5, hh])

ofile = outpath+outfile+outfile_suffix+'.png'
print('saving plot to ' + ofile)
plt.savefig(ofile)
ofile = outpath+outfile+outfile_suffix+'.ps'
print('saving plot to ' + ofile)
plt.savefig(ofile)
