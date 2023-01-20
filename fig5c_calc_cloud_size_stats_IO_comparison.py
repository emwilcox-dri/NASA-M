#calc_cloud_size_stats.py
#same as v2 using pandas, but for absolute cloud number

import fig5c_roca_ram_sizes as rr
import xarray as xr
import pandas as pd
import numpy as np
import measures
import misc_tools
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

path = '/project/NASA_measures/deep_v5/'
#path = '/project/NASA_measures/deep_v5_025/'
#path = '/Project/NASA_measures/tropical_BR_v2/'

outpath = '/Users/ewilcox/NASA_measures/test_output/'
outfile = 'size_stats'
outfile2 = 'size_stats_total'
#outfile = 'size_stats_20N-20S_min300_16bins_number'
#outfile2 = 'size_stats_total_20N-20S_min300_16bins_number'
#outfile = 'size_stats_20N-20S_025_min10_20bins_number'
#outfile2 = 'size_stats_total_20N-20S_025_min10_20bins_number'

nsurftypes=2
areamin = 300
#areamin = 10
areamax = 2.e+6
nareabins = 16
#nareabins = 20
#nareabins = 24
cols = ['b-', 'k-']
colsbp = ['blue','black']

#Indian Ocean
flag=1
if (flag):
	filelist=[]
	#	yyyymm_list = [200501,200502,200503,200511,200512,201101,201102,201103,201111,201112]
	yrlist = [2004,2005,2006,2007,2008,2009,2010,2011]
	month_list = [7]
	outfile_suffix = '_IO_Jul_2004-2011'
	for iyr in yrlist:
		for imn in month_list:
			mnstr = str(imn)
			if (imn<10):
				mnstr = '0'+mnstr
			filelist.append("{0}modis_clouds.{1}{2}.nc".format(path,iyr,mnstr))
	latmin=-30
	latmax=30
	lonmin=30
	lonmax=100

dropvars = [
#					"cloud_areas_1km",
					"cloud_classification_1km",
#					"fraction_of_cloud_over_land",
#					"maximum_view_zenith_within_cloud_area",
					"minimum_IR_brightness_temperature",
					"CAPE_averaged_over_cloud_area",
					"change_of_CAPE_averaged_over_cloud_area",
					"maximum_CAPE_within_cloud_area",
					"maximum_change_of_CAPE_within_cloud_area",
					"shear_averaged_over_cloud_area",
					"maximum_shear_within_cloud_area",
					"CAPE_at_minimum_IR_brightness_temperature",
					"shear_at_minimum_IR_brightness_temperature",
					"CAPE_at_minimum_85GHz_PCT",
					"shear_at_minimum_85GHz_PCT",
					"PCT_85GHz_averaged_over_cloud_area",
					"minimum_85GHz_PCT",
					"fraction_of_cloud_area_colder_than_220K",
					"fraction_of_cloud_area_colder_than_210K",
					"fraction_of_cloud_area_colder_than_200K",
#					"latitude_of_minimum_IR_brightness_temperature",
#					"longitude_of_minimum_IR_brightness_temperature",
					"latitude_of_minimum_85GHz_PCT",
					"longitude_of_minimum_85GHz_PCT",
					"wind_direction_850hPa_at_minimum_IR_brightness_temperature",
					"wind_direction_500hPa_at_minimum_IR_brightness_temperature",
					"wind_direction_200hPa_at_minimum_IR_brightness_temperature",
					"IR_brightness_temperature_at_location_of_minimum_85GHz_PCT",
					"PCT_85GHz_at_location_of_minimum_IR_brightness_temperature",
					"aerosol_optical_thickness",
					"AOT_MERRA_averaged_over_cloud_area",
					"AOT_MERRA_at_location_of_minimum_IR_brightness_temperature",
					"particle_effective_radius_averaged_on_IR_bins"
					]

if ( len(filelist)==1 and bool(re.search('\*',filelist[0]))==False ):
	clouds_all = xr.open_dataset(filelist[0], drop_variables=dropvars)
elif (len(filelist)==1 and bool(re.search('\*',filelist[0]))==True ):
	clouds_all = xr.open_mfdataset(filelist[0],
							drop_variables=dropvars,concat_dim="number_of_clouds_1km",combine="nested")
else:
	clouds_all = xr.open_mfdataset(filelist,
							drop_variables=dropvars,concat_dim="number_of_clouds_1km",combine="nested")

clouds = clouds_all.where(
													(clouds_all["latitude_of_minimum_IR_brightness_temperature"] >= latmin)
													& (clouds_all["latitude_of_minimum_IR_brightness_temperature"] <= latmax)
													& (clouds_all["longitude_of_minimum_IR_brightness_temperature"] >= lonmin)
													& (clouds_all["longitude_of_minimum_IR_brightness_temperature"] <= lonmax)
													& (clouds_all["maximum_view_zenith_within_cloud_area"]<65)
													).dropna(dim="number_of_clouds_1km")

clouds_df = clouds.to_dataframe()

toffset = pd.to_datetime("19900101",format="%Y%m%d")
utc_time = toffset+pd.to_timedelta(clouds_df['minutes_since_Jan_1_1990_UTC'],unit='minutes')
clouds_df['UTC time'] = utc_time

hour_offset = np.floor( (7.5 + clouds_df['longitude_of_minimum_IR_brightness_temperature'])/15. )
local_time = utc_time + pd.to_timedelta(hour_offset, unit="hours")
clouds_df["local time"] = local_time

day_b = pd.cut(pd.DatetimeIndex(clouds_df["local time"]).hour, bins=[0,6,18,24],
							include_lowest=True, labels=['night','day','night'], ordered=False)
land_b = pd.cut(clouds_df.fraction_of_cloud_over_land, bins=[0,0.5,1], include_lowest=True,
						labels=['ocean','land'])

clouds_df["log10_cloud_area"] = np.log10(clouds_df["cloud_areas_1km"])
logareabins = np.arange(nareabins)*( np.log10(areamax) - np.log10(areamin) )/(nareabins-1) + np.log10(areamin)
area_b = pd.cut(clouds_df['log10_cloud_area'], bins=logareabins, include_lowest=True)
grouped = clouds_df['log10_cloud_area'].groupby([area_b,land_b])
areas_cnt = grouped.size()

dlogbin = logareabins[2]-logareabins[1]
areaout = np.power( 10, logareabins[0:nareabins-1] + (dlogbin/2.) )

nclouds=[]
total_area=[]
mean_size=[]
median_size=[]

(rr_area, rr_num_jan, rr_num_jul) = rr.get_roca_ram_sizes()
		
f1 = plt.figure()
s1 = f1.add_subplot()

# 0=ocean/day; 1=ocean/night; 2=land/day; 3=land/night
#for land_flag in range(nsurftypes):
#	for day_flag in range(2):
land_flag=0 #plot only the ocean clouds

#varout = areas_cnt[np.arange(nareabins-1)*4 + (land_flag*2+day_flag)]
varout = areas_cnt[np.arange(nareabins-1)*2 + land_flag]

#note that currently hard-coded for the Roca and Ram January
#s1.plot( areaout,varout/np.sum(varout),cols[0], rr_area, rr_num_jan/np.sum(rr_num_jan), cols[1] )
#note that currently hard-coded for the Roca and Ram June
s1.plot( areaout,varout/np.sum(varout),cols[0], rr_area, rr_num_jul/np.sum(rr_num_jul), cols[1] )

s1.set_ylabel('number of clouds')
s1.set_xlabel('cloud area (km^2)')
s1.set_xscale('log')
if (areamin<100):
	s1.set_xlim([areamin,areamax])
else:
	s1.set_xlim([100,areamax])
#		s1.set_ylim([1,3.e+5])
#s1.set_ylim([1,1.e+6])
s1.set_ylim([3.e-6,0.8])
s1.set_yscale('log')
    
print("saving to {0}".format(outpath+outfile+outfile_suffix+".png"))    
plt.savefig(outpath+outfile+outfile_suffix+".png")
print("saving to {0}".format(outpath+outfile+outfile_suffix+".ps"))  
plt.savefig(outpath+outfile+outfile_suffix+".ps")



