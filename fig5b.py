#calc_cloud_size_stats.py
#same as v2 using pandas, but for absolute cloud number

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
areamax = 1.e+6
nareabins = 16
#nareabins = 20
#nareabins = 24
cols = ['b-', 'k-']
colsbp = ['blue','black']

filelist = [path+'modis_clouds.2005*.nc']
latmin=-20
latmax=20
lonmin=-180
lonmax=180
outfile_suffix = '_20N-20S'

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
grouped = clouds_df['log10_cloud_area'].groupby([area_b,land_b,day_b])
areas_cnt = grouped.size()

dlogbin = logareabins[2]-logareabins[1]
areaout = np.power( 10, logareabins[0:nareabins-1] + (dlogbin/2.) )

varout = areas_cnt[np.arange(nareabins-1)*4 + 3]
#for i in range(nareabins-1):
#	print("{0}, {1}".format(areaout[i],varout[i]))

f1 = plt.figure()
s1 = f1.add_subplot()
tmp = np.loadtxt("data_for_fig5b.csv",delimiter=",",dtype=str)
area = tmp[1:,0].astype(np.float)
num = tmp[1:,1].astype(np.float)
print(area)
print(num)
s1.plot(areaout,varout,'k-', area,num,'k--')

s1.set_ylabel('number of clouds')
s1.set_xlabel('cloud area (km^2)')
s1.set_xscale('log')
s1.set_xlim([10,1.e+6])
#		s1.set_ylim([1,3.e+5])
s1.set_ylim([1,2.e+6])
s1.set_yscale('log')
print("saving to {0}".format(outpath+outfile2+outfile_suffix+".png"))    
plt.savefig(outpath+outfile2+outfile_suffix+".png")
print("saving to {0}".format(outpath+outfile2+outfile_suffix+".ps"))  
plt.savefig(outpath+outfile2+outfile_suffix+".ps")



