################################
# like plot_structure_stats_v2.py but modified for the nc4 data files.

import measures
import misc_tools
import xarray as xr
import pandas as pd
import numpy as np
import re
#from numpy import ma
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

path = '/project/NASA_measures/deep_v6/nc4/'
file_root = 'MODIS_Aqua_L2_DC.V001'

outpath = '/Users/ewilcox/NASA_measures/test_output/'
outfile = 'structure_stats'

#mcs_flag = 0 #0=all clouds, 1=MCS
ocean_only_flag = 0 #0=all clouds, 1=ocean clouds only

if (ocean_only_flag):
	nsurftypes=1
else:
	nsurftypes=2
areamin = 400
areamax = 1.e+6
nareabins = 14
cols = ['b-', 'k-']
colsbp = ['blue','black']

#global tropics
flag=1
if (flag):
	filelist = [path+file_root+'.2005*.nc']
	latmin=-30
	latmax=30
	lonmin=-180
	lonmax=180
	outfile_suffix = "_30S-30N_global_2005"
		
#Warm Pool
flag=0
if (flag):
	filelist = [path+'modis_clouds.2005*.nc']
	latmin=-20
	latmax=20
	lonmin=110
	lonmax=180
	outfile_suffix = "_WP_ocean_2005"

#Amazon
flag=0
if (flag):
	latmin=-10
	latmax=0
	lonmin=-65
	lonmax=-50

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
#					,"fraction_colder_than_220K"
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
#					,"min_IR"
#					,"time"
					,"particle_effective_radius"
					,"PCT_at_min_IR"
#					,"PCT"
					,"shear_at_min_PCT"
					,"shear_at_min_IR"
					,"shear"
					,"wind_direction_200hPa_at_min_IR"
					,"wind_direction_500hPa_at_min_IR"
					,"wind_direction_850hPa_at_min_IR"
					]

# Must screen out -9999s. Can use mask_and_scale=True in open_dataset here with the 
# final product files if the fill_value is defined in the header. For now will filter out
# with where below.
if ( len(filelist)==1 and bool(re.search('\*',filelist[0]))==False ):
	clouds_all = xr.open_dataset(filelist[0], drop_variables=dropvars)
elif (len(filelist)==1 and bool(re.search('\*',filelist[0]))==True ):
	clouds_all = xr.open_mfdataset(filelist[0],
							drop_variables=dropvars,concat_dim="number_of_clouds",combine="nested")
else:
	clouds_all = xr.open_mfdataset(filelist,
							drop_variables=dropvars,concat_dim="number_of_clouds",combine="nested")
print("# of clouds_all = {0}".format(clouds_all["cloud_areas"].values.shape[0]))

clouds = clouds_all[["latitude_of_min_IR"
										,"longitude_of_min_IR"
										,"max_view_zenith"
										,"min_IR","cloud_areas"
										,"time"
										,"fraction_colder_than_220K"
										,"fraction_over_land"]]
print("# of IR clouds = {0}".format(clouds["cloud_areas"].values.shape[0]))

clouds = clouds.where(
										(clouds_all["latitude_of_min_IR"] >= latmin)
										& (clouds_all["latitude_of_min_IR"] <= latmax)
										& (clouds_all["longitude_of_min_IR"] >= lonmin)
										& (clouds_all["longitude_of_min_IR"] <= lonmax)
										& (clouds_all["max_view_zenith"]<65)
										).dropna(dim="number_of_clouds")
print("# of IR clouds after QA step = {0}".format(clouds["cloud_areas"].values.shape[0]))

clouds_df = clouds.to_dataframe()

hour_offset = np.floor( (7.5 + clouds_df["longitude_of_min_IR"])/15. )
local_time = clouds_df["time"] + pd.to_timedelta(hour_offset, unit="hours")
clouds_df["local time"] = local_time


clouds_df["log10_cloud_area"] = np.log10(clouds_df["cloud_areas"])
logareabins = np.arange(nareabins)*( np.log10(areamax) - np.log10(areamin) )/(nareabins-1) + np.log10(areamin)
clouds_df["day-night"] = pd.cut(pd.DatetimeIndex(clouds_df["local time"]).hour, bins=[0,6,18,24],
						include_lowest=True, labels=['night','day','night'], ordered=False)
clouds_df["land-ocean"] = pd.cut(clouds_df.fraction_over_land, bins=[0,0.5,1],
						include_lowest=True, labels=['ocean','land'])
clouds_df["areabins"] = pd.cut(clouds_df['log10_cloud_area'],
						bins=logareabins, include_lowest=True)

dlogbin = logareabins[2]-logareabins[1]
areaout = np.power( 10, logareabins[0:nareabins-1] + (dlogbin/2.) )
		
f1 = plt.figure()
s1 = f1.add_subplot(221)

minIR_ocean_day = clouds_df[["min_IR","areabins","land-ocean","day-night"]].where((clouds_df["land-ocean"]=="ocean") & (clouds_df["day-night"]=="day"))
s1.plot(areaout,minIR_ocean_day.groupby("areabins").mean(),'b-')
minIR_ocean_day.boxplot(column="min_IR",by="areabins",ax=s1,positions=areaout,whis=0,sym='',widths=areaout/10.,color='blue')
minIR_land_day = clouds_df[["min_IR","areabins","land-ocean","day-night"]].where((clouds_df["land-ocean"]=="land") & (clouds_df["day-night"]=="day"))
s1.plot(areaout,minIR_land_day.groupby("areabins").mean(),'k-')
minIR_land_day.boxplot(column="min_IR",by="areabins",ax=s1,positions=areaout,whis=0,sym='',widths=areaout/10.,color='black')
minIR_ocean_night = clouds_df[["min_IR","areabins","land-ocean","day-night"]].where((clouds_df["land-ocean"]=="ocean") & (clouds_df["day-night"]=="night"))
s1.plot(areaout,minIR_ocean_night.groupby("areabins").mean(),'b--')
minIR_land_night = clouds_df[["min_IR","areabins","land-ocean","day-night"]].where((clouds_df["land-ocean"]=="land") & (clouds_df["day-night"]=="night"))
s1.plot(areaout,minIR_land_night.groupby("areabins").mean(),'k--')




clouds_amsr = clouds_all[[
													"latitude_of_min_IR"
													,"longitude_of_min_IR"
													,"max_view_zenith"
													,"PCT"
													,"cloud_areas"
													,"time"
													,"fraction_over_land"
													]]
print("# of AMSR clouds = {0}".format(clouds_amsr["cloud_areas"].values.shape[0]))
										
clouds_amsr = clouds_amsr.where(
													(clouds_all["latitude_of_min_IR"] >= latmin)
													& (clouds_all["latitude_of_min_IR"] <= latmax)
													& (clouds_all["longitude_of_min_IR"] >= lonmin)
													& (clouds_all["longitude_of_min_IR"] <= lonmax)
#													& (clouds_all["max_view_zenith"]<53)
													).dropna(dim="number_of_clouds")
print("# of AMSR clouds after QA step = {0}".format(clouds_amsr["cloud_areas"].values.shape[0]))
clouds_amsr_df = clouds_amsr.to_dataframe()
hour_offset = np.floor( (7.5 + clouds_amsr_df['longitude_of_min_IR'])/15. )
local_time = clouds_amsr_df["time"] + pd.to_timedelta(hour_offset, unit="hours")
clouds_amsr_df["local time"] = local_time
clouds_amsr_df["log10_cloud_area"] = np.log10(clouds_amsr_df["cloud_areas"])
clouds_amsr_df["day-night"] = pd.cut(pd.DatetimeIndex(clouds_amsr_df["local time"]).hour, bins=[0,6,18,24],
						include_lowest=True, labels=['night','day','night'], ordered=False)
clouds_amsr_df["land-ocean"] = pd.cut(clouds_amsr_df.fraction_over_land, bins=[0,0.5,1],
						include_lowest=True, labels=['ocean','land'])
clouds_amsr_df["areabins"] = pd.cut(clouds_amsr_df['log10_cloud_area'],
						bins=logareabins, include_lowest=True)


s2 = f1.add_subplot(222)

PCT_ocean_day = clouds_amsr_df[["PCT","areabins","land-ocean","day-night"]].where((clouds_amsr_df["land-ocean"]=="ocean") & (clouds_amsr_df["day-night"]=="day"))
s2.plot(areaout,PCT_ocean_day.groupby("areabins").mean(),'b-')
PCT_ocean_day.boxplot(column="PCT",by="areabins",ax=s2,positions=areaout,whis=0,sym='',widths=areaout/10.,color='blue')
PCT_land_day = clouds_amsr_df[["PCT","areabins","land-ocean","day-night"]].where((clouds_amsr_df["land-ocean"]=="land") & (clouds_amsr_df["day-night"]=="day"))
s2.plot(areaout,PCT_land_day.groupby("areabins").mean(),'k-')
PCT_land_day.boxplot(column="PCT",by="areabins",ax=s2,positions=areaout,whis=0,sym='',widths=areaout/10.,color='black')
PCT_ocean_night = clouds_amsr_df[["PCT","areabins","land-ocean","day-night"]].where((clouds_amsr_df["land-ocean"]=="ocean") & (clouds_amsr_df["day-night"]=="night"))
s2.plot(areaout,PCT_ocean_night.groupby("areabins").mean(),'b--')
PCT_land_night = clouds_amsr_df[["PCT","areabins","land-ocean","day-night"]].where((clouds_amsr_df["land-ocean"]=="land") & (clouds_amsr_df["day-night"]=="night"))
s2.plot(areaout,PCT_land_night.groupby("areabins").mean(),'k--')


		
l = 0.1
bt = 0.08
w = 0.35
h = 0.4	
s1.set_ylabel(r'min IR Tb (K)')
s1.set_xlabel(r'cloud area (km^2)')
s1.set_xscale(r'log')
s1.set_xlim([areamin,areamax])
s1.set_ylim([180,260])
s1.set_position([l,bt+0.5,w,h]) #[left, bottom, width, height]
s1.set_title('')
s2.set_ylabel(r'ave. 89 GHz PCT (K)')
s2.set_xlabel(r'cloud area (km^2)')
s2.set_xscale(r'log')
s2.set_xlim([areamin,areamax])
s2.set_ylim([260,300])
s2.set_position([l+0.5,bt+0.5,w,h]) #[left, bottom, width, height]
s2.set_title('')
plt.suptitle('')

print("saving to "+outpath+outfile+outfile_suffix+".png")
plt.savefig(outpath+outfile+outfile_suffix+".png")
print("saving to "+outpath+outfile+outfile_suffix+".ps")
plt.savefig(outpath+outfile+outfile_suffix+".ps")










