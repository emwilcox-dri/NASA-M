#measures_setup.py

from pyhdf import SD
import numpy as np
import subprocess
import time_tools3 as tt



def get_data_paths(yr,month,proj_root):
	monstr = str(month)
	if (month<10):
		monstr = "0" + monstr
	yrstr = str(yr)

	modpath = proj_root + "MODIS_" + yrstr + "/MYD06_L2/" + monstr + "/"
	mod02path = proj_root + "MODIS_" + yrstr + "/MYD021KM/" + monstr + "/"
	mod03path = proj_root + "MODIS_" + yrstr + "/MYD03/" + monstr + "/"
	mod08path = proj_root + "MODIS_" + yrstr + "/MYD08_D3/"
	if ( (yr>=2001) & (yr<=2010) ):
		merrapath = "/discover/nobackup/projects/gmao/merra2/data/products/d5124_m2_jan00/Y" + yrstr + "/M" + monstr + "/"
	if (yr>2010):
		merrapath = "/discover/nobackup/projects/gmao/merra2/data/products/d5124_m2_jan10/Y" + yrstr + "/M" + monstr + "/"
	amsrpath = proj_root + "AMSR_E_L2A/"
	
	return (modpath, mod02path, mod03path, merrapath, amsrpath, mod08path)



########################
# updated 5/24/2021 to access MODIS and MERRA2 data in the ccs data holdings
# mounted on discover
# Required shifting the modis file structure from mm/dd to mm/doy
# MODIS_root: /nfs3m/css/curated01/modis/data/Collection6.1/
# merra_root: /css/merra2/
# MOD06_L2: MODIS_root/L2/MYD06_L2/yyyy/doy/
# MOD03: MODIS_root/L1/MYD03/yyyy/doy/
# MOD02: MODIS_root/L1/MOD021KM/yyyy/doy/
# MOD08: MODIS_root/L3/MYD08_D3/yyyy/doy/
# AMSR-E: unchanged, proj_root/AMSR_E_L2A/yyyymm/
# MERRA2: merra_root/d5124_m2_jan00/Yyyyy/Mmm/ or merra_root/d5124_m2_jan10/Yyyyy/Mmm/
########################
def get_data_paths_css(yr,month,day,proj_root,modis_root,merra_root):
	monstr = str(month)
	if (month<10):
		monstr = "0" + monstr
	yrstr = str(yr)
	doy = tt.JDNumberToDOY(tt.mdyhmsToJDNumber(month,day,yr,6,0,0))
	doystr = str(doy)
	if (doy<100):
		doystr = "0" + doystr
	if (doy<10):
		doystr = "0" + doystr

	modpath = modis_root + "L2" + "/MYD06_L2/" + yrstr + "/" + doystr + "/"
	mod02path = modis_root + "L1" + "/MYD021KM/" + yrstr + "/" + doystr + "/"
	mod03path = modis_root + "L1" + "/MYD03/" + yrstr + "/" + doystr + "/"
	mod08path = modis_root + "L3" + "/MYD08_D3/" + yrstr + "/" + doystr + "/"
	if ( (yr>=2001) & (yr<=2010) ):
		merrapath = merra_root + "d5124_m2_jan00/Y" + yrstr + "/M" + monstr + "/"
	if (yr>2010):
		merrapath = merra_root + "d5124_m2_jan10/Y" + yrstr + "/M" + monstr + "/"
	amsrpath = proj_root + "AMSR_E_L2A/"
	
	return (modpath, mod02path, mod03path, merrapath, amsrpath, mod08path)



def write_filelists(yr, month, listpath_root, modpath, mod02path,
										mod03path, mod08path, amsrpath):
	monstr = str(month)
	if (month<10):
		monstr = "0" + monstr
	yrstr = str(yr)

	cmd = "mkdir " + listpath_root
	tmp = subprocess.run(cmd,shell=True)
	modlistroot = listpath_root + "MODIS_" + yrstr
	cmd = "mkdir " + modlistroot
	tmp = subprocess.run(cmd,shell=True)
	modlistpath = modlistroot + "/" + monstr
	cmd = "mkdir " + modlistpath
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}*.hdf -printf \"%f\\n\" > {1}/myd06list".format(modpath,modlistpath)
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}*.hdf -printf \"%f\\n\" > {1}/myd02list".format(mod02path,modlistpath)
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}*.hdf -printf \"%f\\n\" > {1}/myd03list".format(mod03path,modlistpath)
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}*.hdf -printf \"%f\\n\" > {1}/myd08_d3filelist.txt".format(mod08path,modlistroot)
	tmp = subprocess.run(cmd,shell=True)
	
	amsrlistroot = listpath_root + "AMSR_E_L2A"
	cmd = "mkdir " + amsrlistroot
	tmp = subprocess.run(cmd,shell=True)
	amsrlistpath = amsrlistroot + "/" + yrstr + monstr
	cmd = "mkdir " + amsrlistpath
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}{1}{2}/*.hdf -printf \"%f\\n\" > {3}/AMSRlist".format(amsrpath,yrstr,monstr,amsrlistpath)
	tmp = subprocess.run(cmd,shell=True)



########################
# updated 5/24/2021 to access MODIS and MERRA2 data in the ccs data holdings
# mounted on discover
########################
def write_filelists_css(yr, month, day, listpath_root, modpath, mod02path,
										mod03path, mod08path, amsrpath):
	monstr = str(month)
	if (month<10):
		monstr = "0" + monstr
	yrstr = str(yr)
	doy = tt.JDNumberToDOY(tt.mdyhmsToJDNumber(month,day,yr,6,0,0))
	doystr = str(doy)
	if (doy<100):
		doystr = "0" + doystr
	if (doy<10):
		doystr = "0" + doystr

	cmd = "mkdir " + listpath_root
	tmp = subprocess.run(cmd,shell=True)
	modlistroot = listpath_root + "MODIS_" + yrstr
	cmd = "mkdir " + modlistroot
	tmp = subprocess.run(cmd,shell=True)
	modlistpath = modlistroot + "/" + doystr
	cmd = "mkdir " + modlistpath
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}*.hdf -printf \"%f\\n\" > {1}/myd06list".format(modpath,modlistpath)
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}*.hdf -printf \"%f\\n\" > {1}/myd02list".format(mod02path,modlistpath)
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}*.hdf -printf \"%f\\n\" > {1}/myd03list".format(mod03path,modlistpath)
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}*.hdf -printf \"%f\\n\" > {1}/myd08_d3filelist.txt".format(mod08path,modlistpath)
	tmp = subprocess.run(cmd,shell=True)
	
	amsrlistroot = listpath_root + "AMSR_E_L2A"
	cmd = "mkdir " + amsrlistroot
	tmp = subprocess.run(cmd,shell=True)
	amsrlistpath = amsrlistroot + "/" + yrstr + monstr
	cmd = "mkdir " + amsrlistpath
	tmp = subprocess.run(cmd,shell=True)
	cmd = "find {0}{1}{2}/*.hdf -printf \"%f\\n\" > {3}/AMSRlist".format(amsrpath,yrstr,monstr,amsrlistpath)
	tmp = subprocess.run(cmd,shell=True)
	
	
	

def get_grid_bounds(modpath,mod06file):
	sdid = SD.SD(modpath+mod06file)
	sdsid = sdid.select("Longitude")
	lon = sdsid.get()
	sdsid = sdid.select("Latitude")
	lat = sdsid.get()
	sdid.end()

	latmin = np.floor(np.min(lat[np.where(lat>=-90)]))
	latmax = np.ceil(np.max(lat))
	lonmin = np.floor(np.min(lon[np.where(lon>=-180)]))
	lonmax = np.ceil(np.max(lon))
	
	nx = len(lon[0,:])
	na = len(lon[:,0])
	
	if ( (lonmin<-179) & (lonmax>179) ):
		for i in range(nx):
			for j in range(na):
				if (lon[j,i]<0):
					lon[j,i] = 360.+lon[j,i]
				
		lonmin = np.floor(np.min(lon[np.where(lon>=-180)]))
		lonmax = np.ceil(np.max(lon))
	
	return (latmin,latmax,lonmin,lonmax)




def get_granule_midpoint(modpath,mod06file):
	sdid = SD.SD(modpath+mod06file)
	sdsid = sdid.select("Longitude")
	lon = sdsid.get()
	sdsid = sdid.select("Latitude")
	lat = sdsid.get()
	sdid.end()

	ix = int(len(lon[0,:])/2)
	ia = int(len(lon[:,0])/2)
	
	lon_mid = lon[ia,ix]
	lat_mid = lat[ia,ix]
	
	return (lon_mid,lat_mid)








