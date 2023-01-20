/*****************
NOTES:
v12. updated the CAPE calculation to take the mixing ratio straight from the MERRA2
	files after converting from kg/kg to g/kg. Also merra2.c (and .h) updated to use
	the correct value of dlon/dlat, and also to use correct midpoints on grid in get_merra_indices_4D.
	Finally, also updated to use only pressure levels with lower pressure than surface
	pressure. Having bad values in the profile (for pressure levels below the surface
	pressure level) was causing the CAPE routine to return zero for those profiles.
	

v11. partial change to output the cloud mask with indices for each cloud in a netcdf file
	was finished in v12.

v10. updated for locating data archived on CSS directories on Discover.
	MODIS filelists now organized by YYYY/DOY, rather than YYYY/MM

v9. updated for MERRA2 aerosols
merra2 ==> /css/merra2
modis ==> /nfs3m/css/curated01/modis/data/Collection6.1
	MYD03 in L1/
	MYD08 in L2/
	MYD08_D3 in L3/

v6. Previous version was leading to missing CAPE/shear values for many clouds. I think
the problem was that it was looping over the MERRA grid, but that the gridded satellite
data are of finer resolution, so many cloud pixels were not being assigned CAPE/shear
values. This version is updated now to compute CAPE/shear at each cloudy grid cell on the
common grid. Also now outputting the CAPE/shear at the minimum IR brightness temperature
location, in addition to the max and the cloud-averaged value.

v5. Updated with date/times for each cloud in the output and removing the 5km dimension
and variables to aid in automated concatenation of the netcdf output files. Also changed
to output the IR brightness temp bins in their own file, since they are the same
for all output in a single run. Also changing the number of clouds to an unlimited
dimension since it is unlimited.

v4. Updated to include: (1) the bounds of the grid as passed arguments;
(2) revised code to use pressure level MERRA2 data (as opposed to z-level MERRA); and
(3) bug fixes to the effective radius, IR temperature profiles.

No AMSR-E 85GHz feedhorn A data after Sept. 2004

The nomenclature regarding Tb31 and Tb32 are now mixed up quite badly here.
Tb32 used to be for channel 32 brightness temperature, now it mostly refers to the downstream
results from the 1km channel 31 brightness temperature, while Tb31 mostly refers to the downstream
results from the 5km channel 31 brightness temperature. This should be fixed.

longitude characteristics:
MOD06_L2 ==> -180 to 180
MOD08_L3 ==> 0 to 360 (converted to -180 to 180 to match with the others)
AMSR-E ==> -180 to 180
MERRA2 ==> -180 to 180

v4 now handles the special case for granules that span the international dateline

v6:
Corrects the gaps in the CAPE/shear values

v8:
added the AOD values
*****************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include <hdf.h>
#include <mfhdf.h>
#include "misc_tools.h"
#include "MOD06.h"
#include "MOD08.h"
#include "MOD02.h"
#include "amsr.h"
#include "grid.h"
#include "time_tools.h"
#include "das.h"
#include "thermo.h"
#include "clouds.h"
#include "measures.h"
/******
#include "misc_tools.h"
#include "MOD06.h"
#include "MOD08.h"
#include "amsr.h"
#include "grid.h"
#include "das.h"
#include "thermo.h"
#include "clouds.h"
*****/
#include "mask.h"
//#include "mask.h"
#ifdef MERRA_FLAG
#include "merra2.h"
#endif




/************** Begin user defined constants **********/

//Define the width of lat/lon grid cells
#define GRID_DLON 0.125
#define GRID_DLAT 0.125
#ifdef COARSE_GRID_FLAG
//These are for statistics on a coarse grid
#define CGRID_DLON 5
#define CGRID_DLAT 5
#endif

//#define N_LEVELS_CAPE 36 //When using Nv files with 72 levels
#define N_LEVELS_CAPE 25 //When using Np files with 42 levels

//central US
/****
#define LATMIN 30
#define LATMAX 40
#define LONMIN -103
#define LONMAX -85
****/

/****
#ifdef SUMMER_MONSOON
#define LATMIN 0
#define LATMAX 35 
#define LONMIN 65
#define LONMAX 95
#endif

#ifdef WINTER_MONSOON
#define LATMIN -10
#define LATMAX 15 
#define LONMIN 50
#define LONMAX 100
#endif

#define LATMIN -5
#define LATMAX 20 
#define LONMIN 120
#define LONMAX 145
****/

//N. Africa
/****
#define LATMIN 4
#define LATMAX 12
#define LONMIN -20
#define LONMAX 40
****/

#define IR_BTMIN 180.
#define IR_BTMAX 270.
#define NIR_BTBINS 45
#define MAX_MODIS_ZENITH_ANGLE 58. //note that presently this is not being applied to 1km IR brightness temperatures

const char amsr_monthly_dir_flag = 1; //1=AMSR-E tbs divided into monthly directories
//const char *amsrPath_root = "/Data/AMSR-E/tbs/";
//const char *mod08PathRoot = "/Volumes/shores/MODIS/MYD08_D3/col_51/";
const char *mod08filelist = "myd08_d3filelist.txt";
#ifdef DAO_FLAG
const char *daoPath = "/Volumes/shores/NCEP/dao_ops/";
const char *dao_fileroot = "tsyn3d_mis_p_ncep";
#endif

//const char *maskFile =
//  "/Users/ewilcox/misc_tools/land_water_masks/land_ocean_masks_xdeg/land_ocean_mask2_qd.asc";
const char *maskFilename = "land_ocean_mask2_qd.asc";

#ifdef MERRA_FLAG
//const char *merra_fileroot = "MERRA300.prod.assim.inst3_3d_asm_Cp.";
//const char *merra_fileroot = "MERRA300.prod.assim.inst6_3d_ana_Nv.";
//const char *merra_fileroot = "MERRA301.prod.assim.inst6_3d_ana_Nv.";
//const char *merra_fileroot = "MERRA2_300.inst6_3d_ana_Np.";
const char *merra_fileroot = "inst6_3d_ana_Np.";
const char *merra_aer_fileroot = "inst3_2d_gas_Nx.";
#endif

/************** End of user defined constants *********/





const char *aqua = "MYD";
const char *terra = "MOD";
const char *outfile_str = "modis_clouds";
const char *MOD02list = "mod02list";
const char *MOD03list = "mod03list";
const char *MOD04list = "mod04list";
const char *MYD02list = "myd02list";
const char *MYD03list = "myd03list";
const char *MYD04list = "myd04list";

int main(int argc, char *argv[])
{
  FILE *fpo, *fp;
  char MOD06file[150], modPath[150], amsrfile[150], mod08File[50];
  char MOD03file[150], MOD02file[150], MOD04filetmp[150];
  char amsrPath[150], outfile[150], tmpfile[60], dir_list[100];
  char dao_opsfile[150], outPath[150], aodPath[150], maskpath[150];
  char maskfile[200], mod02Path[150], mod03Path[150], sattest[4];
  char yyyyddd[8], ddd[4], hhmm[5], hh[3], mm[3], mn[3];
  char dd[3], yyyy[5], mmdd[5], ccc[4], *s;
  char gridfile[150], *noneStr = "none", cmd[100], amsrPath_root[50];
  char merraPath[150], ismodfile, verstr[3], mod08PathRoot[50];
  char listpath_root[150], amsrlistPath[150], modlistPath[150];
  char mod08listPath[150], gstr[10];
  double jd, jd0;
  long min_cent, *cloud_datetime=NULL;
  char continue_flag=1, dateline_flag=0;
  unsigned char *land_mask=NULL;
  signed char *cldtop_phase=NULL;
  short *land_ocean_qd=NULL, amsr_fileflag=1;
  short *class31=NULL, *class32=NULL;
  short *cloud31_class=NULL, *cloud32_class=NULL;
  int mon, day, minute, hour, year, doy, gran, status, inttmp;
  int dao_lat_offset, dao_lon_offset, tstep, it, ilo, ila, ilev;
  int Tb31_nhistbins, pct_nhistbins, amsr_vers, mod_collection;
  int ncid, tb31dimid, tb32dimid, area31id, area32id, class31id, class32id;
  int landfracid, maxmodzenid, minIRtbid, aveCAPEid, daveCAPEid, maxCAPEid;
  int maxdCAPEid, aveshearid, maxshearid, pctaveid, minPCTid, frac220id;
  int frac210id, frac200id, latminTBid, lonminTBid, latminPCTid, lonminPCTid;
  int IRminPCTid, PCTminIRid, IRbinid, derbinid, irbindimid, derbindims[2];
  int capeminTBid, shearminTBid, capeminpctid, shearminpctid, w850minTBid;
  int w500minTBid, w200minTBid, AODminIRid,AODmerraid,AODmerraminIRid;
  int iband, nmod02bands, timeid, is_merrafile, is_merra_aerfile;
  unsigned int *Tb31_hist=NULL, *pct_hist=NULL;
  unsigned long *Tb31_histcnt=NULL, *pct_histcnt=NULL;
  float Tb31_binmin, Tb31_binmax, Tb31_binsize, pct_binmin;
  float pct_binmax, pct_binsize;
  float *Tb31_histbins=NULL, *pct_histbins=NULL;
  long ngridlons, ngridlats, na_5km, na_1km, ncloud31=0, ncloud32=0, ncloud;
  long i, j, ia, ix, ib, ic, iloc, kloc, jloc, baseloc, toploc, yyyymmdd;
  long *cloud31=NULL, *cloud32=NULL, *cloud=NULL, *ncloud_pix;
  long na_amsr, nx_amsr, *cloud_pix=NULL, *cloudpix_5km=NULL;
  float *land_ocean=NULL, *modflttmp=NULL, *aod_cg=NULL, *aod_mod=NULL;
  float *lon_5km=NULL, *lon_1km=NULL, *Tbs=NULL, *Tb_31=NULL, *Tb_32=NULL;
  float32 *lat_5km=NULL, *lat_1km=NULL;
  float *grid_lons=NULL, *grid_lats=NULL, *Tb31_g=NULL, *Tb32_g=NULL;
  float *tmp_tb=NULL, *mod08lons=NULL, *mod08lats=NULL;
  float *areas=NULL, dlon, dlat, *tmpstats=NULL;
  float *lon_amsr_a=NULL, *lat_amsr_a=NULL, *lon_amsr_b=NULL;
  float *lat_amsr_b=NULL, *tb89_Ha=NULL, *tb89_Va=NULL;
  float *tb89_Hb=NULL, *tb89_Vb=NULL, *Tb89V_g=NULL, *Tb89H_g=NULL;
  float *cloud31_areas=NULL, *cloud32_areas=NULL, *IRrad=NULL, *IRrad_uncert=NULL;
  float *Tb89V_ave=NULL, *Tb89H_ave=NULL, *minTb89V=NULL, *minTb89H=NULL;
  float *Tb31_ave=NULL, *Tb32_ave=NULL, *minTb31=NULL, *minTb32=NULL;
  float *minpct=NULL, *lat_minTb32=NULL, *lon_minTb32=NULL;
  float *lat_minpct=NULL, *lon_minpct=NULL, *Tb31_1km=NULL;
  float *mask_lon=NULL,*mask_lat=NULL;
  float *lmask_a=NULL,*lmask_b=NULL, *landfrac=NULL, *pct=NULL, *pct_ave=NULL;
  float *cape_ave=NULL, *dcape_g=NULL, *dcape_ave, *max_cape=NULL;
  float *max_dcape=NULL, *Tb31_220_flag=NULL, *Tb31_200_flag=NULL;
  float *Tb31_220frac=NULL, *Tb31_200frac=NULL;
  float *Tb31_210frac=NULL, *Tb31_210_flag=NULL;
  float *lon_ncep=NULL, *lat_ncep=NULL, *lev_ncep=NULL, *uwnd=NULL, *vwnd=NULL;
  float *temp_ncep=NULL, *rh_ncep=NULL, *theta_e=NULL, *theta_e_g=NULL;
  float *lons_dao=NULL, *lats_dao=NULL, *levs_dao=NULL;
  float *temp_dao=NULL, *rh_dao=NULL, *z_dao=NULL, *cape_g=NULL, cape;
  float *T_prof=NULL, *q_prof=NULL, *z_prof=NULL, Plcl, *RH_prof=NULL;
  float *P_prof=NULL, *Tv_prof=NULL, *delP_prof=NULL, *P_prof_tmp=NULL;
  float *RH_prof_tmp=NULL, *T_prof_tmp=NULL,*z_prof_tmp=NULL;
  float *reff=NULL, *tmpvar=NULL, *tmpvar1km=NULL;
  double delt_amsr;
  float *gam_prof=NULL, *gamma=NULL, *Tvp=NULL, *Tve=NULL, *Tvp_prof=NULL;
  float *Tve_prof=NULL;
  float *cape_cg=NULL, *dcape_cg=NULL, *shear_g=NULL, *shear_cg=NULL;
  float *Tb240Kfrac_cg=NULL, *crf_lw_cg=NULL, *crf_sw_cg=NULL;
  float *tmp_flag=NULL,dlon_cgrid, dlat_cgrid, *cgrid_lons=NULL;
  float *cgrid_lats=NULL, *latstmp=NULL, *lonstmp=NULL, *areas_cgrid=NULL;
  float *shear_ave=NULL, *max_shear=NULL, *max_modzen=NULL;
  float *mod_zenith=NULL, *mod_azimuth=NULL;
  float *modzen_g=NULL, *der_bin=NULL, *bt_bin=NULL, *der=NULL, *bt=NULL;
  float *Tb31_minpct=NULL, *pct_minTb32=NULL;
  float *cape_minTb32=NULL, *shear_minTb32=NULL, *cape_minpct=NULL, *shear_minpct=NULL;
  float *wdir850_minTb32=NULL, *wdir500_minTb32=NULL, *wdir200_minTb32=NULL;
  float *aod_minTb32=NULL, *aod_merra_minTb32=NULL, *aod_merra_ave=NULL;
  long ncgridlons, ncgridlats;
#ifdef MERRA_FLAG
	char merra_file[50], merra_aer_file[50];
	double *lon_merra=NULL, *lat_merra=NULL, *P_merra=NULL;
	int *t_merra=NULL;
	float *temp_merra=NULL, *q_merra=NULL, *uwnd_merra=NULL, *vwnd_merra=NULL;
	float *delP_merra=NULL, *P0_merra=NULL, merra_fill=0, *aod_merra=NULL;
	float *levs_merra=NULL, ptop, flttmp, *cape_merra=NULL, *shear_merra=NULL;
	float *Ttest=NULL, *RHtest=NULL, *plcltest=NULL;
	long testloc, capecnt=0, ilon, ilat;
#endif //MERRA_FLAG
  float lat_cent, lon_cent, minlo, maxlo, minla, maxla;
  float viewed_area, max_modzenith_thresh;
  float latmin, lonmin, latmax, lonmax, lonmintmp, lontest;
  float *wdir850=NULL, *wdir500=NULL, *wdir200=NULL, *aod_merra_g=NULL;
  size_t sp, cp, cp31, cp32, cpbin, spder[2], cpder[2];
  int finegriddims[2], cld31id, cld32id, londimid, latdimid, longridid;
  int latgridid, capgid, flag, profilev, irgid;

//printf("in the program\n");
//fflush(NULL);

  strcpy(MOD06file,argv[1]);
  strcpy(modPath,argv[2]);
  strcpy(mod02Path,argv[3]);
  strcpy(mod03Path,argv[4]);
  strcpy(merraPath,argv[5]);
  strcpy(amsrPath_root,argv[6]);
  strcpy(mod08PathRoot,argv[7]);
  strcpy(maskpath,argv[8]);
  strcpy(listpath_root,argv[9]);
  strcpy(outPath,argv[10]);
  strcpy(gstr,argv[11]);
  latmin = atof(gstr);
  strcpy(gstr,argv[12]);
  latmax = atof(gstr);
  strcpy(gstr,argv[13]);
  lonmin = atof(gstr);
  strcpy(gstr,argv[14]);
  lonmax = atof(gstr);

//printf("%f %f %f %f\n",latmin,latmax,lonmin,lonmax);
//fflush(NULL);
  
  if (lonmax>180) dateline_flag=1;

  strncpy(yyyy,MOD06file+10,4);
  yyyy[4] = '\0';
  year = atoi(yyyy);
  strncpy(ddd,MOD06file+14,3);
  ddd[3] = '\0';
  strncpy(yyyyddd,MOD06file+10,7);
  yyyyddd[7] = '\0';
  doy = atoi(ddd);
  doy_to_mmdd(doy,year,mmdd);
  strncpy(mn,mmdd,2);
  mn[2] = '\0';
  mon = atoi(mn);
  strncpy(dd,mmdd+2,2);
  dd[2] = '\0';
  day = atoi(dd);
  strncpy(hhmm,MOD06file+18,4);
  hhmm[4] = '\0';
  strncpy(hh,MOD06file+18,2);
  hh[2] = '\0';
  hour = atoi(hh);
  strncpy(mm,hhmm+2,2);
  mm[2] = '\0';
  minute = atoi(mm);
  yyyymmdd = year*10000+mon*100+day;
  strncpy(ccc,MOD06file+23,3);
  ccc[3] = '\0';
  mod_collection = atoi(ccc);
  
  max_modzenith_thresh = MAX_MODIS_ZENITH_ANGLE;

/***** specify the grid *****/
  dlon = GRID_DLON;
  dlat = GRID_DLAT;
  get_grid_specs(lonmin,latmin,lonmax,latmax,dlon,dlat,&ngridlons,&ngridlats);
  grid_lons = (float *) malloc(ngridlons*sizeof(float));
  grid_lats = (float *) malloc(ngridlats*sizeof(float));
  areas = (float *) malloc(ngridlats*sizeof(float));
  get_grid(lonmin,latmin,dlon,dlat,ngridlons,ngridlats,
           grid_lons,grid_lats,areas);
  sprintf(gridfile,"%sgrid.bin",outPath);
#ifdef VERBOSE
  fp = fopen(gridfile,"w");
  if (fp != NULL) {
    fwrite(&ngridlons,sizeof(long),1,fp);
    fwrite(&ngridlats,sizeof(long),1,fp);
    fwrite(grid_lons,sizeof(float),ngridlons,fp);
    fwrite(grid_lats,sizeof(float),ngridlats,fp);
    fclose(fp);
    fp = NULL;
  }
#endif
  printf("lonmin=%f lonmax=%f\n",lonmin,lonmax);
  printf("grid created: dlon=%f, dlat=%f, ngridlons=%ld, ngridlats=%ld\n",dlon,dlat,ngridlons,ngridlats);



/********* Get the MOD03 and MOD02 filenames ***********/
	strncpy(sattest,MOD06file,3);
	sattest[3] = '\0';
	sprintf(modlistPath,"%sMODIS_%s/%s/",listpath_root,yyyy,ddd);
  if ( (inttmp = strcmp(sattest,aqua)) == 0 )
		getFileNames(modlistPath,modlistPath,yyyyddd,hhmm,MYD03list,MYD04list,
										MYD02list,MOD03file,MOD04filetmp,MOD02file);
  else
    getFileNames(modlistPath,modlistPath,yyyyddd,hhmm,MOD03list,MOD04list,
                    MOD02list,MOD03file,MOD04filetmp,MOD02file);
  if (strcmp(MOD03file,noneStr)==0) {
    printf("no MOD03 file.\n");
    exit(0);
  }
  if (strcmp(MOD02file,noneStr)==0) {
    printf("no MOD02 L1B file.\n");
  }
//  printf("%s\n",MOD03file);
//  printf("%s\n",MOD02file);



#ifdef COARSE_GRID_FLAG
/******* specify the coarse grid ******/
  dlon_cgrid = CGRID_DLON;
  dlat_cgrid = CGRID_DLAT;
  get_grid_specs(lonmin,latmin,lonmax,latmax,dlon_cgrid,dlat_cgrid,
                 &ncgridlons,&ncgridlats);
  cgrid_lons = (float *) malloc(ncgridlons*sizeof(float));
  cgrid_lats = (float *) malloc(ncgridlats*sizeof(float));
  areas_cgrid = (float *) malloc(ncgridlats*sizeof(float));
  get_grid(lonmin,latmin,dlon_cgrid,dlat_cgrid,ncgridlons,ncgridlats,
           cgrid_lons,cgrid_lats,areas_cgrid);
  sprintf(gridfile,"%scoarse_grid.bin",outPath);
#ifdef VERBOSE
  fp = fopen(gridfile,"w");
  if (fp != NULL) {
    fwrite(&ncgridlons,sizeof(long),1,fp);
    fwrite(&ncgridlats,sizeof(long),1,fp);
    fwrite(cgrid_lons,sizeof(float),ncgridlons,fp);
    fwrite(cgrid_lats,sizeof(float),ncgridlats,fp);
    fclose(fp);
    fp = NULL;
  }
#endif
#endif



/********* allocate gridded fields *******/
  Tb31_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  Tb32_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  cloud31 = (long *) calloc(ngridlons*ngridlats,sizeof(long));
  cloud32 = (long *) calloc(ngridlons*ngridlats,sizeof(long));
  cloud = (long *) calloc(ngridlons*ngridlats,sizeof(long));
  class31 = (short *) calloc(ngridlons*ngridlats,sizeof(short));
  class32 = (short *) calloc(ngridlons*ngridlats,sizeof(short));
  Tb89V_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  Tb89H_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  theta_e_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  pct = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  modzen_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  land_ocean = (float *) malloc(ngridlons*ngridlats*sizeof(float));
/*  get_grid(lonmin,latmin,dlon,dlat,ngridlons,ngridlats,grid_lons,
           grid_lats,areas);*/
#ifdef DAO_OPS
#endif //end if DAO_OPS
#ifdef COARSE_GRID_FLAG
  cape_cg = (float *) malloc(ncgridlons*ncgridlats*sizeof(float));
  dcape_cg = (float *) malloc(ncgridlons*ncgridlats*sizeof(float));
  shear_cg = (float *) malloc(ncgridlons*ncgridlats*sizeof(float));
  Tb240Kfrac_cg = (float *) malloc(ncgridlons*ncgridlats*sizeof(float));
  crf_lw_cg = (float *) malloc(ncgridlons*ncgridlats*sizeof(float));
  crf_sw_cg = (float *) malloc(ncgridlons*ncgridlats*sizeof(float));
  aod_cg = (float *) malloc(ncgridlons*ncgridlats*sizeof(float));
#endif
#ifdef MERRA_FLAG
  cape_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  dcape_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  shear_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  wdir850 = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  wdir500 = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  wdir200 = (float *) malloc(ngridlons*ngridlats*sizeof(float));
	aod_merra_g = (float *) malloc(ngridlons*ngridlats*sizeof(float));
#endif



/********* get the land-ocean mask ********/
/**** longitudes arranged -180 to 180 ******/
  land_ocean_qd = (short *) malloc(180*360*4*4*sizeof(short));
  mask_lon = (float *) malloc(1440*sizeof(float));
  mask_lat = (float *) malloc(720*sizeof(float));
printf("reading the land-ocean mask ... ");
  sprintf(maskfile,"%s%s",maskpath,maskFilename);
  read_land_ocean_qd(maskfile,land_ocean_qd,mask_lon,mask_lat);
//printf("done");
//fflush(NULL);
  for (i=0; i<ngridlats; i++)
  for (j=0; j<ngridlons; j++) {
    iloc = (i*ngridlons)+j;
    if ( dateline_flag && (grid_lons[j]>180) ) ilo = (int) floor((grid_lons[j]-180)*4);
    else ilo = (int) floor((grid_lons[j]-(-180))*4);
    ila = (int) floor((90-grid_lats[i])*4);
    kloc = (ila*1440)+ilo;
    land_ocean[iloc] = land_ocean_qd[kloc];
  }
  FREE(land_ocean_qd);
  FREE(mask_lon);
  FREE(mask_lat);
printf("done\n");



/******** get the MODIS clouds ***********/
  lat_5km = (float32 *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(float32));
  lon_5km = (float *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(float));
  lat_1km = (float32 *) malloc(N_ACROSS_1KM*N_ALONG_1KM_MAX*sizeof(float32));
  lon_1km = (float *) malloc(N_ACROSS_1KM*N_ALONG_1KM_MAX*sizeof(float));
  Tbs = (float *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*N_IR_BANDS*sizeof(float));
  Tb_31 = (float *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(float));
  Tb_32 = (float *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(float));
  IRrad = (float *) malloc(N_ACROSS_L1B_1KM_MAX*N_ALONG_L1B_1KM_MAX*N_L1B_1KM_BANDS*sizeof(float));
  IRrad_uncert = (float *) malloc(N_ACROSS_L1B_1KM_MAX*N_ALONG_L1B_1KM_MAX*N_L1B_1KM_BANDS*sizeof(float));
  Tb31_1km = (float *) malloc(N_ACROSS_L1B_1KM_MAX*N_ALONG_L1B_1KM_MAX*sizeof(float));
  land_mask = (unsigned char *) malloc(N_ACROSS_1KM*N_ALONG_1KM_MAX*sizeof(unsigned char));
  mod_zenith =
       (float *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(float));
  mod_azimuth =
       (float *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(float));
  tmp_flag = (float *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(float));
	cldtop_phase = (signed char *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*
							 sizeof(signed char));
  reff = (float *) malloc(N_ACROSS_1KM*N_ALONG_1KM_MAX*sizeof(float));
  tmpvar1km = (float *) malloc(N_ACROSS_1KM*N_ALONG_1KM_MAX*sizeof(float));
  tmpvar = (float *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(float));
  readMOD06_IR_tbs(modPath,MOD06file,&na_5km,lat_5km,lon_5km,Tbs);
  readMOD06_sensor_angles(modPath,MOD06file,mod_zenith,
                          mod_azimuth);
  readMOD06(modPath,mod03Path,MOD06file,MOD03file,mod_collection,&na_1km,
            &na_5km,lat_5km,lon_5km,lat_1km,lon_1km,
            tmpvar,reff,tmpvar1km,tmpvar,tmpvar,cldtop_phase,tmpvar1km,
            land_mask);
  FREE(mod_azimuth);
  FREE(tmpvar1km);
  FREE(tmpvar);
	nmod02bands = readMOD021KM_IR(mod02Path,mod03Path,MOD02file,MOD03file,&na_1km,
               lat_1km,lon_1km,IRrad,IRrad_uncert,land_mask);
  printf("number of bands in MOD021KM file = %d\n",nmod02bands);
  if (nmod02bands==16) iband = 10;
  else if (nmod02bands==1) iband = 0;
  else printf("ERROR: number of bands in MOD021KM file not 1 or 16!!\n");
//using band index 10 for channel 31 in the IRrad array (need to check this)
//updated now since MOD021KM files on Discover have only band 31 data in them.
/****************
Nans correspond to small negative values of the radiance and values of the raw emiss
variable of 2000-3000 (valid range 0-32767). Maybe need to screen further based on the
EV_1KM_Emissive_Uncert_Indexes
*****************/

	if (dateline_flag){
		for (ix=0; ix<N_ACROSS_5KM; ix++)
		for (ia=0; ia<na_5km; ia++) {
    	iloc = (ia*N_ACROSS_5KM)+ix;
			if (lon_5km[iloc]<0) lon_5km[iloc] = 360.+lon_5km[iloc];
		}
		for (ix=0; ix<N_ACROSS_1KM; ix++)
		for (ia=0; ia<na_1km; ia++) {
    	iloc = (ia*N_ACROSS_1KM)+ix;
			if (lon_1km[iloc]<0) lon_1km[iloc] = 360.+lon_1km[iloc];
		}
	}

  for (ix=0; ix<N_ACROSS_L1B_1KM_MAX; ix++)
  for (ia=0; ia<na_1km; ia++) {
  	kloc = (ia*N_ACROSS_L1B_1KM_MAX) + ix;
//  	iloc = (ia*N_ACROSS_L1B_1KM_MAX*N_IR_BANDS) + (ix*N_IR_BANDS) + 10;
  	iloc = (iband*N_ACROSS_L1B_1KM_MAX*na_1km) + (ia*N_ACROSS_L1B_1KM_MAX) + ix;
  	if (IRrad[iloc] != MOD02_BADFLAG) Tb31_1km[kloc] = modis_bright(IRrad[iloc],31);
  	else Tb31_1km[kloc] = MOD02_BADFLAG;
//  	if (isnan(Tb31_1km[kloc])==1) printf("[%f %f]  ",IRrad[iloc],Tb31_1km[kloc]);
//  	printf("%f ",Tb31_1km[kloc]);
  }
  
  for (ix=0; ix<N_ACROSS_5KM; ix++)
  for (ia=0; ia<na_5km; ia++) {
    kloc = (ia*N_ACROSS_5KM)+ix;
    if (mod_zenith[kloc]<=max_modzenith_thresh) {
      ib=1;
      iloc = (ib*N_ACROSS_5KM*na_5km) + (ia*N_ACROSS_5KM)+ix;
      Tb_31[kloc] = Tbs[iloc];
      ib=2;
      iloc = (ib*N_ACROSS_5KM*na_5km) + (ia*N_ACROSS_5KM)+ix;
      Tb_32[kloc] = Tbs[iloc];
    } else {
      Tb_31[kloc] = MOD06_BADFLAG;
      Tb_32[kloc] = MOD06_BADFLAG;
    }
  }
  
/****** temporary test Tbs *********
//  lat_cent = latmin + (dlat*ngridlats/2);
//  lon_cent = lonmin + (dlon*ngridlons/2);
	minlo = min(lon_5km,N_ACROSS_5KM*na_5km);
	maxlo = max(lon_5km,N_ACROSS_5KM*na_5km);
	minla = min(lat_5km,N_ACROSS_5KM*na_5km);
	maxla = max(lat_5km,N_ACROSS_5KM*na_5km);
  
	lat_cent = ((maxla-minla)/2.) + minla;
	lon_cent = ((maxlo-minlo)/2.) + minlo;
//  printf("%f %f %f %f lat_cent=%f  lon_cent=%f\n",minlo,maxlo,minla,maxla,lat_cent,lon_cent);
  create_test_image(Tb_31,lon_5km,lat_5km,na_5km,lon_cent,lat_cent);
***********************************/
  
  bin_field(Tb_31,lon_5km,lat_5km,MOD06_BADFLAG,Tb31_g,'f',
            N_ACROSS_5KM*na_5km,lonmin,latmin,ngridlons,
            ngridlats,dlon,dlat,grid_lons,grid_lats);
/*  bin_field(Tb_32,lon_5km,lat_5km,MOD06_BADFLAG,Tb32_g,'f',
            N_ACROSS_5KM*na_5km,lonmin,latmin,ngridlons,
            ngridlats,dlon,dlat,grid_lons,grid_lats);*/
  bin_field(Tb31_1km,lon_1km,lat_1km,MOD02_BADFLAG,Tb32_g,'f',
            N_ACROSS_L1B_1KM_MAX*na_1km,lonmin,latmin,ngridlons,
            ngridlats,dlon,dlat,grid_lons,grid_lats);
  bin_field(mod_zenith,lon_5km,lat_5km,MOD06_BADFLAG,modzen_g,'f',
            N_ACROSS_5KM*na_5km,lonmin,latmin,ngridlons,
            ngridlats,dlon,dlat,grid_lons,grid_lats);
            
  for (ix=0; ix<N_ACROSS_5KM; ix++)
  for (ia=0; ia<na_5km; ia++) {
    iloc = (ia*N_ACROSS_5KM)+ix;
    if (Tb_31[iloc] < 240) tmp_flag[iloc] = 1.;
    else tmp_flag[iloc] = 0.;
  }
#ifdef COARSE_GRID_FLAG
  bin_field(tmp_flag,lon_5km,lat_5km,MOD06_BADFLAG,Tb240Kfrac_cg,'f',
            N_ACROSS_5KM*na_5km,lonmin,latmin,ncgridlons,
            ncgridlats,dlon_cgrid,dlat_cgrid,cgrid_lons,cgrid_lats);
            
#endif
  FREE(tmp_flag);
  FREE(lat_5km);
  FREE(lon_5km);
  FREE(IRrad);
  FREE(IRrad_uncert);
  FREE(land_mask);
  FREE(Tbs);
  FREE(Tb_31);
  FREE(Tb_32);
  FREE(mod_zenith);
  
/***** calculate total viewed area of granule ******/
  viewed_area=0.;
  for (i=0; i<ngridlats; i++)
    for (j=0; j<ngridlons; j++) {
      iloc = (i*ngridlons)+j;
      if (Tb32_g[iloc]!=MOD06_BADFLAG) viewed_area += areas[i];
  }
  

#ifdef COARSE_GRID_FLAG

/************************ get MODIS AOD *****************************/
/**** MODIS daily AOD [lons:0.5-359.5 lats:-89.5-89.5 ****/
/****/
  modflttmp = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
  mod08lons = (float *) malloc(N_MOD08_LONS*sizeof(float));
  mod08lats = (float *) malloc(N_MOD08_LATS*sizeof(float));
  aod_mod = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
/****
  sprintf(aodPath,"%s%d/%d/",mod08PathRoot,year,doy);
  sprintf(cmd,"ls %s",aodPath);
  if ( (fpo = popen(cmd,"r")) != NULL) {
    while ( (s=fgets(dir_list,100,fpo)) != NULL ) {
      inttmp = pclose(fpo);
      sscanf(dir_list,"%s",mod08File);
****/
			sprintf(mod08listPath,"%sMODIS_%s/%s/",listpath_root,yyyy,ddd);
      ismodfile = get_MOD08_filename(mod08listPath,mod08filelist,year,doy,
                                     mod08File);
      if (ismodfile) read_mod08_aerosol_cloud(mod08PathRoot,mod08File,mod08lons,
                              mod08lats,modflttmp,modflttmp,modflttmp,modflttmp,
                               aod_mod,modflttmp,modflttmp,61);
      lonstmp = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
      latstmp = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
//printf("\n\n");
      for (i=0; i<N_MOD08_LATS; i++)
      for (j=0; j<N_MOD08_LONS; j++) {
        iloc = (i*N_MOD08_LONS)+j;
        latstmp[iloc] = mod08lats[i];
        if ( (!dateline_flag) && (mod08lons[j]>180) ) lonstmp[iloc] = mod08lons[j]-360.;
        else lonstmp[iloc] = mod08lons[j];
//printf("%f ",aod_mod[iloc]);
      }
      bin_field(aod_mod,lonstmp,latstmp,MOD08_BADFLAG,aod_cg,'f',
            N_MOD08_LONS*N_MOD08_LATS,lonmin,latmin,ncgridlons,
            ncgridlats,dlon_cgrid,dlat_cgrid,cgrid_lons,cgrid_lats);
      FREE(lonstmp);
      FREE(latstmp);
/****
    }
  } else printf("cannot list directory\n");
****/
  FREE(modflttmp);
  FREE(mod08lons);
  FREE(mod08lats);
  FREE(aod_mod);
#endif //end if COARSE_GRID_FLAG
  



/********* get AMSR-E brightness temperatures **********/
/***** longitudes are -179.99968 to 179.99968 ******/
  if (amsr_monthly_dir_flag) sprintf(amsrPath,"%s%s%s/",amsrPath_root,yyyy,mn);
  else sprintf(amsrPath,"%s",amsrPath_root);
  if (amsr_monthly_dir_flag) sprintf(amsrlistPath,"%sAMSR_E_L2A/%s%s/",listpath_root,yyyy,mn);
  else sprintf(amsrlistPath,"%s",listpath_root);
//printf("%s %d %d %d %d %d\n",amsrPath,year,mon,day,hour,minute);
  delt_amsr = getAMSR_fileName(amsrlistPath,year,mon,day,hour,minute,tmpfile,38);
  delt_amsr *= (24.*60.);
//printf("%s\n",tmpfile);
//fflush(NULL);
  if ( (strcmp(tmpfile,noneStr) != 0) || (delt_amsr>60.) ) {
    printf("%lf minutes between AMSR and MODIS files. ",delt_amsr);
    strncpy(verstr,tmpfile+35,2);
    amsr_vers = atoi(verstr);
    printf("AMSR file version=%d\n",amsr_vers);
//fflush(NULL);
    
    lon_amsr_a = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*
                                  sizeof(float));
    lat_amsr_a = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*
                                  sizeof(float));
    lon_amsr_b = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*
                                  sizeof(float));
    lat_amsr_b = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*
                                  sizeof(float));
    tb89_Ha = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*sizeof(float));
    tb89_Va = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*sizeof(float));
    tb89_Hb = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*sizeof(float));
    tb89_Vb = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*sizeof(float));
    lmask_a = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*sizeof(float));
    lmask_b = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*sizeof(float));
    strcpy(amsrfile,amsrPath);
    strncat(amsrfile,tmpfile,NCHAR_FNAME_AMSR_TBS);
    printf("%s\n",amsrfile);
    fflush(NULL);
//printf("%d ",tb89_Vb);
    read_AMSR_pixel_tbs(amsrfile,&na_amsr,&nx_amsr,lon_amsr_a,
                        lat_amsr_a,lon_amsr_b,lat_amsr_b,
                        tb89_Ha,tb89_Va,tb89_Hb,tb89_Vb,lmask_a,lmask_b,amsr_vers);
//printf("%d\n",tb89_Vb);
//printf("%ld %ld\n",na_amsr,nx_amsr);
//for (i=0; i<na_amsr*nx_amsr; i++) printf("[%f %f] ",lon_amsr_b[i],lat_amsr_b[i]);
		if (dateline_flag) {
			for (i=0; i<na_amsr*nx_amsr; i++) {
				if (lon_amsr_a[i]<0) lon_amsr_a[i] = 360.+lon_amsr_a[i];
				if (lon_amsr_b[i]<0) lon_amsr_b[i] = 360.+lon_amsr_b[i];
			}
		}
//}
    bin_field(tb89_Hb,lon_amsr_b,lat_amsr_b,AMSR_BADFLAG,Tb89H_g,'f',
              na_amsr*nx_amsr,lonmin,latmin,ngridlons,
              ngridlats,dlon,dlat,grid_lons,grid_lats);
    bin_field(tb89_Vb,lon_amsr_b,lat_amsr_b,AMSR_BADFLAG,Tb89V_g,'f',
              na_amsr*nx_amsr,lonmin,latmin,ngridlons,
              ngridlats,dlon,dlat,grid_lons,grid_lats);
/*****
for (i=0; i<ngridlats; i++)
for (j=0; j<ngridlons; j++) {
iloc = (i*ngridlons)+j;
printf("[%f %f] ",Tb89H_g[iloc],Tb89V_g[iloc]);
}
*****/
/****
    bin_field(lmask_b,lon_amsr_b,lat_amsr_b,AMSR_BADFLAG,land_ocean,'f',
              na_amsr*nx_amsr,lonmin,latmin,ngridlons,
              ngridlats,dlon,dlat,grid_lons,grid_lats);
****/
    for (i=0; i<ngridlats; i++)
      for (j=0; j<ngridlons; j++) {
        iloc = (i*ngridlons)+j;
//        if (land_ocean[iloc] != AMSR_BADFLAG) land_ocean[iloc] /= 100.;
        if ( (Tb89H_g[iloc]!=AMSR_BADFLAG) &&
             (Tb89V_g[iloc]!=AMSR_BADFLAG) ) {
          pct[iloc] = (1.818*Tb89V_g[iloc]) - (0.818*Tb89H_g[iloc]);
        } else pct[iloc] = AMSR_BADFLAG;
    }
    FREE(lon_amsr_a);
    FREE(lat_amsr_a);
    FREE(lon_amsr_b);
    FREE(lat_amsr_b);
    FREE(tb89_Ha);
    FREE(tb89_Va);
    FREE(tb89_Hb);
    FREE(tb89_Vb);
    FREE(lmask_a);
    FREE(lmask_b);
  } else {
    printf("no AMSR-E file found.\n");
    amsr_fileflag = 0;
    for (i=0; i<ngridlats; i++)
      for (j=0; j<ngridlons; j++) {
        iloc = (i*ngridlons)+j;
        pct[iloc] = AMSR_BADFLAG;
    }
  }




#ifdef DAO_FLAG
/************ legacy ********/
#endif //end if DAO_FLAG





#ifdef DAS_FLAG
/***************
reverse the bad_flag sign so that bad flagged grid cells don't
appear to DAS as really cold cloud tops
This has to be done here because the sign of the IR
brightness temperatures is flipped inside go_das
***************/
  for (i=0; i<ngridlats; i++)
    for (j=0; j<ngridlons; j++) {
      iloc = (i*ngridlons)+j;
      if (Tb31_g[iloc] == MOD06_BADFLAG) Tb31_g[iloc]=-MOD06_BADFLAG;
      if (Tb32_g[iloc] == MOD06_BADFLAG) Tb32_g[iloc]=-MOD06_BADFLAG;
  }
#endif //end if DAS_FLAG



/********* find the clouds *********/
/***



Need to flip the IR bad_flag sign from -9999. to 9999!!!!!



****/
  tmp_tb = (float *) malloc(ngridlons*ngridlats*sizeof(float));
  
#ifdef DAS_FLAG
/***** Sign of the Tbs flipped at beginning of go_das ****/
//printf("doing DAS on Tb_31\n");
#ifdef VERBOSE
	go_das(ngridlats,ngridlons,Tb31_g,class31,cloud31,&ncloud31,0,yyyyddd,hhmm,outPath);
#endif //end if VERBOSE
#ifndef VERBOSE
  go_das(ngridlats,ngridlons,Tb31_g,class31,cloud31,&ncloud31,0);
#endif //end ifn VERBOSE
#endif //end if DAS_FLAG
#ifndef DAS_FLAG
  for (i=0; i<ngridlats; i++)
    for (j=0; j<ngridlons; j++) {
      iloc = (i*ngridlons)+j;
      if (Tb31_g[iloc] != MOD06_BADFLAG) tmp_tb[iloc] = -Tb31_g[iloc];
      else tmp_tb[iloc] = MOD06_BADFLAG;
  }
  label_region(ngridlats,ngridlons,tmp_tb,cloud31,-240,&ncloud31);
#endif //end if not DAO_FLAG

#ifdef DAS_FLAG
/***** Sign of the Tbs already flipped at beginning of go_das ****/
//printf("doing DAS on Tb_32\n");
#ifndef VERBOSE
  go_das(ngridlats,ngridlons,Tb32_g,class32,cloud32,&ncloud32,0);
#endif //end ifn VERBOSE
#endif //end if DAS_FLAG
#ifndef DAS_FLAG
    for (i=0; i<ngridlats; i++)
    for (j=0; j<ngridlons; j++) {
      iloc = (i*ngridlons)+j;
      if (Tb32_g[iloc] != MOD06_BADFLAG) tmp_tb[iloc] = -Tb32_g[iloc];
      else tmp_tb[iloc] = MOD06_BADFLAG;
  }
  label_region(ngridlats,ngridlons,tmp_tb,cloud32,-240,&ncloud32);
#endif //end if not DAO_FLAG
  FREE(tmp_tb);

printf("ncloud31=%ld\n",ncloud31);
printf("ncloud32=%ld\n",ncloud32);
fflush(NULL);



	if (ncloud31==0) continue_flag=0;
	if (!continue_flag) {
		printf("no clouds found. Aborting.\n");
		exit(0);
	}

#ifdef DAS_FLAG
/***************
Reverse the badflag values back to my standard value
***************/
  for (i=0; i<ngridlats; i++)
    for (j=0; j<ngridlons; j++) {
      iloc = (i*ngridlons)+j;
      if (Tb31_g[iloc] == -MOD06_BADFLAG) Tb31_g[iloc]=MOD06_BADFLAG;
      if (Tb32_g[iloc] == -MOD06_BADFLAG) Tb32_g[iloc]=MOD06_BADFLAG;
  }
#endif //end if not DAS_FLAG




/******* Get the reff vs IR_Tb profile *********
1)Need to decide on number of bt bins to average the results over
2)Allocate array of size [nclouds,nbtbins]
3)create loop over the clouds creating a new vector of der and bt of the values for the given cloud each time through
4)call get_reff_profile
get_reff_profile(float *der,float *bt,int N,float btmin,
                      float btmax,float nbtbins,float *der_mean,
                      float *bins)
**/
	cloud_pix = (long *) malloc(N_ACROSS_L1B_1KM_MAX*N_ALONG_L1B_1KM_MAX*sizeof(long));
	cloudpix_5km = (long *) malloc(N_ACROSS_5KM*N_ALONG_5KM_MAX*sizeof(long));
	ncloud_pix = (long *) calloc(ncloud32,sizeof(long));
//	labelCloudPixels(lon_1km,lat_1km,N_ACROSS_L1B_1KM_MAX*na_1km,ngridlons,ngridlats,
//									dlon,dlat,lonmin,latmin,cloud31,ncloud_pix,cloud_pix);
	labelCloudPixels(lon_1km,lat_1km,N_ACROSS_L1B_1KM_MAX*na_1km,ngridlons,ngridlats,
									dlon,dlat,lonmin,latmin,cloud32,ncloud_pix,cloud_pix,grid_lons,grid_lats);
									
	der_bin = (float *) calloc(ncloud32*NIR_BTBINS,sizeof(float));
	bt_bin = (float *) malloc(NIR_BTBINS*sizeof(float));
	
	for (ic=0; ic<ncloud32; ic++) {
//		printf("\n[%ld %ld %ld] ",ic,ncloud_pix[ic],ncloud32); //ncloud_pix=0 for most clouds
		fflush(NULL);
		der = (float *) malloc(ncloud_pix[ic]*sizeof(float));
		bt = (float *) malloc(ncloud_pix[ic]*sizeof(float));
		i=0;
		for (ix=0; ix<N_ACROSS_1KM; ix++) //Seg fault somewhere in this loop
		for (ia=0; ia<na_1km; ia++) {
  		iloc = (ia*N_ACROSS_1KM) + ix;
  		if (cloud_pix[iloc]==ic+1) {
			if (i<ncloud_pix[ic]) {
  				der[i] = reff[iloc];
  				bt[i] = Tb31_1km[iloc];
			}
//			if (i==ncloud_pix[ic]-1) printf("[%ld %ld] ",i,ncloud_pix[ic]);
//			if (ic==80) {
//  			printf("%ld %ld ",i,iloc);
//  			fflush(NULL);
//			}
//  			printf("%f %f ",reff[iloc],Tb31_1km[iloc]);
//  			fflush(NULL);
  			i++;
  		}
  	}
	 	get_reff_profile(der,bt,ncloud_pix[ic],IR_BTMIN,IR_BTMAX,NIR_BTBINS,der_bin+(ic*NIR_BTBINS),bt_bin);
		FREE(der);
		FREE(bt);
	}

//test that bt_bin is getting properly populated

	FREE(cloud_pix);
	FREE(cloudpix_5km);
//	FREE(ncloud_pix);
	FREE(reff);
	FREE(cldtop_phase);
  FREE(lat_1km);
  FREE(lon_1km);
  FREE(Tb31_1km);
	
	
	
	
	
	ncloud = ncloud32; //which set of clouds to use for cloud statistics?
  for (i=0; i<ngridlats; i++)
	for (j=0; j<ngridlons; j++) {
		iloc = (i*ngridlons)+j;
		cloud[iloc] = cloud32[iloc];
	}		
	
//	exit(0);
	
/********* calc the cloud stats *********/
  cloud31_areas = (float *) calloc(ncloud31,sizeof(float));
  cloud32_areas = (float *) calloc(ncloud32,sizeof(float));
  cloud31_class = (short *) calloc(ncloud31,sizeof(short));
  cloud32_class = (short *) calloc(ncloud32,sizeof(short));
  minTb32 = (float *) calloc(ncloud,sizeof(float));
  Tb89V_ave = (float *) calloc(ncloud,sizeof(float));
  Tb89H_ave = (float *) calloc(ncloud,sizeof(float));
  minTb89V = (float *) calloc(ncloud,sizeof(float));
  minTb89H = (float *) calloc(ncloud,sizeof(float));
  landfrac = (float *) calloc(ncloud,sizeof(float));
  pct_ave = (float *) calloc(ncloud,sizeof(float));
  minpct = (float *) calloc(ncloud,sizeof(float));
  cape_ave = (float *) calloc(ncloud,sizeof(float));
  dcape_ave = (float *) calloc(ncloud,sizeof(float));
  max_cape = (float *) calloc(ncloud,sizeof(float));
  max_dcape = (float *) calloc(ncloud,sizeof(float));
  shear_ave = (float *) calloc(ncloud,sizeof(float));
  max_shear = (float *) calloc(ncloud,sizeof(float));
  max_modzen = (float *) calloc(ncloud,sizeof(float));
  lat_minTb32 = (float *) calloc(ncloud,sizeof(float));
  lon_minTb32 = (float *) calloc(ncloud,sizeof(float));
  lat_minpct = (float *) calloc(ncloud,sizeof(float));
  lon_minpct = (float *) calloc(ncloud,sizeof(float));
  Tb31_220frac = (float *) calloc(ncloud,sizeof(float));
  Tb31_210frac = (float *) calloc(ncloud,sizeof(float));
  Tb31_200frac = (float *) calloc(ncloud,sizeof(float));
  Tb31_minpct = (float *) calloc(ncloud,sizeof(float));
  pct_minTb32 = (float *) calloc(ncloud,sizeof(float));
  cloud_datetime = (long *) calloc(ncloud,sizeof(long));
  cape_minTb32 = (float *) calloc(ncloud,sizeof(float));
  shear_minTb32 = (float *) calloc(ncloud,sizeof(float));
  cape_minpct = (float *) calloc(ncloud,sizeof(float));
  shear_minpct = (float *) calloc(ncloud,sizeof(float));
  wdir850_minTb32 = (float *) calloc(ncloud,sizeof(float));
  wdir500_minTb32 = (float *) calloc(ncloud,sizeof(float));
  wdir200_minTb32 = (float *) calloc(ncloud,sizeof(float));
  aod_minTb32 = (float *) calloc(ncloud,sizeof(float)); //note that this is 5-deg MODIS, so not convective-scale
  aod_merra_minTb32 = (float *) calloc(ncloud,sizeof(float)); //This at MERRA2 native resolution
  aod_merra_ave = (float *) calloc(ncloud,sizeof(float));
  
  mdyhms_to_juldate(mon,day,year,hour,minute,0.,&jd);
  mdyhms_to_juldate(1,1,1990,0,0,0.,&jd0);
  min_cent = (long) ((jd-jd0)*24.*60.);
  for (i=0; i<ncloud; i++) {
  	cloud_datetime[i] = min_cent;
//  	cloud_datetime[i] = malloc(17*sizeof(char));
//  	sprintf(cloud_datetime[i],"%d-%d-%d %d:%d",mon,day,year,hour,minute);
  }

  calc_cloudarea_pix(cloud31,areas,ngridlats,ngridlons,
                     ncloud31,cloud31_areas);
  calc_cloudarea_pix(cloud32,areas,ngridlats,ngridlons,
                     ncloud32,cloud32_areas);
  calc_mincloudstat_bin(cloud,ngridlons,ngridlats,ncloud,Tb32_g,
                        MOD06_BADFLAG,minTb32);
  if (amsr_fileflag) {
		calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,Tb89V_g,
											 AMSR_BADFLAG,Tb89V_ave);
		calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,Tb89H_g,
											 AMSR_BADFLAG,Tb89H_ave);
		calc_mincloudstat_bin(cloud,ngridlons,ngridlats,ncloud,Tb89V_g,
											 AMSR_BADFLAG,minTb89V);
		calc_mincloudstat_bin(cloud,ngridlons,ngridlats,ncloud,Tb89H_g,
											 AMSR_BADFLAG,minTb89H);
		calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,land_ocean,
											 AMSR_BADFLAG,landfrac);
		calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,pct,
											 AMSR_BADFLAG,pct_ave);
		calc_mincloudstat_bin(cloud,ngridlons,ngridlats,ncloud,pct,
											 AMSR_BADFLAG,minpct);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,pct,
															Tb32_g,AMSR_BADFLAG,pct_minTb32);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,Tb32_g,
															pct,AMSR_BADFLAG,Tb31_minpct);
		for (i=0; i<ncloud; i++) {
			lat_minpct[i] = AMSR_BADFLAG;
			lon_minpct[i] = AMSR_BADFLAG;
			cape_minpct[i] = AMSR_BADFLAG;
			shear_minpct[i] = AMSR_BADFLAG;
		}
		find_stat_location_bin(cloud,minpct,pct,grid_lons,grid_lats,
									 ngridlons,ngridlats,lon_minpct,lat_minpct);
  } else {
    for (i=0; i<ncloud; i++) {
      Tb89V_ave[i] = AMSR_BADFLAG;
      Tb89H_ave[i] = AMSR_BADFLAG;
      minTb89V[i] = AMSR_BADFLAG;
      minTb89H[i] = AMSR_BADFLAG;
      landfrac[i] = AMSR_BADFLAG;
      pct_ave[i] = AMSR_BADFLAG;
      minpct[i] = AMSR_BADFLAG;
      lon_minpct[i] = AMSR_BADFLAG;
      lat_minpct[i] = AMSR_BADFLAG;
      Tb31_minpct[i] = AMSR_BADFLAG;
      pct_minTb32[i] = AMSR_BADFLAG;
    }
  }
  
/**** This block was inside the if (amsr_fileflag) block just above up through v8, but ****/
/****	probably shouldn't have been ****/ 
	Tb31_220_flag = (float *) malloc(ngridlons*ngridlats*sizeof(float));
	Tb31_210_flag = (float *) malloc(ngridlons*ngridlats*sizeof(float));
	Tb31_200_flag = (float *) malloc(ngridlons*ngridlats*sizeof(float));

	for (i=0; i<ngridlats; i++)
	for (j=0; j<ngridlons; j++) {

		iloc = (i*ngridlons)+j;
		if (Tb31_g[iloc]<220 && Tb31_g[iloc]!=-9999.)
			Tb31_220_flag[iloc]=1; 
		else Tb31_220_flag[iloc]=0;
		if (Tb31_g[iloc]<210 && Tb31_g[iloc]!=-9999.)
			Tb31_210_flag[iloc]=1;
		else Tb31_210_flag[iloc]=0;
		if (Tb31_g[iloc]<200 && Tb31_g[iloc]!=-9999.) 
			Tb31_200_flag[iloc]=1;
		else Tb31_200_flag[iloc]=0;

		if (cloud31[iloc] > 0) cloud31_class[cloud31[iloc]-1] = class31[iloc];
		if (cloud32[iloc] > 0) cloud32_class[cloud32[iloc]-1] = class32[iloc];

	}
	calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,Tb31_220_flag,
										 MOD06_BADFLAG,Tb31_220frac);
	calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,Tb31_210_flag,
										 MOD06_BADFLAG,Tb31_210frac);
	calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,Tb31_200_flag,
										 MOD06_BADFLAG,Tb31_200frac);
	FREE(Tb31_220_flag);
	FREE(Tb31_210_flag);
	FREE(Tb31_200_flag);
/**** end of the block that was misplaced up to v8 ****/
	
  calc_maxcloudstat_bin(cloud,ngridlons,ngridlats,ncloud,modzen_g,
                     -9999.,max_modzen);
  for (i=0; i<ncloud; i++) {
    lat_minTb32[i] = AMSR_BADFLAG;
    lon_minTb32[i] = AMSR_BADFLAG;
  }
  find_stat_location_bin(cloud,minTb32,Tb32_g,grid_lons,grid_lats,
                 ngridlons,ngridlats,lon_minTb32,lat_minTb32);

#ifdef COARSE_GRID_FLAG
/**************
Here to assign the coarse gridded AOD corresponding to the lat/lon of the
minimum IR Tb in the cloud. At the moment, the routine assumes that the two variables
being matched up are on the same grid. Need to generalize this.
	calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,aod_cg,
															Tb32_g,MOD08_BADFLAG,aod_minTb32);
**************/
	for (i=0; i<ncloud; i++) {
		get_grid_indices(lat_minTb32[i],lon_minTb32[i],ncgridlons,ncgridlats,
                      dlat_cgrid,dlon_cgrid,latmin,lonmin,&ila,&ilo);
    if ( (ila != -9999.) && (ilo != -9999.) ) {
	iloc = (long) (ila*ncgridlons)+ilo;
    	aod_minTb32[i] = aod_cg[iloc];
    } else aod_minTb32[i] = -9999.;
	}
#endif
#ifndef COARSE_GRID_FLAG
	for (i=0; i<ncloud; i++) aod_minTb32[i] = -9999.;
#endif
                 
                 
                 
  Tb31_binmin=180;
  Tb31_binmax=300;
  Tb31_binsize=1;
  pct_binmin=100;
  pct_binmax=300;
  pct_binsize=1;
  Tb31_nhistbins = get_histogram_nbins(Tb31_binmin,Tb31_binmax,Tb31_binsize);
  Tb31_hist = (unsigned int *)
              calloc(ncloud*Tb31_nhistbins,sizeof(unsigned int));
  Tb31_histcnt = (unsigned long *) calloc(ncloud,sizeof(unsigned long));
  Tb31_histbins = (float *) calloc(Tb31_nhistbins,sizeof(float));
  get_histogram_stats(Tb31_binmin,Tb31_binsize,Tb31_nhistbins,Tb31_histbins);
  get_cloudstat_histogram(cloud,Tb31_hist,Tb31_histcnt,Tb32_g,ngridlons,
                          ngridlats,Tb31_nhistbins,Tb31_binmin,Tb31_binmax,
                          Tb31_binsize);
  pct_nhistbins = get_histogram_nbins(pct_binmin,pct_binmax,pct_binsize);
  pct_hist = (unsigned int *)
              calloc(ncloud*pct_nhistbins,sizeof(unsigned int));
  pct_histcnt = (unsigned long *) calloc(ncloud,sizeof(unsigned long));
  pct_histbins = (float *) calloc(pct_nhistbins,sizeof(float));
  get_histogram_stats(pct_binmin,pct_binsize,pct_nhistbins,pct_histbins);
  get_cloudstat_histogram(cloud,pct_hist,pct_histcnt,pct,ngridlons,
                          ngridlats,pct_nhistbins,pct_binmin,pct_binmax,
                          pct_binsize);
#ifdef VERBOSE
  sprintf(outfile,"%sTb_histogram_bins.bin",outPath);
  fp = fopen(outfile,"w");
  if (fp != NULL) {
    fwrite(&Tb31_nhistbins,sizeof(int),1,fp);
    fwrite(Tb31_histbins,sizeof(float),Tb31_nhistbins,fp);
    fwrite(&pct_nhistbins,sizeof(int),1,fp);
    fwrite(pct_histbins,sizeof(float),pct_nhistbins,fp);
    fclose(fp);
    fp = NULL;
  }
  sprintf(outfile,"%sTb_histograms.bin",outPath);
  fp = fopen(outfile,"a");
  if (fp != NULL) {
    fwrite(&doy,sizeof(int),1,fp);
    fwrite(&hour,sizeof(int),1,fp);
    fwrite(&minute,sizeof(int),1,fp);
    fwrite(&ncloud,sizeof(long),1,fp);
    fwrite(Tb31_hist,sizeof(unsigned int),Tb31_nhistbins*ncloud,fp);
    fwrite(Tb31_histcnt,sizeof(unsigned int),ncloud,fp);
    fwrite(pct_hist,sizeof(unsigned int),pct_nhistbins*ncloud,fp);
    fwrite(pct_histcnt,sizeof(unsigned int),ncloud,fp);
    fclose(fp);
    fp = NULL;
  }
#endif
  FREE(Tb31_histbins);
  FREE(Tb31_hist);
  FREE(Tb31_histcnt);
  FREE(pct_histbins);
  FREE(pct_hist);
  FREE(pct_histcnt);
  

/********************
Flip from 0-360 to -180 to 180 in lon_minTb32 and lon_minpct if the
dateline_flag is true
********************/
	if (dateline_flag) {
    for (i=0; i<ncloud; i++) {
    	if (lon_minTb32[i]>180) lon_minTb32[i] = lon_minTb32[i]-360.;
    	if (lon_minpct[i]>180) lon_minpct[i] = lon_minpct[i]-360.;
    }
  }
  
  
//  exit(0);




/****************************************
CAPE calculation moved here to speed up code to just calculate where
needed.
MERRA2 longitudes are -180 to 180
****************************************/
#ifdef MERRA_FLAG
  lon_merra = (double *) malloc(N_MERRA_NLONS*sizeof(double));
  lat_merra = (double *) malloc(N_MERRA_NLATS*sizeof(double));
  P_merra = (double *) malloc(N_MERRA_NLEVS*sizeof(double));
  levs_merra = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
  t_merra = (int *) malloc(N_MERRA_NTIMES*sizeof(int));
  temp_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*
  													N_MERRA_NTIMES*sizeof(float));
  q_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*
  													N_MERRA_NTIMES*sizeof(float));
  uwnd_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*
  													N_MERRA_NTIMES*sizeof(float));
  vwnd_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*
  													N_MERRA_NTIMES*sizeof(float));
  aod_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*
  													N_MERRA_INST3_NTIMES*sizeof(float));
//  delP_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*
//  													N_MERRA_NTIMES*sizeof(float));
  P0_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*
  													N_MERRA_NTIMES*sizeof(float));
	if ( (year>=2001) && (year<=2010) ) {
		sprintf(merra_file,"MERRA2_300.%s%ld.nc4",merra_fileroot,yyyymmdd);
		sprintf(merra_aer_file,"MERRA2_300.%s%ld.nc4",merra_aer_fileroot,yyyymmdd);
	}
  if (year>2010) {
  	sprintf(merra_file,"MERRA2_400.%s%ld.nc4",merra_fileroot,yyyymmdd);
		sprintf(merra_aer_file,"MERRA2_400.%s%ld.nc4",merra_aer_fileroot,yyyymmdd);
  }
  printf("reading %s%s\n",merraPath,merra_file);
  fflush(NULL);
  is_merrafile = read_merra_float(merraPath,merra_file,"lon",lon_merra,&merra_fill);
  is_merrafile = read_merra_float(merraPath,merra_file,"lat",lat_merra,&merra_fill);
  is_merrafile = read_merra_float(merraPath,merra_file,"lev",P_merra,&merra_fill);
  is_merrafile = read_merra_float(merraPath,merra_file,"time",t_merra,&merra_fill);
  is_merrafile = read_merra_float(merraPath,merra_file,"T",temp_merra,&merra_fill);
  if ( is_merrafile && (merra_fill!=0) ) {
//	printf("merra_fill=%f\n",merra_fill);
//	fflush(NULL);
  	for (i=0; i<N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*N_MERRA_NTIMES; i++)
  		if (temp_merra[i]==merra_fill) temp_merra[i]=MERRA_BADFLAG;
  }
  read_merra_float(merraPath,merra_file,"QV",q_merra,&merra_fill);
  if ( is_merrafile && (merra_fill!=0) ) {
//        printf("merra_fill=%f\n",merra_fill);
//        fflush(NULL);
  	for (i=0; i<N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*N_MERRA_NTIMES; i++)
  		if (q_merra[i]==merra_fill) q_merra[i]=MERRA_BADFLAG;
  }
  read_merra_float(merraPath,merra_file,"U",uwnd_merra,&merra_fill);
  if ( is_merrafile && (merra_fill!=0) ) {
//        printf("merra_fill=%f\n",merra_fill);
//        fflush(NULL);
  	for (i=0; i<N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*N_MERRA_NTIMES; i++)
  		if (uwnd_merra[i]==merra_fill) uwnd_merra[i]=MERRA_BADFLAG;
  }
  read_merra_float(merraPath,merra_file,"V",vwnd_merra,&merra_fill);
  if ( is_merrafile && (merra_fill!=0) ) {
//        printf("merra_fill=%f\n",merra_fill);
//        fflush(NULL);
  	for (i=0; i<N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NLEVS*N_MERRA_NTIMES; i++)
  		if (vwnd_merra[i]==merra_fill) vwnd_merra[i]=MERRA_BADFLAG;
  }
  read_merra_float(merraPath,merra_file,"PS",P0_merra,&merra_fill);
  if ( is_merrafile && (merra_fill!=0) ) {
  	for (i=0; i<N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_NTIMES; i++)
  		if (P0_merra[i]==merra_fill) P0_merra[i]=MERRA_BADFLAG;
  }
  printf("reading %s%s\n",merraPath,merra_aer_file);
  fflush(NULL);
  is_merra_aerfile = read_merra_float(merraPath,merra_aer_file,"AODANA",
  																		aod_merra,&merra_fill);
  if ( is_merra_aerfile && (merra_fill!=0) ) {
  	for (i=0; i<N_MERRA_NLONS*N_MERRA_NLATS*N_MERRA_INST3_NTIMES; i++)
  		if (aod_merra[i]==merra_fill) aod_merra[i]=MERRA_BADFLAG;
  }
//  read_merra_float(merraPath,merra_file,"DELP",delP_merra);
  printf("done reading MERRA\n");
  fflush(NULL);

//exit(0);
	
	if (is_merrafile) {
	
		T_prof = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
		T_prof_tmp = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
		q_prof = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
		RH_prof = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
		RH_prof_tmp = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
	//  z_prof = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
	//  z_prof_tmp = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
		P_prof = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
		P_prof_tmp = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
		Tv_prof = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
	//  delP_prof = (float *) malloc(N_MERRA_NLEVS*sizeof(float));
	//  sprintf(merra_file,"%s%ld.hdf",merra_fileroot,yyyymmdd);
	
		for (ilev=0; ilev<N_MERRA_NLEVS; ilev++) levs_merra[ilev] = P_merra[ilev]*100.;

		cape_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*
															N_MERRA_NTIMES*sizeof(float));
		shear_merra = (float *) malloc(N_MERRA_NLONS*N_MERRA_NLATS*
															N_MERRA_NTIMES*sizeof(float));
		for (i=0; i<N_MERRA_NLATS; i++)
		for (j=0; j<N_MERRA_NLONS; j++) 
		for (it=0; it<N_MERRA_NTIMES; it++) {
			iloc = (it*N_MERRA_NLATS*N_MERRA_NLONS) + (i*N_MERRA_NLONS) + j;
			cape_merra[iloc] = MERRA_BADFLAG;
			shear_merra[iloc] = MERRA_BADFLAG;
		}
	
	/*****  
		for (i=0; i<N_MERRA_NLATS; i++) printf("%f ",lat_merra[i]);
		printf("\n\n");
		for (j=0; j<N_MERRA_NLONS; j++) printf("%f ",lon_merra[j]);
	****/
flag=2;
		for (i=0; i<N_MERRA_NLATS; i++)
		for (j=0; j<N_MERRA_NLONS; j++) {
		
	//		if ( (dateline_flag) && (lon_merra[j]<0) ) ilon = (int) floor((360.+lon_merra[j]-lonmin)/dlon);
	//    else ilon = (int) floor((lon_merra[j]-lonmin)/dlon);
	//    ilat = (int) floor((lat_merra[i]-latmin)/dlat);
	//    kloc = (ilat*ngridlons)+ilon;
	//    printf("[%f %f %ld %ld] ",lat_merra[i],lon_merra[j],ilat,ilon);
	
	//    if ( (ilon>=0) && (ilon<ngridlons) && (ilat>=0) && (ilat<ngridlats) &&
	//    			(Tb31_g[kloc] != MOD06_BADFLAG) ) {
					
			if ( (dateline_flag) && (lon_merra[j]<0) ) lontest = lon_merra[j]+360.;
			else lontest = lon_merra[j];
			if ( (lontest>=lonmin) && (lontest<=lonmax) &&
						(lat_merra[i]>=latmin) && (lat_merra[i]<=latmax) ) {
			
	//    printf("[%f %f %f %f] ",lat_merra[i],lon_merra[j],grid_lats[ilat],grid_lons[ilon]);

	//    if (lon_merra[j]<lonmin+2 && lat_merra[i]<latmin+2)
	//    	printf("[%ld %ld %f]\n",ilon,ilat,Tb31_g[kloc]);
	/**
				capecnt++;
				if ( (capecnt % 100) == 0) {
					printf("%ld [%f %f]\n",capecnt,lat_merra[i],lon_merra[j]);
					fflush(NULL);
				}
	**/
	/****    
				baseloc = get_merra_indices_4D(lat_merra[i],lon_merra[j],850,hour,N_MERRA_NLONS,
																			N_MERRA_NLATS,N_MERRA_NLEVS,N_MERRA_NTIMES,P_merra,
																			 &ila,&ilo,&ilev,&it);
				toploc = get_merra_indices_4D(lat_merra[i],lon_merra[j],200,hour,N_MERRA_NLONS,
																			N_MERRA_NLATS,N_MERRA_NLEVS,N_MERRA_NTIMES,P_merra,
																			&ila,&ilo,&ilev,&it);
	*****/
				it = (int) floor((float) hour/(24./N_MERRA_NTIMES));
				baseloc = (it*N_MERRA_NLEVS*N_MERRA_NLATS*N_MERRA_NLONS) + 
									(6*N_MERRA_NLATS*N_MERRA_NLONS) + (i*N_MERRA_NLONS) + j;
				toploc = (it*N_MERRA_NLEVS*N_MERRA_NLATS*N_MERRA_NLONS) + 
									(22*N_MERRA_NLATS*N_MERRA_NLONS) + (i*N_MERRA_NLONS) + j;
																		
				jloc = (it*N_MERRA_NLATS*N_MERRA_NLONS) +
									(i*N_MERRA_NLONS) + j;
			
				shear_merra[jloc] = 
						calc_shear(uwnd_merra[baseloc],uwnd_merra[toploc],
										vwnd_merra[baseloc],vwnd_merra[toploc]);
	
				tstep=0;
				while (tstep <=1) {
					jloc = ((it-tstep)*N_MERRA_NLATS*N_MERRA_NLONS) +
									(i*N_MERRA_NLONS) + j;
					profilev=0;
					for (ilev=0; ilev<=N_LEVELS_CAPE; ilev++) {
						kloc = ((it-tstep)*N_MERRA_NLEVS*N_MERRA_NLATS*N_MERRA_NLONS) +
									 (ilev*N_MERRA_NLATS*N_MERRA_NLONS) +
									 (i*N_MERRA_NLONS) + j;
						if (P0_merra[jloc] >= P_merra[ilev]*100.) {
							T_prof[profilev] = temp_merra[kloc];
							q_prof[profilev] = q_merra[kloc];
							P_prof[profilev] = P_merra[ilev];
							profilev++;
						} else {
							
						}
	//					delP_prof[ilev] = delP_merra[kloc];
					}
	//				calc_merra_native_P_profile(delP_prof,P_prof);
//					for (ilev=0; ilev<N_LEVELS_CAPE; ilev++) {
//						P_prof[ilev] = P_merra[ilev+1]*100.;
//						RH_prof[ilev] = calc_e_from_q(P_prof[ilev],q_prof[ilev]) / calc_es(T_prof[ilev]);
//						Tv_prof[ilev] = calc_Tv(T_prof[ilev],RH_prof[ilev],P_prof[ilev]);
//					}
	//				calc_merra_native_z_profile(delP_prof,Tv_prof,P_prof,z_prof);
	/*********
	For MERRA inst6_3d_ana_Nv data is arranged from TOA to surface, so need to flip it for the
	CAPE calculation.
	For MERRA2 inst6_3d_ana_Np data is arranged surface to TOA, so no need to flip
					for (ilev=0; ilev<N_MERRA_NLEVS; ilev++) {
						T_prof_tmp[N_MERRA_NLEVS-1-ilev] = T_prof[ilev];
						RH_prof_tmp[N_MERRA_NLEVS-1-ilev] = RH_prof[ilev];
						P_prof_tmp[N_MERRA_NLEVS-1-ilev] = P_merra[ilev]*100.;
	//					z_prof_tmp[N_MERRA_NLEVS-1-ilev] = z_prof[ilev];
					}
	*****/
	/*****
	for (ilev=0; ilev<N_MERRA_NLEVS; ilev++) {
	printf("%f %f %f %f\n",z_prof[ilev],P_prof[ilev],RH_prof[ilev],T_prof[ilev]);
	printf("%f %f %f %f\n",z_prof_tmp[ilev],P_prof_tmp[ilev],RH_prof_tmp[ilev],T_prof_tmp[ilev]);
	}
	*****/
	//				cape_merra[jloc] =
	//					calc_CAPE_emanuel(N_LEVELS_CAPE,T_prof_tmp,P_prof_tmp,RH_prof_tmp);
	//				cape_merra[jloc] =
	//					calc_CAPE_emanuel(N_LEVELS_CAPE,T_prof,P_prof,RH_prof);
					cape_merra[jloc] =
						calc_CAPE_emanuel(N_LEVELS_CAPE,T_prof,P_prof,q_prof); //now taking mixing ratio in kg/kg instead or RH and converting to g/kg in the subroutine
					
/*					if (flag<1 && cape_merra[jloc]>0) {
						printf("CAPE=0\n");
						printf("T\n");
						for (ilev=0; ilev<N_LEVELS_CAPE-1; ilev++) printf("%f ",T_prof[ilev]);
						printf("\n");
						printf("P\n");
						for (ilev=0; ilev<N_LEVELS_CAPE-1; ilev++) printf("%f ",P_prof[ilev]);
						printf("\n");
						printf("q\n");
						for (ilev=0; ilev<N_LEVELS_CAPE-1; ilev++) printf("%f ",q_prof[ilev]);
						printf("\n");
						flag++;
					}
					if (flag==1 && cape_merra[jloc]==0.) {
						printf("CAPE=0\n");
						printf("T\n");
						for (ilev=0; ilev<N_LEVELS_CAPE-1; ilev++) printf("%f ",T_prof[ilev]);
						printf("\n");
						printf("P\n");
						for (ilev=0; ilev<N_LEVELS_CAPE-1; ilev++) printf("%f ",P_prof[ilev]);
						printf("\n");
						printf("q\n");
						for (ilev=0; ilev<N_LEVELS_CAPE-1; ilev++) printf("%f ",q_prof[ilev]);
						printf("\n");
						flag++;
					} */
					
	//				if ( (tstep==0) && (cape_merra[jloc] <= 0) ) {
					if ( (tstep==0) ) {
	//					printf("CAPE = %f\n",cape_merra[jloc]);
	//					for (ilev=0; ilev<N_LEVELS_CAPE; ilev++) printf("%d %f %f %f\n",ilev,P_prof[ilev],T_prof[ilev],RH_prof[ilev]);
					}
					
	//exit(0);
					tstep++;
					if (it==0) {
						tstep++; //add one more to end the loop
		//      	printf("no earlier time step in the MERRA file\n");
					}
				}
			}
		}

	//exit(0);
	
		for (i=0; i<ngridlats; i++)
		for (j=0; j<ngridlons; j++) {
			iloc = (i*ngridlons)+j;
	//printf("%ld ",iloc);
	fflush(NULL);
			if ( (dateline_flag) && (grid_lons[j]>180) ) lontest = grid_lons[i]-360;
			else lontest = grid_lons[j];
			testloc = get_merra_indices_4D(grid_lats[i],lontest,850,hour,N_MERRA_NLONS,
																		N_MERRA_NLATS,N_MERRA_NLEVS,N_MERRA_NTIMES,P_merra,
																		 &ila,&ilo,&ilev,&it);
			jloc = (it*N_MERRA_NLATS*N_MERRA_NLONS) +
								(ila*N_MERRA_NLONS) + ilo;
			shear_g[iloc] = shear_merra[jloc];
		
			jloc = (it*N_MERRA_NLEVS*N_MERRA_NLATS*N_MERRA_NLONS) + 
									(6*N_MERRA_NLATS*N_MERRA_NLONS) + (ila*N_MERRA_NLONS) + ilo;
			wdir850[iloc] = 270. - (atan2(vwnd_merra[jloc],uwnd_merra[jloc]) * (180./3.1415926));
			if (wdir850[iloc]>=360) wdir850[iloc] = wdir850[iloc] - 360.;
			jloc = (it*N_MERRA_NLEVS*N_MERRA_NLATS*N_MERRA_NLONS) + 
									(16*N_MERRA_NLATS*N_MERRA_NLONS) + (ila*N_MERRA_NLONS) + ilo;
			wdir500[iloc] = 270. - (atan2(vwnd_merra[jloc],uwnd_merra[jloc]) * (180./3.1415926));
			if (wdir500[iloc]>=360) wdir500[iloc] = wdir500[iloc] - 360.;
			jloc = (it*N_MERRA_NLEVS*N_MERRA_NLATS*N_MERRA_NLONS) + 
									(22*N_MERRA_NLATS*N_MERRA_NLONS) + (ila*N_MERRA_NLONS) + ilo;
			wdir200[iloc] = 270. - (atan2(vwnd_merra[jloc],uwnd_merra[jloc]) * (180./3.1415926));
			if (wdir200[iloc]>=360) wdir200[iloc] = wdir200[iloc] - 360.;
		
			tstep=0;
			while (tstep <=1) {
				jloc = ((it-tstep)*N_MERRA_NLATS*N_MERRA_NLONS) +
								(ila*N_MERRA_NLONS) + ilo;
				if (tstep==0) cape_g[iloc] = cape_merra[jloc];
				if (tstep==1) dcape_g[iloc] = (cape_g[iloc] - cape_merra[jloc])/6.;
				tstep++;
			}
		
			testloc = get_merra_indices_4D(grid_lats[i],lontest,850,hour,N_MERRA_NLONS,
																		N_MERRA_NLATS,N_MERRA_NLEVS,N_MERRA_INST3_NTIMES,P_merra,
																		 &ila,&ilo,&ilev,&it);
			jloc = (it*N_MERRA_NLATS*N_MERRA_NLONS) +
								(ila*N_MERRA_NLONS) + ilo;
			aod_merra_g[iloc] = aod_merra[jloc];
		
	//  printf("[%ld %ld %f %f] ",i,j,cape_g[iloc],shear_g[iloc]);
	//  printf("%f ",cape_g[iloc]);
		}
	
		FREE(cape_merra);
		FREE(shear_merra);

	printf("Done with CAPE\n");
	fflush(NULL);
	
	#ifdef COARSE_GRID_FLAG
		latstmp = (float *) malloc(ngridlons*ngridlats*sizeof(float));
		lonstmp = (float *) malloc(ngridlons*ngridlats*sizeof(float));
		for (i=0; i<ngridlats; i++)
		for (j=0; j<ngridlons; j++) {
			iloc = (i*ngridlons)+j;
			latstmp[iloc] = grid_lats[i];
			lonstmp[iloc] = grid_lons[j];
		}
		bin_field(cape_g,lonstmp,latstmp,-9999.,cape_cg,'f',
							ngridlons*ngridlats,lonmin,latmin,ncgridlons,
							ncgridlats,dlon_cgrid,dlat_cgrid,cgrid_lons,cgrid_lats);
		bin_field(dcape_g,lonstmp,latstmp,-9999.,dcape_cg,'f',
							ngridlons*ngridlats,lonmin,latmin,ncgridlons,
							ncgridlats,dlon_cgrid,dlat_cgrid,cgrid_lons,cgrid_lats);
		bin_field(shear_g,lonstmp,latstmp,-9999.,shear_cg,'f',
							ngridlons*ngridlats,lonmin,latmin,ncgridlons,
							ncgridlats,dlon_cgrid,dlat_cgrid,cgrid_lons,cgrid_lats);
		FREE(latstmp);
		FREE(lonstmp);
	#endif //end if COARSE_GRID_FLAG
	
	//  FREE(delP_merra);
	  FREE(P0_merra);
		FREE(T_prof);
		FREE(T_prof_tmp);
		FREE(RH_prof);
		FREE(RH_prof_tmp);
		FREE(P_prof);
		FREE(P_prof_tmp);
		FREE(q_prof);
	//  FREE(z_prof);
	//  FREE(z_prof_tmp);
		FREE(Tv_prof);
	//  FREE(delP_prof);


		calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,cape_g,
											 AMSR_BADFLAG,cape_ave);
		calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,dcape_g,
											 AMSR_BADFLAG,dcape_ave);
		calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,shear_g,
											 -9999.,shear_ave);
		calc_cloudstat_bin(cloud,ngridlons,ngridlats,ncloud,aod_merra_g,
											 -9999.,aod_merra_ave);
		calc_maxcloudstat_bin(cloud,ngridlons,ngridlats,ncloud,cape_g,
											 AMSR_BADFLAG,max_cape);
		calc_maxcloudstat_bin(cloud,ngridlons,ngridlats,ncloud,dcape_g,
											 AMSR_BADFLAG,max_dcape);
		calc_maxcloudstat_bin(cloud,ngridlons,ngridlats,ncloud,shear_g,
											 -9999.,max_shear);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,cape_g,
																Tb32_g,AMSR_BADFLAG,cape_minTb32);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,shear_g,
																Tb32_g,AMSR_BADFLAG,shear_minTb32);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,cape_g,
																pct,AMSR_BADFLAG,cape_minpct);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,shear_g,
																pct,AMSR_BADFLAG,shear_minpct);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,wdir850,
																Tb32_g,AMSR_BADFLAG,wdir850_minTb32);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,wdir500,
																Tb32_g,AMSR_BADFLAG,wdir500_minTb32);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,wdir200,
																Tb32_g,AMSR_BADFLAG,wdir200_minTb32);
		calc_stat_at_mincloudstat(cloud,ngridlons,ngridlats,ncloud,aod_merra_g,
																Tb32_g,AMSR_BADFLAG,aod_merra_minTb32);
	} else { //end if (is_merrafile)
	
		for (i=0; i<ngridlats; i++)
		for (j=0; j<ngridlons; j++) {
			iloc = (i*ngridlons)+j;
			cape_g[iloc] = -9999.;
			dcape_g[iloc] = -9999.;
			shear_g[iloc] = -9999.;
		}
		for (i=0;i<ncloud;i++) {
			cape_ave[i] = -9999.;
			dcape_ave[i] = -9999.;
			shear_ave[i] = -9999.;
			max_cape[i] = -9999.;
			max_dcape[i] = -9999.;
			max_shear[i] = -9999.;
			cape_minTb32[i] = -9999.;
			shear_minTb32[i] = -9999.;
			cape_minpct[i] = -9999.;
			shear_minpct[i] = -9999.;
			wdir850_minTb32[i] = -9999.;
			wdir500_minTb32[i] = -9999.;
			wdir200_minTb32[i] = -9999.;
			aod_merra_ave[i] = -9999.;
			aod_merra_minTb32[i] = -9999.;
		}
	} //end else ~is_merrafile
	
	FREE(lon_merra);
	FREE(lat_merra);
	FREE(P_merra);
	FREE(levs_merra);
	FREE(t_merra);
	FREE(temp_merra);
	FREE(q_merra);
	FREE(uwnd_merra);
	FREE(vwnd_merra);
	FREE(aod_merra);

#endif //MERRA_FLAG
#ifndef MERRA_FLAG
  for (i=0; i<ngridlats; i++)
  for (j=0; j<ngridlons; j++) {
    iloc = (i*ngridlons)+j;
    cape_g[iloc] = -9999.;
    dcape_g[iloc] = -9999.;
    shear_g[iloc] = -9999.;
  }
  for (i=0;i<ncloud;i++) {
  	cape_ave[i] = -9999.;
  	dcape_ave[i] = -9999.;
  	shear_ave[i] = -9999.;
  	max_cape[i] = -9999.;
  	max_dcape[i] = -9999.;
  	max_shear[i] = -9999.;
  	cape_minTb32[i] = -9999.;
  	shear_minTb32[i] = -9999.;
  	cape_minpct[i] = -9999.;
  	shear_minpct[i] = -9999.;
  	wdir850_minTb32[i] = -9999.;
  	wdir500_minTb32[i] = -9999.;
  	wdir200_minTb32[i] = -9999.;
  	aod_merra_ave[i] = -9999.;
  	aod_merra_minTb32[i] = -9999.;
  }
#endif //NOT MERRA_FLAG


/*	for (i=0;i<ncloud;i++) {
		printf("%f ",max_cape[i]);
	}
	printf("\n");*/
  
//  exit(0);

	if (ncloud>0) { //only write the output netcdf files if there were clouds found
		sprintf(outfile,"%sIR_brightness_temperature_bins.nc",outPath);
		if ( (status = nc_create(outfile,NC_NOCLOBBER,&ncid)) == NC_NOERR ) {
			nc_def_dim(ncid,"number_of_IR_temp_bins_for_profile",NIR_BTBINS,&irbindimid);
			nc_def_var(ncid,"IR_brightness_temperature_bins_for_profile",NC_FLOAT,1,&irbindimid,&IRbinid);
			nc_enddef(ncid);
			sp = 0;
			cpbin = (size_t) NIR_BTBINS;
			status = nc_put_vara_float(ncid,IRbinid,&sp,&cpbin,bt_bin);
			nc_close(ncid);
		} else if (status!=-35) printf("Error opening file to write IR brightness temperature bins.\n");
	

		sprintf(outfile,"%s%s.%s.%s.nc",outPath,outfile_str,yyyyddd,hhmm);
		if ( (status = nc_create(outfile,NC_NOCLOBBER,&ncid)) == NC_NOERR ) {
	//  	nc_def_dim(ncid,"number_of_clouds_5km",ncloud31,&tb31dimid);
			nc_def_dim(ncid,"number_of_clouds_1km",NC_UNLIMITED,&tb32dimid);
			nc_def_dim(ncid,"number_of_IR_temp_bins_for_profile",NIR_BTBINS,&irbindimid);
	//  	nc_def_var(ncid,"cloud_areas_5km",NC_FLOAT,1,&tb31dimid,&area31id);
			nc_def_var(ncid,"minutes_since_Jan_1_1990_UTC",NC_LONG,1,&tb32dimid,&timeid);
			nc_def_var(ncid,"cloud_areas_1km",NC_FLOAT,1,&tb32dimid,&area32id);
	//  	nc_def_var(ncid,"cloud_classification_5km",NC_SHORT,1,&tb31dimid,&class31id);
			nc_def_var(ncid,"cloud_classification_1km",NC_SHORT,1,&tb32dimid,&class32id);
			nc_def_var(ncid,"fraction_of_cloud_over_land",NC_FLOAT,1,&tb32dimid,&landfracid);
			nc_def_var(ncid,"maximum_view_zenith_within_cloud_area",NC_FLOAT,1,&tb32dimid,&maxmodzenid);
			nc_def_var(ncid,"minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&minIRtbid);
			nc_def_var(ncid,"CAPE_averaged_over_cloud_area",NC_FLOAT,1,&tb32dimid,&aveCAPEid);
			nc_def_var(ncid,"change_of_CAPE_averaged_over_cloud_area",NC_FLOAT,1,&tb32dimid,&daveCAPEid);
			nc_def_var(ncid,"maximum_CAPE_within_cloud_area",NC_FLOAT,1,&tb32dimid,&maxCAPEid);
			nc_def_var(ncid,"maximum_change_of_CAPE_within_cloud_area",NC_FLOAT,1,&tb32dimid,&maxdCAPEid);
			nc_def_var(ncid,"shear_averaged_over_cloud_area",NC_FLOAT,1,&tb32dimid,&aveshearid);
			nc_def_var(ncid,"maximum_shear_within_cloud_area",NC_FLOAT,1,&tb32dimid,&maxshearid);
			nc_def_var(ncid,"CAPE_at_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&capeminTBid);
			nc_def_var(ncid,"shear_at_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&shearminTBid);
			nc_def_var(ncid,"CAPE_at_minimum_85GHz_PCT",NC_FLOAT,1,&tb32dimid,&capeminpctid);
			nc_def_var(ncid,"shear_at_minimum_85GHz_PCT",NC_FLOAT,1,&tb32dimid,&shearminpctid);
			nc_def_var(ncid,"PCT_85GHz_averaged_over_cloud_area",NC_FLOAT,1,&tb32dimid,&pctaveid);
			nc_def_var(ncid,"minimum_85GHz_PCT",NC_FLOAT,1,&tb32dimid,&minPCTid);
			nc_def_var(ncid,"fraction_of_cloud_area_colder_than_220K",NC_FLOAT,1,&tb32dimid,&frac220id);
			nc_def_var(ncid,"fraction_of_cloud_area_colder_than_210K",NC_FLOAT,1,&tb32dimid,&frac210id);
			nc_def_var(ncid,"fraction_of_cloud_area_colder_than_200K",NC_FLOAT,1,&tb32dimid,&frac200id);
			nc_def_var(ncid,"latitude_of_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&latminTBid);
			nc_def_var(ncid,"longitude_of_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&lonminTBid);
			nc_def_var(ncid,"latitude_of_minimum_85GHz_PCT",NC_FLOAT,1,&tb32dimid,&latminPCTid);
			nc_def_var(ncid,"longitude_of_minimum_85GHz_PCT",NC_FLOAT,1,&tb32dimid,&lonminPCTid);
			nc_def_var(ncid,"wind_direction_850hPa_at_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&w850minTBid);
			nc_def_var(ncid,"wind_direction_500hPa_at_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&w500minTBid);
			nc_def_var(ncid,"wind_direction_200hPa_at_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&w200minTBid);
			nc_def_var(ncid,"IR_brightness_temperature_at_location_of_minimum_85GHz_PCT",NC_FLOAT,1,&tb32dimid,&IRminPCTid);
			nc_def_var(ncid,"PCT_85GHz_at_location_of_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&PCTminIRid);
			nc_def_var(ncid,"aerosol_optical_thickness",NC_FLOAT,1,&tb32dimid,&AODminIRid);
			nc_def_var(ncid,"AOT_MERRA_averaged_over_cloud_area",NC_FLOAT,1,&tb32dimid,&AODmerraid);
			nc_def_var(ncid,"AOT_MERRA_at_location_of_minimum_IR_brightness_temperature",NC_FLOAT,1,&tb32dimid,&AODmerraminIRid);
			derbindims[0] = tb32dimid;
			derbindims[1] = irbindimid;
			status = nc_def_var(ncid,"particle_effective_radius_averaged_on_IR_bins",NC_FLOAT,2,derbindims,&derbinid);
	//  	printf("status=%d\n",status);
			nc_enddef(ncid);
			sp = 0;
			cp = (size_t) ncloud;
			cp31 = (size_t) ncloud31;
			cp32 = (size_t) ncloud32;
			cpbin = (size_t) NIR_BTBINS;
			spder[0] = 0;
			spder[1] = 0;
			cpder[0] = ncloud32;
			cpder[1] = NIR_BTBINS;
	//		status = nc_put_vara_float(ncid,area31id,&sp,&cp31,cloud31_areas);
			status = nc_put_vara_long(ncid,timeid,&sp,&cp32,cloud_datetime);
			status = nc_put_vara_float(ncid,area32id,&sp,&cp32,cloud32_areas);
	//		status = nc_put_vara_short(ncid,class31id,&sp,&cp31,cloud31_class);
			status = nc_put_vara_short(ncid,class32id,&sp,&cp32,cloud32_class);
			status = nc_put_vara_float(ncid,landfracid,&sp,&cp,landfrac);
			status = nc_put_vara_float(ncid,maxmodzenid,&sp,&cp,max_modzen);
			status = nc_put_vara_float(ncid,minIRtbid,&sp,&cp,minTb32);
			status = nc_put_vara_float(ncid,aveCAPEid,&sp,&cp,cape_ave);
			status = nc_put_vara_float(ncid,daveCAPEid,&sp,&cp,dcape_ave);
			status = nc_put_vara_float(ncid,maxCAPEid,&sp,&cp,max_cape);
			status = nc_put_vara_float(ncid,maxdCAPEid,&sp,&cp,max_dcape);
			status = nc_put_vara_float(ncid,aveshearid,&sp,&cp,shear_ave);
			status = nc_put_vara_float(ncid,maxshearid,&sp,&cp,max_shear);
			status = nc_put_vara_float(ncid,capeminTBid,&sp,&cp,cape_minTb32);
			status = nc_put_vara_float(ncid,shearminTBid,&sp,&cp,shear_minTb32);
			status = nc_put_vara_float(ncid,capeminpctid,&sp,&cp,cape_minpct);
			status = nc_put_vara_float(ncid,shearminpctid,&sp,&cp,shear_minpct);
			status = nc_put_vara_float(ncid,pctaveid,&sp,&cp,pct_ave);
			status = nc_put_vara_float(ncid,minPCTid,&sp,&cp,minpct);
			status = nc_put_vara_float(ncid,frac220id,&sp,&cp,Tb31_220frac);
			status = nc_put_vara_float(ncid,frac210id,&sp,&cp,Tb31_210frac);
			status = nc_put_vara_float(ncid,frac200id,&sp,&cp,Tb31_200frac);
			status = nc_put_vara_float(ncid,latminTBid,&sp,&cp,lat_minTb32);
			status = nc_put_vara_float(ncid,lonminTBid,&sp,&cp,lon_minTb32);
			status = nc_put_vara_float(ncid,latminPCTid,&sp,&cp,lat_minpct);
			status = nc_put_vara_float(ncid,lonminPCTid,&sp,&cp,lon_minpct);
			status = nc_put_vara_float(ncid,w850minTBid,&sp,&cp,wdir850_minTb32);
			status = nc_put_vara_float(ncid,w500minTBid,&sp,&cp,wdir500_minTb32);
			status = nc_put_vara_float(ncid,w200minTBid,&sp,&cp,wdir200_minTb32);
			status = nc_put_vara_float(ncid,IRminPCTid,&sp,&cp,Tb31_minpct);
			status = nc_put_vara_float(ncid,PCTminIRid,&sp,&cp,pct_minTb32);
			status = nc_put_vara_float(ncid,AODminIRid,&sp,&cp,aod_minTb32);
			status = nc_put_vara_float(ncid,AODmerraid,&sp,&cp,aod_merra_ave);
			status = nc_put_vara_float(ncid,AODmerraminIRid,&sp,&cp,aod_merra_minTb32);
			status = nc_put_vara_float(ncid,derbinid,spder,cpder,der_bin);
	//  	printf("status=%d\n",status);
			nc_close(ncid);

	//    fwrite(&viewed_area,sizeof(float),1,fp);
	//    fwrite(Tb89V_ave,sizeof(float),ncloud,fp);
	//    fwrite(Tb89H_ave,sizeof(float),ncloud,fp);
	//    fwrite(minTb89V,sizeof(float),ncloud,fp);
	//    fwrite(minTb89H,sizeof(float),ncloud,fp);
		} else if (status==-35) printf("Error: output file already exists\n");
		else printf("Error creating output file, status=%d\n",status);
	} //end if (cloud>0)


#ifdef FINE_GRID_FLAG
#ifdef VERBOSE

  sprintf(outfile,"%scloud_labels.%s.%s.nc",outPath,yyyyddd,hhmm);
	if ( (status = nc_create(outfile,NC_NOCLOBBER,&ncid)) == NC_NOERR ) {
		nc_def_dim(ncid,"longitude",ngridlons,&londimid);
		nc_def_dim(ncid,"latitude",ngridlats,&latdimid);
//			nc_def_var(ncid,"minutes_since_Jan_1_1990_UTC",NC_LONG,1,&tb32dimid,&timeid);
//		finegriddims[0] = londimid;
//		finegriddims[1] = latdimid;
		finegriddims[0] = latdimid;
		finegriddims[1] = londimid;
		nc_def_var(ncid,"longitude",NC_FLOAT,1,&londimid,&longridid);
		nc_def_var(ncid,"latitude",NC_FLOAT,1,&latdimid,&latgridid);
		nc_def_var(ncid,"cloud31",NC_LONG,2,finegriddims,&cld31id);
		nc_def_var(ncid,"cloud32",NC_LONG,2,finegriddims,&cld32id);
		nc_def_var(ncid,"CAPE",NC_FLOAT,2,finegriddims,&capgid);
		nc_def_var(ncid,"Tb_IR",NC_FLOAT,2,finegriddims,&irgid);
		nc_enddef(ncid);
//			status = nc_put_vara_long(ncid,timeid,);
		spder[0] = 0;
		spder[1] = 0;
//		cpder[0] = ngridlons;
//		cpder[1] = ngridlats;
		cpder[0] = ngridlats;
		cpder[1] = ngridlons;
//		status = nc_put_vara_float(ncid,longridid,spder,cpder,grid_lons);
//		status = nc_put_vara_float(ncid,latgridid,spder+1,cpder+1,grid_lats);
		status = nc_put_vara_float(ncid,longridid,spder+1,cpder+1,grid_lons);
		status = nc_put_vara_float(ncid,latgridid,spder,cpder,grid_lats);
		status = nc_put_vara_long(ncid,cld31id,spder,cpder,cloud31);
		status = nc_put_vara_long(ncid,cld32id,spder,cpder,cloud32);
		status = nc_put_vara_float(ncid,capgid,spder,cpder,cape_g);
		status = nc_put_vara_float(ncid,irgid,spder,cpder,Tb31_g);
		nc_close(ncid);
	}
/****************
  sprintf(outfile,"%sfine_stats.bin",outPath);
  fp = fopen(outfile,"a");
  if (fp != NULL) {
    fwrite(land_ocean,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(Tb31_g,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(Tb32_g,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(Tb89V_g,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(Tb89H_g,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(pct,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(cloud31,sizeof(long),ngridlons*ngridlats,fp);
    fwrite(cloud32,sizeof(long),ngridlons*ngridlats,fp);
    fwrite(cape_g,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(modzen_g,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(dcape_g,sizeof(float),ngridlons*ngridlats,fp);
    fwrite(shear_g,sizeof(float),ngridlons*ngridlats,fp);
//    fwrite(gamma,sizeof(float),ngridlons*ngridlats*N_DAO_LEVS,fp);
//    fwrite(Tvp,sizeof(float),ngridlons*ngridlats*N_DAO_LEVS,fp);
//    fwrite(Tve,sizeof(float),ngridlons*ngridlats*N_DAO_LEVS,fp);
    fclose(fp);
    fp = NULL;
  }
*****************/
#endif
#endif



#ifdef COARSE_GRID_FLAG
#ifndef MERRA_FLAG
  for (i=0; i<ncgridlats; i++)
  for (j=0; j<ncgridlons; j++) {
    iloc = (i*ncgridlons)+j;
    cape_cg[iloc] = -9999.;
    dcape_cg[iloc] = -9999.;
    shear_cg[iloc] = -9999.;
  }
#endif //end if not MERRA_FLAG
#ifdef VERBOSE
  sprintf(outfile,"%scoarse_stats.bin",outPath);
  fp = fopen(outfile,"a");
  if (fp != NULL) {
    fwrite(cape_cg,sizeof(float),ncgridlons*ncgridlats,fp);
    fwrite(dcape_cg,sizeof(float),ncgridlons*ncgridlats,fp);
    fwrite(shear_cg,sizeof(float),ncgridlons*ncgridlats,fp);
    fwrite(Tb240Kfrac_cg,sizeof(float),ncgridlons*ncgridlats,fp);
    fwrite(crf_lw_cg,sizeof(float),ncgridlons*ncgridlats,fp);
    fwrite(crf_sw_cg,sizeof(float),ncgridlons*ncgridlats,fp);
    fwrite(aod_cg,sizeof(float),ncgridlons*ncgridlats,fp);
    fclose(fp);
    fp = NULL;
  }
#endif
#endif

  FREE(grid_lons);
  FREE(grid_lats);
  FREE(areas);
  FREE(Tb31_g);
  FREE(Tb32_g);
  FREE(cloud31);
  FREE(cloud32);
  FREE(cloud);
  FREE(class31);
  FREE(class32);
  FREE(cloud31_areas);
  FREE(cloud32_areas);
  FREE(cloud31_class);
  FREE(cloud32_class);
  FREE(Tb89H_g);
  FREE(Tb89V_g);
  FREE(minTb32);
  FREE(Tb89V_ave);
  FREE(Tb89H_ave);
  FREE(minTb89V);
  FREE(minTb89H);
  FREE(land_ocean);
  FREE(landfrac);
  FREE(pct);
  FREE(pct_ave);
  FREE(minpct);
  FREE(lat_minTb32);
  FREE(lon_minTb32);
  FREE(lat_minpct);
  FREE(lon_minpct);
  FREE(Tb31_minpct);
  FREE(pct_minTb32);
  FREE(Tb31_220frac);
  FREE(Tb31_210frac);
  FREE(Tb31_200frac);
  FREE(cloud_datetime);
  									
	FREE(der_bin);
	FREE(bt_bin);

  FREE(cape_ave);
  FREE(dcape_ave);
  FREE(max_cape);
  FREE(max_dcape);
  FREE(shear_ave);
  FREE(max_shear);
  FREE(cape_minTb32);
  FREE(shear_minTb32);
  FREE(cape_minpct);
  FREE(shear_minpct);
  FREE(wdir850_minTb32);
  FREE(wdir500_minTb32);
  FREE(wdir200_minTb32);
  FREE(aod_minTb32);
  FREE(aod_merra_minTb32);
  FREE(aod_merra_ave);
  FREE(gamma);
  FREE(Tvp);
  FREE(Tve);
  FREE(modzen_g);
#ifdef COARSE_GRID_FLAG
	FREE(cgrid_lons);
	FREE(cgrid_lats);
	FREE(areas_cgrid);
  FREE(cape_cg);
  FREE(dcape_cg);
  FREE(shear_cg);
  FREE(Tb240Kfrac_cg);
  FREE(crf_lw_cg);
  FREE(crf_sw_cg);
  FREE(aod_cg);
#endif //end if COARSE_GRID_FLAG
#ifdef MERRA_FLAG
  FREE(cape_g);
  FREE(dcape_g);
  FREE(shear_g);
  FREE(wdir850);
  FREE(wdir500);
  FREE(wdir200);
	FREE(aod_merra_g);
#endif //end if MERRA_FLAG

}
