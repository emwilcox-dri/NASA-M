#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "time_tools.h"
#include "amsr.h"


void read_L2Ocean(char *filename,float *lat,float *lon,long *na_amsr,
                  long *nx_amsr,float *clw)
{
  char sds_name[128];
  char fname[150];
  short *short_tmp=NULL;
  int sd_id, sds_id, sds_index, attr_index;
  int i, status, start[3], edge[3];
  int rank, dims[3], type, attributes;
  long ix, ia, iloc;
  float scale;

  for (i=0; i<3; i++) start[i] = edge[i] = dims[i] = 0;

  printf("reading %s\n",filename);
  sd_id = SDstart(filename, DFACC_READ);

  sprintf(sds_name,"Longitude");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  for (i=0; i<rank; i++) {
    start[i] = 0;
    edge[i] = dims[i];
  }
  *na_amsr = dims[0];
  *nx_amsr = dims[1];
  status = SDreaddata(sds_id,start,NULL,edge,lon);
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Latitude");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,lat);
  status = SDendaccess(sds_id);

  sprintf(sds_name,"High_res_cloud");
  short_tmp = (signed short *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR*sizeof(signed short));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,short_tmp);
  attr_index = SDfindattr(sds_id,"Scale");
  status = SDreadattr(sds_id,attr_index,&scale);
  status = SDendaccess(sds_id);

  for (ix=0; ix<N_ACROSS_AMSR; ix++)
    for (ia=0; ia<(*na_amsr); ia++) {
      iloc = (ia*N_ACROSS_AMSR)+ix;
      clw[iloc] = ((float) short_tmp[iloc])*scale;
  }
  FREE(short_tmp);

  status = SDend(sd_id);
}



void read_L2Ocean_QA(char *filename,unsigned char *qa)
{
  char sds_name[128];
  unsigned char *qatmp;
  int i, sd_id, sds_id, sds_index;
  int status, start[3], edge[3], rank, dims[3], type, attributes;
  long na, nx, ix, ia, iloc, qaloc, nbytes;

  for (i=0; i<3; i++) start[i] = edge[i] = dims[i] = 0;

  sd_id = SDstart(filename, DFACC_READ);

  sprintf(sds_name,"Ocean_products_quality_flag");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  for (i=0; i<rank; i++) {
    start[i] = 0;
    edge[i] = dims[i];
  }
  na = dims[0];
  nx = dims[1];
  qatmp = (unsigned char *) calloc((N_ACROSS_AMSR*N_ALONG_AMSR*
          N_QAFLAG_BYTES_AMSR),sizeof(unsigned char));
  status = SDreaddata(sds_id,start,NULL,edge,qatmp);
  status = SDendaccess(sds_id);

  status = SDend(sd_id);

  for (ix=0; ix<nx; ix++)
    for (ia=0; ia<na; ia++) {
      iloc = (ia*nx*N_QAFLAG_BYTES_AMSR)+
               (ix*N_QAFLAG_BYTES_AMSR);
      qaloc = (ia*nx*N_QAFLAG_TYPES_AMSR)+
               (ix*N_QAFLAG_TYPES_AMSR);
      qa[qaloc] = (qatmp[iloc] & 3);
      qaloc++;
      qa[qaloc] = (qatmp[iloc] & 4) >> 2;                     
      qa[qaloc] += (qatmp[iloc] & 8) >> 3;                    
      qaloc++;                                                
      qa[qaloc] = (qatmp[iloc] & 16) >> 4;                    
      qaloc++;                                                
      qa[qaloc] = (qatmp[iloc] & 32) >> 5;                    
      qaloc++;                                                
      qa[qaloc] = (qatmp[iloc] & 64) >> 6;                    
      qaloc++;                                                
      qa[qaloc] = (qatmp[iloc] & 128) >> 7;                   
      qaloc++;
      iloc++;                                                 
      qa[qaloc] = (qatmp[iloc] & 3);                          
      qaloc++;                                                
      qa[qaloc] = (qatmp[iloc] & 12) >> 2;                    
      qaloc++;                                                
      qa[qaloc] = (qatmp[iloc] & 48) >> 4;                    
      qaloc++;                                                
      qa[qaloc] = (qatmp[iloc] & 192) >> 6;                   
  }                                                           
  free(qatmp);                                                
                                                              
}



double getAMSR_fileName(char *amsrPath,int year,int mon,int day,int hour,
                      int min,char *amsrFile,int ichar_date)
{
  FILE *fp=NULL;
  char filelist[150];
  char tmpfile[60];
  char yyyy[5],mm[3],dd[3],hh[3];
  char *noneStr = "none";
  int endflag=0;
  double mod_jdate, amsr_jdate, del_t, tmp;
  double delt_out=-1;

  mod_jdate = amsr_jdate = del_t = tmp = 0.;

// each file has about 50 minutes of data
  mdyhms_to_juldate(1,1,1999,10,0,0.,&del_t);
//  mdyhms_to_juldate(1,1,1999,9,15,0,&tmp);
  mdyhms_to_juldate(1,1,1999,9,5,0.,&tmp);
  del_t -= tmp;
//printf("del_t=%lf\n",del_t);
  strcpy(amsrFile,noneStr);

  mdyhms_to_juldate(mon,day,year,hour,min,0.,&mod_jdate);
  sprintf(filelist,"%sAMSRlist",amsrPath);
//printf("%s\n",filelist);
//fflush(NULL);
  fp = fopen(filelist,"r");
  while (!endflag) {
    if ( fgets(tmpfile,60,fp) != NULL ) {
//printf("%s\n",tmpfile);
      strncpy(yyyy,tmpfile+ichar_date,4);
      yyyy[4] = '\0';
      year = atoi(yyyy);
      strncpy(mm,tmpfile+ichar_date+4,2);
      mm[2] = '\0';
      mon = atoi(mm);
      strncpy(dd,tmpfile+ichar_date+6,2);
      dd[2] = '\0';
      day = atoi(dd);
      strncpy(hh,tmpfile+ichar_date+8,2);
      hh[2] = '\0';
      hour = atoi(hh);
      strncpy(mm,tmpfile+ichar_date+10,2);
      mm[2] = '\0';
      min = atoi(mm);
      mdyhms_to_juldate(mon,day,year,hour,min,0.,&amsr_jdate);
      if ( (mod_jdate > amsr_jdate) &&
           ((mod_jdate-amsr_jdate) < del_t) ) {
//        strncpy(amsrFile,tmpfile,38);
        sprintf(amsrFile,"%s",tmpfile);
        delt_out = mod_jdate-amsr_jdate;
      }
    } else endflag = 1;
  }
  fclose(fp);
  return delt_out;
}






/**********************************
header from read_amsr_day_v5.pro
;output products (lon,lat,asc/dsc)
time =fltarr(1440,720,2)
sst  =fltarr(1440,720,2)
wind =fltarr(1440,720,2)
vapor=fltarr(1440,720,2)
cloud=fltarr(1440,720,2)
rain =fltarr(1440,720,2)

 The routine returns:
   time, sst, wind, vapor, cloud, rain real arrays sized (1440,720,2)
   time  is the mean gmt time in fractional hours of the observations within that grid cell
   sst   is the sea surface temperature in degree Celcius, valid range=[-3.0,34.5]
   wind  is the 10 meter surface wind speed in m/s,  valid range=[0.,50.]
   vapor is the columnar atmospheric water vapor in mm,  valid range=[0.,75.]
   cloud is the liquid cloud water in mm, valid range = [0.,2.5]
   rain  is the derived radiometer rain rate in mm/hr,  valid range = [0.,24.95]

 Longitude  is 0.25*(xdim+1)-0.125    East longitude
 Latitude   is 0.25*(ydim+1)-90.125



v5 multipliers to change binary data to real data
xscale=[0.1,0.15,.2,.3,.01,.1]
offset=[0,-3,0,0,0,0]
**********************************/
void read_v5_daily_gridded(char *filename,float *lat,float *lon,
                           float *time,float *lwp,float *rain,
                           float *sst)
{
  FILE *fp;
//  char tmpvar[N_DAILY_PASS_AMSR][N_DAILY_VARS_AMSR][N_DAILY_LATS_AMSR][N_DAILY_LONS_AMSR];
//  char tmp;
  char *tmpvar;
  unsigned char tmp;
  int ilo,ila,iasc;
  long iloc, kloc;
  float lwp_scale=0.01, lwp_offset=0.;
  float time_scale=0.1, time_offset=0.;
  float rain_scale=0.1, rain_offset=0.;
  float sst_scale=0.15, sst_offset=-3.;

  tmpvar = (char *) malloc(N_DAILY_PASS_AMSR*N_DAILY_VARS_V5_AMSR*
              N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR*sizeof(char));
  fp = fopen(filename,"r");
  if (fp != NULL) {
    fread(tmpvar,sizeof(char),N_DAILY_LONS_AMSR*N_DAILY_LATS_AMSR*
      N_DAILY_VARS_V5_AMSR*N_DAILY_PASS_AMSR,fp);
    fclose(fp);
  }

  for (ilo=0; ilo<N_DAILY_LONS_AMSR; ilo++)
  for (ila=0; ila<N_DAILY_LATS_AMSR; ila++)
  for (iasc=0; iasc<N_DAILY_PASS_AMSR; iasc++) {

    iloc = (ilo*N_DAILY_LATS_AMSR*N_DAILY_PASS_AMSR) +
           (ila*N_DAILY_PASS_AMSR) + iasc;

    kloc = (iasc*N_DAILY_VARS_V5_AMSR*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (4*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) + 
           (ila*N_DAILY_LONS_AMSR) + ilo;
//    tmp = tmpvar[iasc][4][ila][ilo];
    tmp = tmpvar[kloc];
    if (tmp<=250) lwp[iloc] = ((float) tmp)*lwp_scale + lwp_offset;
    else lwp[iloc] = AMSR_BADFLAG;

    kloc = (iasc*N_DAILY_VARS_V5_AMSR*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (0*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) + 
           (ila*N_DAILY_LONS_AMSR) + ilo;
//    tmp = tmpvar[iasc][0][ila][ilo];
    tmp = tmpvar[kloc];
    if (tmp<=250) time[iloc] =((float) tmp)*time_scale + time_offset;
    else time[iloc] = AMSR_BADFLAG;

    kloc = (iasc*N_DAILY_VARS_V5_AMSR*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (5*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (ila*N_DAILY_LONS_AMSR) + ilo;
//    tmp = tmpvar[iasc][5][ila][ilo];
    tmp = tmpvar[kloc];
    if (tmp<=250) rain[iloc] = ((float) tmp)*rain_scale + rain_offset;
    else rain[iloc] = AMSR_BADFLAG;

    kloc = (iasc*N_DAILY_VARS_V5_AMSR*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (1*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (ila*N_DAILY_LONS_AMSR) + ilo;
//    tmp = tmpvar[iasc][1][ila][ilo];
    tmp = tmpvar[kloc];
    if (tmp<=250) sst[iloc] = ((float) tmp)*sst_scale + sst_offset;
    else sst[iloc] = AMSR_BADFLAG;
  }

  for (ilo=0; ilo<N_DAILY_LONS_AMSR; ilo++) lon[ilo] = 0.25*(ilo+1)-0.125;
  for (ila=0; ila<N_DAILY_LATS_AMSR; ila++) lat[ila] = 0.25*(ila+1)-90.125;
}





/*********************************
ver is version: 5 or 7

v5 multipliers to change binary data to real data
xscale=[0.1,0.15,.2,.3,.01,.1]
offset=[0,-3,0,0,0,0]

v7 multipliers
xscale=[0.1,0.15,0.2,0.2,0.3,0.01,0.1]
offset=[0,-3,0,0,0,0.05,0]
**********************************/
void read_daily_gridded(char *filename,float *lat,float *lon,
                           float *time,float *lwp,float *rain,
                           float *sst,float *cwv,int ver)
{
  FILE *fp;
//  char tmpvar[N_DAILY_PASS_AMSR][N_DAILY_VARS_AMSR][N_DAILY_LATS_AMSR][N_DAILY_LONS_AMSR];
//  char tmp;
  char *tmpvar;
  unsigned char tmp;
  int ilo,ila,iasc,ilwp,irain, itime=0, isst=1,icwv;
  long iloc, kloc, nvars;
  float lwp_scale, lwp_offset;
  float time_scale, time_offset;
  float rain_scale, rain_offset;
  float sst_scale, sst_offset;
  float cwv_scale, cwv_offset;
  float v5_scale[6] = {0.1, 0.15, .2, .3, .01, .1};
  float v5_offset[6] = {0, -3, 0, 0, 0, 0};
  float v7_scale[7] = {0.1, 0.15, 0.2, 0.2, 0.3, 0.01, 0.1};
  float v7_offset[7] = {0, -3, 0, 0, 0, -0.05, 0};

	if (ver==5) {
	  icwv=3;
		ilwp=4;
		irain=5;
		lwp_scale = v5_scale[ilwp];
		lwp_offset = v5_offset[ilwp];
		rain_scale = v5_scale[irain];
		rain_offset = v5_offset[irain];
		time_scale = v5_scale[itime];
		time_offset = v5_offset[itime];
		sst_scale = v5_scale[isst];
		sst_offset = v5_offset[isst];
		nvars = N_DAILY_VARS_V5_AMSR;
	} else {
	  icwv=4;
		ilwp=5;
		irain=6;
		lwp_scale = v7_scale[ilwp];
		lwp_offset = v7_offset[ilwp];
		rain_scale = v7_scale[irain];
		rain_offset = v7_offset[irain];
		time_scale = v7_scale[itime];
		time_offset = v7_offset[itime];
		sst_scale = v7_scale[isst];
		sst_offset = v7_offset[isst];
		cwv_scale = v7_scale[icwv];
		cwv_offset = v7_offset[icwv];
		nvars = N_DAILY_VARS_V7_AMSR;
	}
 
	tmpvar = (char *) malloc(N_DAILY_PASS_AMSR*nvars*
              N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR*sizeof(char));

  fp = fopen(filename,"r");
  if (fp != NULL) {
    fread(tmpvar,sizeof(char),N_DAILY_LONS_AMSR*N_DAILY_LATS_AMSR*
      nvars*N_DAILY_PASS_AMSR,fp);
    fclose(fp);
  }

  for (ilo=0; ilo<N_DAILY_LONS_AMSR; ilo++)
  for (ila=0; ila<N_DAILY_LATS_AMSR; ila++)
  for (iasc=0; iasc<N_DAILY_PASS_AMSR; iasc++) {

    iloc = (ilo*N_DAILY_LATS_AMSR*N_DAILY_PASS_AMSR) +
           (ila*N_DAILY_PASS_AMSR) + iasc;

    kloc = (iasc*nvars*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (ilwp*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) + 
           (ila*N_DAILY_LONS_AMSR) + ilo;
    tmp = tmpvar[kloc];
    if (tmp<=250) lwp[iloc] = ((float) tmp)*lwp_scale + lwp_offset;
    else lwp[iloc] = AMSR_BADFLAG;

    kloc = (iasc*nvars*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (itime*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) + 
           (ila*N_DAILY_LONS_AMSR) + ilo;
    tmp = tmpvar[kloc];
    if (tmp<=250) time[iloc] =((float) tmp)*time_scale + time_offset;
    else time[iloc] = AMSR_BADFLAG;

    kloc = (iasc*nvars*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (irain*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (ila*N_DAILY_LONS_AMSR) + ilo;
    tmp = tmpvar[kloc];
    if (tmp<=250) rain[iloc] = ((float) tmp)*rain_scale + rain_offset;
    else rain[iloc] = AMSR_BADFLAG;

    kloc = (iasc*nvars*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (isst*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (ila*N_DAILY_LONS_AMSR) + ilo;
    tmp = tmpvar[kloc];
    if (tmp<=250) sst[iloc] = ((float) tmp)*sst_scale + sst_offset;
    else sst[iloc] = AMSR_BADFLAG;

    kloc = (iasc*nvars*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (icwv*N_DAILY_LATS_AMSR*N_DAILY_LONS_AMSR) +
           (ila*N_DAILY_LONS_AMSR) + ilo;
    tmp = tmpvar[kloc];
    if (tmp<=250) cwv[iloc] = ((float) tmp)*cwv_scale + cwv_offset;
    else cwv[iloc] = AMSR_BADFLAG;
  }
  FREE(tmpvar);

  for (ilo=0; ilo<N_DAILY_LONS_AMSR; ilo++) lon[ilo] = 0.25*(ilo+1)-0.125;
  for (ila=0; ila<N_DAILY_LATS_AMSR; ila++) lat[ila] = 0.25*(ila+1)-90.125;
}






void read_AMSR_pixel_tbs(char *filename,long *na_amsr,long *nx_amsr,
                    float *lon_a,float *lat_a,float *lon_b,
                    float *lat_b,float *tb89_Ha,float *tb89_Va,
                    float *tb89_Hb,float *tb89_Vb,float *landfrac_a,
                    float *landfrac_b,int version)
{
  char sds_name[128];
  unsigned char *uchartmp=NULL;
  short *shrttmp=NULL;
  int sd_id, sds_id, sds_index, attr_index;
  int i, status, start[2], edge[2], alatid, alonid, blatid, blonid;
  int rank, dims[2], type, attributes;
  long ix, ia, iloc;
  float *flttmp=NULL, scale, offset;
  
  if (version==9) {
    alatid = 2;
    alonid = 3;
    blatid = 4;
    blonid = 5;
  }
  else if ( (version==10) || (version==11) ) {
    alatid = 68;
    alonid = 69;
    blatid = 80;
    blonid = 81;
  }
  else if ( (version==12) || (version==13) ) {
    alatid = 69;
    alonid = 70;
    blatid = 81;
    blonid = 82;
  } else printf("ERROR: need to specify HDF Variable indices for AMSR-E high-res lat/lons.");

  printf("reading %s\n",filename);
  sd_id = SDstart(filename, DFACC_READ);

  sds_id = SDselect(sd_id,alonid);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
//  printf("%d %d\n",sds_id,status);
  for (i=0; i<rank; i++) {
    start[i] = 0;
    edge[i] = dims[i];
  }
  *na_amsr = dims[0];
  *nx_amsr = dims[1];
  status = SDreaddata(sds_id,start,NULL,edge,lon_a);
  status = SDendaccess(sds_id);
  sds_id = SDselect(sd_id,alatid);
  status = SDreaddata(sds_id,start,NULL,edge,lat_a);
  status = SDendaccess(sds_id);
  sds_id = SDselect(sd_id,blonid);
  status = SDreaddata(sds_id,start,NULL,edge,lon_b);
  status = SDendaccess(sds_id);
  sds_id = SDselect(sd_id,blatid);
  status = SDreaddata(sds_id,start,NULL,edge,lat_b);

  sprintf(sds_name,"89.0V_Res.5A_TB_(not-resampled)");
  shrttmp = (short *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*sizeof(short));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,shrttmp);
  attr_index = SDfindattr(sds_id,"SCALE FACTOR");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"OFFSET");
  status = SDreadattr(sds_id,attr_index,&offset);
  for (ix=0; ix<(*nx_amsr); ix++)
    for (ia=0; ia<(*na_amsr); ia++) {
      iloc = (ia*(*nx_amsr))+ix;
      tb89_Va[iloc] = ((float) shrttmp[iloc]*scale)+offset;
  }
  status = SDendaccess(sds_id);

  sprintf(sds_name,"89.0H_Res.5A_TB_(not-resampled)");  
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,shrttmp);
  attr_index = SDfindattr(sds_id,"SCALE FACTOR");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"OFFSET");
  status = SDreadattr(sds_id,attr_index,&offset);
  for (ix=0; ix<(*nx_amsr); ix++)
    for (ia=0; ia<(*na_amsr); ia++) {
      iloc = (ia*(*nx_amsr))+ix;
      tb89_Ha[iloc] = ((float) shrttmp[iloc]*scale)+offset;
  }
  status = SDendaccess(sds_id);

  sprintf(sds_name,"89.0V_Res.5B_TB_(not-resampled)");
  flttmp = (float *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*sizeof(float));
  sds_index = SDnametoindex(sd_id,sds_name);
//printf("%d ",sds_index);
  sds_id = SDselect(sd_id,sds_index);
//printf("%d ",sds_id);
  status = SDreaddata(sds_id,start,NULL,edge,shrttmp);
//printf("%d ",status);
  attr_index = SDfindattr(sds_id,"SCALE FACTOR");
//printf("%d ",attr_index);
  status = SDreadattr(sds_id,attr_index,&scale);
//printf("%d ",status);
  attr_index = SDfindattr(sds_id,"OFFSET");
//printf("%d ",attr_index);
  status = SDreadattr(sds_id,attr_index,&offset);
//printf("%d ",status);
  for (ix=0; ix<(*nx_amsr); ix++)
    for (ia=0; ia<(*na_amsr); ia++) {
      iloc = (ia*(*nx_amsr))+ix;
      tb89_Vb[iloc] = ((float) shrttmp[iloc]*scale)+offset;
//printf("%f ",tb89_Vb[iloc]);
  }
  status = SDendaccess(sds_id);
//printf("%d ",status);

  sprintf(sds_name,"89.0H_Res.5B_TB_(not-resampled)");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,shrttmp);
  attr_index = SDfindattr(sds_id,"SCALE FACTOR");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"OFFSET");
  status = SDreadattr(sds_id,attr_index,&offset);
  for (ix=0; ix<(*nx_amsr); ix++)
    for (ia=0; ia<(*na_amsr); ia++) {
      iloc = (ia*(*nx_amsr))+ix;
      tb89_Hb[iloc] = ((float) shrttmp[iloc]*scale)+offset;
  }
  status = SDendaccess(sds_id);

  FREE(shrttmp);

  uchartmp = (unsigned char *) malloc(N_ALONG_AMSR*N_ACROSS_AMSR_5KM*
                                    sizeof(unsigned char));
  sprintf(sds_name,"Res5A_Surf");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,uchartmp);
  attr_index = SDfindattr(sds_id,"SCALE FACTOR");
  status = SDreadattr(sds_id,attr_index,&scale);
  for (ix=0; ix<(*nx_amsr); ix++)
    for (ia=0; ia<(*na_amsr); ia++) {
      iloc = (ia*(*nx_amsr))+ix;
      landfrac_a[iloc] = (float) uchartmp[iloc]*scale;
  }
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Res5B_Surf");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,uchartmp);
  attr_index = SDfindattr(sds_id,"SCALE FACTOR");
  status = SDreadattr(sds_id,attr_index,&scale);
  for (ix=0; ix<(*nx_amsr); ix++)
    for (ia=0; ia<(*na_amsr); ia++) {
      iloc = (ia*(*nx_amsr))+ix;
      landfrac_b[iloc] = (float) uchartmp[iloc]*scale;
  }
  status = SDendaccess(sds_id);

  FREE(uchartmp);

  status = SDend(sd_id);
//printf("%d\n",status);
  
//  printf("done.\n");
/*
for (ix=0; ix<(*nx_amsr); ix++)
for (ia=0; ia<(*na_amsr); ia++) {
iloc = (ia*(*nx_amsr))+ix;
printf("%f ",tb89_Vb[iloc]);
}
*/
  fflush(NULL);

}


