#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "time_tools.h"
#include "grid.h"
#include "MOD06.h"


const char *reff_name_C004 = "Effective_Particle_Radius";
const char *reff_name_C005 = "Cloud_Effective_Radius";
const char *reff_1621_name_C005 = "Cloud_Effective_Radius_1621";
const char *lwp_name_C004 = "Water_Path";
const char *lwp_name_C005 = "Cloud_Water_Path";
const char *lwp_1621_name_C005 = "Cloud_Water_Path_1621";

void readMOD06(char *modPath,char *mod03path,char *MOD06file,char *MOD03file,
               int collection,long *na_1km,long *na_5km,
               float *lat_5km,float *lon_5km,float *lat_1km,
               float *lon_1km,float *fraction,float *reff,
               float *cld_tau,float *cldtop_temp,float *cldtop_pres,
               signed char *cldtop_phase,float *lwp,unsigned char *land_mask)
{
  char sds_name[128];
  char fname[150];
  signed char *tmpfrac, frac_fill;
  short *tmptemp, *tmppres, *tmpreff, *tmptau, *tmplwp;
  signed short reff_fill, tmp_fill, phs_fill, tau_fill;
  signed short lwp_fill, pres_fill;
  int sd_id, sds_id, sds_index, attr_index;
  int status, start[2], edge[2];
  int rank, dims[2], type, attributes;
  long ix, ia, iloc;
  double reff_scale, frac_scale, tmp_offset, tmp_scale;
  double phs_scale, tau_scale, pres_scale;

  sprintf(fname,"%s%s",modPath,MOD06file);
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  start[0] = start[1] = 0;
  edge[1] = N_ACROSS_5KM;

  sprintf(sds_name,"Latitude");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  if ( (dims[0]>N_ALONG_5KM_MAX) || (dims[1]>N_ACROSS_5KM) ) {
    printf("MODIS 5km dimensions = [%d %d] - too large\n",dims[0],dims[1]);
    status = SDendaccess(sds_id);
    status = SDend(sd_id);
    exit(0);
  }
  *na_5km = dims[0];
  edge[0] = *na_5km;
  status = SDreaddata(sds_id,start,NULL,edge,lat_5km);
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Longitude");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,lon_5km);
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Cloud_Fraction");
  tmpfrac = (signed char *) calloc((N_ACROSS_5KM*N_ALONG_5KM_MAX),
            sizeof(signed char));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmpfrac);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&frac_scale);
  status = SDgetfillvalue(sds_id,&frac_fill);
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Cloud_Top_Temperature");
  tmptemp = (short *) calloc((N_ACROSS_5KM*N_ALONG_5KM_MAX),sizeof(short));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmptemp);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&tmp_scale);
  status = SDgetfillvalue(sds_id,&tmp_fill);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&tmp_offset);
  status = SDendaccess(sds_id);

//  sprintf(sds_name,"Cloud_Top_Pressure_From_Ratios"); //note this is rank=3
  sprintf(sds_name,"Cloud_Top_Pressure_Infrared");
  tmppres = (short *) calloc((N_ACROSS_5KM*N_ALONG_5KM_MAX),sizeof(short));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmppres);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&pres_scale);
  status = SDgetfillvalue(sds_id,&tmp_fill);
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Cloud_Phase_Infrared");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,cldtop_phase);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&phs_scale);
  status = SDgetfillvalue(sds_id,&phs_fill);
  status = SDendaccess(sds_id);

  for (ix=0; ix<N_ACROSS_5KM; ix++)                                
    for (ia=0; ia<(*na_5km); ia++) {
      iloc = (ia*N_ACROSS_5KM)+ix;
      if (tmpfrac[iloc] != frac_fill) fraction[iloc] = 
	      ((float) tmpfrac[iloc]) * ((float) frac_scale);
      else fraction[iloc] = MOD06_BADFLAG;
      if (tmptemp[iloc] != tmp_fill)
	cldtop_temp[iloc] = ( ((float) tmptemp[iloc]) - tmp_offset )*tmp_scale;
      else cldtop_temp[iloc] = MOD06_BADFLAG;
      if (tmppres[iloc] != pres_fill)
        cldtop_pres[iloc] = ((float) tmppres[iloc])*pres_scale;
      else cldtop_pres[iloc] = MOD06_BADFLAG;
  }
  FREE(tmpfrac);
  FREE(tmptemp);
  FREE(tmppres);

  edge[1] = N_ACROSS_1KM;
/***
  if (collection==4) sprintf(sds_name,reff_name_C004);
  if (collection==5) sprintf(sds_name,reff_name_C005);
***/
  if (collection==4) strcpy(sds_name,reff_name_C004);
//  if ( (collection==5) || (collection==51) ) strcpy(sds_name,reff_name_C005);
  else strcpy(sds_name,reff_name_C005);
//printf("collection = %d\n",collection);
//printf("%s ",sds_name);
  tmpreff = (short *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX),sizeof(short));
  sds_index = SDnametoindex(sd_id,sds_name);
//printf("status %d ",sds_index);
  sds_id = SDselect(sd_id,sds_index);
//printf("%d ",sds_id);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
//printf("%d ",status);
  if ( (dims[0]>N_ALONG_1KM_MAX) || (dims[1]>N_ACROSS_1KM) ) {
    printf("MODIS 1km dimensions = [%d %d] - too large\n",dims[0],dims[1]);
    status = SDendaccess(sds_id);
    status = SDend(sd_id);
    exit(0);
  }
  *na_1km = dims[0];
  edge[0] = *na_1km;
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&reff_scale);
  status = SDgetfillvalue(sds_id,&reff_fill);
  status = SDreaddata(sds_id,start,NULL,edge,tmpreff);
//printf("%d\n",status);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<(*na_1km); ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
//      printf("%d ",tmpreff[iloc]);
/*      iloc = (ix*(*na_1km))+ia; */
      if (tmpreff[iloc] != reff_fill) reff[iloc] =
	            ((float) tmpreff[iloc])*reff_scale;
      else reff[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);
  FREE(tmpreff);

  sprintf(sds_name,"Cloud_Optical_Thickness");
  tmptau = (short *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX),sizeof(short));
  sds_index = SDnametoindex(sd_id,sds_name);
//printf("%s ",sds_name);
//printf("status %d ",sds_index);
  sds_id = SDselect(sd_id,sds_index);
//printf("%d ",sds_id);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
//printf("%d ",status);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&tau_scale);
  status = SDgetfillvalue(sds_id,&tau_fill);
  status = SDreaddata(sds_id,start,NULL,edge,tmptau);
//printf("%d\n",status);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<(*na_1km); ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      if (tmptau[iloc] != tau_fill) cld_tau[iloc] =
                    ((float) tmptau[iloc])*tau_scale;
      else cld_tau[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);
  FREE(tmptau);

/***
  if (collection==4) sprintf(sds_name,lwp_name_C004);
  if (collection==5) sprintf(sds_name,lwp_name_C005);
***/
  if (collection==4) strcpy(sds_name,lwp_name_C004);
  if (collection==5) strcpy(sds_name,lwp_name_C005);
  tmplwp = (short *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX),sizeof(short));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  status = SDgetfillvalue(sds_id,&lwp_fill);
  status = SDreaddata(sds_id,start,NULL,edge,tmplwp);
//printf("%d ",status);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<(*na_1km); ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      if (tmplwp[iloc] != lwp_fill) lwp[iloc] =
                    ((float) tmplwp[iloc]);
      else lwp[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);
  FREE(tmplwp);


  readMOD03(mod03path,MOD03file,na_1km,lat_1km,lon_1km,land_mask);

  status = SDend(sd_id);

}




void read_1621_and_reff_diff(char *modPath,char *MOD06file,long *na_1km,
                             float *reff_1621,float *tau_1621,
                             float *lwp_1621,float *reff_diff)
{
  char sds_name[128];
  char fname[150];
  short *tmp;
  signed short fill;
  int sd_id, sds_id, sds_index, attr_index;
  int status, start[3], edge[3];
  int rank, dims[3], type, attributes;
  long ix, ia, iloc;
  double scale;

  sprintf(fname,"%s%s",modPath,MOD06file);
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  start[0] = start[1] = start[2] = 0;
  edge[1] = N_ACROSS_1KM;

//  sprintf(sds_name,reff_1621_name_C005);
  strcpy(sds_name,reff_1621_name_C005);
  tmp = (short *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX*2),sizeof(short));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  if ( (dims[0]>N_ALONG_1KM_MAX) || (dims[1]>N_ACROSS_1KM) ) {
    printf("MODIS 1km dimensions = [%d %d] - too large\n",dims[0],dims[1]);
    status = SDendaccess(sds_id);
    status = SDend(sd_id);
    exit(0);
  }
  *na_1km = dims[0];
  edge[0] = *na_1km;
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  status = SDgetfillvalue(sds_id,&fill);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<(*na_1km); ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      if (tmp[iloc] != fill) reff_1621[iloc] =
                    ((float) tmp[iloc])*scale;
      else reff_1621[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Cloud_Optical_Thickness_1621");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  status = SDgetfillvalue(sds_id,&fill);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<(*na_1km); ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      if (tmp[iloc] != fill) tau_1621[iloc] =
                    ((float) tmp[iloc])*scale;
      else tau_1621[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);

//  sprintf(sds_name,lwp_1621_name_C005);
  strcpy(sds_name,lwp_1621_name_C005);
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  status = SDgetfillvalue(sds_id,&fill);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<(*na_1km); ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      if (tmp[iloc] != fill) lwp_1621[iloc] =
                    ((float) tmp[iloc]);
      else lwp_1621[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);

  edge[0]=2;
  edge[1]=*na_1km;
  edge[2]=N_ACROSS_1KM;
  sprintf(sds_name,"Effective_Radius_Difference");
    sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  status = SDgetfillvalue(sds_id,&fill);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<(*na_1km); ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      if (tmp[iloc] != fill) reff_diff[iloc] =
                    ((float) tmp[iloc]);
      else reff_diff[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);
  FREE(tmp);

  status = SDend(sd_id);
}






void readMOD03(char *modPath,char *MOD03file,long *na_1km,
               float *lat_1km,float *lon_1km,unsigned char *land_mask)
{
  char sds_name[128];
  char fname[150];
  int sd_id, sds_id, sds_index;
  int status, start[2], edge[2];
  int rank, dims[2], type, attributes;

  sprintf(fname,"%s%s",modPath,MOD03file);
/*  strcpy(fname,modPath);
  strncat(fname,MOD03file,NCHAR_FNAME_MOD03); */
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  start[0] = start[1] = 0;
  edge[1] = N_ACROSS_1KM;

  sprintf(sds_name,"Latitude");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
//  printf("%d %d %d %d\n",sd_id,sds_index,sds_id,status);
  *na_1km = dims[0];
  edge[0] = *na_1km;
  status = SDreaddata(sds_id,start,NULL,edge,lat_1km);
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Longitude");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,lon_1km);
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Land/SeaMask");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,land_mask);
  status = SDendaccess(sds_id);

  status = SDend(sd_id);
}





void readMOD03_sensor_angles(char *modPath,char *MOD03file,
               float *sensor_zenith,float *sensor_azimuth)
{
  char sds_name[128];
  char fname[150];
  short *tmpvar=NULL;
  signed short fill;
  int sd_id, sds_id, sds_index, attr_index;
  int status, start[2], edge[2];
  int rank, dims[2], type, attributes;
  long ix, ia, iloc, na_1km;
  double scale;

  strcpy(fname,modPath);
  strncat(fname,MOD03file,NCHAR_FNAME_MOD03);
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  start[0] = start[1] = 0;
  edge[1] = N_ACROSS_1KM;

  tmpvar = (short *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX),sizeof(short));
  sprintf(sds_name,"SensorZenith");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
  na_1km = dims[0];
  edge[0] = na_1km;
  status = SDreaddata(sds_id,start,NULL,edge,tmpvar);
  status = SDgetfillvalue(sds_id,&fill);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<na_1km; ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      if (tmpvar[iloc] != fill) sensor_zenith[iloc] =
                    ((float) tmpvar[iloc])*scale;
      else sensor_zenith[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);

  sprintf(sds_name,"SensorAzimuth");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmpvar);
  status = SDgetfillvalue(sds_id,&fill);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<na_1km; ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      if (tmpvar[iloc] != fill) sensor_azimuth[iloc] =
                    ((float) tmpvar[iloc])*scale;
      else sensor_zenith[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);

  FREE(tmpvar)
  status = SDend(sd_id);
}





void readMOD06_sensor_angles(char *modPath,char *MOD06file,
               float *sensor_zenith,float *sensor_azimuth)
{
  char sds_name[128];
  char fname[150];
  short *tmpvar=NULL;
  signed short fill;
  int sd_id, sds_id, sds_index, attr_index;
  int status, start[2], edge[2];
  int rank, dims[2], type, attributes;
  long ix, ia, iloc, na_5km;
  double scale;

  strcpy(fname,modPath);
  strncat(fname,MOD06file,NCHAR_FNAME_MOD06);
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  start[0] = start[1] = 0;
  edge[1] = N_ACROSS_5KM;

  tmpvar = (short *) calloc((N_ACROSS_5KM*N_ALONG_5KM_MAX),sizeof(short));
  sprintf(sds_name,"Sensor_Zenith");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
//  printf("%d %d %d\n",sds_index,sds_id,status);
  na_5km = dims[0];
  edge[0] = na_5km;
  status = SDreaddata(sds_id,start,NULL,edge,tmpvar);
  status = SDgetfillvalue(sds_id,&fill);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  for (ix=0; ix<N_ACROSS_5KM; ix++)
    for (ia=0; ia<na_5km; ia++) {
      iloc = (ia*N_ACROSS_5KM)+ix;
      if (tmpvar[iloc] != fill) sensor_zenith[iloc] =
                    ((float) tmpvar[iloc])*scale;
      else sensor_zenith[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);

  sprintf(sds_name,"Sensor_Azimuth");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmpvar);
  status = SDgetfillvalue(sds_id,&fill);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  for (ix=0; ix<N_ACROSS_5KM; ix++)
    for (ia=0; ia<na_5km; ia++) {
      iloc = (ia*N_ACROSS_5KM)+ix;
      if (tmpvar[iloc] != fill) sensor_azimuth[iloc] =
                    ((float) tmpvar[iloc])*scale;
      else sensor_zenith[iloc] = MOD06_BADFLAG;
  }
  status = SDendaccess(sds_id);

  FREE(tmpvar)
  status = SDend(sd_id);
}






void readMOD06_1km_cldmask(char *modPath,char *MOD06file,long na_1km,
                           signed char *cloud_mask)
{
  char sds_name[128];
  char fname[150];
  signed char *masktmp;
  int sd_id, sds_id, sds_index;
  int status, start[3], edge[3];
  long ix, ia, iloc, maskloc;

  sprintf(fname,"%s%s",modPath,MOD06file);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  start[0] = start[1] = start[2] = 0;
  edge[0] = na_1km;
  edge[1] = N_ACROSS_1KM;
  edge[2] = N_MASK_BYTES;
  sprintf(sds_name,"Cloud_Mask_1km");
  masktmp = (signed char *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX*
              N_MASK_BYTES),sizeof(signed char));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,masktmp);

  for (ix=0; ix<N_ACROSS_1KM; ix++)                                
    for (ia=0; ia<na_1km; ia++) {
      iloc = (ia*N_ACROSS_1KM*N_MASK_BYTES)+
               (ix*N_MASK_BYTES);
      maskloc = (ia*N_ACROSS_1KM*N_MASK_TYPES) + 
	        (ix*N_MASK_TYPES);
      cloud_mask[maskloc] = (masktmp[iloc] & 1);
      maskloc++;
      cloud_mask[maskloc] = (masktmp[iloc] & 6) >> 1;
      maskloc++;
      cloud_mask[maskloc] = (masktmp[iloc] & 8) >> 3;
      maskloc++;
      cloud_mask[maskloc] = (masktmp[iloc] & 16) >> 4;
      maskloc++;
      cloud_mask[maskloc] = (masktmp[iloc] & 32) >> 5;
      maskloc++;
      cloud_mask[maskloc] = (masktmp[iloc] & 192) >> 6;
      iloc++;
      maskloc++;
      cloud_mask[maskloc] = (masktmp[iloc] & 1);
      maskloc++;
      cloud_mask[maskloc] = (masktmp[iloc] & 2) >> 1;
      maskloc++;
      cloud_mask[maskloc] = (masktmp[iloc] & 4) >> 2;
  }
  FREE(masktmp);
  status = SDend(sd_id);

}




void readMOD06_QA_1km(char *modPath,char *MOD06file,long na_1km,
                           signed char *qa_1km)
{
  char sds_name[128];
  char fname[150];
  signed char *qatmp;
  int sd_id, sds_id, sds_index;
  int status, start[3], edge[3];
  long ix, ia, iloc, qaloc;

  sprintf(fname,"%s%s",modPath,MOD06file);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  start[0] = start[1] = start[2] = 0;
  edge[0] = na_1km;
  edge[1] = N_ACROSS_1KM;
  edge[2] = N_QA1KM_BYTES;
  sprintf(sds_name,"Quality_Assurance_1km");
  qatmp = (signed char *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX*
            N_QA1KM_BYTES),sizeof(signed char));
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,qatmp);

  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<na_1km; ia++) {
      iloc = (ia*N_ACROSS_1KM*N_QA1KM_BYTES)+
               (ix*N_QA1KM_BYTES);
      qaloc = (ia*N_ACROSS_1KM*N_QA1KM_TYPES) +
                (ix*N_QA1KM_TYPES);
      qa_1km[qaloc] = (qatmp[iloc] & 1);
      qaloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 6) >> 1;
      qaloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 24) >> 3;
      qaloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 32) >> 5;
      qaloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 192) >> 6;
      qaloc++;
      iloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 1);
      qaloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 6) >> 1;
      qaloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 56) >> 3;
      iloc++;
      qaloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 7);
      qaloc++;
      qa_1km[qaloc] = (qatmp[iloc] & 8) >> 3;
  }
  FREE(qatmp);
  status = SDend(sd_id);

}




void match_1km_2_5km(float *lon_1km,float *lat_1km,float *lon_5km,
                     float *lat_5km,long na_1km,long na_5km,
                     long *index_1to5km,float *dist_1to5km)
{
  long ix, ia, iloc, ia_5km, ix_5km; 
  float dx, dy, deg2rad;

  deg2rad = PI/180.;

  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<na_1km; ia++) {
      ia_5km = floor(ia/5.);
      ix_5km = floor(ix/5.);
      iloc = (ia*N_ACROSS_1KM)+ix;
      if ( (ia_5km<na_5km) && (ix_5km<N_ACROSS_5KM) ) {
        index_1to5km[iloc] = (ia_5km*N_ACROSS_5KM) + ix_5km;
        dx = RE*(lon_1km[iloc]-lon_5km[index_1to5km[iloc]])*deg2rad*
             cos(lat_1km[iloc]*deg2rad);                              
        dy = RE*(lat_1km[iloc]-lat_5km[index_1to5km[iloc]])*deg2rad;  
        dist_1to5km[iloc] = sqrt((dx*dx) + (dy*dy));                  
      } else {
        index_1to5km[iloc] = MOD06_BADFLAG;
        dist_1to5km[iloc] = MOD06_BADFLAG;
      }
  }
}










void screen_reff(float *lon_1km,float *lat_1km,long na_1km,
                 float *reff_raw,float *tau,signed char *cldtop_phase,
		 signed char *cloud_mask,long *index_1to5km,
                 signed char *qa_1km,long *npix,long *nglnt_pix,
                 long *ncld_pix,long *nliq_pix,long *nscrn_pix,float *reff)
{
  long ix, ia, ixstep, iastep, iloc, maskloc, qaloc;
  long ixstart, ixend, iastart, iaend, steploc, ixtmp;
  long np, nc, nl, ns, ng;
  float dx, dy, deg2rad, dist;

  deg2rad = PI/180.;

/* remove the last 4 pixels of each scan since they have*/
/* no corresponding 5km pixel*/
  for (ix=N_ACROSS_1KM-4; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<na_1km; ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      reff_raw[iloc] = MOD06_BADFLAG;
  }

  np = *npix;
  nc = *ncld_pix;
  nl = *nliq_pix;
  ns = *nscrn_pix;
  ng = *nglnt_pix;

  ixstart = NXSTEP;
  iastart = NASTEP;
  ixend = N_ACROSS_1KM-5;
  ixtmp = N_ACROSS_1KM-NXSTEP-1;
  if (ixtmp < ixend) ixend=ixtmp;
  iaend = na_1km-NASTEP-1;

  for (ix=0; ix<ixstart; ix++)
    for (ia=0; ia<na_1km; ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      reff[iloc] = MOD06_SCRNFLAG;
  }
  for (ix=ixend; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<na_1km; ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      reff[iloc] = MOD06_SCRNFLAG;
  }
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=0; ia<iastart; ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      reff[iloc] = MOD06_SCRNFLAG;
  }
  for (ix=0; ix<N_ACROSS_1KM; ix++)
    for (ia=iaend; ia<na_1km; ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      reff[iloc] = MOD06_SCRNFLAG;
  }

  for (ix=ixstart; ix<ixend; ix++)
    for (ia=iastart; ia<iaend; ia++) {
      iloc = (ia*N_ACROSS_1KM)+ix;
      np++;

/* exclude the sun glint pixels */
      maskloc = (ia*N_ACROSS_1KM*N_MASK_TYPES) +
                (ix*N_MASK_TYPES)+3;
      qaloc = (ia*N_ACROSS_1KM*N_QA1KM_TYPES) +
                (ix*N_QA1KM_TYPES)+8;
      if (cloud_mask[maskloc] == 1) { /* no glint? */
/* limits on unambiguous retrieval - Nakajima and King (1990) */ 
      if ( (reff_raw[iloc]>=6) && (tau[iloc]>=4) &&
           (index_1to5km[iloc]!=MOD06_BADFLAG) ) {
	nc++;
/* is it liquid? */
/*        if ( (cldtop_phase[index_1to5km[iloc]] == 1) ||
             (cldtop_phase[index_1to5km[iloc]] == 5) ) { */
        if (qa_1km[qaloc] == 2) {
        
	reff[iloc]=reff_raw[iloc];
	nl++;
	ns++;

/* screen out pixels without a cloudy pixel with DIST_THRESH */
	for (ixstep=-NXSTEP; ixstep<=NXSTEP; ixstep++)
          for (iastep=-NASTEP; iastep<=NASTEP; iastep++) {          
            steploc = ((ia+iastep)*N_ACROSS_1KM) + (ix+ixstep);
            dx = RE*(lon_1km[iloc]-lon_1km[steploc])*deg2rad *      
                 cos(lat_1km[iloc]*deg2rad);                        
            dy = RE*(lat_1km[iloc]-lat_1km[steploc])*deg2rad;       
            dist = sqrt((dx*dx)+(dy*dy));                           
            if ( (dist<DIST_THRESH) && (reff_raw[steploc] <= 0.) ) {
              reff[iloc] = MOD06_SCRNFLAG;                          
	      ns--;
              ixstep=NXSTEP+1;                                      
              iastep=NASTEP+1;                                      
            }
        }
	} else reff[iloc] = MOD06_SCRNFLAG; /*ice*/
      } else reff[iloc] = MOD06_SCRNFLAG;   /*re or tau too low*/
      } else {
        reff[iloc] = MOD06_SCRNFLAG;        /*sun glint*/
        ng++;
      }

/* label the likely clear pixels */
    maskloc = (ia*N_ACROSS_1KM*N_MASK_TYPES) +
                (ix*N_MASK_TYPES)+1;
    if ( (cloud_mask[maskloc] > 1) && (tau[iloc]<4) )
      reff[iloc] = MOD06_CLRFLAG;

  }
  *npix = np;
  *ncld_pix = nc;
  *nliq_pix = nl;
  *nscrn_pix = ns;
  *nglnt_pix = ng;
}



void getMODneighbors(float lon,float lat,float dist_thresh,
                     float *lon_mod,float *lat_mod,long na,long nx,
                     long *ineighbors,long *cnt,long ixmin, long ixmax,
                     long iamin,long iamax)
{
  long ix, ia, iloc, count, ixpnt, iapnt;
  float dx, dy, deg2rad, dist, mindist=10000.;

  deg2rad = PI/180.;
  count = 0;  //counts the number of neighbors

/*** first skip through the image to get close to the point ***/
  for (ix=ixmin; ix<ixmax; ix+=100)
    for (ia=iamin; ia<iamax; ia+=100) {
      iloc = (ia*nx)+ix;
      dx = RE*(lon-lon_mod[iloc])*deg2rad*cos(lat*deg2rad);
      dy = RE*(lat-lat_mod[iloc])*deg2rad;
      dist = sqrt((dx*dx)+(dy*dy));
      if (dist<mindist) {
        mindist=dist;
        ixpnt = ix;
        iapnt = ia;
      }
  }
  if (mindist>100) {
    printf("no MODIS sample found within 100 km of the AMSR sample\n");
  } else {
    ixmin=ixpnt-100;
    ixmax=ixpnt+100;
    iamin=iapnt-100;
    iamax=iapnt+100;
    for (ix=ixmin; ix<ixmax; ix++)
    for (ia=iamin; ia<iamax; ia++) {
      iloc = (ia*nx)+ix;
      dx = RE*(lon-lon_mod[iloc])*deg2rad*cos(lat*deg2rad);
      dy = RE*(lat-lat_mod[iloc])*deg2rad;
      dist = sqrt((dx*dx)+(dy*dy));
      if (dist < dist_thresh) {
        ineighbors[count] = iloc;
        count++;
      }
    }
  }
  *cnt = count;
}


/********************
void getMODneighbors(float lon,float lat,float dist_thresh,
                     float *lon_mod,float *lat_mod,long na,long nx,
                     long *ineighbors,long *cnt,long ixmin, long ixmax,
                     long iamin,long iamax)
{
  int jump_flag;
  long ix, ia, iloc, count, loopcnt;
  float dx, dy, deg2rad, dist;

  deg2rad = PI/180.;

  count = 0;  //counts the number of neighbors
  loopcnt = 0; //counts the number of test points since the last neighbor was found
  for (ix=ixmin; ix<ixmax; ix++) {
if (ix>ixmin) {
printf("%d %d %d\n",count,loopcnt,jump_flag);
exit(0);
}
    for (ia=iamin; ia<iamax; ia++) {
      iloc = (ia*nx)+ix;
      dx = RE*(lon-lon_mod[iloc])*deg2rad*cos(lat*deg2rad);
      dy = RE*(lat-lat_mod[iloc])*deg2rad;
      dist = sqrt((dx*dx)+(dy*dy));
      jump_flag=0; //checks for case of jumping to a neighbor
//if (iloc > 1255900 && iloc <1256000) printf("%f %d %d\n",dist,count,loopcnt);
      if ( (dist > 100.) && ((count==0) || ((count>0 && loopcnt>10000))) ) {
        ia+=5;
        jump_flag++;
      }
      if (dist > 1000.) {
        ia+=95;
        jump_flag++;
      }
      loopcnt++;
    
//if (ia==(na/2) && ix==(nx/2)) printf("%f\n",dist);
//printf("%f ",dist);
      if (dist < dist_thresh) {
        if (!jump_flag) {
          ineighbors[count] = iloc;
          count++;
          loopcnt=0;
//printf("[%f %f %f %d]\n",lat_mod[iloc],lon_mod[iloc],dist,count);
        } else {
          if (jump_flag==2) {
            ia-=95;
            jump_flag--;
          }
          if (jump_flag==1) {
            ia-=5;
          }
        }
      }
    }
    if (jump_flag>0
  }
  *cnt = count;
}
*******************/



/***********************
* note that the MOD02 files were subsetted in the initial version using MOD02
* as a result the filenames were longer in terms of # of characters.
* now modified to use standard product files with a 45 character filename.
***********************/
void getFileNames(char *modPath, char *mod02Path, char *yyyyddd,char *hhmm,
									const char *MOD03list,const char *MOD04list,const char *MOD02list,
                  char *MOD03file,char *MOD04file,char *MOD02file)
{
  FILE *fp=NULL;
  char filelist[150];
  char tmpfile[80];
  char *noneStr = "none";
  char tmp_yyyyddd[8], tmp_hhmm[5];
  int endflag=0;

  strcpy(MOD03file,noneStr);
  strcpy(MOD04file,noneStr);
  strcpy(MOD02file,noneStr);
  sprintf(filelist,"%s%s",modPath,MOD03list);
  fp = fopen(filelist,"r");
  while (!endflag) {
    if ( fgets(tmpfile,50,fp) != NULL ) {
      strncpy(tmp_yyyyddd,tmpfile+7,7);
      tmp_yyyyddd[7] = '\0';
      strncpy(tmp_hhmm,tmpfile+15,4);
      tmp_hhmm[4] = '\0';
      if ( (strcmp(tmp_yyyyddd,yyyyddd)==0) && (strcmp(tmp_hhmm,hhmm)==0) ) {
        strncpy(MOD03file,tmpfile,41);
        MOD03file[41] = '\0';
        endflag=1;
      }
    } else endflag=1;
  }
  fclose(fp);
//  if ( strcmp(MOD03file,noneStr) != 0 )

  endflag = 0;
  sprintf(filelist,"%s%s",modPath,MOD04list);
  if ( (fp = fopen(filelist,"r")) ) {
		while (!endflag) {
			if ( fgets(tmpfile,50,fp) != NULL ) {
				strncpy(tmp_yyyyddd,tmpfile+10,7);
				strncpy(tmp_hhmm,tmpfile+18,4);
				if ( (strcmp(tmp_yyyyddd,yyyyddd)==0) && (strcmp(tmp_hhmm,hhmm)==0) ) {
					strncpy(MOD04file,tmpfile,44);
        	MOD04file[44] = '\0';
					endflag=1;
				}
			} else endflag=1;
		}
  	fclose(fp);
	}

  endflag = 0;
  sprintf(filelist,"%s%s",mod02Path,MOD02list);
  fp = fopen(filelist,"r");
  while (!endflag) {
    if ( fgets(tmpfile,80,fp) != NULL ) {
      strncpy(tmp_yyyyddd,tmpfile+10,7);
      strncpy(tmp_hhmm,tmpfile+18,4);
      if ( (strcmp(tmp_yyyyddd,yyyyddd)==0) && (strcmp(tmp_hhmm,hhmm)==0) ) {
        strncpy(MOD02file,tmpfile,62);
//        MOD02file[62] = '\0';
        MOD02file[44] = '\0';
        endflag=1;
      }
    } else endflag=1;
  }
  fclose(fp);
}



void readMOD06_IR_tbs(char *modPath,char *MOD06file,long *na_5km,
                      float32 *lat_5km,float *lon_5km,float *Tbs)
{
  FILE *f;
  char sds_name[128];
  char fname[150];
  short *tmpvar=NULL;
  signed short fill;
  int sd_id, sds_id, sds_index, attr_index;
  int status, start[3], edge[3];
  int rank, dims[3], type, attributes;
  long ix, ia, ib, iloc;
  double scale, offset;
  int32 i, e; 
	const char *str;
	float32 *flttmp;
  
  for (ix=0; ix<3; ix++) {
    start[ix]=0;
    edge[ix]=0;
  }

  sprintf(fname,"%s%s",modPath,MOD06file);
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
//	printf("%d ",sd_id);
  start[0] = start[1] = start[2] = 0;
  edge[1] = N_ACROSS_5KM;

  sprintf(sds_name,"Latitude");
  sds_index = SDnametoindex(sd_id,sds_name);
//	printf("reading lats: %d ",sds_index);
  sds_id = SDselect(sd_id,sds_index);
//	printf("%d ",sds_id);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
//	printf("%d ",status);
  if ( (dims[0]>N_ALONG_5KM_MAX) || (dims[1]>N_ACROSS_5KM) ) {
    printf("MODIS 5km dimensions = [%d %d] - too large\n",dims[0],dims[1]);
    status = SDendaccess(sds_id);
    status = SDend(sd_id);
    exit(0);
  }
  *na_5km = dims[0];
  edge[0] = *na_5km;
//printf("lat_5km = %ld\n",lat_5km);
  status = SDreaddata(sds_id,start,NULL,edge,lat_5km);
//  status = SDreaddata(sds_id,start,NULL,edge,flttmp);
//printf("\n%d %d %d %d %d %d\n",start[0],start[1],start[2],edge[0],edge[1],edge[2]);
//printf("%d ",status);
if (status==FAIL) {
//printf("into error handling ... ");
/****
i = 0; 
//while ((e = HEvalue(i)) != DFE_NONE) {
while (i<5) {
  e = HEvalue(i);
  str = HEstring(e); 
	printf("%s\n",str);
			i++;
}
printf("end\n");
f = fopen("/opt/eric/code/scripts/error.txt","a");
	HEprint(f, 0);
fclose(f);
***/
}
  status = SDendaccess(sds_id);
//printf("%d\n",status);

  sprintf(sds_name,"Longitude");
  sds_index = SDnametoindex(sd_id,sds_name);
//printf("reading lons: %d ",sds_index);
  sds_id = SDselect(sd_id,sds_index);
//printf("%d ",sds_id);
  status = SDreaddata(sds_id,start,NULL,edge,lon_5km);
//printf("%d ",status);
  status = SDendaccess(sds_id);
//printf("%d\n",status);
  
/*  for (ix=0; ix<N_ACROSS_5KM; ix++)
  for (ia=0; ia<(*na_5km); ia++) {
    iloc = (ia*N_ACROSS_5KM) + ix;
    printf("[%f %f]  ",lat_5km[iloc],lon_5km[iloc]);
  }*/

  tmpvar = (short *) calloc((N_ACROSS_5KM*N_ALONG_5KM_MAX*N_IR_BANDS),
                            sizeof(short));
  sprintf(sds_name,"Brightness_Temperature");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  edge[0] = N_IR_BANDS;
  edge[2] = N_ACROSS_5KM;
  edge[1] = *na_5km;
  status = SDreaddata(sds_id,start,NULL,edge,tmpvar);
  status = SDgetfillvalue(sds_id,&fill);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  status = SDendaccess(sds_id);
  for (ix=0; ix<N_ACROSS_5KM; ix++)
  for (ia=0; ia<(*na_5km); ia++)
  for (ib=0; ib<N_IR_BANDS; ib++) {
    iloc = (ia*N_ACROSS_5KM*N_IR_BANDS) + (ix*N_IR_BANDS) + ib;
    if (tmpvar[iloc] != fill) Tbs[iloc] = scale *
          (((float) tmpvar[iloc])-offset);
    else Tbs[iloc] = -9999.;
  }


  status = SDend(sd_id);

  FREE(tmpvar);
}




void readMOD06_lwpuncertainty(char *modPath,char *MOD06file,long na_1km,
                              float *lwp_err,float *tau_err)
{
  char sds_name[128];
  char fname[150];
  short *tmpvar=NULL;
  signed short fill;
  int sd_id, sds_id, sds_index, attr_index;
  int status, start[2], edge[2];
  long ix, ia, iloc;
  double scale, offset;

  sprintf(fname,"%s%s",modPath,MOD06file);
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  start[0] = start[1] = 0;
  edge[0] = na_1km;
  edge[1] = N_ACROSS_1KM;

  tmpvar = (short *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX),sizeof(short));

  sprintf(sds_name,"Cloud_Water_Path_Uncertainty");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmpvar);
  status = SDgetfillvalue(sds_id,&fill);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  status = SDendaccess(sds_id);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
  for (ia=0; ia<na_1km; ia++) {
    iloc = (ia*N_ACROSS_1KM)+ix;
    if (tmpvar[iloc] != fill) lwp_err[iloc] = scale *
          (((float) tmpvar[iloc])-offset);
    else lwp_err[iloc] = -9999.;
  }
  FREE(tmpvar);

  tmpvar = (short *) calloc((N_ACROSS_1KM*N_ALONG_1KM_MAX),sizeof(short));
  sprintf(sds_name,"Cloud_Optical_Thickness_Uncertainty");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmpvar);
  status = SDgetfillvalue(sds_id,&fill);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  status = SDendaccess(sds_id);
  for (ix=0; ix<N_ACROSS_1KM; ix++)
  for (ia=0; ia<na_1km; ia++) {
    iloc = (ia*N_ACROSS_1KM)+ix;
    if (tmpvar[iloc] != fill) tau_err[iloc] = scale *
          (((float) tmpvar[iloc])-offset);
    else tau_err[iloc] = -9999.;
//printf("%f ",tau_err[iloc]);
  }
  FREE(tmpvar);

  status = SDend(sd_id);
}





char get_MOD06_filename(const char *MOD06_path,const char *list_filename,
                        int year,int doy,char *MOD06_file)
{
  FILE *fp=NULL;
  char filelist[150], *s;
  char tmpfile[43];
  char yyyy[5], ddd[4];
  char outfile=0;
  int yr, day, endflag=0;
  
  sprintf(filelist,"%s%s",MOD06_path,list_filename);
  fp = fopen(filelist,"r");
  while (!feof(fp)) {
    if ( fgets(tmpfile,41,fp) != NULL ) {
      strncpy(yyyy,tmpfile+10,4);
      yyyy[4] = '\0';
      yr = atoi(yyyy);
      strncpy(ddd,tmpfile+14,3);
      ddd[3] = '\0';
      day = atoi(ddd);
      if ( (yr==year) && (day==doy) ) {
        strncpy(MOD06_file,tmpfile,39);
        endflag=1;
      }
    }
  }
  fclose(fp);
  return endflag;
}








/************************************************
returns fraction of the day between the specified time and the modis granual
*************************************************/
float get_nearest_MOD06_gran(double jdate,char *MOD06_path,const char *MOD06list,
                          const char *MYD06list,const char *MOD04list,
                          const char *MYD04list,const char *MOD03list,
                          const char *MYD03list,char *MOD06file,char *MOD03file)
{
  FILE *fp=NULL;
  char terra_flag=-1;
  char filelist[150], MOD04file[50], MOD02file[50], *mod02path="none";
  char tmpfile[43];
  char hh[3],mm[3], yyyyddd[8], hhmm[5];
  char yyyy[5], dd[3], ddd[4], mmdd[5];
  int mon,day,yr,hr,min,doy;
  long yyyymmdd;
  float sec;
  double jdate_test, dt=9999, dt_test;
  
  const char MYD02list[150],MOD02list[150];
  
// find the nearest Terra overpass
  sprintf(filelist,"%s%s",MOD06_path,MOD06list);
  fp = fopen(filelist,"r");
  while (!feof(fp)) {
    if ( fgets(tmpfile,41,fp) != NULL ) {
      strncpy(yyyy,tmpfile+10,4);
      yyyy[4] = '\0';
      yr = atoi(yyyy);
      strncpy(ddd,tmpfile+14,3);
      ddd[3] = '\0';
      doy = atoi(ddd);
      doy_to_mmdd(doy,yr,mmdd);
      strncpy(mm,mmdd,2);
      mm[2] = '\0';
      mon = atoi(mm);
      strncpy(dd,mmdd,2);
      dd[2] = '\0';
      day = atoi(dd);
      strncpy(hh,MOD06file+18,2);
      hh[2] = '\0';
      hr = atoi(hh);
      strncpy(mm,MOD06file+20,2);
      mm[2] = '\0';
      min = atoi(mm);
      
      mdyhms_to_juldate(mon,day,yr,hr,min,0,&jdate_test);
      dt_test = fabs(jdate_test-jdate);
      
      if (dt_test < dt) {
        dt = dt_test;
        strcpy(MOD06file,tmpfile);
        terra_flag=1;
        strcpy(yyyyddd,yyyy);
        strcat(yyyyddd,ddd);
        strcpy(hhmm,hh);
        strcat(hhmm,mm);
      }
    }
  }
  fclose(fp);
  
// now see if there is an Aqua pass that is closer
  sprintf(filelist,"%s%s",MOD06_path,MYD06list);
  fp = fopen(filelist,"r");
  while (!feof(fp)) {
    if ( fgets(tmpfile,41,fp) != NULL ) {
      strncpy(yyyy,tmpfile+10,4);
      yyyy[4] = '\0';
      yr = atoi(yyyy);
      strncpy(ddd,tmpfile+14,3);
      ddd[3] = '\0';
      doy = atoi(ddd);
      doy_to_mmdd(doy,yr,mmdd);
      strncpy(mm,mmdd,2);
      mm[2] = '\0';
      mon = atoi(mm);
      strncpy(dd,mmdd,2);
      dd[2] = '\0';
      day = atoi(dd);
      strncpy(hh,MOD06file+18,2);
      hh[2] = '\0';
      hr = atoi(hh);
      strncpy(mm,MOD06file+20,2);
      mm[2] = '\0';
      min = atoi(mm);
      
      mdyhms_to_juldate(mon,day,yr,hr,min,0,&jdate_test);
      dt_test = fabs(jdate_test-jdate);
      
      if (dt_test < dt) {
        dt = dt_test;
        strcpy(MOD06file,tmpfile);
        terra_flag=0;
        strcpy(yyyyddd,yyyy);
        strcat(yyyyddd,ddd);
        strcpy(hhmm,hh);
        strcat(hhmm,mm);
      }
    }
  }
  fclose(fp);
  
  if (terra_flag)
    getFileNames(MOD06_path,mod02path,yyyyddd,hhmm,MOD03list,MOD04list,MYD02list,
                 MOD03file,MOD04file,MOD02file);
  else
    getFileNames(MOD06_path,mod02path,yyyyddd,hhmm,MYD03list,MYD04list,MYD02list,
                 MOD03file,MOD04file,MOD02file);
    
  
  return (float) dt;
}






void create_test_image(float *ir_img,float *lons,float *lats,
											long na_5km,float lon_cent,float lat_cent)
{
	long ix,ia,iloc;
	float dist;
	
  for (ix=0; ix<N_ACROSS_5KM; ix++)
    for (ia=0; ia<na_5km; ia++) {
      iloc = (ia*N_ACROSS_5KM)+ix;
			dist = calc_great_circle_distance(lats[iloc],lons[iloc],lat_cent,lon_cent);
//			printf("%f %f %f %f %f\n",lats[iloc],lat_cent,lons[iloc],lon_cent,dist);
			ir_img[iloc] = 180. + (dist/5.);
	}
	
}







