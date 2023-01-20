#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "MOD08.h"
#include "time_tools.h"

const char *reff_nameC004 = "Cloud_Effective_Radius_Liquid_Mean_Mean";
const char *reff_nameC051 = "Cloud_Effective_Radius_Liquid_Mean";
const char *cldfrac_nameC004 = "Cloud_Fraction_Mean_Mean";
const char *cldfrac_nameC051 = "Cloud_Fraction_Mean";
const char *cldtmp_nameC004 = "Cloud_Top_Temperature_Mean_Mean";
const char *cldtmp_nameC051 = "Cloud_Top_Temperature_Mean";
const char *cldtau_nameC004 = "Cloud_Optical_Thickness_Liquid_Mean_Mean";
const char *cldtau_nameC051 = "Cloud_Optical_Thickness_Liquid_Mean";
const char *cldpres_nameC051 = "Cloud_Top_Pressure_Mean";
const char *lwp_nameC004 = "Cloud_Water_Path_Liquid_Mean_Mean";
const char *lwp_nameC051 = "Cloud_Water_Path_Liquid_Mean";
const char *aodDB_nameC051 = "Deep_Blue_Aerosol_Optical_Depth_Land_Mean";
const char *aod_nameC004 = "Optical_Depth_Land_And_Ocean_Mean_Mean";
const char *aod_nameC051 = "Optical_Depth_Land_And_Ocean_Mean";
const char *aod_nameC061 = "Aerosol_Optical_Depth_Land_Ocean_Mean";

/*****************************
* shifts data so that lons:0.5-359.5 lats:-89.5-89.5
*****************************/
void read_mod08_aerosol_cloud(const char *modPath,char *MOD08file,float *lons,
                              float *lats,float *reff,
                              float *cld_tau,float *cldtop_temp,float *lwp,
                              float *aod,float *cld_frac,float *cldtop_pres,
                              int collection)
{
  char sds_name[128];
  char fname[150];
  short fill, *tmp=NULL;
  int ilo, ila, i;
  int sd_id, sds_id, sds_index, attr_index;
  int status, start[2], edge[2];
  int rank, dims[2], type, attributes;
  long iloc, cloc;
  float *flttmp=NULL;
  double scale, offset;

  sprintf(fname,"%s%s",modPath,MOD08file);
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
//  printf("sd_id=%d",sd_id);
  start[0] = start[1] = 0;
  edge[1] = N_MOD08_LONS;
  edge[0] = N_MOD08_LATS;

  flttmp = (float *) calloc(N_MOD08_LONS,sizeof(float));
  sprintf(sds_name,"XDim");
  sds_index = SDnametoindex(sd_id,sds_name);
//  printf("  sds_index=%d",sds_index);
  sds_id = SDselect(sd_id,sds_index);
//  printf("  sds_id=%d",sds_id);
  status = SDreaddata(sds_id,start+1,NULL,edge+1,flttmp);
//  printf("  %d",status);
  status = SDendaccess(sds_id);
//  printf("  %d\n",status);
  for (ilo=0; ilo<180; ilo++) lons[ilo] = flttmp[ilo+180];
  for (ilo=180; ilo<360; ilo++) lons[ilo] = 360+flttmp[ilo-180];
  FREE(flttmp);

  flttmp = (float *) calloc(N_MOD08_LATS,sizeof(float));
  sprintf(sds_name,"YDim");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,flttmp);
  status = SDendaccess(sds_id);
  for (ila=0; ila<N_MOD08_LATS; ila++) lats[ila] = flttmp[N_MOD08_LATS-1-ila];
  FREE(flttmp);

  tmp = (short *) calloc(N_MOD08_LATS*N_MOD08_LONS,sizeof(short));
	if (collection==4) strcpy(sds_name,reff_nameC004);
	if (collection==51) strcpy(sds_name,reff_nameC051);
	if (collection==61) strcpy(sds_name,reff_nameC051);
//  sprintf(sds_name,"Cloud_Effective_Radius_Liquid_Mean");
//  sprintf(sds_name,"Cloud_Effective_Radius_Liquid_Mean_Mean");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  attr_index = SDfindattr(sds_id,"_FillValue");
  status = SDreadattr(sds_id,attr_index,&fill);
  status = SDendaccess(sds_id);
  for (ilo=0; ilo<N_MOD08_LONS; ilo++)
  for (ila=0; ila<N_MOD08_LATS; ila++) {
    iloc = (ila*N_MOD08_LONS)+ilo;
    if (ilo<180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo+180);
    if (ilo>=180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo-180);
    if (tmp[cloc] != fill) reff[iloc] = scale * (((float) tmp[cloc])-offset);
    else reff[iloc] = MOD08_BADFLAG;
  }
  FREE(tmp);

  tmp = (short *) calloc(N_MOD08_LATS*N_MOD08_LONS,sizeof(short));
	if (collection==4) strcpy(sds_name,cldfrac_nameC004);
	if (collection==51) strcpy(sds_name,cldfrac_nameC051);
	if (collection==61) strcpy(sds_name,cldfrac_nameC051);
//  sprintf(sds_name,"Cloud_Fraction_Mean");  
//  sprintf(sds_name,"Cloud_Fraction_Mean_Mean");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  attr_index = SDfindattr(sds_id,"_FillValue");
  status = SDreadattr(sds_id,attr_index,&fill);
  status = SDendaccess(sds_id);
  for (ilo=0; ilo<N_MOD08_LONS; ilo++)
  for (ila=0; ila<N_MOD08_LATS; ila++) {
    iloc = (ila*N_MOD08_LONS)+ilo;
    if (ilo<180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo+180);
    if (ilo>=180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo-180);
    if (tmp[cloc] != fill) cld_frac[iloc] = scale * (((float) tmp[cloc])-offset);
    else cld_frac[iloc] = MOD08_BADFLAG;
  }
  FREE(tmp);

  tmp = (short *) calloc(N_MOD08_LATS*N_MOD08_LONS,sizeof(short));
	if (collection==4) strcpy(sds_name,cldtmp_nameC004);
	if (collection==51) strcpy(sds_name,cldtmp_nameC051);
	if (collection==61) strcpy(sds_name,cldtmp_nameC051);
//  sprintf(sds_name,"Cloud_Top_Temperature_Mean");
//  sprintf(sds_name,"Cloud_Top_Temperature_Mean_Mean");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  attr_index = SDfindattr(sds_id,"_FillValue");
  status = SDreadattr(sds_id,attr_index,&fill);
  status = SDendaccess(sds_id);
  for (ilo=0; ilo<N_MOD08_LONS; ilo++)
  for (ila=0; ila<N_MOD08_LATS; ila++) {
    iloc = (ila*N_MOD08_LONS)+ilo;
    if (ilo<180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo+180);
    if (ilo>=180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo-180);
    if (tmp[cloc] != fill) cldtop_temp[iloc] = scale * (((float) tmp[cloc])-offset);
    else cldtop_temp[iloc] = MOD08_BADFLAG;
  }
  FREE(tmp);

  tmp = (short *) calloc(N_MOD08_LATS*N_MOD08_LONS,sizeof(short));
//  sprintf(sds_name,"Cloud_Optical_Thickness_Liquid_Mean");
//  sprintf(sds_name,"Cloud_Optical_Thickness_Liquid_Mean_Mean");
	if (collection==4) strcpy(sds_name,cldtau_nameC004);
	if (collection==51) strcpy(sds_name,cldtau_nameC051);
	if (collection==61) strcpy(sds_name,cldtau_nameC051);
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  attr_index = SDfindattr(sds_id,"_FillValue");
  status = SDreadattr(sds_id,attr_index,&fill);
  status = SDendaccess(sds_id);
  for (ilo=0; ilo<N_MOD08_LONS; ilo++)
  for (ila=0; ila<N_MOD08_LATS; ila++) {
    iloc = (ila*N_MOD08_LONS)+ilo;
    if (ilo<180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo+180);
    if (ilo>=180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo-180);
    if (tmp[cloc] != fill) cld_tau[iloc] = scale * (((float) tmp[cloc])-offset);
    else cld_tau[iloc] = MOD08_BADFLAG;
  }
  FREE(tmp);

  tmp = (short *) calloc(N_MOD08_LATS*N_MOD08_LONS,sizeof(short));
//  sprintf(sds_name,"Cloud_Water_Path_Liquid_Mean");
//  sprintf(sds_name,"Cloud_Water_Path_Liquid_Mean_Mean");
	if (collection==4) strcpy(sds_name,lwp_nameC004);
	if (collection==51) strcpy(sds_name,lwp_nameC051);
	if (collection==61) strcpy(sds_name,lwp_nameC051);
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  attr_index = SDfindattr(sds_id,"_FillValue");
  status = SDreadattr(sds_id,attr_index,&fill);
  status = SDendaccess(sds_id);
  for (ilo=0; ilo<N_MOD08_LONS; ilo++)
  for (ila=0; ila<N_MOD08_LATS; ila++) {
    iloc = (ila*N_MOD08_LONS)+ilo;
    if (ilo<180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo+180);
    if (ilo>=180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo-180);
    if (tmp[cloc] != fill) lwp[iloc] = scale * (((float) tmp[cloc])-offset);
    else lwp[iloc] = MOD08_BADFLAG;
  }
  FREE(tmp);

  tmp = (short *) calloc(N_MOD08_LATS*N_MOD08_LONS,sizeof(short));
#ifdef DEEPBLUE_FLAG
//  sprintf(sds_name,"Deep_Blue_Aerosol_Optical_Depth_Land_Mean");
	if (collection==51) strcpy(sds_name,aodDB_nameC051);
	if (collection==61) strcpy(sds_name,aodDB_nameC051);
#else
//  sprintf(sds_name,"Optical_Depth_Land_And_Ocean_Mean");
//  sprintf(sds_name,"Optical_Depth_Land_And_Ocean_Mean_Mean");
	if (collection==4) strcpy(sds_name,aod_nameC004);
	if (collection==51) strcpy(sds_name,aod_nameC051);
	if (collection==61) strcpy(sds_name,aod_nameC061);
#endif //DEEPBLUE_FLAG
  sds_index = SDnametoindex(sd_id,sds_name);
//printf("%d ",sds_index);
  sds_id = SDselect(sd_id,sds_index);
//printf("%d ",sds_id);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
//printf("%d\n",status);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  attr_index = SDfindattr(sds_id,"_FillValue");
  status = SDreadattr(sds_id,attr_index,&fill);
  status = SDendaccess(sds_id);
  for (ilo=0; ilo<N_MOD08_LONS; ilo++)
  for (ila=0; ila<N_MOD08_LATS; ila++) {
    iloc = (ila*N_MOD08_LONS)+ilo;
    if (ilo<180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo+180);
    if (ilo>=180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo-180);
    if (tmp[cloc] != fill) aod[iloc] = scale * (((float) tmp[cloc])-offset);
    else aod[iloc] = MOD08_BADFLAG;
  }
  FREE(tmp);

  tmp = (short *) calloc(N_MOD08_LATS*N_MOD08_LONS,sizeof(short));
	if (collection==51 || collection==61) strcpy(sds_name,cldpres_nameC051);
	else printf("no cloud pressure variable name for collections other than 5.1 and 6.1\n");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDreaddata(sds_id,start,NULL,edge,tmp);
  attr_index = SDfindattr(sds_id,"scale_factor");
  status = SDreadattr(sds_id,attr_index,&scale);
  attr_index = SDfindattr(sds_id,"add_offset");
  status = SDreadattr(sds_id,attr_index,&offset);
  attr_index = SDfindattr(sds_id,"_FillValue");
  status = SDreadattr(sds_id,attr_index,&fill);
  status = SDendaccess(sds_id);
  for (ilo=0; ilo<N_MOD08_LONS; ilo++)
  for (ila=0; ila<N_MOD08_LATS; ila++) {
    iloc = (ila*N_MOD08_LONS)+ilo;
    if (ilo<180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo+180);
    if (ilo>=180) cloc = ((N_MOD08_LATS-1-ila)*N_MOD08_LONS)+(ilo-180);
    if (tmp[cloc] != fill) cldtop_pres[iloc] = scale * (((float) tmp[cloc])-offset);
    else cldtop_pres[iloc] = MOD08_BADFLAG;
  }
  FREE(tmp);

  status = SDend(sd_id);
}






char get_MOD08_filename(const char *MOD08_path,const char *list_filename,
                        int year,int doy,char *MOD08_file)
{
  FILE *fp=NULL;
  char filelist[150], *s;
  char tmpfile[43];
  char yyyy[5], ddd[4];
  char outfile=0;
  int yr, day, endflag=0;
  
  sprintf(filelist,"%s%s",MOD08_path,list_filename);
  fp = fopen(filelist,"r");
//  while ((fgets(tmpfile,40,fp) != EOF) || !endflag ) {
//    if ( fgets(tmpfile,40,fp) != NULL ) {
//  printf("%s\n",tmpfile);
//printf("s=%c\n",s);
  while (!feof(fp) && !endflag) {
    if ( fgets(tmpfile,41,fp) != NULL ) {
      strncpy(yyyy,tmpfile+10,4);
      yyyy[4] = '\0';
      yr = atoi(yyyy);
      strncpy(ddd,tmpfile+14,3);
      ddd[3] = '\0';
      day = atoi(ddd);
      if ( (yr==year) && (day==doy) ) {
        strncpy(MOD08_file,tmpfile,39);
//        sprintf(MOD08_file,"%s",tmpfile);
        endflag=1;
      }
    }
  }
  fclose(fp);
  return endflag;
}





int calc_monthly_climatology(const char *MOD08_path,const char *list_filename,
                            int month,int year_start,int year_end,char *var_name,
                            float *lons,float *lats,float *climatology)
{
  char MOD08_file[50], isfile;
  int doy, month_cnt=0, var_num=-1;
  long iyr, yyyymmdd, i, j, iloc;
  float *cnt=NULL, *cldfrac=NULL;
  float *aod=NULL, *reff=NULL, *cldtau=NULL, *cldtmp=NULL, *lwp=NULL, *cldP=NULL;
  
  if (strcmp(var_name,"aod")==0) var_num=0;
  if (strcmp(var_name,"reff")==0) var_num=1;
  if (strcmp(var_name,"cld_tau")==0) var_num=2;
  if (strcmp(var_name,"lwp")==0) var_num=3;
  if (strcmp(var_name,"cld_temp")==0) var_num=4;
  if (strcmp(var_name,"cld_frac")==0) var_num=5;

  aod = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
  reff = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
  cldtau = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
  cldtmp = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
  cldP = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
  cldfrac = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
  lwp = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));
  cnt = (float *) malloc(N_MOD08_LONS*N_MOD08_LATS*sizeof(float));

  for (iyr=year_start; iyr<=year_end; iyr++) {
    yyyymmdd = (iyr*10000) + (month*100) + 1;
    doy = yyyymmdd_to_doy(yyyymmdd);
    isfile = get_MOD08_filename(MOD08_path,list_filename,(int) iyr,doy,MOD08_file);
    if (isfile) {
      month_cnt++;
      read_mod08_aerosol_cloud(MOD08_path,MOD08_file,lons,lats,reff,
                               cldtau,cldtmp,lwp,aod,cldfrac,cldP,51);
      for (i=0;i<N_MOD08_LONS;i++)
      for (j=0;j<N_MOD08_LATS;j++) {
        iloc = (i*N_MOD08_LATS)+j;
        
        switch (var_num)
        {
          case 0:
            if (aod[iloc] != MOD08_BADFLAG) {
              climatology[iloc] += aod[iloc];
              cnt[iloc]++;
            }
            break;
          case 1:
            if (reff[iloc] != MOD08_BADFLAG) {
              climatology[iloc] += reff[iloc];
              cnt[iloc]++;
            }
            break;
          case 2:
            if (cldtau[iloc] != MOD08_BADFLAG) {
              climatology[iloc] += cldtau[iloc];
              cnt[iloc]++;
            }
            break;
          case 3:
            if (lwp[iloc] != MOD08_BADFLAG) {
              climatology[iloc] += lwp[iloc];
              cnt[iloc]++;
            }
            break;
          case 4:
            if (cldtmp[iloc] != MOD08_BADFLAG) {
              climatology[iloc] += cldtmp[iloc];
              cnt[iloc]++;
            }
            break;
          case 5:
            if (cldfrac[iloc] != MOD08_BADFLAG) {
              climatology[iloc] += cldfrac[iloc];
              cnt[iloc]++;
            }
            break;
        }
      }
    }
  }
  
  for (i=0;i<N_MOD08_LONS;i++)
  for (j=0;j<N_MOD08_LATS;j++) {
    iloc = (i*N_MOD08_LATS)+j;
    if (cnt[iloc]>0) climatology[iloc] /= cnt[iloc];
    else climatology[iloc] = -9999.;
  }
  
  FREE(aod);
  FREE(reff);
  FREE(cldtau);
  FREE(cldtmp);
  FREE(cldP);
  FREE(cldfrac);
  FREE(lwp);
  
  return month_cnt;
}








