#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include "merra2.h"
#include "thermo.h"

/************
Don't need this one because in netCDF-4 files the dimension values are stored as DATASETs.
*************/
/*************
void read_merra_dim(char *path,char *file,char *dim_name,size_t *dim)
{
  char fname[150];
  int sd_id, sds_id, sds_index, dim_id;
  int status;

	sprintf(fname,"%s%s",path,file);
	status = nc__open(fname,NC_NOWRITE,NC_SIZEHINT_DEFAULT,&sd_id);
	if (status == NC_NOERR) {
		status = nc_inq_dimid(sd_id,dim_name,&dim_id);
		nc_close(sd_id);
	} else printf("error opening %s\n",fname);

Old MERRA HDF code below here.
	sd_id = SDstart(fname, DFACC_READ);
	if (sd_id != -1) {
		sds_index = SDnametoindex(sd_id,dim_name);
		sds_id = SDselect(sd_id,sds_index);
		dim_id = SDgetdimid(sds_id,0);
		status = SDgetdimscale(dim_id,dim);
		status = SDendaccess(sds_id);
		status = SDend(sd_id);
	} else printf("error opening %s\n",fname);
}
****************/





//void read_merra_float(char *path,char *file,char *var_name,float *var)
int read_merra_float(char *path,char *file,char *var_name,void *var,float *fill)
{
  char fname[150];
  char tmp[65];
  int sd_id, sds_id, sds_index, dim_id, var_id;
  int status, start[5], edge[5], is_scale, is_offset, is_fill;
  int rank, dims[10], type, attributes, is_file;
  long i;
  long iloc;
  float add_offset, scale;

	sprintf(fname,"%s%s",path,file);
//	printf("reading %s\n",fname);
//	status = nc__open(fname,NC_NOWRITE,NC_SIZEHINT_DEFAULT,&sd_id);
	status = nc_open(fname,NC_NOWRITE,&sd_id);
//	printf("%d ",status);
	if (status == NC_NOERR) {
		is_file=1;
		status = nc_inq_varid(sd_id,var_name,&var_id);
		status = nc_get_var(sd_id,var_id,var);
//		printf("%d ",status);
		is_fill = nc_get_att(sd_id,var_id,"_FillValue",fill);
		is_offset = nc_get_att(sd_id,var_id,"add_offset",&add_offset);
		is_scale = nc_get_att(sd_id,var_id,"scale_factor",&scale);
		nc_close(sd_id);
		
		if (is_fill!=0) *fill=0;
		
	} else {
		printf("error opening %s\n",fname);
		is_file=0;
	}
	return is_file;
}



/**********
added this, but didn't use it because the 4D one returns the component indices, so I
can just exclude the index for the 4th dimension that I'm not interested in
**********/
long get_merra_indices_3D(float lat,float lon,int hour,int nlons,int nlats,int nts,
												int *ila,int *ilo,int *it)
{
	int i;
	long iloc;
	float lonmin, latmin;
	
//time index
	*it = (int) floor((float) hour/(24./nts));

//calculate the longitude index
	lonmin = -180;
	if (lon>180) lon-=360;
	*ilo = (int) floor((lon-lonmin)/(360./(float) nlons));
	if (*ilo>=nlons) *ilo=nlons;
	
//calculate the latitude index
	latmin = -90;
	*ila = (int) floor((lat-latmin)/(180./(float) nlats));
	if (*ila>=nlats) *ila=nlats;
	
	
	iloc = (*it*nlats*nlons) + (*ila*nlons) + *ilo;

	return iloc;
}




long get_merra_indices_4D(float lat,float lon,float pres,int hour,
												int nlons,int nlats,int nlevs,int nts,
												double *P_merra,int *ila,int *ilo,int *iz,
												int *it)
{
	int i;
	long iloc;
	float lonmin, latmin, dlat, dlon;
	
//time index
	*it = (int) floor((float) hour/(24./nts));

//find the pressure level index
	i=0;
	*iz=0;
	while (i<nlevs) {
		if ( fabs(pres-P_merra[i])<10. ) {
			*iz=i;
			i=nlevs;
		}
		i++;
	}

//calculate the longitude index
	lonmin = -180;
	if (lon>180) lon-=360;
	*ilo = (int) floor((lon-lonmin+(N_MERRA_DLON/2.))/N_MERRA_DLON);
	if (*ilo>=nlons) *ilo=nlons;
	
//calculate the latitude index
	latmin = -90;
	*ila = (int) floor((lat-latmin+(N_MERRA_DLAT/2.))/N_MERRA_DLAT);
	if (*ila>=nlats) *ila=nlats;
	
	
	iloc = (*it*nlevs*nlats*nlons) + (*iz*nlats*nlons) +
				(*ila*nlons) + *ilo;

	return iloc;
}








												
void calc_merra_native_P_profile(float *delP,float *P)
{
	int i;
	float ptop;
	
	ptop = P_TOP;
	for (i=0; i<N_MERRA_NLEVS; i++) {
		P[i] = ptop + (delP[i]/100.)/2.;
		ptop = ptop+(delP[i]/100.);
	}
}








/*************************
* NOTE: this is height above surface, not sea level
* Calculation is the hypsometric equation, however sign opposite
* because delP is positive even though pressure decreases with height.
*************************/
void calc_merra_native_z_profile(float *delP,float *Tv,float *P,float *z)
{
	int i;
	float delz, zbot=0;
	
	i=N_MERRA_NLEVS-1;
	while (i>=0) {
		delz = ((R*Tv[i])/(G*P[i]*100.))*delP[i];
//printf("%d %f %f %f %f\n",i,P[i]*100.,Tv[i],delP[i],delz);
		z[i] = zbot + (delz/2);
		zbot = zbot + delz;
		i--;
	}
}




/*****************************************
float calc_CAPE_emanuel(int nlevs,float *T,float *P,float *rh)
{
	FILE *fp=NULL, *fpo=NULL;
	char cmd[100], tmp[65], cape_str[12];
#ifdef SUMMER_MONSOON
	char *tmpfile="/Users/ewilcox/tmpprofile_sm.txt";
#else
	char *tmpfile="/Users/ewilcox/tmpprofile.txt";
#endif
	char *capefile="/Users/ewilcox/calcsound.out";
	int inttmp, i, ncolumns=64;
	float mr, cape=-9999., flttmp;
	
	fp = fopen(tmpfile,"w");
	fprintf(fp,"  N=%3d\n",nlevs);
	fprintf(fp,"   Pressure (mb)     Temperature (C)     Mixing Ratio (g/kg)\n");
	fprintf(fp,"   -------------     ---------------     -------------------\n");
	for (i=0; i<nlevs; i++) {
		mr = calc_mr(T[i],P[i],rh[i])*1000.;
		fprintf(fp,"     %9.3f        %9.3f            %9.3f\n",P[i]/100.,T[i]-T_ABS,mr);
	}
	fclose(fp);

//	system("calcsound");

#ifdef SUMMER_MONSOON
//	printf("reading %s\n",tmpfile);
	fflush(NULL);
	if ( (fpo = popen("calcsound_sm","r")) != NULL) {
#else
	if ( (fpo = popen("calcsound","r")) != NULL) {
#endif
		fgets(cape_str,12,fpo);
		cape = atof(cape_str);
		inttmp = pclose(fpo);
	} else printf("error running calcsound\n");

	return cape;	
}
*****************************************/



