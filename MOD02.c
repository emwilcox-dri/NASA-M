#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf.h>
#include <mfhdf.h>
#include "MOD06.h"
#include "MOD02.h"

/**********
* uncertainty(%) = specified_uncertainty*exp(uncertainty_index/scaling_factor)
**********/
int readMOD021KM_IR(char *mod02Path,char *mod03Path,char *MOD02file,char *MOD03file,long *na_1km,
               float *lat_1km,float *lon_1km,float *IR_radiances,float *IR_raduncert,
               unsigned char *land_mask)
{
	char fname[150];
	char sds_name[128];
	unsigned char fill1, *IRuncert=NULL;
	unsigned short fill2, *IRemiss=NULL;
	int sd_id, sds_id, sds_index, attr_index;
  int status, start[3], edge[3], nbands;
  int rank, dims[3], type, attributes, ix, ia, ib, na;
  long iloc;
  float scales[16], offsets[16];
  
  IRemiss = (unsigned short *) malloc(N_ACROSS_L1B_1KM_MAX*N_ALONG_L1B_1KM_MAX*N_L1B_1KM_BANDS*sizeof(unsigned short));
  IRuncert = (unsigned char *) malloc(N_ACROSS_L1B_1KM_MAX*N_ALONG_L1B_1KM_MAX*N_L1B_1KM_BANDS*sizeof(unsigned char));

	sprintf(fname,"%s%s",mod02Path,MOD02file);
  printf("reading %s\n",fname);
  sd_id = SDstart(fname, DFACC_READ);
//printf("%d ",sd_id);
  if (sd_id==-1) {
  	printf("cannot open %s\n",fname);
  	exit(0);
  }
  
  sprintf(sds_name,"EV_1KM_Emissive");
  sds_index = SDnametoindex(sd_id,sds_name);
//printf("%d ",sds_index);
  sds_id = SDselect(sd_id,sds_index);
//printf("%d ",sds_id);
  status = SDgetinfo(sds_id,sds_name,&rank,dims,&type,&attributes);
//printf("%d ",status);
  if ( (dims[1]>N_ALONG_L1B_1KM_MAX) || (dims[2]>N_ACROSS_L1B_1KM_MAX) ) {
  	printf("MODIS L1B 1km dimensions = [%d %d] - too large\n",dims[1],dims[2]);
    status = SDendaccess(sds_id);
    status = SDend(sd_id);
    exit(0);
  }
  start[0] = start[1] = start[2] = 0;
  edge[0] = dims[0];
  edge[1] = dims[1];
  edge[2] = dims[2];
  nbands = dims[0];
  status = SDgetfillvalue(sds_id,&fill2);
  status = SDreaddata(sds_id,start,NULL,edge,IRemiss);
  attr_index = SDfindattr(sds_id,"radiance_scales");
  status = SDreadattr(sds_id,attr_index,scales);
  attr_index = SDfindattr(sds_id,"radiance_offsets");
  status = SDreadattr(sds_id,attr_index,offsets);
  status = SDendaccess(sds_id);
  
  na=dims[1];
  for (ix=0; ix<N_ACROSS_L1B_1KM_MAX; ix++)
  for (ia=0; ia<na; ia++)
  for (ib=0; ib<nbands; ib++) {
//    iloc = (ia*N_ACROSS_L1B_1KM_MAX*N_IR_BANDS) + (ix*N_IR_BANDS) + ib;
    iloc = (ib*N_ACROSS_L1B_1KM_MAX*na) + (ia*N_ACROSS_L1B_1KM_MAX) + ix;
    if (IRemiss[iloc] != fill2) IR_radiances[iloc] = scales[ib]*((float) IRemiss[iloc]-offsets[ib]);
    else IR_radiances[iloc] = MOD02_BADFLAG;
  }
  
  sprintf(sds_name,"EV_1KM_Emissive_Uncert_Indexes");
  sds_index = SDnametoindex(sd_id,sds_name);
  sds_id = SDselect(sd_id,sds_index);
  status = SDgetfillvalue(sds_id,&fill1);
  status = SDreaddata(sds_id,start,NULL,edge,IRuncert);
  attr_index = SDfindattr(sds_id,"scaling_factor");
  status = SDreadattr(sds_id,attr_index,scales);
  attr_index = SDfindattr(sds_id,"specified_uncertainty");
  status = SDreadattr(sds_id,attr_index,offsets);
  status = SDendaccess(sds_id);
  
//  for (ib=0; ib<N_L1B_1KM_BANDS; ib++) printf("[%f %f] ",scales[ib],offsets[ib]);
  
  for (ix=0; ix<N_ACROSS_L1B_1KM_MAX; ix++)
  for (ia=0; ia<na; ia++)
  for (ib=0; ib<nbands; ib++) {
//    iloc = (ia*N_ACROSS_L1B_1KM_MAX*N_IR_BANDS) + (ix*N_IR_BANDS) + ib;
		iloc = (ib*N_ACROSS_L1B_1KM_MAX*na) + (ia*N_ACROSS_L1B_1KM_MAX) + ix;
		if (IRuncert[iloc] != fill1) IR_raduncert[iloc] = offsets[ib]*exp(((float) IRuncert[iloc])/scales[ib]);
		else IR_raduncert[iloc] = MOD02_BADFLAG;
//    if ( (ib==10) && (IR_radiances[iloc]<0) ) printf("[%d %f %f]  ",IRemiss[iloc],IR_radiances[iloc],IR_raduncert[iloc]);
//    if (IR_radiances[iloc]<0) printf("[%d %f %d]  ",IRemiss[iloc],IR_radiances[iloc],IRuncert[iloc]);
  }
  
  readMOD03(mod03Path,MOD03file,na_1km,lat_1km,lon_1km,land_mask);

  status = SDend(sd_id);

	free(IRemiss);
	free(IRuncert);

	return nbands;
}




/***********************
Compute brightness temperature for an EOS-AM MODIS infrared band. (Works for Aqua?)
 Adapted from modis_bright.pro, received from Tianle Yuan in Jan 2019
 originally written by Liam Gumley, SSEC, U. Wisc.
 and adapted for Matlab by Shaima Nasiri in 2000

VIS/NIR and IR channel number breakdown:  
  8,9,10,11,12,13lo,13hi,14lo,14hi,15,16,17,18,19,26
  20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36

Note: channel 26 is absent with zeros as placeholders
************************/
float modis_bright(float rad, int band)
{
	float bt;
//- Effective central wavenumber (inverse centimenters)
	float cwn[17] = {
		2.641775E+03, 2.505277E+03, 2.518028E+03, 2.465428E+03,
		2.235815E+03, 2.200346E+03, 0., 1.477967E+03,
		1.362737E+03, 1.173190E+03, 1.027715E+03, 9.080884E+02,
		8.315399E+02, 7.483394E+02, 7.308963E+02, 7.188681E+02,
		7.045367E+02};

//- Temperature correction slope (no units)
	float tcs[17] = {
		9.993411E-01, 9.998646E-01, 9.998584E-01, 9.998682E-01,
		9.998819E-01, 9.998845E-01, 0., 9.994877E-01,
		9.994918E-01, 9.995495E-01, 9.997398E-01, 9.995608E-01,
		9.997256E-01, 9.999160E-01, 9.999167E-01, 9.999191E-01,
		9.999281E-01};

//- Temperature correction intercept (Kelvin)
	float tci[17] = {
		4.770532E-01, 9.262664E-02, 9.757996E-02, 8.929242E-02,
		7.310901E-02, 7.060415E-02, 0., 2.204921E-01,
		2.046087E-01, 1.599191E-01, 8.253401E-02, 1.302699E-01,
		7.181833E-02, 1.972608E-02, 1.913568E-02, 1.817817E-02,
		1.583042E-02};
		
	bt = ( inv_planck_wvl(1.0e+4/cwn[band-20], rad) - tci[band-20] )
					/ tcs[band-20];
	
	return bt;
}



/***********************
Compute brightness temperature for an EOS-AM MODIS infrared band. (Works for Aqua?)
 Adapted from modis_bright.pro, received from Tianle Yuan in Jan 2019
 originally written by Liam Gumley, SSEC, U. Wisc.
 and adapted for Matlab by Shaima Nasiri in 2000
 
	//input:
//wvl:       in micron
//intensity: in W/m^2/micron
************************/
float inv_planck_wvl(float wvl, float intensity)
{
	float bt, a1, a2, x1, x2;

//h = 6.62606957E-34  ;in J*s
//k = 1.3806488E-23   ;in J/K
//c = 2.99792458E14   ;in micron/s

	a1 = 1.19104283E-04;    //a1 = 2.0*h*c*c
	a2 = 1.43877705E+04;    //a2 = h*c/k
	//x1 = 2.0*h*c*c/wvl/wvl/wvl/wvl/wvl;
	x1 = a1/wvl/wvl/wvl/wvl/wvl;
	x2 = (float) log(x1/(intensity*1.0E-12) + 1.0);
	//bt = h*c/k/wvl/x2; 
	bt = a2/wvl/x2;
	
	return bt;
}



