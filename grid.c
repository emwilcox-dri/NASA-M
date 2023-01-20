
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"

#define EARTH_RADIUS  6370 /*in kilometers*/
#define PI            3.1415926



void get_grid_specs(float lonmin,float latmin,float lonmax,float latmax,
              float dlon,float dlat,long *nlons,long *nlats)
{
  float dellon, dellat;

  dellon = lonmax - lonmin;
  dellat = latmax - latmin;
  *nlons = floor((double) dellon/dlon);
  *nlats = floor((double) dellat/dlat);
}




void get_grid(float lonmin,float latmin,float dlon,float dlat,
              long nlons,long nlats,float *lons,float *lats,
              float *areas)
{
  int i;
  float lat1,lat2;
  float area;

  for (i=0; i<nlons; i++) lons[i] = (i*dlon)+lonmin+(dlon/2);
  for (i=0; i<nlats; i++) {
      lats[i] = (i*dlat)+latmin+(dlat/2);
//      areas[i] = EARTH_RADIUS * (dlon*PI/180.) *
//                 EARTH_RADIUS * cos(lats[i]*PI/180.) * (dlat*PI/180.);
      lat1 = (lats[i] - (dlat/2.))*(PI/180.);
      lat2 = (lats[i] + (dlat/2.))*(PI/180.);
      areas[i] = (PI/180.)*EARTH_RADIUS*EARTH_RADIUS *
                fabs(sin(lat1)-sin(lat2))*dlon;
  }
}




void get_grid_indices(float lat,float lon,long nlons,long nlats,
                      float dlat,float dlon,float latmin,
                      float lonmin,int *ilat,int *ilon)
{
  *ilat = (int) floor((lat-latmin)/dlat);
  *ilon = (int) floor((lon-lonmin)/dlat);
  if ( (*ilat<0) || (*ilat>=nlats) ) *ilat=-9999;
  if ( (*ilon<0) || (*ilon>=nlons) ) *ilon=-9999;
}




void bin_field_tmp(void *var,float *lon,float *lat,float bad_flag,
              void *varbin,char vartype,long npix,float lonmin,float latmin,
              long nlons,long nlats,float dlon,float dlat,
              float *lons,float *lats)
{
  long i,j,ilon,ilat,iloc,cnt=0;
  long *varcnt=NULL, tmpcnt=0;
  float flttmp;
  double *varsum=NULL;
  varcnt = (long *) malloc(nlons*nlats*sizeof(long));
  if (varcnt == NULL) printf("varcnt not allocated\n");
  varsum = (double *) malloc(nlons*nlats*sizeof(double));
  if (varsum == NULL) printf("varsum not allocated\n");
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
        iloc = (i*nlons)+j;
        varsum[iloc] = 0;
        varcnt[iloc] = 0;
  }

printf("npix=%ld bad_flag=%f lonmin=%f latmin=%f\n",
          npix, bad_flag, lonmin,latmin);
          
      
  switch (vartype) {


  case 'f': {    
  for (i=0; i<npix; i++) {
if (*((float *) var+i) != bad_flag) tmpcnt++;
//if (cnt<100) printf("[%f %f %f] ",lon[i],lat[i],*((float *) var+i));
//cnt++;
    ilon = (long) floor((lon[i]-lonmin)/dlon);
    ilat = (long) floor((lat[i]-latmin)/dlat);
if (lon[i]>=lonmin && lat[i]>=latmin)
//printf("[%f %f %ld %ld] ",lon[i],lat[i],ilon,ilat);
    if ( (*((float *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((float *) var+i);
//printf("[%f %f] ",*((float *) var+i),varsum[iloc]);
      ++varcnt[iloc];
    }
//else printf("[%ld %ld %f] ",ilon,ilat,*((float *) var+i));
  }
  
  printf("goodpix=%ld\n",tmpcnt);
  
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((float *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((float *) varbin+iloc) = bad_flag;
      if (cnt<100) printf("%f ",*((float *) varbin+iloc));
      cnt++;

  }
  break;
  }




  case 'd': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((double *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
if (cnt<100) printf("[%f %f %lf] ",lon[i],lat[i],*((double *) var+i));
cnt++;
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((double *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((double *) varbin+iloc) = varsum[iloc]/varcnt[iloc];      else *((double *) varbin+iloc) = bad_flag;
  }
  break;
  }

  }
  FREE(varcnt);
  FREE(varsum);
}









void bin_field(void *var,float *lon,float *lat,float bad_flag,
              void *varbin,char vartype,long npix,float lonmin,float latmin,
              long nlons,long nlats,float dlon,float dlat,
              float *lons,float *lats)
{
  long i,j,ilon,ilat,iloc;
  long *varcnt=NULL;
  double *varsum=NULL;

  varcnt = (long *) malloc(nlons*nlats*sizeof(long));
  if (varcnt == NULL) printf("varcnt not allocated\n");
  varsum = (double *) malloc(nlons*nlats*sizeof(double));
  if (varsum == NULL) printf("varsum not allocated\n");
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
        iloc = (i*nlons)+j;
        varsum[iloc] = 0;
        varcnt[iloc] = 0;
  }

  switch (vartype) {


  case 'f': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((float *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {

      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((float *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((float *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((float *) varbin+iloc) = bad_flag;

  }
  break;
  }


  case 't': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);

    if ( (*((float *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((float *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((float *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((float *) varbin+iloc) = bad_flag;
  }
  break;
  }


  case 'd': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((double *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((double *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((double *) varbin+iloc) = varsum[iloc]/varcnt[iloc];      else *((double *) varbin+iloc) = bad_flag;
  }
  break;
  }


  case 's': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((short *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((short *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((short *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((short *) varbin+iloc) = bad_flag;
  }
  break;
  }

  }
  FREE(varcnt);
  FREE(varsum);
}




void bin_field_cnt(void *var,float *lon,float *lat,float bad_flag,
              void *varbin,long *varcnt,char vartype,long npix,
              float lonmin,float latmin,long nlons,long nlats,
              float dlon,float dlat,float *lons,float *lats)
{
  long i,j,ilon,ilat,iloc;
  double *varsum=NULL;

  varsum = (double *) malloc(nlons*nlats*sizeof(double));
  if (varsum == NULL) printf("varsum not allocated\n");
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
        iloc = (i*nlons)+j;
        varsum[iloc] = 0;
        varcnt[iloc] = 0;
  }

  switch (vartype) {

  case 'f': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((float *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {

      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((float *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((float *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((float *) varbin+iloc) = bad_flag;

  }
  break;
  }


  case 't': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);

    if ( (*((float *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((float *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((float *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((float *) varbin+iloc) = bad_flag;
  }
  break;
  }


  case 'd': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((double *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((double *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((double *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((double *) varbin+iloc) = bad_flag;
  }
  break;
  }


  case 's': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((short *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((short *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
//      if (varcnt[iloc] > 0) *((short *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      if (varcnt[iloc] > 0) *((short *) varbin+iloc) = (short) round(varsum[iloc]/varcnt[iloc]);
      else *((short *) varbin+iloc) = bad_flag;
//      if (varcnt[iloc] > 0) printf("[%lf %ld %d] ",varsum[iloc],varcnt[iloc],*((short *) varbin+iloc) );
  }
  break;
  }

  }
  FREE(varsum);
}




void bin_field_ind(void *var,float *lon,float *lat,float bad_flag,
              void *varbin,long *ind,char vartype,long npix,
              float lonmin,float latmin,long nlons,long nlats,
              float dlon,float dlat,float *lons,float *lats)
{
  long i,j,ilon,ilat,iloc;
  long *varcnt=NULL;
  double *varsum=NULL;

  varcnt = (long *) malloc(nlons*nlats*sizeof(long));
  if (varcnt == NULL) printf("varcnt not allocated\n");
  varsum = (double *) malloc(nlons*nlats*sizeof(double));
  if (varsum == NULL) printf("varsum not allocated\n");
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
        iloc = (i*nlons)+j;
        varsum[iloc] = 0;
        varcnt[iloc] = 0;
  }
  for (i=0; i<npix; i++) ind[i] = (long) bad_flag;

  switch (vartype) {

  case 'f': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((float *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {

      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((float *) var+i);
      ind[i] = iloc;
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((float *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((float *) varbin+iloc) = bad_flag;

  }
  break;
  }



  case 't': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);

    if ( (*((float *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((float *) var+i);
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((float *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((float *) varbin+iloc) = bad_flag;
  }
  break;
  }


  case 'd': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((double *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((double *) var+i);
      ind[i] = iloc;
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((double *) varbin+iloc) = varsum[iloc]/varcnt[iloc];      else *((double *) varbin+iloc) = bad_flag;
  }
  break;
  }


  case 's': {
  for (i=0; i<npix; i++) {
    ilon = (int) floor((lon[i]-lonmin)/dlon);
    ilat = (int) floor((lat[i]-latmin)/dlat);
    if ( (*((short *) var+i) != bad_flag) &&
         (ilon >= 0) && (ilat >= 0) &&
         (ilon < nlons) && (ilat < nlats) ) {
      iloc = (ilat*nlons)+ilon;
      varsum[iloc] += *((short *) var+i);
      ind[i] = iloc;
      ++varcnt[iloc];
    }
  }
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      if (varcnt[iloc] > 0) *((short *) varbin+iloc) = varsum[iloc]/varcnt[iloc];
      else *((short *) varbin+iloc) = bad_flag;
  }
  break;
  }

  }
  FREE(varcnt);
  FREE(varsum);
}









/************
* Calculates the great circle distance between two point.
* lats and lons are input in degrees.
* distance returned in kilometers.
************/
float calc_great_circle_distance(float lat1,float lon1,float lat2,float lon2)
{
  float dist, dlon, dlat;
  
  dlon = fabs(lon2-lon1) * (PI/180.);
  dlat = fabs(lat2-lat1) * (PI/180.);
  lat1 = lat1*(PI/180.);
  lat2 = lat2*(PI/180.);
  lon1 = lon1*(PI/180.);
  lon2 = lon2*(PI/180.);
  dist = EARTH_RADIUS*2.*asin( sqrt( (sin(dlat/2.)*sin(dlat/2.)) +
                          (cos(lat1)*cos(lat2)*sin(dlon/2.)*sin(dlon/2.)) ) );
/*  printf("%f %f %f %f %f %f %f\n",dlon*180./PI,dlat*180./PI,dist,
  				(sin(dlat/2.)*sin(dlat/2.)),
  				(cos(lat1)*cos(lat2)*sin(dlon/2.)*sin(dlon/2.)),
  				sqrt( (sin(dlat/2.)*sin(dlat/2.)) +
                          (cos(lat1)*cos(lat2)*sin(dlon/2.)*sin(dlon/2.)) ),
          acos( sqrt( (sin(dlat/2.)*sin(dlat/2.)) +
                          (cos(lat1)*cos(lat2)*sin(dlon/2.)*sin(dlon/2.)) ) ) ); */
  
  return dist;
}





