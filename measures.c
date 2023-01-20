#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "measures.h"



void get_reff_profile(float *der,float *bt,int N,float btmin,
                      float btmax,int nbtbins,float *der_mean,
                      float *bins)
{
  int i, ibt;
  float dbt, *der_cnt=NULL;
	
  dbt = (btmax-btmin)/nbtbins;
  der_cnt = (float *) calloc(nbtbins,sizeof(float));
  
  for (i=0; i<nbtbins; i++) bins[i] = (i*dbt) + btmin + (dbt/2.);
  
  for (i=0; i<N; i++) {
    ibt = (int) floor((bt[i]-btmin)/dbt);
    if ( (ibt>=0) && (ibt<nbtbins) && (der[i]!=-9999) ) {
      der_mean[ibt] += der[i];
      der_cnt[ibt]++;
    }
  }
  
  for (i=0; i<nbtbins; i++) {
    if (der_cnt[i]>0) der_mean[i]/=der_cnt[i];
    else der_mean[i] = -9999;
  }
  free(der_cnt);
}




float glaciation_temp_diff(float *der,float *bt,int N)
{
  float gt = -9999;
  
  return gt;
}




float glaciation_temp_maxval(float *der,float *bt,int N)
{
  int i;
  float max_der=0;
  float gt = -9999;
  
  for (i=0; i<N; i++) {
    if (der[i]>max_der) gt = bt[i];
  }
  
  return gt;
}




void boxcar_smooth(float *array,int nbox,float *smoothed_array)
{
}




