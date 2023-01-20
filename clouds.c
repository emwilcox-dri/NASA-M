#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "clouds.h"

/*************************
calc_cloudarea_pix() takes imager detected clouds and calculates the
area of each cloud.
This routine is identical to calc_cloudarea_bin, except allows for pixels
whose size vary with latitude.
*************************/
void calc_cloudarea_pix(long *cloud,float *pixarea,long nlats,
                        long nlons,long nclouds,float *cloudarea)
{
  long c, i, j, iloc;
  long lab;
  double area;

  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++)
      {
        iloc = (i*nlons)+j;
        lab = cloud[iloc];
        if (lab > 0) cloudarea[lab-1] += pixarea[i-1];
      }
  for (c=0; c<nclouds; c++)
    if (cloudarea[c] == 0) cloudarea[c] = CLOUDS_BADFLAG;
}





void calc_cloudstat_bin(long *cloud,long nlons,long nlats,long nclouds,
                        float *var,float bad_flag,float *cloudstat)

{
  long c, i, j, iloc;

  long *totcnt=NULL;

  totcnt = (long *) calloc(nclouds,sizeof(long));
  if (totcnt == NULL) printf("ERROR: cannot allocate totcnt in calc_cloudstat_bin\n");
//  for (i=0; i<nlons; i++)
//    for (j=0; j<nlats; j++) {
//      iloc = (i*nlats)+j;
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      c = cloud[iloc];
      if ( (c > 0) && (var[iloc] != bad_flag) ) {
        totcnt[c-1]++;
        cloudstat[c-1] += var[iloc];
      }
  }
  for (c=0; c<nclouds; c++)
    if (totcnt[c] > 0) cloudstat[c] /= totcnt[c];
    else cloudstat[c] = bad_flag;
  FREE(totcnt);
}





/***********************
calc_maxcloudstat_bin() same as above, but does the max value
***********************/
void calc_maxcloudstat_bin(long *cloud,long nlons,long nlats,long nclouds,
                             float *var,float bad_flag,float *maxcloudstat)
{
  long c, i, j, iloc;
  long *totcnt=NULL;

  totcnt = (long *) calloc(nclouds,sizeof(long));
  if (totcnt == NULL) printf("ERROR: cannot allocate totcnt in calc_cloudstat_bin\n");
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
        iloc = (i*nlons)+j;
        c = cloud[iloc];
        if ( (c > 0) && (var[iloc] != bad_flag) ) {
            totcnt[c-1]++;
            if (totcnt[c-1] == 1) maxcloudstat[c-1] = var[iloc];
            else if (var[iloc] > maxcloudstat[c-1])
                maxcloudstat[c-1] = var[iloc];
        }
  }
 for (c=0; c<nclouds; c++)
    if (totcnt[c] == 0) maxcloudstat[c] = bad_flag;
  FREE(totcnt);
}






/***********************
calc_mincloudstat_bin() same as above, but does the max value
***********************/
void calc_mincloudstat_bin(long *cloud,long nlons,long nlats,long nclouds,
                             float *var,float bad_flag,float *mincloudstat)
{
  long c, i, j, iloc;
  long *totcnt;

  totcnt = (long *) calloc(nclouds,sizeof(long));
  if (totcnt == NULL) printf("ERROR: cannot allocate totcnt in calc_mincloudstat_bin\n");
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
        iloc = (i*nlons)+j;
        c = cloud[iloc];
        if ( (c > 0) && (var[iloc] != bad_flag) ) {
            totcnt[c-1]++;
            if (totcnt[c-1] == 1) mincloudstat[c-1] = var[iloc];
            else if (var[iloc] < mincloudstat[c-1])
                mincloudstat[c-1] = var[iloc];
        }
  }
 for (c=0; c<nclouds; c++)
    if (totcnt[c] == 0) mincloudstat[c] = bad_flag;
  FREE(totcnt);
}



/***********************
calc_stat_at_mincloudstat() finds the value of a quantity at the location
of the minimum value of another quantity
***********************/
void calc_stat_at_mincloudstat(long *cloud,long nlons,long nlats,long nclouds,
                             float *var,float *minvar,float bad_flag,
                             float *stat_at_mincloudstat)
{
  long c, i, j, iloc;
  long *totcnt;
  float *mincloudstat;

  totcnt = (long *) calloc(nclouds,sizeof(long));
  mincloudstat = (float *) calloc(nclouds,sizeof(float));
  if (totcnt == NULL) printf("ERROR: cannot allocate totcnt in calc_mincloudstat_bin\n");
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
        iloc = (i*nlons)+j;
        c = cloud[iloc];
        if ( (c > 0) && (minvar[iloc] != bad_flag) ) {
            totcnt[c-1]++;
            if (totcnt[c-1] == 1) {
              mincloudstat[c-1] = minvar[iloc];
              stat_at_mincloudstat[c-1] = var[iloc];
            } else {
              if (var[iloc] < mincloudstat[c-1]) {
                mincloudstat[c-1] = minvar[iloc];
                stat_at_mincloudstat[c-1] = var[iloc];
              }
            }
        }
  }
 for (c=0; c<nclouds; c++)
    if (totcnt[c] == 0) stat_at_mincloudstat[c] = bad_flag;
  FREE(totcnt);
  FREE(mincloudstat);
}



/*********************
Find the lat and lon of the grid cell containing a particular stat for a
particular cloud
********************/
void find_stat_location_bin(long *cloud,float *stat,float *var,
                            float *grid_lons,float *grid_lats,
                            long nlons,long nlats,float *stat_lon,
                            float *stat_lat)
{
  long c,i,j,iloc;

  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      c = cloud[iloc];
      if ( (c>0) && (stat[c-1]!=-9999.) && (var[iloc]==stat[c-1]) ) {
        stat_lon[c-1] = grid_lons[j];
        stat_lat[c-1] = grid_lats[i];
      }
    }
}



int get_histogram_nbins(float binmin,float binmax,float binsize)
{
  int nbins;

  nbins = (int) floor( (binmax-binmin)/binsize );
  
  return nbins;
}


void get_histogram_stats(float binmin,float binsize,int nbins,float *bins)
{
  long i;

  for (i=0; i<nbins; i++) bins[i] = (i*binsize) + binmin + (binsize/2.);
}



void get_cloudstat_histogram(long *cloud,unsigned int *histogram,
                             unsigned long *hist_count,float *var,
                             long nlons,long nlats,int nbins,float binmin,
                             float binmax,float binsize)
{
  int ibin;
  long c,i,j,iloc,hloc;
  
  for (i=0; i<nlats; i++)
    for (j=0; j<nlons; j++) {
      iloc = (i*nlons)+j;
      c = cloud[iloc];
      if ( (c>0) && (var[iloc]>=binmin) && (var[iloc]<binmax) ) {
        ibin = (int) floor( (var[iloc]-binmin)/binsize );
//        histogram[c-1][ibin]++;
        hloc = ((c-1)*nbins)+ibin;
        histogram[hloc]++;
        hist_count[c-1]++;
      }
  }
}




/*******************
* label pixels of a high-resolution dataset with the cloud index from a coarser-grid
* set of clouds found with DAS
* e.g. label the 1km MODIS L1B brightness temperatures with the corresponding cloud index
* from a gridded brightness temperature passed through DAS
********************/
void labelCloudPixels(float *lon,float *lat,long npix, long nlons,long nlats,float dlon,
											float dlat,float lonmin,float latmin,long *cloud,
											long *ncloud_pix,long *cloud_pix,float *grid_lons,float *grid_lats)
{
	long i, iloc, j;
	int ilat, ilon;

	for (i=0; i<npix; i++) cloud_pix[i]=0;

	for (i=0; i<npix; i++) {
//		printf("%ld ",i);
//		fflush(NULL);
//		printf("%f %f %ld %ld %f %f %f %f ",lat[i],lon[i],nlons,nlats,dlat,dlon,latmin,lonmin);
//		fflush(NULL);
		get_grid_indices(lat[i],lon[i],nlons,nlats,dlat,dlon,latmin,lonmin,&ilat,&ilon);
//		printf("[%d %d]\n",ilat,ilon);
//		fflush(NULL);
		if ( (ilat!=-9999) && (ilon!=-9999) ) {
			iloc = (ilat*nlons)+ilon;
//			if (cloud[iloc]>0) printf("[%f %f %f %f %ld]  ",lat[i],grid_lats[ilat],lon[i],grid_lons[ilon],cloud[iloc]);
			cloud_pix[i] = cloud[iloc];
			if (cloud[iloc]>0) ncloud_pix[cloud[iloc]-1]++;
//if (cloud[iloc]==83) printf("%ld ",ncloud_pix[cloud[iloc]-1]);
		}
	}
/*j=0;
for (i=0; i<npix; i++) {
if (cloud_pix[i]==83) j++;
}
printf("# of cloud_pix=83 = %ld\n",j);*/

}






void get_reff_irbt_profile(float *der,float *bt,int N,float btmin,
                      float btmax,float nbtbins,float *der_mean,
                      float *der_dev_mean,float *bins)
{
  int i, ibt;
  float dbt, *der_cnt=NULL;
  
  dbt = (btmax-btmin)/nbtbins;
  der_cnt = (float *) malloc(nbtbins*sizeof(float));
  
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
  
  for (i=0; i<N; i++) {
    ibt = (int) floor((bt[i]-btmin)/dbt);
  	if ( (ibt>=0) && (ibt<nbtbins) && (der[i]!=-9999) ) {
  		der_dev_mean[ibt] += fabs(der[i]-der_mean[i]);
//  		finish calculation of mean deviation
		}
	}
  
  free(der_cnt);
}




