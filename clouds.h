#define CLOUDS_BADFLAG -9999.0
#define N_MAX_HIST_BINS 250

#define FREE(x) if(x) { free(x); x = NULL; }

void calc_cloudarea_pix(long *cloud,float *pixarea,long nlats,
                        long nlons,long nclouds,float *cloudarea);

void calc_cloudstat_bin(long *cloud,long nlons,long nlats,long nclouds,
                        float *var,float bad_flag,float *cloudstat);

void calc_maxcloudstat_bin(long *cloud,long nlons,long nlats,long nclouds,
                             float *var,float bad_flag,float *maxcloudstat);

void calc_mincloudstat_bin(long *cloud,long nlons,long nlats,long nclouds,
                             float *var,float bad_flag,float *mincloudstat);

void calc_stat_at_mincloudstat(long *cloud,long nlons,long nlats,long nclouds,
                               float *var,float *minvar,float bad_flag,
                               float *stat_at_mincloudstat);

void find_stat_location_bin(long *cloud,float *stat,float *var,
                            float *grid_lons,float *grid_lats,
                            long nlons,long nlats,float *stat_lon,
                            float *stat_lat);

int get_histogram_nbins(float binmin,float binmax,float binsize);

void get_histogram_stats(float binmin,float binsize,int nbins,float *bins);

void get_cloudstat_histogram(long *cloud,unsigned int *histogram,
                             unsigned long *max_hist,float *var,
                             long nlons,long nlats,int nbins,float binmin,
                             float binmax,float binsize);

void labelCloudPixels(float *lon,float *lat,long npix, long nlons,long nlats,float dlon,
											float dlat,float lonmin,float latmin,long *cloud,
											long *ncloud_pix,long *cloud_pix,float *grid_lons,float *grid_lats);


