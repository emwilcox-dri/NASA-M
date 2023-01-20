/* general tools for specifying and averaging over a grid */

#define FREE(x) if(x) { free(x); x = NULL; }

void get_grid_specs(float lonmin,float latmin,float lonmax,float latmax,
              float dlon,float dlat,long *nlons,long *nlats);

void get_grid(float lonmin,float latmin,float dlon,float dlat,
              long nlons,long nlats,float *lons,float *lats,float *areas);

void get_grid_indices(float lat,float lon,long nlons,long nlats,
                      float dlat,float dlon,float latmin,
                      float lonmin,int *ilat,int *ilon);

void bin_field_tmp(void *var,float *lon,float *lat,float bad_flag,
              void *varbin,char vartype,long npix,float lonmin,float latmin,
              long nlons,long nlats,float dlon,float dlat,
              float *lons,float *lats);

void bin_field(void *var,float *lon,float *lat,float bad_flag,
              void *varbin,char vartype,long npix,float lonmin,float latmin,
              long nlons,long nlats,float dlon,float dlat,
              float *lons,float *lats);

void bin_field_cnt(void *var,float *lon,float *lat,float bad_flag,
              void *varbin,long *varcnt,char vartype,long npix,
              float lonmin,float latmin,long nlons,long nlats,
              float dlon,float dlat,float *lons,float *lats);

float calc_great_circle_distance(float lat1,float lon1,float lat2,float lon2);