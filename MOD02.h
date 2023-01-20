#define N_L1B_1KM_BANDS			16
#define N_ALONG_L1B_1KM_MAX	2050
#define N_ACROSS_L1B_1KM_MAX	1354
#define MOD02_BADFLAG   -9999.

int readMOD021KM_IR(char *mod02Path,char *mod03Path,char *MOD02file,char *MOD03file,long *na_1km,
               float *lat_1km,float *lon_1km,float *IRemiss,float *IR_radiances,
               unsigned char *land_mask);
               
float modis_bright(float rad, int band);

float inv_planck_wvl(float wv1, float intensity);