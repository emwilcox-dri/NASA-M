void get_reff_profile(float *der,float *bt,int N,float btmin,
                      float btmax,int nbtbins,float *der_mean,
                      float *bins);
                      
float glaciation_temp_diff(float *der,float *bt,int N);

float glaciation_temp_maxval(float *der,float *bt,int N);

void boxcar_smooth(float *array,int nbox,float *smoothed_array);