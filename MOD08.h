#define MOD08_BADFLAG   -9999.
#define N_MOD08_LONS      360
#define N_MOD08_LATS      180
#define NCHAR_FNAME_MOD08  39
//#define NCHAR_FNAME_MOD08  40 /*includes trailing character*/

#define FREE(x) if(x) { free(x); x = NULL; }

void read_mod08_aerosol_cloud(const char *modPath,char *MOD08file,float *lons,
                              float *lats,float *reff,
                              float *cld_tau,float *cldtop_temp,float *lwp,
                              float *aod,float *cld_frac,float *cldtop_pres,
                              int collection);

char get_MOD08_filename(const char *MOD08_path,const char *list_filename,
                        int year,int doy,char *MOD08_file);
                        
int calc_monthly_climatology(const char *MOD08_path,const char *list_filename,
                            int month,int year_start,int year_end,char *var_name,
                            float *lons,float *lats,float *climatology);