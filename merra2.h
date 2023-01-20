#define MERRA_BADFLAG   -9999.
#define P_TOP		0.01 //pressure (in hPa) at the top of the model atmosphere

// dimensions for coarse grid, time-averaged fields
#define N_MERRA_CLONS      288   //-179.375 to 179.375
#define N_MERRA_CLATS      144   //-89.375 to 89.375
#define N_MERRA_CTIMES       8		//minutes since beginning of day (0 to 1260)
#define N_MERRA_CLEVS       42		//pressure levels in hPa (1000 to 0.1)

// dimensions for native grid, instantaneous fields
#define N_MERRA_NLONS      576   //-180 to 179.375
#define N_MERRA_NLATS      361   //-90 to 90
#define N_MERRA_NTIMES       4		//minutes since beginning of day (0 to 1260)
#define N_MERRA_NLEVS       42		//pressure levels in hPa (1000 to 0.1) Changed to 42 from 72 when changed to the NP files
#define N_MERRA_DLAT			 0.5
#define N_MERRA_DLON		 0.625

#define N_MERRA_INST3_NTIMES 8		// while met files are 4 times daily, the aerosol files are 8 times daily

#define FREE(x) if(x) { free(x); x = NULL; }

//void read_merra_dim(char *path,char *file,char *dim_name,float64 *dim);

//void read_merra_float(char *path,char *file,char *var_name,float *var);
int read_merra_float(char *path,char *file,char *var_name,void *var,float *fill);

long get_merra_indices_3D(float lat,float lon,int hour,int nlons,int nlats,int nts,
												int *ila,int *ilo,int *it);

long get_merra_indices_4D(float lat,float lon,float pres,int hour,
		int nlons,int nlats,int nlevs,int nts,
		double *P_merra,int *ila,int *ilo,int *iz,
		int *it);
												
void calc_merra_native_P_profile(float *delP,float *P);

void calc_merra_native_z_profile(float *delP,float *Tv,float *P,float *z);

/********
void read_merra_lons(char *path,char *file,char *var_name,float *var);

void read_merra_lats(char *path,char *file,char *var_name,float *var);

void read_merra_levs(char *path,char *file,char *var_name,float *var);

void read_merra_times(char *path,char *file,char *var_name,float *var);

void read_merra_float(char *path,char *file,char *var_name,float *var);
*********/
