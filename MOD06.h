#define MOD06_BADFLAG   -9999.
#define MOD06_CLRFLAG   -7777.
#define MOD06_SCRNFLAG  -8888.
#define N_ACROSS_1KM     1354
#define N_ACROSS_5KM      270
#define N_ALONG_1KM_MAX  2050
#define N_ALONG_5KM_MAX   410
#define N_MASK_TYPES        9
#define N_MASK_BYTES        2
#define N_QA1KM_TYPES      10
#define N_QA1KM_BYTES       5
#define NCHAR_FNAME_MOD06  44
#define NCHAR_FNAME_MOD03  41
#define N_IR_BANDS          7

/***********
const char *reff_name_C004 = "Effective_Particle_Radius";
const char *reff_name_C005 = "Cloud_Effective_Radius";
const char *reff_1621_name_C005 = "Cloud_Effective_Radius_1621";
const char *lwp_name_C004 = "Water_Path";
const char *lwp_name_C005 = "Cloud_Water_Path";
const char *lwp_1621_name_C005 = "Cloud_Water_Path_1621";
***********/

      /*  cloud mask indices:                                                 */
      /*  cloud_mask(*,*,0) = cloud determined(1)/undetermined(0)             */
      /*                 1  = cloudy(0)/66% conf. clear(1)/95% conf. clear(2) */
      /*                      /99% conf. clear(3)                             */
      /*                 2  = night(0)/day(1)                                 */
      /*                 3  = sun glint yes(0)/no(1)                          */
      /*                 4  = surface snow or ice yes(0)/no(1)                */
      /*                 5  = surface ocean(0)/coast(1)/desert(2)/land(3)     */
      /*                 6  = heavy aerosol yes(0)/no(1)                      */
      /*                 7  = thin cirrus yes(0)/no(1)                        */
      /*                 8  = shadow yes(0)/no(1)                             */

/* QA 1km indices:                                                            */
/*     qa_1km(*,*,0) = cloud tau QA not useful(0)/useful(1)                   */
/*                1  = cloud tau conf. 0-3                                    */
/*                2  = cloud tau bounds t<=100(0)/100<t<150(1)/               */
/*                     t>150(2)/surf. reflectance too large(3)                */
/*                3  = reff QA not useful(0)/useful(1)                        */
/*                4  = reff conf. 0-3                                         */
/*                5  = LWP QA not useful(0)/useful(1)                         */
/*                6  = LWP conf. 0-3                                          */
/*                7  = cloud phase (SWIR) alg. not run(0)/                    */
/*                       clear(1)/water(2)/ice(3)/                            */
/*                       mixed or undetermined(4)                             */
/*                8  = cloud phase used cloud mask undetermined(0)/           */
/*                       decision tree 'stop'(1)/liquid(2)/ice(3)             */
/*                       undetermined(4)/mixed(5)                             */
/*                9  = retrieval outcome not attempted or unsuccessful(0)/    */
/*                       successful(1)                                        */

#define PI                  3.1415926
#define RE               6378.1

#define DIST_THRESH    4.
#define NASTEP        15   /*number of steps along track */
#define NXSTEP         5   /*number of steps across track*/
/*#define DIST_THRESH    5.*/
/*#define NASTEP        18*/
/*#define NXSTEP         6*/
/*#define DIST_THRESH    6.*/
/*#define NASTEP        21*/
/*#define NXSTEP         7*/
/*#define DIST_THRESH    8.*/
/*#define NASTEP        30*/
/*#define NXSTEP        10*/
/*#define DIST_THRESH    10.*/
/*#define NASTEP        36*/
/*#define NXSTEP        12*/

#define FREE(x) if(x) { free(x); x = NULL; }

/*  land_mask values                      */
/*      0:  shallow ocean                 */
/*      1:  land                          */
/*      2:  ocean/lake shorline           */
/*      3:  shallow inland water          */
/*      4:  ephemeral water               */
/*      5:  deep inland water             */
/*      6:  moderate or continental ocean */
/*      7:  deep ocean                    */

/* cloud thermodynamic phase values       */
/*      clear = 0                         */
/*      water = 1 or 5                    */
/*      ice = 2 or 4                      */
/*      mixed = 3                         */
/*      uncertain = 6                     */

/* IR Brightness_Temperature variable contains following 7 bands */
/*    29 = 8.4 - 8.7 um                                          */
/*    31 = 10.78 - 11.28                                         */
/*    32 = 11.77 - 12.27                                         */
/*    33 = 13.185 - 13.485                                       */
/*    34 = 13.485 - 13.785                                       */
/*    35 = 13.785 - 14.085                                       */
/*    36 = 14.085 - 14.385                                       */

/* data arrays are arranged iloc=(ia*N_ACROSS_5KM)+ix */
/* where ia is the along-track index and ix is the across-track index */
void readMOD06(char *modPath,char *mod03path,char *MOD06file,char *MOD03file,
               int collection,long *na_1km,long *na_5km,
               float *lat_5km,float *lon_5km,float *lat_1km,
               float *lon_1km,float *fraction,float *reff,
               float *cld_tau,float *cldtop_temp,float *cldtop_pres,
               signed char *cldtop_phase,float *lwp,unsigned char *land_mask);

void read_1621_and_reff_diff(char *modPath,char *MOD06file,long *na_1km,
                             float *reff_1621,float *tau_1621,
                             float *lwp_1621,float *reff_diff);

void readMOD03(char *modPath,char *MOD03file,long *na_1km,
               float *lat_1km,float *lon_1km,unsigned char *land_mask);

void readMOD03_sensor_angles(char *modPath,char *MOD03file,
               float *sensor_zenith,float *sensor_azimuth);

void readMOD06_sensor_angles(char *modPath,char *MOD06file,
               float *sensor_zenith,float *sensor_azimuth);

void readMOD06_1km_cldmask(char *modPath,char *MOD06file,long na_1km,
                           signed char *cloud_mask);

void readMOD06_QA_1km(char *modPath,char *MOD06file,long na_1km,
                           signed char *qa_1km);

void match_1km_2_5km(float *lon_1km,float *lat_1km,float *lon_5km,
                     float *lat_5km,long na_1km,long na_5km,
		     long *index_1to5km,float *dist_1to5km);

void screen_reff(float *lon_1km,float *lat_1km,long na_1km,
                 float *reff_raw,float *tau,signed char *cldtop_phase,
                 signed char *cloud_mask,long *index_1to5km,
                 signed char *qa_1km,long *npix,long *nglnt_pix,
                 long *ncld_pix,long *nliq_pix,long *nscrn_pix,float *reff);

void getMODneighbors(float lon,float lat,float dist_thresh,
                     float *lon_mod,float *lat_mod,long na,long nx,
                     long *ineighbors,long *cnt,long ixmin, long ixmax,
                     long iamin,long iamax);

void getFileNames(char *modPath, char *mod02Path, char *yyyyddd,char *hhmm,
									const char *MOD03list,const char *MOD04list,const char *MOD02list,
                  char *MOD03file,char *MOD04file,char *MOD02file);

void readMOD06_IR_tbs(char *modPath,char *MOD06file,long *na_5km,
//                    float32 *lat_5km,float *lon_5km,float *Tbs);
                    float *lat_5km,float *lon_5km,float *Tbs);

void readMOD06_lwpuncertainty(char *modPath,char *MOD06file,long na_1km,
                              float *lwp_err,float *tau_err);
                              
char get_MOD06_filename(const char *MOD06_path,const char *list_filename,
                        int year,int doy,char *MOD06_file);
                        
float get_nearest_MOD06_gran(double jdate,char *MOD06_path,const char *MOD06list,
                          const char *MYD06list,const char *MOD04list,
                          const char *MYD04list,const char *MOD03list,
                          const char *MYD03list,char *MOD06file,char *MOD03file);

void create_test_image(float *ir_img,float *lons,float *lats,
											long na_5km,float lon_cent,float lat_cent);

