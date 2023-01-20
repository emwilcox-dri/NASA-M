#define AMSR_BADFLAG   -9999
#define N_ALONG_AMSR    2010
/*#define N_ACROSS_AMSR    196*/
#define N_ACROSS_AMSR   250 
#define N_ACROSS_AMSR_5KM 486
#define N_QAFLAG_BYTES_AMSR 6
#define N_QAFLAG_TYPES_AMSR 10
#define NCHAR_FNAME_AMSR  38
#define NCHAR_FNAME_AMSR_TBS 56
#define N_DAILY_LONS_AMSR    1440
#define N_DAILY_LATS_AMSR     720
#define N_DAILY_VARS_V5_AMSR       6
#define N_DAILY_VARS_V7_AMSR       7
#define N_DAILY_PASS_AMSR       2


/************ Quality flag definitions ********
qa_flag(*,*,0): 0=ocean; 1=coast; 2=land
qa_flag(*,*,1): 0=no ice (climate); 1=ice (climate); 2=no ice (Tbs); 3=ice (Tbs)qa_flag(*,*,2): 0=6.9GHz good; 1=6.9GHz bad
qa_flag(*,*,3): 0=10.7GHz good; 1=10.7GHz bad
qa_flag(*,*,4): 0=18.7GHz good; 1=18.7GHz bad
qa_flag(*,*,5): 0=36.5GHz good; 1=36.5GHz bad
qa_flag(*,*,6): very low res retrieval 0=good; 1=bad; 2=incomplete
qa_flag(*,*,7): low res retrieval 0=good; 1=bad; 2=incomplete
qa_flag(*,*,8): medium res retrieval 0=good; 1=bad; 2=incomplete
qa_flag(*,*,9): high res retrieval 0=good; 1=bad; 2=incomplete
**********************************************/


#define FREE(x) if(x) { free(x); x = NULL; }


void read_L2Ocean(char *filename,float *lat,float *lon,long *na_amsr,
                  long *nx_amsr,float *clw);

void read_L2Ocean_QA(char *filename,unsigned char *qa);

double getAMSR_fileName(char *amsrPath,int year,int mon,int day,int hour,
                      int min,char *amsrFile,int ichar_date);


void tmp(char *filename,float *lat,float *lon,
                           float *time,float *lwp,float *rain,
                           float *sst);

void read_v5_daily_gridded(char *filename,float *lat,float *lon,
                           float *time,float *lwp,float *rain,
                           float *sst);
                           
void read_daily_gridded(char *filename,float *lat,float *lon,
                           float *time,float *lwp,float *rain,
                           float *sst,float *cwv,int ver);

void read_AMSR_pixel_tbs(char *filename,long *na_amsr,long *nx_amsr,
                    float *lon_a,float *lat_a,float *lon_b,
                    float *lat_b,float *tb89_Ha,float *tb89_Va,
                    float *tb89_Hb,float *tb89_Vb,float *landfrac_a,
                    float *landfrac_b,int version);

