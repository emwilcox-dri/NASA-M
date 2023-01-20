//#define  SOLAR_CONST 1369
#define  SOLAR_CONST 1361
#define  PI	     3.1415926

float calc_dailymean_dec(long yymmdd,float lat,float lon);

float calc_dailymean_dec_y2k(long yyyymmdd,float lat,float lon);

void calc_diurnalannualave_flux(long yy,long nlats,float *lats,float *ave_sol);

void sunae1(double jd, float lat ,float lon, float *soldec, float *solza);
