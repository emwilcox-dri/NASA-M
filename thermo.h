#define P0    100000.     // in Pascals
#define R        287.     // in J/deg/kg
#define CP      1004.     // in J/deg/kg
#define CV       719.     // in J/deg/kg
#define EPS        0.622
#define L    2500000.     // in J/kg
#define G          9.8    // in m/s/s
#define T_ABS    273.15   // in K (degC to K ==> K = T_ABS+degC)
#define R_e  6370000.     // radius of earth in meters
#define PI         3.1415926

float calc_theta(float T,float P);
float calc_es(float T);
float calc_q_from_RH_T(float RH,float T,float P);
float calc_theta_e(float T, float P);
float calc_h(float T,float z,float q);
float calc_e_from_q(float P,float q);
float calc_e_from_rh(float T,float rh);
float calc_Tv(float T,float q,float P);
float calc_CAPE(float *T,float *rh,int nlevs,float *lev,float *z,float *P_lcl);
float calc_Td(float T,float rh);
float calc_lcl(float T,float rh);
float calc_gamma_m(float T,float P,float rh);
float calc_gamma_m_debug(long i,long j,float T,float P,float rh);
float calc_mr(float T,float P,float rh);
float calc_CAPE_debug(long i,long j,float *T,float *rh,int nlevs,
                float *lev,float *z,
                float *P_lcl,float *gam_prof,float *Tvp_prof,
                float *Tve_prof);
float calc_shear(float u_base, float u_top,float v_base,float v_top);
/*float calc_T_conv_horiz(float u,float v,float *lons,float *lats,float *T) */

float calc_CAPE_emanuel(int nlevs,float *T,float *P,float *rh);

void calcsound_(int *, float *, float *, float *, float *);
