#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermo.h"

/**********************************
NOTE: All pressures must be passed to these routines in Pascals, 
NOT millibar or hPa
ALSO: All relative humidities are fractions, NOT percents
**********************************/


float calc_theta(float T,float P)
{
//  float P0 = 100000.; // in Pascals
//  float R = 287.;    // in J/deg/kg
//  float Cp = 1004.;  // in J/deg/kg
  float theta;
  
  theta = T * pow( P0/P, R/CP  );

  return theta;
}


/******
calculates the saturation vapor pressure for a specified T.
uses the Goff and Gratch, 1946 equation
According to wikipedia, e_s is in hPa in the equation for log_e below,
which must be true since the reference pressure at the end of the equation
(1013.246) is clearly in hPa.
******/
float calc_es(float T)
{
  float log_e, e_s, a, b;

  a = pow( 10, 11.344*(1-(T/373.16)) );
  b = pow( 10, -3.49149*((373.16/T)-1) );
  log_e = (-7.90298*((373.16/T)-1)) 
          + (5.02808*log10(373.16/T))
          - (1.3816*pow(10,-7)*(a-1))
          + (8.1328*pow(10,-3)*(b-1))
          + log10(1013.246);
  e_s = pow(10,log_e);
  e_s *= 100.; //convert from hPa to Pa

  return e_s;
}


float calc_q_from_RH_T(float RH,float T,float P)
{
//  float eps=0.622;
  float es, e, w, q;

  es = calc_es(T);
  e = es*RH;
  w = (EPS*e)/(P-e);
  q = w/(1+w);

  return q;
}


float calc_theta_e(float T, float P)
{
//  float eps=0.622;
//  float L = 2500000.;  // in J/kg
//  float Cp = 1004.;  // in J/deg/kg
  float es, ws, theta, theta_e;

  es = calc_es(T);
  ws = (EPS*es)/(P-es);
  theta = calc_theta(T,P);
  theta_e = theta * exp( (L*ws)/(CP*T) );

  return theta_e;
}


float calc_h(float T,float z,float q)
{
  float h;

  h = CP*T + G*z + L*q;

  return h;
}


/*********
Note: my derivation results in a different equation than in Ram's notes.
2/11/2012 - This doesn't seem to work right. Need to check the equation.
*********/
float calc_e_from_q(float P,float q)
{
  float e;

//  e = (P*(1-q)) / (q*(EPS-1) + 1);
  e = (P*q)/(EPS+q);

  return e;
}


float calc_e_from_rh(float T,float rh)
{
  float es, e;

  es = calc_es(T);
//  e = rh*es/EPS; //where did I get the EPS from!? - typo in Ram's class notes!
  e = rh*es;

  return e;
}


float calc_Tv(float T,float rh,float P)
{
  float e, Tv;

  e = calc_e_from_rh(T,rh);
  Tv = T / (1-(e/P)*(1-EPS));

  return Tv;
}


/*************
* this code assumes that the first element of the profiles
* is the surface value.
* NOTE: P must be in pascals
*************/
/********
float calc_CAPE(float *T,float *rh,int nlevs,float *P,float *z,float *P_lcl)
{
  int ilev, ilev_lcl, ilev_neut;
  float z_lcl, zlev_lcl, zdiff, zdifftmp, cape;
  float gamma, Tp, Tvp, Tve;

  cape = 0;

//    find LCL 
  z_lcl = calc_lcl(T[0],rh[0]);
  zdiff = fabs(z[0]-z_lcl);
  ilev_lcl = 0;
  *P_lcl = P[ilev_lcl];
  zlev_lcl = z[ilev_lcl];
  for (ilev=1; ilev<nlevs; ilev++) {
//printf("[%d %f %f %f %f]  ",ilev,z[ilev],z_lcl,zdifftmp,*P_lcl);
    zdifftmp = fabs(z[ilev]-z_lcl);
//fflush(NULL);
    if (zdifftmp<zdiff) {
      ilev_lcl = ilev;
      *P_lcl = P[ilev_lcl];
      zlev_lcl = z[ilev_lcl];
    }
  }
//  printf("%d ",ilev_lcl);
//  exit(0);
//  fflush(NULL);

// This isn't strictly true - Parcel temp doesn't have to equal environmental tempat the LCL
  if (ilev_lcl>0) {
//    Tp = T[0]-((G/CP)*(zlev_lcl-z[0]));
    Tp = T[ilev_lcl];
  } else Tp = T[0];
  
  
  ilev=ilev_lcl;
  Tvp = calc_Tv(Tp,1.,*P_lcl);
  Tve = calc_Tv(T[ilev],rh[ilev],P[ilev]);
  while (Tvp > Tve) {
    gamma = calc_gamma_m(T[ilev],P[ilev],rh[ilev]);
    Tp -= (gamma*(z[ilev+1]-z[ilev]));
    ilev++;
    cape += R * (Tvp - Tve) * (log(P[ilev-1])-log(P[ilev]));
    Tvp = calc_Tv(Tp,1.,P[ilev]);
    Tve = calc_Tv(T[ilev],rh[ilev],P[ilev]);
  }

  return cape;
}
***********/



float calc_CAPE(float *T,float *rh,int nlevs,float *P,float *z,float *P_lcl)
{
  int ilev, ilev_lcl, ilev_nb;
  float z_lcl, zlev_lcl, zdiff, zdifftmp, cape;
  float gamma, Tvp, Tve, cape_p, cape_n;
  float q0, ep, esp, gam;
  float *Tp, *rhp;

  cape_p = 0;
  cape_n = 0;
  cape = 0;
  
  Tp = (float *) malloc(nlevs*sizeof(float));
  rhp = (float *) malloc(nlevs*sizeof(float));

/**** find LCL *****/
  z_lcl = calc_lcl(T[0],rh[0]);
  zdiff = fabs(z[0]-z_lcl);
  ilev_lcl = 0;
  *P_lcl = P[ilev_lcl];
  zlev_lcl = z[ilev_lcl];
  for (ilev=1; ilev<nlevs; ilev++) {
    zdifftmp = fabs(z[ilev]-z_lcl);
    if (zdifftmp<zdiff) {
      ilev_lcl = ilev;
      *P_lcl = P[ilev_lcl];
      zlev_lcl = z[ilev_lcl];
    }
  }
  
  printf("old ilev_lcl=%d\n",ilev_lcl);
  
  
  ilev_lcl=-1;
  ilev_nb=-1;
//  for (ilev=0; ilev<=ilev_lcl; ilev++) {
  q0 = calc_q_from_RH_T(rh[0],T[0],P[0]);
  Tp[0] = T[0];
  rhp[0] = rh[0];
  for (ilev=1; ilev<=nlevs-1; ilev++) {
    Tve = calc_Tv(T[ilev],rh[ilev],P[ilev]);
    if (ilev_lcl==-1) {
      Tp[ilev] = Tp[ilev-1] - ((G/CP)*(z[ilev]-z[ilev-1]));
      gam = G/CP;
      ep = q0*P[ilev]/EPS;
      esp = calc_es(Tp[ilev]);
      rhp[ilev] = ep/esp;
      if (ilev==0) cape_n -= R * (Tvp - Tve) * (log(P0)-log(P[ilev]));
      else cape_n -= R * (Tvp - Tve) * (log(P[ilev-1])-log(P[ilev]));
      if (rhp[ilev]>=1.) {
        rhp[ilev]=1;
        ilev_lcl=ilev;
      }
      Tvp = calc_Tv(Tp[ilev],rhp[ilev],P[ilev]);
    } else {
      gamma = calc_gamma_m(Tp[ilev-1],P[ilev-1],1.);
      gam=gamma;
      Tp[ilev] = Tp[ilev-1] - (gamma*(z[ilev]-z[ilev-1]));
      rhp[ilev]=1.;
      Tvp = calc_Tv(Tp[ilev],rhp[ilev],P[ilev]);
      if (Tvp<Tve) {
        ilev_nb = ilev;
//        ilev=nlevs;
      }
      cape_p += R * (Tvp - Tve) * (log(P[ilev-1])-log(P[ilev]));
    }
    printf("%d %f %f %f %f %f %f %f\n",ilev,Tp[ilev],rhp[ilev],T[ilev],rh[ilev],gam,Tvp,Tve);
  }
  
  printf("ilev_lcl=%d ilev_nb=%d cape_n=%f cape_p=%f cape=%f\n",
      ilev_lcl,ilev_nb,cape_n,cape_p,cape_p-cape_n);

  free(Tp);
  free(rhp);
  
  return cape;
}


/*************
Approximate: From Lawrence 2005
*************/
float calc_Td(float T,float rh)
{
  float Td;

  Td = (T-T_ABS) - ((100.-(rh*100.))/5.)*pow((T/300.),2) -
       0.00135*pow(((rh*100.)-84.),2) + 0.35;

  return Td;
}



/*************
Approximate: From Lawrence 2005
Returns height of lcl in meters.
*************/
float calc_lcl(float T,float rh)
{
  float z_lcl;

  z_lcl = (20.+((T-T_ABS)/5.))*(100.-(rh*100.));

  return z_lcl;
}



/**********
pseudoadiabatic lapse rate
**********/
float calc_gamma_m(float T,float P,float rh)
{
  float w, gamma_m;

  w = calc_mr(T,P,rh);
  gamma_m = G * ( (1+w)*(1 + (L*w/(R*T))) ) /
             ( CP + (w*CV) + (L*L*w*(EPS+w)/(R*R*R)) );
  if (gamma_m > G/CP) gamma_m = G/CP; //not to exceed the dry adiabat

  return gamma_m;
}


float calc_gamma_m_debug(long i,long j,float T,float P,float rh)
{
  float w, gamma_m;

  w = calc_mr(T,P,rh);
  gamma_m = G * ( (1+w)*(1 + (L*w/(R*T))) ) /
             ( CP + (w*CV) + (L*L*w*(EPS+w)/(R*R*R)) );
if (i==50 && j==200) {
printf("%f %f %f\n",G,(1+w)*(1 + (L*w/(R*T))),( CP + (w*CV) + (L*L*w*(EPS+w)/(R*R*R)) ) );
}

  return gamma_m;
}



/**********
water vapor mixing ratio
**********/
float calc_mr(float T,float P,float rh)
{
  float mr, e;

  e = calc_e_from_rh(T,rh);
  mr = (EPS*e)/(P-e);
//  if (mr > G/CP) mr = G/CP; //not to exceed the dry adiabat - what!?!?!?

  return mr;
}




float calc_CAPE_debug(long i,long j,float *T,float *rh,int nlevs,
                float *P,float *z,
                float *lcl_out,float *gam_prof,float *Tvp_prof,
                float *Tve_prof)
{
	char cloud_base=0,cloud_top=0,stable=1;
  int ilev, ilev_lcl, ilev_neut;
  float z_lcl, zlev_lcl, zdiff, zdifftmp, P_lcl, cape;
  float gamma, Tp, Tvp, Tve;

  cape = 0;
  for (ilev=0; ilev<nlevs; ilev++) {
    Tvp_prof[ilev]=-9999.;
    Tve_prof[ilev]=-9999.;
    gam_prof[ilev]=-9999.;
  }

/**** find LCL *****/
  z_lcl = calc_lcl(T[0],rh[0]);
  zdiff = fabs(z[0]-z_lcl);
  ilev_lcl = 0;
  P_lcl = P[ilev_lcl];
  zlev_lcl = z[ilev_lcl];
  for (ilev=1; ilev<nlevs; ilev++) {
    zdifftmp = fabs(z[ilev]-z_lcl);
    if (zdifftmp<zdiff) {
      ilev_lcl = ilev;
      P_lcl = P[ilev_lcl];
      zlev_lcl = z[ilev_lcl];
    }
  }
  *lcl_out = P_lcl;

  if (ilev_lcl>0) {
//    Tp = T[0]-((G/CP)*(zlev_lcl-z[0]));
    Tp = T[ilev_lcl];
  } else Tp = T[0];
//if (i==50 && j==200) {
printf("\n\n");
printf("calc T at LCL: ");
printf("gamma_d=%f T0=%f del_z=%f, Tp=%f\n",
      G/CP,T[0],zlev_lcl-z[0],Tp);
//}
  ilev=ilev_lcl;
  Tvp = calc_Tv(Tp,1.,P_lcl);
  Tve = calc_Tv(T[ilev],rh[ilev],P[ilev]);
  Tvp_prof[ilev]=Tvp;
  Tve_prof[ilev]=Tve;
//if (i==50 && j==200)
printf("Tvp at LCL=%f Tve at LCL=%f\n",Tvp,Tve);
  while (Tvp > Tve) {
//if (i==50 && j==200)
	printf("T=%f P=%f rh=%f mr=%f e=%f gam_m=%f\n",
                     T[ilev],P[ilev],rh[ilev],
                     calc_mr(T[ilev],P[ilev],rh[ilev]),
                     calc_e_from_rh(T[ilev],rh[ilev]),
                     calc_gamma_m(T[ilev],P[ilev],rh[ilev]));
//    gamma = calc_gamma_m_debug(i,j,T[ilev],P[ilev],rh[ilev]);
    gamma = calc_gamma_m(T[ilev],P[ilev],rh[ilev]);
    gam_prof[ilev] = gamma;
    Tp -= (gamma*(z[ilev+1]-z[ilev]));
    ilev++;
    cape += R * (Tvp - Tve) * (log(P[ilev-1])-log(P[ilev]));
    Tvp = calc_Tv(Tp,1.,P[ilev]);
    Tve = calc_Tv(T[ilev],rh[ilev],P[ilev]);
    Tvp_prof[ilev]=Tvp;
    Tve_prof[ilev]=Tve;
//if (i==50 && j==200)
printf("gamma_m=%f del_Z=%f Tvp=%f Tve=%f\n",gamma,z[ilev+1]-z[ilev],Tvp,Tve);
  }
//if (i==50 && j==200)
	printf("cloud_base=%f cld_top=%f CAPE=%f\n\n\n",
                     P_lcl,P[ilev],cape);

  return cape;
}







float calc_shear(float u_base, float u_top,float v_base,float v_top)
{
  float shear;
  shear = sqrt( pow((u_base-u_top),2) + pow((v_base-v_top),2) );

  return shear;
}




/********************
float calc_T_advect_horiz calculates the horizontal temperature convergence
negative velocity dot grad temperature
********************/
/******
float calc_T_conv_horiz(float u,float v,float *lons,float *lats,float *T)
{
  float T_conv

  dx = R_e * (lons[2]-lons[0]) * cos(lats[1]*PI/180.) * (PI/180.);
  dy = R_e * (lats[2]-lats[0]) * (PI/180.);
  T_conv = -(u * 
}
******/




/*****************************
* NOTE: was originally taking RH and converting to mixing ration in g/kg
* but for MERRA2 now just taking the mixing ratio from the MERRA2 files in kg/kg (here in *rh)
* and converting to g/kg
*****************************/
float calc_CAPE_emanuel(int nlevs,float *T,float *P,float *rh)
{
  int i;
  float cape, Tc[80], Pp[80], mr[80];
  
//printf("calc_CAPE_emanuel in thermo.c\n");
  for (i=0; i<nlevs; i++) {
    Pp[i] = P[i];
    Tc[i] = T[i]-T_ABS;
//    mr[i] = calc_mr(T[i],P[i],rh[i])*1000.;
		mr[i] = rh[i]*1000.;
//printf("%f %f %f\n",Pp[i],Tc[i],mr[i]);
  }
  
  calcsound_(&nlevs, Pp, Tc, mr, &cape);
  
  return cape;
}




