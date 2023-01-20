#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <limits.h>
#include "time_tools.h"
#include "solar.h"


/*********
calc_dailymean_dec takes the date specified as yymmdd and
returns a floating point scaler containing the daily-mean
solar declination for that day as calculated by the
sunae1 routine below.
*********/
float calc_dailymean_dec(long yymmdd,float lat,float lon)
{
  long cnt=0;
  float dec, solza;
  float mean_dec=0;
  double jd;
  double i;
  int mon,day,yr,hour,min;
  float sec;
  
  yymmdd_to_juldate(yymmdd, &jd);
  for (i=jd; i<jd+.5; i+=.001)
    {
      juldate_to_mdyhms(i,&mon,&day,&yr,&hour,&min,&sec);
      sunae1(i, lat, lon, &dec, &solza);
/*      printf("%lf %d %d %d %d %d %f %lf\n",i,mon,day,yr,hour,min,sec,dec);*/
      mean_dec += dec;
      cnt++;
    }
  mean_dec /= cnt;
  return mean_dec;
}


float calc_dailymean_dec_y2k(long yyyymmdd,float lat,float lon)
{
  long cnt=0;
  float dec, solza;
  float mean_dec=0;
  double jd;
  double i;
  int mon,day,yr,hour,min;
  float sec;

  yyyymmdd_to_juldate(yyyymmdd, &jd);
  for (i=jd; i<jd+.5; i+=.001)
    {
      juldate_to_mdyhms(i,&mon,&day,&yr,&hour,&min,&sec);
      sunae1(i, lat, lon, &dec, &solza);
/*      printf("%lf %d %d %d %d %d %f %lf\n",i,mon,day,yr,hour,min,sec,dec);*/
      mean_dec += dec;
      cnt++;
    }
  mean_dec /= cnt;
  return mean_dec;
}


/*********
calc_diurnalannualave_flux calculates a vector of the diurnal-annual averaged
solar flux for the year and list of latitudes supplied.  The year is non-
Y2K compliant, e.g. 1998 is specified as yy=98.
*********/
void calc_diurnalannualave_flux(long yy,long nlats,float *lats,float *ave_sol)
{
  long yymmdd;
  int ilat, i;
  float lat, dec;
  float x;
  float meandist_sqr;
  double juldate;
  
  for (ilat=0; ilat<nlats; ilat++) {
    yymmdd = (yy*10000)+101; 
    lat = lats[ilat]*PI/180.;
    for (i=0; i<365; i++) {
      yymmdd_to_juldate(yymmdd,&juldate);
      x = 2*PI*(juldate-1.)/365.;
      meandist_sqr = 1.000110 + (0.034221*cos(x)) + 
                      (0.001280*sin(x)) + (0.000719*cos(2*x)) + 
                      (0.000077*sin(2*x));
      dec = calc_dailymean_dec(yymmdd,lat,0);
      ave_sol[ilat] += (SOLAR_CONST*meandist_sqr/PI) *
              (acos(-tan(lat)*tan(dec))*sin(lat)*sin(dec) +
              cos(lat)*cos(dec)*sin(acos(-tan(lat)*tan(dec))));
      yymmdd = get_next_day(yymmdd);
    }
    ave_sol[ilat] /= 365.;
  }
}

void sunae1(double jd, float lat ,float lon, float *soldec, float *solza)

{
/*    This a C version of the routine sunae1 
found in stratus_insat.f by Sandrine Bony
translated by Remy Roca
     Purpose:
     Calculates azimuth and elevation of sun

updated 5/99 by Eric Wilcox to accept julian date (jd) as a
parameter and pass solar zenith angle (instead of elevation
angle) and solar declination back out.  Also changed so that
julian date, lat and lon are passed by value, not by reference.
NOTE: soldec given in radians
      solza given in degrees
*/

/*float  pi,twopi,rad,delta,leap,jd,time,mnlong,mnanom,eclong;*/
float  pi,twopi,rad,delta,leap,time,mnlong,mnanom,eclong;
float  oblqec, num, den, ra, dec, gmst, lmst, ha;
float  refrac,el;
int mon,day,yr,tmp,min;
float sec,hour;


/*
	year=date[0];
	day=date[1];
	hour=date[2];
*/
/* calculate hour from jd */
      juldate_to_mdyhms(jd, &mon, &day, &yr, &tmp, &min, &sec);
      hour = (float) tmp + (float) min/60;
      jd -= 2400000;
      jd += .5;
      
	solza[0]=0.;

      pi=4.* atan(1.);
      twopi=2.* pi;
      rad=pi/180.;
/*     get the current Julian date*/
/*
      delta=year-1949.;
      leap=delta/4;
      jd=32916.5+delta*365.+leap+day+hour/24.;
*/

/*     calculate ecliptic coordinates*/
      time=jd-51545.0;
/*     force mean longitude between 0 and 360 degs*/
      mnlong=280.460+(0.9856474 * time);
      mnlong=fmodf(mnlong,360.);

      if (mnlong < 0.) 
		mnlong=mnlong+360.;

/*     mean anomaly in radians between 0, 2*pi */

      mnanom=357.528+(0.9856003 * time);
      mnanom=fmodf(mnanom,360.);
      if(mnanom < 0.) mnanom=mnanom+360.;
      mnanom=mnanom* rad;

/*     compute ecliptic longitude and obliquity of ecliptic */

      eclong=mnlong+1.915*sin(mnanom)+0.020*sin(2.*mnanom);
      eclong=fmodf(eclong,360.);
      if(eclong < 0.0) eclong=eclong+360.;
      oblqec=23.429 -(0.0000004 * time);
      eclong=eclong* rad;
      oblqec=oblqec* rad;

/*     calculate right ascention and declination */
      num=cos(oblqec)*sin(eclong);
      den=cos(eclong);
      ra=atan(num/den);

/*     force ra between 0 and 2*pi */
      if(den < 0.) 
         ra=ra+pi;
      if(den < 0. && num< 0.)
         ra=ra+twopi;

/*     dec in radians */
      dec=asin(sin(oblqec)*sin(eclong));

/*     calculate Greenwich mean sidereal time in hours*/
      gmst=6.697375+0.0657098242*time+hour;
/*     hour not changed to sidereal sine "time" includes the fractional day*/
      gmst=fmodf(gmst,24.);
      if(gmst < 0.)
	 gmst=gmst+24.;

/*     calculate local mean sidereal time in radians */
      lmst=gmst+lon/15.;
      lmst=fmodf(lmst,24.);
      if(lmst<0.)
	 lmst=lmst+24.;
      lmst=lmst*15.* rad;

/*     calculate hour angle in radians between -pi, pi*/
      ha = lmst - ra;
      if(ha+pi < 0)
	 ha=ha+twopi;
      if(ha > pi)
	 ha=ha-twopi;
      lat=lat * rad;
/*C     calculate azimuth and elevation*/

      el=asin(sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(ha));

/*C     this puts azimuth between 0 and 2*pi radians
     calculate refraction correction for US stand. atm.*/

      el=el/rad;
      if(el>-0.56)
         refrac=3.51561 * ( 0.1594 + (0.0196 * el) + 0.00002 * (el * el) )/(1. + (0.505 * el)+0.0845 * (el * el));
      if (el<=-0.56)
         refrac=0.56;

      	el=el+refrac;
/*	solza[0]=el;*/
	solza[0]=90-el; /*pass the zenith angle, not elevation angle*/
	*soldec = dec;
}



/*****************
void main()
{
int mon,day,yr,hour,min;
float sec;
float lat, lon;
float el[1];
float date[3];
double jd, tmp, tmp2;
int i;

mon = 4;
day = 6;
yr = 1999;
sec = 0;

lat=5.;
lon=74.;

el[0]=-1.;

for (i=0; i<96; i++)
{
  tmp = modf(((double) i)/4,&tmp2);
  hour = (int) tmp2;
  min=0;
  if (floor(tmp*4) == 0) min=0;
  if (floor(tmp*4) == 1) min=15;
  if (floor(tmp*4) == 2) min=30;
  if (floor(tmp*4) == 3) min=45;
  mdyhms_to_juldate(mon,day,yr,hour,min,sec,&jd);
  sunae1(jd, lat, lon, el);
}

}
*********************/
