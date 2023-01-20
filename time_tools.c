
/*
  Utility time routines used by the fil layer and user
  
  NOTE: EMW 5/5/99
    EPS time is a 2 element long array where the first element
    is (julian date + 1/2 day) and the second element is the number
    of thousandths of a second since midnight.
    e.g.:
      jan. 13, 1998, 14:24:00
	Julian date = 2450827.1
	EPS time = [2450827, 51840000] equivalent to 2450827.6
*/
/******************
    EPS time is a 2 element long array where the first element
    is (julian date + 1/2 day) and the second element is the number
    of thousandths of a second since midnight.
    e.g.:
      jan. 13, 1998, 14:24:00
	Julian date = 2450827.1
	EPS time = [2450827, 51840000] equivalent to 2450827.6
******************/
#include <math.h>
#ifdef unix
#include <memory.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "time_tools.h"
/*#include "epsystem.h"*/

#define GREGORIAN (15+31L*(10+12L*1582))
#define JULGREG   2299161


static char SCCSid[]="@(#)fil_time.c	2.5    1/21/92";

static char ShortMonths[12][4] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
  };

static char LongMonths[12][10] = {
  "January", "February", "March",
  "April",   "May",      "June",
  "July",    "August",   "September",
  "October", "November", "December"
  };

static char ShortWeekDays[7][4] = {
  "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"
  };

static char LongWeekDays[7][10] = {
  "Sunday", "Monday", "Tuesday", "Wednesday",
  "Thursday", "Friday", "Saturday"
  };

void 
ep_time_to_mdyhms(time, mon, day, yr, hour, min, sec)
     long *time;
     int *mon, *day, *yr, *hour, *min;
     float *sec;
{
/*
 * convert eps time format to mdy hms
 */
  long ja, jalpha, jb, jc, jd, je;
  
/*  printf("%ld %ld %ld\n",*time, time[0], time[1]); */
  if(time[0] >= JULGREG) {
    jalpha=((float) (time[0]-1867216)-0.25)/36524.25;
    ja=time[0]+1+jalpha-(long)(0.25*jalpha);
  } else
    ja=time[0];
  jb=ja+1524;
  jc=6680.0+((float)(jb-2439870)-122.1)/365.25;
  jd=365*jc+(0.25*jc);
  je=(jb-jd)/30.6001;
  *day=jb-jd-(int)(30.6001*je);
  *mon=je-1;
  if(*mon > 12) *mon -= 12;
  *yr=jc-4715;
  if(*mon > 2) --(*yr);
  if(*yr <=0) --(*yr);

  ja = time[1]/1000;
  *hour = ja/3600;
  *min = (ja - (*hour)*3600)/60;
  *sec = (float)(time[1] - ((*hour)*3600 + (*min)*60)*1000)/1000.0;
}


void 
mdyhms_to_ep_time(int mon,int day,int yr,int hour,int min,
    double sec,long *time)
{
/*
 * convert mdy hms to eps time
 */
  long jul, ja, jy, jm;

  if(yr < 0) ++yr;
  if(mon > 2) {
    jy = yr;
    jm = mon+1;
  } else {
    jy = yr-1;
    jm = mon+13;
  }
  jul = (long)(floor(365.25*jy)+floor(30.6001*jm)+day+1720995);
  if(day+31L*(mon+12L*yr) >= GREGORIAN) {
    ja=0.01*jy;
    jul += 2-ja+(int)(0.25*ja);
  }
  time[0]=jul;
  time[1]=(hour*3600+min*60)*1000+(sec*1000);
}


int 
eps_CountLetters(s)
     char **s;
{
/*
 * CountLetters -- count the run of letters in s matching 
 * the first letter in s.  Advance s to point to the next
 * letter that does not match the first character.  Return
 * the count.  This is a helper function for eptimetostr.
 */
  int n;
  char c;
  
  for(n=0, c = **s; **s == c; n++, (*s)++);
  return n;
}

void 
ep_time_to_str(time, frmt, str)
     long *time;
     char *frmt, **str;
{
/*
 * convert eps time format to a string using a format string
 * (see EPS manual)
 */
  int Count, fract;
  char *sptr, *fptr;
  char tmp[256];
  int mon, day, yr, hr, min, dayofweek;
  float sec;

  memset(tmp, 0, 256);
  sptr = tmp;

  fptr = frmt;

  (void) ep_time_to_mdyhms(time,&mon,&day,&yr,&hr,&min,&sec);
  dayofweek = (time[0]+1)%7;

  while(*fptr) {
    switch(*fptr) {
    case 'D':				/* day */
      Count = eps_CountLetters(&fptr);
      if(Count >= 4) {
	strcat(sptr, LongWeekDays[dayofweek]);
      } else if (Count == 3) {
	strcat(sptr, ShortWeekDays[dayofweek]);
      } else {
	if((Count == 2) && (day < 10))
	  strcat(sptr++, "0");
	sprintf(sptr,"%d",day);
      }
      while(*sptr)
	sptr++;
      break;
    case 'M':				/* month */
      Count = eps_CountLetters(&fptr);
      if(Count >= 4) {
	strcat(sptr, LongMonths[mon-1]);
      } else if (Count == 3) {
	strcat(sptr, ShortMonths[mon-1]);
      } else {
	if((Count == 2) && (mon < 10))
	  strcat(sptr++, "0");
	sprintf(sptr,"%d",mon);
      }
      while(*sptr)
	sptr++;
      break;
    case 'Y':				/* year */
      Count = eps_CountLetters(&fptr);
      if(Count >= 4) {
	sprintf(sptr,"%d",yr);
      } else {
	if((yr%100) < 10)
	  strcat(sptr++, "0");
	sprintf(sptr,"%d",(yr%100));
      }
      while(*sptr)
	sptr++;
      break;
    case 'h':				/* hour */
      Count = eps_CountLetters(&fptr);
      if(Count >= 2 && (hr < 10))
	strcat(sptr++, "0");
      sprintf(sptr,"%d",hr);
      while(*sptr)
	sptr++;
      break;
    case 'm':				/* minute */
      Count = eps_CountLetters(&fptr);
      if(Count >= 2 && (min < 10))
	strcat(sptr++, "0");
      sprintf(sptr,"%d",min);
      while(*sptr)
	sptr++;
      break;
    case 's':				/* second */
      Count = eps_CountLetters(&fptr);
      if(Count >= 2 && ((int)sec < 10))
	strcat(sptr++, "0");
      sprintf(sptr,"%d",(int)sec);
      while(*sptr)
	sptr++;
      break;
    case 'f':				/* fraction of second */
      Count = eps_CountLetters(&fptr);
      fract = ((int)sec*1000)%1000;
      if(Count >= 3) {
	sprintf(sptr,"%03d",fract);
      } else if(Count == 2) {
	sprintf(sptr,"%02d",fract/10);
      } else {
	sprintf(sptr,"%01d",fract/10);
      }
      while(*sptr)
	sptr++;
      break;
    case '\\':				/* copy next character */
      fptr++;
    default:
      *sptr = *fptr;
      sptr++;
      fptr++;
      break;
    }
  }
  *sptr=0;

/* commented out because of type error!!!!! */
/*  *str=strdup(tmp); */
}

void 
ep_time_sub(time1, time2, delta)
     long *time1, *time2, *delta;
{
/*
 * substract time1 from time2
 */
  delta[0] = time2[0] - time1[0];
  delta[1] = time2[1] - time1[1];

  if((delta[0] >0) && (delta[1] < 0)) {
/*
 * subtract one day from delta[0] and add one day to delta[1]
 */
    delta[0]--;
    delta[1] = delta[1] + 86400000;
  }
  if((delta[0] <0) && (delta[1] > 0)) {
/*
 * add one day to delta[0] and subract on day from delta[1]
 */
    delta[0]++;
    delta[1] = delta[1] - 86400000;
  }
}

int 
ep_time_intervals(time1, time2, delta)
     long *time1, *time2, *delta;
{
/*
 * integer function that returns the number of delta time 
 * intervals between time1 and time2
 */
  long diff[2];
  int result;
  float num, denom;

  (void) ep_time_sub(time1, time2, diff);
  num=diff[0] + (diff[1])/86400000.0;
  denom=delta[0] + (delta[1])/86400000.0;

  if(denom != 0.0)
    result = num/denom+0.5;
  else
    result = 0;

  return(result);
}

void 
ep_time_to_ydn(time, year, yrday)
     long *time;
     int *year, *yrday;
{
/*
 * converts eps time to year and day number
 */
  int mon, day, yr, hr, min;
  long time2[2];
  float sec;

  (void) ep_time_to_mdyhms(time, &mon, &day, &yr, &hr, &min, &sec);

  mon=day=1;
  hr=min=0;
  sec=0.0;

  mdyhms_to_ep_time(mon, day, yr, hr, min, sec, time2);

  *year=yr;
  *yrday=time[0] - time2[0] + 1;
}

void 
ydn_to_ep_time(year, yrday, time)
     int year, yrday;
     long *time;
{
/*
 * converts year and day number to eps time
 */
  int mon, day, yr, hr, min;
  float sec;

  mon=day=1;
  hr=min=0;
  sec=0.0;
  yr= year;

  mdyhms_to_ep_time(mon, day, yr, hr, min, sec, time);

  time[0]=time[0]+yrday-1;
}

void 
ep_time_add(time1, delta, time2)
     long *time1, *delta, *time2;
{
/*
 * add time1 and delta
 */
  time2[0] = time1[0] + delta[0];
  time2[1] = time1[1] + delta[1];

  if(time2[1] > 86400000) {
    time2[0]++;
    time2[1] = time2[1] - 86400000;
  }
}


void yymmdd_to_mdyhms(long yymmdd,int *mon,int *day,int *year,int *hour,
		      int *min,float *sec)
{
  int yy;
  
  yy = (int) floor(yymmdd/10000);
  *mon = (int) floor((yymmdd-(yy*10000))/100);
  *mon = (int) floor((yymmdd-(yy*10000))/100);
  *day = (int) yymmdd-((yy*10000)+(*mon*100));
  *year = 1900 + yy;
  *hour = 0;
  *min = 0;
  *sec = 0;
}


void yyyymmdd_to_mdyhms(long yyyymmdd,int *mon,int *day,int *year,
                        int *hour,int *min,float *sec)
{
  int yy;
 
  yy = (int) floor(yyyymmdd/10000);
  *mon = (int) floor((yyyymmdd-(yy*10000))/100);
  *mon = (int) floor((yyyymmdd-(yy*10000))/100);
  *day = (int) yyyymmdd-((yy*10000)+(*mon*100));
  *year = yy;
  *hour = 0;
  *min = 0;
  *sec = 0;
}


/*****************
yymmdd_to_juldate() converts yymmdd to julian date at midnight
*****************/
void yymmdd_to_juldate(long yymmdd, double *jdate)
{
  int year, mon, day, hour, min;
  float sec;
  
  yymmdd_to_mdyhms(yymmdd,&mon,&day,&year,&hour,&min,&sec);
  mdyhms_to_juldate(mon,day,year,hour,min,sec,jdate);
  
}



void yyyymmdd_to_juldate(long yyyymmdd, double *jdate)
{
  int year, mon, day, hour, min;
  float sec;

  yyyymmdd_to_mdyhms(yyyymmdd,&mon,&day,&year,&hour,&min,&sec);
  mdyhms_to_juldate(mon,day,year,hour,min,sec,jdate);

}


void juldate_to_eps(double jdate, long epsdate[2])
{
  double dayfrac, day;

  dayfrac = modf(jdate,&day);
  dayfrac += .5;
  if (dayfrac > 1)
    {
      dayfrac--;
      day++;
    }
  dayfrac = floor(dayfrac*86400000);
  epsdate[0] = (long) day;
  epsdate[1] = (long) dayfrac;
}



void eps_to_juldate(long epsdate[2], double *jdate)
{
  long day;
  double dayfrac;

  day = epsdate[0];
  dayfrac = (double) epsdate[1]/86400000;
  dayfrac -= .5;
  if (dayfrac < 1)
    {
      dayfrac++;
      day--;
    }
  *jdate = (double) day+dayfrac;
}



void juldate_to_mdyhms(double jdate,int *mon,int *day,int *yr,int *hour,
		      int *min,float *sec)
{
  long epsdate[2];

  juldate_to_eps(jdate,epsdate);
  ep_time_to_mdyhms(epsdate, mon, day, yr, hour, min, sec);
/* added emw after getting hour=24 */
  if (*hour==24) {
    *hour=0;
    if (*day==num_of_days_in_month(*mon,*yr)) {
      *day=1;
      if (*mon==12) {
        *mon=1;
        *yr=*yr+1;
      } else *mon=*mon+1;
    } else *day=*day+1;
  }
}




void mdyhms_to_juldate(int mon,int day,int yr,int hour,
		      int min,float sec,double *jdate)
{
  long epsdate[2];
  double tmpsec;
  
  tmpsec = (double) sec;
  mdyhms_to_ep_time(mon, day, yr, hour, min, tmpsec, epsdate);
  eps_to_juldate(epsdate,jdate);
}



void juldate_to_yymmdd(double jdate,long *yymmdd)
{
  int mon,day,yr,hour,min;
  long lontmp;
  float sec;
  
  juldate_to_mdyhms(jdate,&mon,&day,&yr,&hour,&min,&sec);
  yr -= 1900;
  lontmp = (long) yr*10000;
  lontmp += (long) ((mon*100) + day);
  *yymmdd = lontmp;
}




void juldate_to_yyyymmdd(double jdate,long *yyyymmdd)
{
  int mon,day,yr,hour,min;
  long lontmp;
  float sec;
 
  juldate_to_mdyhms(jdate,&mon,&day,&yr,&hour,&min,&sec);
  lontmp = (long) yr*10000;
  lontmp += (long) ((mon*100) + day);
  *yyyymmdd = lontmp;
}



void get_month_name(int mon,char *name)
{
  switch (mon)
  {
    case 1:
      {
	sprintf(name,"Jan");
	break;
      }
    case 2:
      {
	sprintf(name,"Feb");
	break;
      }
    case 3:
      {
	sprintf(name,"Mar");
	break;
      }
    case 4:
      {
	sprintf(name,"Apr");
	break;
      }
    case 5:
      {
	sprintf(name,"May");
	break;
      }
    case 6:
      {
	sprintf(name,"Jun");
	break;
      }
    case 7:
      {
	sprintf(name,"Jul");
	break;
      }
    case 8:
      {
	sprintf(name,"Aug");
	break;
      }
    case 9:
      {
	sprintf(name,"Sep");
	break;
      }
    case 10:
      {
	sprintf(name,"Oct");
	break;
      }
    case 11:
      {
	sprintf(name,"Nov");
	break;
      }
    case 12:
      {
	sprintf(name,"Dec");
	break;
      }
  }  
}


long get_next_day(long yymmdd)
{
  int days_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31};
  int year, mon, day, hour, min;
  long next_day;
  float sec;
  
  yymmdd_to_mdyhms(yymmdd,&mon,&day,&year,&hour,&min,&sec);
  year -= 1900;
/***
  year = (int) floor(yymmdd/10000);
  mon = (int) floor((yymmdd-(year*10000))/100);
  day = (int) yymmdd-((year*10000)+(mon*100));
***/
  day++;
  if (day > days_in_month[mon-1])
    {
      if (mon == 12)
        {
          year++;
          mon = 1;
          day = 1;
        }
      else
        {
          mon++;
          day = 1;
        }
    }
  next_day = (long) year*10000+mon*100+day;
  return next_day;
}



long get_next_day_y2k(long yyyymmdd)
{
  int days_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31};
  int year, mon, day, hour, min;
  long next_day;
  float sec;
 
  yyyymmdd_to_mdyhms(yyyymmdd,&mon,&day,&year,&hour,&min,&sec);
  day++;
  if (day > days_in_month[mon-1])
    {
      if (mon == 12)
        {
          year++;
          mon = 1;
          day = 1;
        }
      else
        {
          mon++;
          day = 1;
        }
    }
  next_day = (long) year*10000+mon*100+day;
  return next_day;
}




long get_previous_day(long yymmdd)
{
  int days_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31};
  int year, mon, day, hour, min;
  long previous_day;
  float sec;
  
  yymmdd_to_mdyhms(yymmdd,&mon,&day,&year,&hour,&min,&sec);
/***
  year = (int) floor(yymmdd/10000);
  mon = (int) floor((yymmdd-(year*10000))/100);
  day = (int) yymmdd-((year*10000)+(mon*100));
***/
  day--;
  if (day < 1)
    {
      if (mon == 1)
        {
          year--;
          mon = 12;
          day = 31;
        }
      else
        {
          mon--;
          day = days_in_month[mon-1];
        }
    }
  if (year > 100) year -= 1900;
  previous_day = (long) year*10000+mon*100+day;
  return previous_day;
}



/******************************************************************/
/*      
 *      Returns number of days in a month for a given year; the
 *      months are indexed 1-12; the year is specified in full 
 *      (e.g. 1970). 
 *                   
 */
int     num_of_days_in_month(int month, int year)

{
switch (month)
        {
        case 1: return (nod_jan);
                break;
        case 2: if (isly(year) == 1)
                        return (nod_febl);
                else
                        return (nod_feb);
                break;
        case 3: return nod_mar;
                break;
        case 4: return nod_apr;
                break;
        case 5: return nod_may;
                break;
        case 6: return nod_jun;
                break;
        case 7: return nod_jul;
                break;
        case 8: return nod_aug;
                break;
        case 9: return nod_sep;
                break;
        case 10: return nod_oct;
                break;
        case 11: return nod_nov;
                break;
        case 12: return nod_dec;
                break;

        }

fprintf(stderr,"Month argument out of range. \n");
exit(EXIT_FAILURE);

}




void doy_to_mmdd(int doy,int yyyy,char *mmdd)
{
  char strzero[2];
  char strtmp[10];
  int mon, dd;

  mon = 0;

  while (doy > 0) {
    mon++;
    dd = doy;
    doy -= num_of_days_in_month(mon,yyyy);
  }
  if (mon<10 && dd<10) sprintf(mmdd,"0%d0%d",mon,dd);
  if (mon<10 && dd>=10) sprintf(mmdd,"0%d%d",mon,dd);
  if (mon>=10 && dd<10) sprintf(mmdd,"%d0%d",mon,dd);
  if (mon>=10 && dd>=10) sprintf(mmdd,"%d%d",mon,dd);
printf("mmdd=%s\n",mmdd);
}




/***************************************
NOTE: This also exists in the time_lib.c library as isly()
***************************************/
int is_leap_year(int year)
{
  int leap_flag=0;

  if ( ((fmod(year,4)==0) && (year/100!=0)) ||
        (fmod(year,400)==0) ) leap_flag=1;

  return leap_flag;
}




/******************
NOTE: The num_of_days_in_month() routine takes care of the leap year
******************/
int yyyymmdd_to_doy(long yyyymmdd)
{
  int leap_flag=0;
  int mon,day,yr,hour,min,doy,i;
  float sec;

  yyyymmdd_to_mdyhms(yyyymmdd,&mon,&day,&yr,&hour,&min,&sec);

  doy=0;
  if (mon>1) for (i=1;i<mon;i++) doy+=num_of_days_in_month(i,yr);
  doy+=day;

  return doy;
}




double yyyymmdd_hhmm_to_juldate(long yyyymmdd,int hhmm)
{
  float hr, mn, fday;
  double jd;
  
  yyyymmdd_to_juldate(yyyymmdd, &jd);
  hr = (float) floor(hhmm/100.);
  mn = ((float) hhmm)-(hr*100);
  fday = ((hr*60.)+mn)/(24.*60.);
  jd += fday;
  
  return jd;
}



