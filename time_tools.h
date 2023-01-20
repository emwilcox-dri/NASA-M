/*
Various routines for switching between time and date formats.
Note that the mdyhms format is NOT Y2K compliant (i.e. the
year 1999 is expressed as 99.
*/
#include <math.h>
#ifdef unix
#include <memory.h>
#endif
#include <stdio.h>
#include <string.h>
/*#include "epsystem.h"*/

#define GREGORIAN (15+31L*(10+12L*1582))
#define JULGREG   2299161

#define nod_jan         31
#define nod_feb         28
#define nod_febl        29
#define nod_mar         31
#define nod_apr         30
#define nod_may         31
#define nod_jun         30
#define nod_jul         31
#define nod_aug         31
#define nod_sep         30
#define nod_oct         31
#define nod_nov         30
#define nod_dec         31

#define isly(y)         (((y % 4) == 0) ? 1 : 0)

/**********
static char SCCSid[]="@(#)fil_time.c    2.5    1/21/92";

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
*******************/

void yymmdd_to_mdyhms(long yymmdd,int *mon,int *day,int *yr,int *hour,
		      int *min,float *sec);

void yyyymmdd_to_mdyhms(long yymmdd,int *mon,int *day,int *yr,
                        int *hour,int *min,float *sec);

void yymmdd_to_juldate(long yymmdd, double *jdate);

void yyyymmdd_to_juldate(long yyyymmdd, double *jdate);

void juldate_to_mdyhms(double jdate,int *mon,int *day,int *yr,int *hour,
                      int *min,float *sec);
                      
void eps_to_juldate(long epsdate[2], double *jdate);

void juldate_to_eps(double jdate, long epsdate[2]);

void mdyhms_to_juldate(int mon,int day,int yr,int hour,
                      int min,float sec,double *jdate);
		      
void juldate_to_yymmdd(double jdate,long *yymmdd);

void juldate_to_yyyymmdd(double jdate,long *yyyymmdd);

void get_month_name(int mon,char *name);

long get_next_day(long yymmdd);

long get_next_day_y2k(long yyyymmdd);

long get_previous_day(long yymmdd);

int  num_of_days_in_month(int month, int year);

void doy_to_mmdd(int doy,int yyyy,char *mmdd);

int is_leap_year(int year);

int yyyymmdd_to_doy(long yyyymmdd);

double yyyymmdd_hhmm_to_juldate(long yyyymmdd,int hhmm);
