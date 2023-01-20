#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mask.h"


void read_land_ocean_qd(const char *filename,short *land_ocean,
                        float *lon,float *lat)
{
  FILE *fp=NULL;
  char endflag=0;
  char tmp[2885];
  char tmpchar[5];
  int lat_cnt=0;
  long i, iloc;

  fp = fopen(filename,"r");
  while ( fgets(tmp,2881,fp) != NULL) {
    for (i=0; i<1440; i++) {
//      strncat(tmpchar,tmp+(i*2),1);
      memcpy(tmpchar,tmp+(i*2),1);
      iloc = (lat_cnt*1440) + i;
      land_ocean[iloc] = atoi(tmpchar);
    }
    lat_cnt++;
  }
  fclose(fp);
  for (i=0; i<1440; i++) lon[i] = ((float) i/4.)-179.875;
  for (i=0; i<720; i++) lat[i] = 89.875-((float) i/4.);
}


//void read_land_ocean_qd(const char *filename,short *land_ocean)
/******
void read_land_ocean_qd(const char *filename,short *land_ocean,
                        float *lon,float *lat)
{
  FILE *fp=NULL;
  char endflag=0;
  char tmp[2885];
  char tmpchar[5];
  int lon_cnt=0;
  long i;

  fp = fopen(filename,"r");
//  while ( fgets(tmp,2880,fp) != NULL ) {
  while (lon_cnt < 720) {
    fgets(tmp,2880,fp);
printf("%d ",lon_cnt);
    for (i=0; i<1440; i++) {
if (lon_cnt==719) printf("%ld ",i);
      strncat(tmpchar,tmp+(i*2),1);
//      land_ocean[i*lon_cnt] = atoi(tmpchar);
    }
printf("%d ",lon_cnt);
    lon_cnt++;
  }
  fclose(fp);
  for (i=0; i<1440; i++) lon[i] = ((float) i/4.)-180.;
  for (i=0; i<720; i++) lat[i] = 90-((float) i/4.);
}
*******/
