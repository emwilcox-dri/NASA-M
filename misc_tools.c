#include <stdio.h>
#include <stdlib.h>
#include "misc_tools.h"




float min(float *var,long npix)
{
	long i;
	float min_var=9999999999.;
	
	for (i=0; i<npix; i++) if (var[i]<min_var) min_var = var[i];
	
	return min_var;
}





float max(float *var,long npix)
{
	long i;
	float max_var=-9999999999.;
	
	for (i=0; i<npix; i++) if (var[i]>max_var) max_var = var[i];
	
	return max_var;
}