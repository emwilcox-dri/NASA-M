/*****
Modified by emw 3/11/2019: in line 164 changed abs to labs
Modified by emw 5/27/99
created header file das.h and moved some of the settings stuff there.

Also added the label_region routine to emulate the label_region routine
in IDL.  Really it just calls cloud_id_coarse with the class
field and classe set to zero.
*****/
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "das.h"

void   go_das(long R,long C,float *ir,short *class,long *cloud,long *ncloud,
	      int verbose)
{
FILE *fp_out;
long r,c,idx,i;
long step,substep;
long lab,*rcloud,*ccloud;
short classe,*class_simple;
float THC,THB;
float *mir;
long ncloudin;
long *category;
float *threshold_Detect;
float *threshold_Spread;
float *deltafine;

int pas_spread;

/* Allocate for memory here */
mir = (float *) calloc(R*C,sizeof(float));
class_simple = (short *) calloc(R*C,sizeof(short));
category = (long *) calloc(STEP,sizeof(long));
threshold_Detect = (float *) calloc(STEP,sizeof(float));
threshold_Spread = (float *) calloc(STEP,sizeof(float));
deltafine = (float *) calloc(STEP,sizeof(float));
rcloud =(long *) calloc(100000,sizeof(long));
if (rcloud==NULL) printf("ERROR: rcloud not allocated\n");
ccloud =(long *) calloc(100000,sizeof(long));
if (ccloud==NULL) printf("ERROR: ccloud not allocated\n");

/* Simplify the class information 
into low and NON LOW */

for(r=0;r<R;r++)
 for(c=0;c<C;c++)
        {
        idx=r*C+c;
	class_simple[idx] = class[idx];
	if(class[idx] >= 3) class_simple[idx] = 2;
        }


/*set the thresholds*/

if (STEP > 0) {
  threshold_Detect[0] = DETECT_0;
  threshold_Spread[0]= SPREAD_0;
  category[0] = CATEGORY_0;
  deltafine[0] = DELTA_SPREAD_0;
}
if (STEP > 1) {
  threshold_Detect[1] = DETECT_1;
  threshold_Spread[1]= SPREAD_1;
  category[1] = CATEGORY_1;
  deltafine[1] = DELTA_SPREAD_1;
}
if (STEP > 2) {
threshold_Detect[2] = DETECT_2;
        threshold_Spread[2]= SPREAD_2;
        category[2] = CATEGORY_2;
        deltafine[2] = DELTA_SPREAD_2;
}
if (STEP > 3) {
threshold_Detect[3] = DETECT_3;
        threshold_Spread[3]= SPREAD_3;
        category[3] = CATEGORY_3;
        deltafine[3] = DELTA_SPREAD_3;
}

*ncloud =0;

/*Negate IR into mir*/
for(r=0;r<R;r++)
 for(c=0;c<C;c++)
	{
	idx=r*C+c;
	mir[idx] =(-1.) *  ir[idx];
	}

if (verbose) printf("ir[0] mir[0] class[10000] class_simple[10000] ncloud in  DAS %f %f %ld %ld %ld\n",ir[0],mir[0],(long)(class[10000]), (long)(class_simple[10000]),*ncloud);

for(step=0;step<STEP;step++)
{
if (verbose)
 {
printf("Doing the Detect and Spread at step # %ld out of %d\n",step+1,STEP);
printf("Detect: %f Spread: %f Classe: %ld Delta K: %f\n",threshold_Detect[step],threshold_Spread[step],category[step],
deltafine[step]);
 }
/*Update the thresholds */
classe = category[step];
THC   = threshold_Detect[step];

/* Detect step */
/* Detect according to TH and classification */
if (verbose) printf("%ld %ld\n",R,C);
cloud_id_coarse(R,C,classe,mir,class_simple,cloud,THC,ncloud);
if (verbose) printf("Ncloud in GODAS %ld \n",*ncloud);

if (verbose) {
/**** write cloud to simple binary file here
Use step and "detect" in the filename for this stage ****/
}

/*extract a point per cloud */
/* from which spread is started */
for(r=0;r<R;r++)
 for(c=0;c<C;c++)
        {
        idx=r*C+c;
	lab=cloud[idx];
	rcloud[lab] = r;
	ccloud[lab] = c;
        }
if (verbose) printf("Filling departure point is OK\n");

/* Spread step */
/* Spread according to TH and classification */
/* resetting mir to positive values for pixels that have been labeled,
after being made very negative in cloud_id_coarse */
/**** commented out emw on 1/11/2011 because I'm not sure what this accomplishes
note: I also commented out the subtraction of 1000 further down
for (r=0;r<R;r++)
        for (c=0;c<C;c++)
                {
                idx     = r*C+c;
                if (cloud[idx] != 0)
                  mir[idx] = mir[idx] + 1000.0;
                }
****/

if (verbose) printf("Negating once again IR is OK\n");

/* Spread in multiple steps for deltafine K cloud separation */

	THB = threshold_Detect[step];
	substep = (long)((threshold_Spread[step] - threshold_Detect[step])/deltafine[step]);
	if (verbose) printf("substep=%ld\n",substep);
fflush(NULL);
	for (i=0;i<substep;i++)
		{
		THB=THB+deltafine[step];
	if (verbose) printf("	Doing the spreading substep at %f ...",THB);
fflush(NULL);

cloud_id_fine(R,C,classe,mir,class_simple,cloud,THB,ncloud,rcloud,ccloud);
  if (verbose) printf("  Done.\n");
fflush(NULL);
		/* resetting cloud to possitive values */
		/* and */
		/* turn mir negative for all labeled clouds */
	        for (r=0;r<R;r++)
                for (c=0;c<C;c++)
                        {
			idx     = r*C+c;
                        cloud[idx] = labs(cloud[idx]);
//			if (cloud[idx] != 0)  mir[idx] = mir[idx] - 1000.0; //changed emw 1/10/2012
//			if (cloud[idx] != 0)  mir[idx] = mir[idx] + 1000.0; //commented both out 1/11/2012 - not sure what the purpose of the adding/subtracting 1000 is
                        }


		}

if (verbose) {
/****** write cloud to simple binary file here
use step and "spread" in the filename *****/
}

/****
if (step == 1) {
  fp_out = fopen("990205.1.bin","w");
  fwrite(cloud,sizeof(long),(R*C),fp_out);
  fclose(fp_out);
}
***/
} /*end of the step loop*/
classe = category[STEP-1];
THC = threshold_Spread[STEP-1];
cloud_id_coarse(R,C,classe,mir,class_simple,cloud,THC,ncloud);
if (verbose) printf("Ncloud in GODAS %ld \n",*ncloud);

FREE(mir);
FREE(class_simple);
FREE(category);
FREE(threshold_Detect);
FREE(threshold_Spread);
FREE(deltafine);
FREE(rcloud);
FREE(ccloud);
}





/* "is_cloud_coarse" checks whether pixel r,c is a cloudy pixel; i.e. its
optical thickness is greater or equal to THRSHLD_C Erwin
Now takes also into accound the classe parameter
Modified Jan 18th by Remy

now only for > threshold - not >= threshold
Modified Jan. 1999 by emw
*/
int     is_cloud_coarse(long R,long C,long r, long c, long classe,
	float *mir,short *class,long *cloud,float THC)
{
int     outcome;
long    idx;

idx     = r*C+c;
outcome = 0;
if (r == R || c == C || r == 0 || c == 0) return (outcome);

//if (mir[idx]>=-252 && mir[idx]<-251) printf("[%ld %ld %ld]  ",r,c,cloud[idx]);
// if (mir[idx]>=THC) 
if (mir[idx]>THC)
        if (cloud[idx] == 0)
//		if(class[idx] == classe)  //commented out emw 11/26/2011
                        outcome = 1;
return  (outcome);
}






/*
"cloud_lb_coarse" labels a whole new cloud if pixel r,c is a cloudy pixel;
because checked cloud pixels are given zero thickness, relabeling is 
impossible. Now takes also into accound the classe parameter.
Modified Jan 18th by Remy
*/
void    cloud_lb_coarse(long R,long C,long r, long c, long cld_label,
			long classe,float *mir, short *class,long *cloud,
			float THC)
{
long    i;

i = r*C+c;

if (r == R || c == C || r == 0 || c == 0) return;

//if (r==48 && c>=89 && c<=101)
//printf("[%ld %ld %f] ",r,c,mir[r*C+c]);

cloud[i] = cld_label;
class[i] = (short) classe; //added emw 11/26/2011

if (is_cloud_coarse(R,C,r+1,c,classe,mir,class,cloud,THC)) 
 cloud_lb_coarse(R,C,r+1,c,cld_label,classe,mir,class,cloud,THC);
if (is_cloud_coarse(R,C,r-1,c,classe,mir,class,cloud,THC))
 cloud_lb_coarse(R,C,r-1,c,cld_label,classe,mir,class,cloud,THC);
if (is_cloud_coarse(R,C,r,c+1,classe,mir,class,cloud,THC))
 cloud_lb_coarse(R,C,r,c+1,cld_label,classe,mir,class,cloud,THC);
if (is_cloud_coarse(R,C,r,c-1,classe,mir,class,cloud,THC))
 cloud_lb_coarse(R,C,r,c-1,cld_label,classe,mir,class,cloud,THC);

/* 8-connected */
/******
if (is_cloud_coarse(R,C,r+1,c+1,classe,mir,class,cloud,THC))
 cloud_lb_coarse(R,C,r+1,c+1,cld_label,classe,mir,class,cloud,THC);
if (is_cloud_coarse(R,C,r-1,c-1,classe,mir,class,cloud,THC))
 cloud_lb_coarse(R,C,r-1,c-1,cld_label,classe,mir,class,cloud,THC);
if (is_cloud_coarse(R,C,r-1,c+1,classe,mir,class,cloud,THC)) 
 cloud_lb_coarse(R,C,r-1,c+1,cld_label,classe,mir,class,cloud,THC);
if (is_cloud_coarse(R,C,r+1,c-1,classe,mir,class,cloud,THC)) 
 cloud_lb_coarse(R,C,r+1,c-1,cld_label,classe,mir,class,cloud,THC);
******/
}





/*
"cloud_id_coarse" checks the whole regeon for cloudy pixels, and if one if found cloud_lb_coarse is called.  Now takes also into accound the classe parameter.
Modified Jan 18th by Remy
*/
void    cloud_id_coarse(long R,long C,long classe,float *mir,short *class,
			long *cloud,float THC,long *ncloud)
{

long    r,c;
long    cld_label;

cld_label       = *ncloud; 
/* Initial cloud label is set at 0 but 1 is the first one used.*/

/**
if (cld_label==1) {
for (r=0;r<R;r++)
for (c=0;c<C;c++) {
  if (mir[r*C+c]>-250 && cloud[r*C+c]==0)
    printf("[%ld %ld %f %ld] ",r,c,mir[r*C+c],cloud[r*C+c]);
}
}
**/

for (r=0;r<R;r++)
        for (c=0;c<C;c++)
                if(is_cloud_coarse(R,C,r,c,classe,mir,class,cloud,THC))
                        {
                        cld_label++;  /* First cloud label is 1 */  
                        cloud_lb_coarse(R,C,r,c,cld_label,classe,mir,class,cloud,THC);
			/* a new cloud is found */
                        }
                        
//changed emw 11/26/2011
//*ncloud   += cld_label;
*ncloud   = cld_label;

//printf("Ncloud in cloud_id_coarse: %ld",*ncloud);
}





/*
"cloud_id_fine" checks the whole region for cloudy pixels (THRSHLD_B<=THRSHLD<THRSHLD_C), and determines to which cloud it belongs.
*/
void    cloud_id_fine(long R,long C,long classe,float *mir,short *class,
		      long *cloud,float THC,long *ncloud,long *rcloud,
		      long *ccloud)
{
long    i,r,c;
long    label;

//printf("ncloud=%ld\n",*ncloud);
for (i=1;i<=*ncloud;i++)
        {
	r=rcloud[i];
	c=ccloud[i];
//printf("[%ld %ld %ld] ",i,r,c);
//printf("[%ld %ld] ",i,*ncloud);
//fflush(NULL);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,r,c,i,0); /*itter = 0 initialy*/
        }

}





/*
"cloud_lb_fine" searches for pixels with THRSHLD_B<=THRSHLD<THRSHLD_C and determinesto which cloud they belong.
*/
/**************
This version is the routine as given to emw, but seems not correct. The following routine is the emw
experimental one.
**************/
//void    cloud_lb_fine_orig(long R, long C, long classe, float *mir,
void    cloud_lb_fine(long R, long C, long classe, float *mir,
		      short *class, long *cloud, float THC, long *ncloud,
		      long k,long l,long cld_label, long itter)
{/* 1 */
int     GOON;
long    idx;

if (k == R || l == C || k == 0 || l == 0) return;


GOON = 0;

//if (k==48 && l>=89 && l<=101)
//printf("[%ld %ld %f] ",k,l,mir[k*C+l]);

idx     = k*C+l;
//printf("[%ld %ld %ld] ",k,l,itter);
//printf("[%ld %ld %f] ",k,l,mir[idx]);
//fflush(NULL);
if (itter < 50)
    { /* 2 */ /* these were added 7/31/2000 because I think they were left out */
        if (cloud[idx] != cld_label)
                {/* 3 */
                if (cloud[idx] != -cld_label)
                        {/* 4 */
                        if (cloud[idx] == 0)
                                {/* 5 *//* has not been checked before in coarse nor fine search */
//                                if (mir[idx] >= THC && class[idx] == classe) original emw
                                if (mir[idx] > THC)
                                        {/* 6 *//* it is a cloudy pixel that needs a cloud center */
                                        cloud[idx] = -cld_label;/* mark it as if belonging to cloud cld_label */
                                        GOON = 1;
                                        itter++; /* only once the edge has been reached. */
                                        }/* 6 */
                                }/* 5 */
                        }/* 4 */
                }/* 3 */
        else
                {/* 3 */ /* first time this coarsly identified pixel is entered */
                cloud[idx] = -cloud[idx];
                GOON = 1;
                }/* 3 */
    }/* 2 */
if (GOON)
        {
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k+1,l  ,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k-1,l  ,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k  ,l+1,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k  ,l-1,cld_label,itter);
        /* 8-connected */
/*******
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k+1,l+1,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k-1,l-1,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k-1,l+1,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k+1,l-1,cld_label,itter);
*******/
        }

}/* 1 */





/*
"cloud_lb_fine" searches for pixels with THRSHLD_B<=THRSHLD<THRSHLD_C and determinesto which cloud they belong.
*/
void    cloud_lb_fine_exp(long R, long C, long classe, float *mir,
//void    cloud_lb_fine(long R, long C, long classe, float *mir,
		      short *class, long *cloud, float THC, long *ncloud,
		      long k,long l,long cld_label, long itter)
{/* 1 */
int     GOON;
long    idx;

if (k == R || l == C || k == 0 || l == 0) return;


GOON = 0;

//printf("[%ld %ld %ld] ",k,l,cld_label);

idx     = k*C+l;
//printf("[%ld %ld %ld] ",k,l,itter);
//printf("[%ld %ld %f] ",k,l,mir[idx]);
//fflush(NULL);
if (itter < 50)
//    { /* 2 */ /* these were added 7/31/2000 because I think they were left out */
        if (cloud[idx] != cld_label)
                {/* 3 */
                if (cloud[idx] != -cld_label)
                        {/* 4 */
                        if (cloud[idx] == 0)
                                {/* 5 *//* has not been checked before in coarse nor fine search */
                                if (mir[idx] >= THC)
                                        {/* 6 *//* it is a cloudy pixel that needs a cloud center */
                                        cloud[idx] = -cld_label;/* mark it as if belonging to cloud cld_label */
                                        class[idx] = -classe;
                                        GOON = 1;
                                        itter++; /* only once the edge has been reached. */
                                        }/* 6 */
                                }/* 5 */
                        }/* 4 */
      					else
      									{/* 4 */
                				cloud[idx] = -cloud[idx];
                				class[idx] = -class[idx];
                				GOON = 1;
                				}/* 4 */
                }/* 3 */
//    }/* 2 */
if (GOON)
        {
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k+1,l  ,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k-1,l  ,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k  ,l+1,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k  ,l-1,cld_label,itter);
        /* 8-connected */
/*******
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k+1,l+1,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k-1,l-1,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k-1,l+1,cld_label,itter);
        cloud_lb_fine(R,C,classe,mir,class,cloud,THC,ncloud,k+1,l-1,cld_label,itter);
*******/
        }

}/* 1 */







/*************	label_region **************
A reimplementation of code above, simplified for just cloud
clustering
*/



void label_region(long R,long C,float *image,
		  long *cloud,float thresh,long *ncloud)
{
  FILE *fp_out=NULL;
  long    r,c;
  long    cld_label;
  long	  idx;
  extern long is_cloudCount, cld_lbCount;

  cld_label = *ncloud;
/* Initial cloud label is set at 0 but 1 is the first one used.*/

  for (r=0;r<R;r++)
    for (c=0;c<C;c++) {
      idx = r*C+c;
      if(is_cloud(R,C,r,c,image,cloud,thresh))
	{
	  cld_label++;  // First cloud label is 1
	  cloud_lb(R,C,r,c,cld_label,image,cloud,thresh); // a new cloud is found		
	}
    }
  *ncloud   += cld_label;
  
}





int is_cloud(long R,long C,long r, long c,
	     float *image,long *cloud,float thresh)
{
  int     outcome;
  long    idx;
  extern long is_cloudCount;

  idx     = r*C+c;
  outcome = 0;
  if (r == (R-1) || c == (C-1) || r == 0 || c == 0) return (outcome);

  if (image[idx]>thresh)
    if (cloud[idx] == 0)
      outcome = 1;
  return  (outcome);
}

/*****
"cloud_lb_coarse" labels a whole new cloud if pixel r,c is a cloudy pixel;
because checked cloud pixels are given zero thickness, relabeling is 
impossible.
*****/
void cloud_lb(long R,long C,long r, long c, long cld_label,
		     float *image,long *cloud,float thresh)
{
  long    i;
  extern long cld_lbCount;

//  cld_lbCount++;
//  printf("cld_lbCnt=%ld ",cld_lbCount);
//  fflush(NULL);
  i = r*C+c;
  
//  if (cld_label == 6) printf("%ld %ld  ",r,c);

  if (r == (R-1) || c == (C-1) || r == 0 || c == 0) return;

//  if (i > (2001.*1751.)) printf("%ld %ld %ld     ",r,c,i);
//  fflush(NULL);
  cloud[i] = cld_label;
  if (is_cloud(R,C,r+1,c,image,cloud,thresh)) 
    cloud_lb(R,C,r+1,c,cld_label,image,cloud,thresh);
  if (is_cloud(R,C,r-1,c,image,cloud,thresh))
    cloud_lb(R,C,r-1,c,cld_label,image,cloud,thresh);
  if (is_cloud(R,C,r,c+1,image,cloud,thresh))
    cloud_lb(R,C,r,c+1,cld_label,image,cloud,thresh);
  if (is_cloud(R,C,r,c-1,image,cloud,thresh))
    cloud_lb(R,C,r,c-1,cld_label,image,cloud,thresh);

// 8-connected
/***********
  if (is_cloud(R,C,r+1,c+1,image,cloud,thresh)) 
    cloud_lb(R,C,r+1,c+1,cld_label,image,cloud,thresh);
  if (is_cloud(R,C,r-1,c-1,image,cloud,thresh))
    cloud_lb(R,C,r-1,c-1,cld_label,image,cloud,thresh);
  if (is_cloud(R,C,r-1,c+1,image,cloud,thresh))
    cloud_lb(R,C,r-1,c+1,cld_label,image,cloud,thresh);
  if (is_cloud(R,C,r+1,c-1,image,cloud,thresh))
    cloud_lb(R,C,r+1,c-1,cld_label,image,cloud,thresh);
************/  
}


/*********
label_region_latlot takes an image and a threshold and does a simple clustering.
This uses the algorithm as described in Mapes and Houze (1993) whereby
the image is processed line by line, starting from the northern edge finding
"line clusters".  Then verically adjacent line clusters are attached when
they share a common column.

the _latlon monicor means that it assumes that the image is ordered such that
all values at a given latitude are adjacent in the array.
*********/
/*********
void label_region_latlon(void *image,long nlons,long nlats,long *cloud,
			 float thresh,long *ncloud)
{
  long i,j,iloc;
  long cloud_count = 1;
  
    // special case for the first line
  if (image[0] > thresh) //special case for the first point
    {
      cloud[0] = cloud_count;
      cloud_count++;
    }
  for (j=1; j<nlons; n++)
    {
    }
    
    // do the rest of the image
  for (i=1; i<nlats; i++)
    {
	  // special case for the first point in the line
      iloc = i*nlons;
      uploc = (i-1)*nlons;
      if (image[iloc] > thresh)
	{
	  if (cloud[uploc]) cloud[iloc] = cloud[uploc]
	  else
	    {
	      cloud[iloc] = cloud_count;
	      cloud_count++;
	    }
	}
      
      for (j=1; j<nlons; n++)
	{
	  iloc = (i*nlons)+nlons;
	  uploc = ((i-1)*nlons)+nlons);
	  lastloc = (i*nlons)+(nlons-1);
	  if (image[iloc] > thresh)
	    if (cloud[lastloc])
	      {
		cloud[iloc] = cloud[lastloc]
	    
	}
    }
  *ncloud = cloud_count;
} 
*********/
