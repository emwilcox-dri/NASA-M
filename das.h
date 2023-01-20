/***********
Created by emw 5/27/99 to make the DAS code more available to the rest
of my code.  This way I can envoke subsets of the DAS routines to do the
label_region stuff.
***********/


/**********
PERFORMING THE LABEL_REGION ROUTINE
just call cloud_lb_coarse with class field set to zero.  Note that it
finds pixels > threshold (NOT >=)
**********/


/***********************
set the number of detect steps, thresholds and number of spread sub-steps
for the go_das routine here.
***/

/****** Remy's Meteosat settings ******/
// set all CATEGORYs to 0 to and set the class array to zeros, to turn
// off the class stuff.

//Also used in first round of terra/aqua runs (Summer 2012)
// /Volumes/shores/terraqua/winter_monsoon ... /summer_monsoon
/****************************
#define STEP	    3
// step 1 
#define DETECT_0	-240
#define SPREAD_0	-260
#define CATEGORY_0	1
#define DELTA_SPREAD_0	-2

// step 2 
#define DETECT_1	-250
#define SPREAD_1	-270
#define CATEGORY_1	2
#define DELTA_SPREAD_1	-2

// step 3 
#define DETECT_2	-260
#define SPREAD_2	-280
#define CATEGORY_2	3
#define DELTA_SPREAD_2	-2

// step 4 
#define DETECT_3	-280
#define SPREAD_3	-300
#define CATEGORY_3	0
#define DELTA_SPREAD_3	-10
****************************/
/****** End Remy's Meteosat settings ******/


//Roca and Ramanathan values
// /Volumes/shores/terraqua/run4/winter_monsoon ... /summer_monsoon
/*******************************/
#ifdef DEEP_CLOUDS
#define STEP	    3
// step 1 
#define DETECT_0	-220
#define SPREAD_0	-240
#define CATEGORY_0	1
#define DELTA_SPREAD_0	-2

// step 2 
#define DETECT_1	-235
#define SPREAD_1	-255
#define CATEGORY_1	2
#define DELTA_SPREAD_1	-2

// step 3 
#define DETECT_2	-255
#define SPREAD_2	-255
#define CATEGORY_2	3
#define DELTA_SPREAD_2	-2

// step 4 
#define DETECT_3	-280
#define SPREAD_3	-300
#define CATEGORY_3	0
#define DELTA_SPREAD_3	-10
#endif //DEEP_CLOUDS
/******************************/

//Boer and Ramanathan values
#ifdef BOER_RAM_CLOUDS
#define STEP        4
// step 1 
#define DETECT_0        -240
#define SPREAD_0        -260
#define CATEGORY_0      1
#define DELTA_SPREAD_0  -6.6

// step 2 
#define DETECT_1        -255
#define SPREAD_1        -275
#define CATEGORY_1      2
#define DELTA_SPREAD_1  -6.6

// step 3 
#define DETECT_2        -270
#define SPREAD_2        -280
#define CATEGORY_2      3
#define DELTA_SPREAD_2  -3.3

// step 4 
#define DETECT_3        -280
#define SPREAD_3        -285
#define CATEGORY_3      0
#define DELTA_SPREAD_3  -3
#endif //BOER_RAM_CLOUDS


#define FREE(x) if(x) { free(x); x = NULL; }

#ifdef VERBOSE
void    go_das(long R,long C,float *ir,short *class,long *cloud,
	       long *ncloud,int verbose,char *yyyyddd,char *hhmm,char *outPath);
#endif //end if VERBOSE
#ifndef VERBOSE
void    go_das(long R,long C,float *ir,short *class,long *cloud,
	       long *ncloud,int verbose);
#endif
	      
int     is_cloud_coarse(long R,long C,long r, long c, long classe,
			float *mir,short *class,long *cloud,float THC);
	
void    cloud_lb_coarse(long R,long C,long r, long c, long cld_label,
			long classe,float *mir, short *class,long *cloud,
			float THC);
			
void    cloud_id_coarse(long R,long C,long classe,float *mir,short *class,
			long *cloud,float THC,long *ncloud);

void    cloud_id_fine(long R,long C,long classe,float *mir,short *class,
		      long *cloud,float THC,long *ncloud,long *rcloud,
		      long *ccloud);
		      
//void    cloud_lb_fine_orig(long R, long C, long classe, float *mir, short *class,
void    cloud_lb_fine(long R, long C, long classe, float *mir, short *class,
		      long *cloud, float THC, long *ncloud, long k,long l,
		      long cld_label, long itter);
		      
void    cloud_lb_fine_exp(long R, long C, long classe, float *mir, short *class,
//void    cloud_lb_fine(long R, long C, long classe, float *mir, short *class,
		      long *cloud, float THC, long *ncloud, long k,long l,
		      long cld_label, long itter);


/******* label_region ******/

void label_region(long R,long C,float *image,
		  long *cloud,float thresh,long *ncloud);
		  
int is_cloud(long R,long C,long r, long c,
	     float *image,long *cloud,float thresh);
	     
void cloud_lb(long R,long C,long r, long c, long cld_label,
		     float *image,long *cloud,float thresh);
		     
