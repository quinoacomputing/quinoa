
 
/* statcoll.h  for ANSI C  */

#ifndef STATCOLL_H
#define STATCOLL_H
 


typedef struct {
   double *V;
   long Dim;
   long NObs;
   char *Desc;
} statcoll_Collector;


statcoll_Collector * statcoll_Create (long N, const char Desc[]);



statcoll_Collector * statcoll_Delete (statcoll_Collector *S);



void statcoll_Init (statcoll_Collector *S, long N);



void statcoll_SetDesc (statcoll_Collector *S, const char Desc[]);



void statcoll_AddObs (statcoll_Collector *S, double x);



void statcoll_Write (statcoll_Collector *S, int k, int p1, int p2, int p3);



double statcoll_Average (statcoll_Collector *S);



double statcoll_Variance (statcoll_Collector *S);



double statcoll_AutoCovar (statcoll_Collector *S, int k);



double statcoll_Covar (statcoll_Collector *S1, statcoll_Collector *S2);

 
#endif
 

