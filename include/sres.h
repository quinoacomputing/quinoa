 
/* sres.h  for ANSI C */
#ifndef SRES_H
#define SRES_H
 
#include "gofw.h"
#include "statcoll.h"


typedef struct {

   statcoll_Collector *sVal1, *pVal1;


   gofw_TestArray sVal2, pVal2;


   char *name;


} sres_Basic;



sres_Basic * sres_CreateBasic (void);



void sres_DeleteBasic (sres_Basic *res);



void sres_InitBasic (sres_Basic *res, long N, char *nam);



void sres_GetNormalSumStat (sres_Basic *res);


typedef struct {

   statcoll_Collector *sVal1;


   double sVal2;


   double pLeft, pRight, pVal2;


   char *name;


} sres_Disc;



sres_Disc * sres_CreateDisc (void);



void sres_DeleteDisc (sres_Disc *res);



void sres_InitDisc (sres_Disc *res, long N, char *nam);


typedef struct {

   double *NbExp;
   long *Count;
   long *Loc;
   long jmin;
   long jmax;


   long degFree;


   statcoll_Collector *sVal1, *pVal1;


   gofw_TestArray sVal2, pVal2;


   char *name;


} sres_Chi2;



sres_Chi2 * sres_CreateChi2 (void);



void sres_DeleteChi2 (sres_Chi2 *res);



void sres_InitChi2 (sres_Chi2 *res, long N, long jmax, char *nam);



void sres_GetChi2SumStat (sres_Chi2 *res);


typedef struct {

   double Lambda, Mu;


   statcoll_Collector *sVal1;


   double sVal2;


   double pLeft, pRight, pVal2;


   char *name;


} sres_Poisson;



sres_Poisson * sres_CreatePoisson (void);



void sres_DeletePoisson (sres_Poisson *res);



void sres_InitPoisson (sres_Poisson *res, long N, double Lambda, char *nam);

 
#endif
 

