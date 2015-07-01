 
/* sspacings.h  for ANSI C */
#ifndef SSPACINGS_H
#define SSPACINGS_H
 
#include "unif01.h"
#include "sres.h"


typedef struct {

   sres_Basic **LogCEMu;
   sres_Basic **LogCAMu;
   sres_Basic **SquareCEMu;
   sres_Basic **SquareCAMu;


   double *LogCESig_sVal, *LogCESig_pVal;
   double *LogCASig_sVal, *LogCASig_pVal;
   double *SquareCESig_sVal, *SquareCESig_pVal;
   double *SquareCASig_sVal, *SquareCASig_pVal;


   int imax;


   char *name;


   statcoll_Collector ** Collectors;
   int smax, step;


} sspacings_Res;



sspacings_Res * sspacings_CreateRes (void);



void sspacings_DeleteRes (sspacings_Res *res);


void sspacings_SumLogsSpacings (unif01_Gen *gen, sspacings_Res *res,
                                long N, long n, int r, int m);



void sspacings_SumSquaresSpacings (unif01_Gen *gen, sspacings_Res *res,
                                   long N, long n, int r, int m);



void sspacings_ScanSpacings (unif01_Gen *gen, sspacings_Res *res,
                             long N, long n, int r, double d);



void sspacings_AllSpacings (unif01_Gen *gen, sspacings_Res *res,
                            long N, long n, int r, int m0, int m1, int d,
                            int LgEps);



void sspacings_AllSpacings2 (unif01_Gen *gen, sspacings_Res *res,
                             long N, long n, int r, int m0, int m1, int d,
                             int LgEps);

 
#endif
 

