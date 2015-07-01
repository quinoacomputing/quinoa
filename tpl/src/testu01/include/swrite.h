 
/* swrite.h for ANSI C */
#ifndef SWRITE_H
#define SWRITE_H
 
#include "gdef.h"
#include "chrono.h"
#include "unif01.h"
#include "sres.h"


extern lebool swrite_Basic;           /* Prints basic results           */
extern lebool swrite_Parameters;      /* Prints details on parameters   */
extern lebool swrite_Collectors;      /* Prints statistical collectors  */
extern lebool swrite_Classes;         /* Prints classes for ChiSquare   */
extern lebool swrite_Counters;        /* Prints counters                */



extern lebool swrite_Host;


extern char swrite_ExperimentName[];



void swrite_SetExperimentName (char Name[]);



void swrite_Head (unif01_Gen *gen, char *TestName, long N, long n, int r);



void swrite_Final (unif01_Gen *gen, chrono_Chrono *Timer);



void swrite_NormalSumTest (long N, sres_Basic *res);



void swrite_AddStrChi (char S[], int len, long d);



void swrite_Chi2SumTest (long N, sres_Chi2 *res);



void swrite_Chi2SumTestb (long N, double sval, double pval, long deg);
 
#endif
 
