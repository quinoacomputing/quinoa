 
/* ftab.h for ANSI C */
#ifndef FTAB_H
#define FTAB_H
 
#include "ffam.h"
#include "unif01.h"

 
typedef enum {
   ftab_NotInit,              /* Uninitialized */
   ftab_pVal1,                /* One-sided p-value */
   ftab_pVal2,                /* Two-sided p-value */
   ftab_pLog10,               /* Logarithm of p-value in base 10 */
   ftab_pLog2,                /* Logarithm of p-value in base 2 */
   ftab_Integer,              /* Integer number */
   ftab_Real,                 /* Real number */
   ftab_String                /* String */
} ftab_FormType;
 

typedef struct {
   double **Mat;
   int *LSize;
   int Nr, Nc;
   int j1, j2, jstep;
   ftab_FormType Form;
   char *Desc;
   char **Strings;
   int Ns;
} ftab_Table;


ftab_Table *ftab_CreateTable (int Nr, int j1, int j2, int jstep,
                              char *Desc, ftab_FormType Form, int Ns);



void ftab_DeleteTable (ftab_Table *T);



void ftab_SetDesc (ftab_Table *T, char *Desc);



void ftab_InitMatrix (ftab_Table *T, double x);



typedef void (*ftab_CalcType) (ffam_Fam *fam, void *res, void *cho,
                               void *par, int LSize, int j,
                               int irow, int icol);



void ftab_MakeTables (ffam_Fam *fam, void *res, void *cho, void *par,
                      ftab_CalcType Calc, 
                      int Nr, int j1, int j2, int jstep);


typedef enum {
   ftab_Plain,                /* To print tables in plain text */
   ftab_Latex                 /* To print tables in Latex format */
} ftab_StyleType;



extern ftab_StyleType ftab_Style;



extern double ftab_Suspectp;



extern int ftab_SuspectLog2p;



void ftab_PrintTable (ftab_Table *T);



void ftab_PrintTable2 (ftab_Table *T1, ftab_Table *T2, lebool ratioF);


 
#endif
 

