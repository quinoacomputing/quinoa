
 
/* fres.h for ANSI C */
#ifndef FRES_H
#define FRES_H
 
#include "gofw.h"
#include "ftab.h"
#include "ffam.h"
#include "bitset.h"



typedef struct {
   ftab_Table *PVal [gofw_NTestTypes];
   bitset_BitSet Active;
   char *name;
} fres_Cont;



fres_Cont *fres_CreateCont (void);



void fres_DeleteCont (fres_Cont *res);

 

void fres_InitCont (ffam_Fam *fam, fres_Cont *res, int N,
                    int Nr, int j1, int j2, int jstep, char *nam);

 

void fres_PrintCont (fres_Cont *res);



void fres_FillTableEntryC (fres_Cont *res, gofw_TestArray pval, int N,
                           int irow, int icol);


typedef struct {
   ftab_Table *PLeft;
   ftab_Table *PRight;
   ftab_Table *PVal2;
   char *name;
} fres_Disc;



fres_Disc * fres_CreateDisc (void);



void fres_DeleteDisc (fres_Disc *res);

 

void fres_InitDisc (ffam_Fam *fam, fres_Disc *res,
                    int Nr, int j1, int j2, int jstep, char *nam);

 

void fres_PrintDisc (fres_Disc *res, lebool LR);



void fres_FillTableEntryD (fres_Disc *res, double pLeft, double pRight,
                           double pVal2, int irow, int icol);


typedef struct {
   ftab_Table *Exp;
   ftab_Table *Obs;
   ftab_Table *PLeft;
   ftab_Table *PRight;
   ftab_Table *PVal2;
   char *name;
} fres_Poisson;



fres_Poisson * fres_CreatePoisson (void);



void fres_DeletePoisson (fres_Poisson *res);

 

void fres_InitPoisson (ffam_Fam *fam, fres_Poisson *res,
                       int Nr, int j1, int j2, int jstep, char *nam);

 

void fres_PrintPoisson (fres_Poisson *res, lebool LR, lebool Ratio);



void fres_FillTableEntryPoisson (fres_Poisson *res, double Exp, double Obs, 
                                 double pLeft, double pRight, double pVal2,
                                 int irow, int icol);

 
#endif
 

