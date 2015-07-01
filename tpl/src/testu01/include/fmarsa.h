 
/* fmarsa.h for ANSI C */
#ifndef FMARSA_H
#define FMARSA_H
 
#include "ffam.h"
#include "fres.h"
#include "fcho.h"


extern long fmarsa_Maxn, fmarsa_MaxL;


typedef struct {
   fres_Cont *GCD;
   fres_Cont *NumIter; 
} fmarsa_Res2;



fmarsa_Res2 * fmarsa_CreateRes2 (void);



void fmarsa_DeleteRes2 (fmarsa_Res2 *res);


fcho_Cho * fmarsa_CreateBirthEC (long N, int t, double EC);



void fmarsa_DeleteBirthEC (fcho_Cho *cho);


void fmarsa_SerialOver1 (void);



void fmarsa_CollisionOver1 (void);



void fmarsa_BirthdayS1 (ffam_Fam *fam, fres_Poisson *res, fcho_Cho2 *cho,
                        long N, int r, int t, int p,
                        int Nr, int j1, int j2, int jstep);



void fmarsa_MatrixR1 (ffam_Fam *fam, fres_Cont *res, fcho_Cho2 *cho,
                      long N, long n, int r, int s, int L,
                      int Nr, int j1, int j2, int jstep);



void fmarsa_GCD1 (ffam_Fam *fam, fmarsa_Res2 *res, fcho_Cho *cho,
                      long N, int r, int s,
                      int Nr, int j1, int j2, int jstep);


 
#endif
 

