 
/* smarsa.h for ANSI C */
#ifndef SMARSA_H
#define SMARSA_H
 
#include "unif01.h"
#include "sres.h"


extern double smarsa_Maxk;


typedef struct {
   sres_Basic *Bas;
   sres_Poisson *Pois;
} smarsa_Res;



smarsa_Res * smarsa_CreateRes (void);



void smarsa_DeleteRes (smarsa_Res *res);


typedef struct {
   sres_Chi2 *GCD;
   sres_Chi2 *NumIter;
} smarsa_Res2;



smarsa_Res2 * smarsa_CreateRes2 (void);



void smarsa_DeleteRes2 (smarsa_Res2 *res);



void smarsa_SerialOver (unif01_Gen *gen, sres_Basic *res,
                        long N, long n, int r, long d, int t);



void smarsa_CollisionOver (unif01_Gen *gen, smarsa_Res *res,
                           long N, long n, int r, long d, int t);



void smarsa_Opso (unif01_Gen *gen, smarsa_Res *res,
                  long N, int r, int p);



void smarsa_CAT (unif01_Gen *gen, sres_Poisson *res,
                 long N, long n, int r, long d, int t, long S[]);



void smarsa_CATBits (unif01_Gen *gen, sres_Poisson *res, long N, long n,
                     int r, int s, int L, unsigned long Key);



void smarsa_BirthdaySpacings (unif01_Gen *gen, sres_Poisson *res,
                              long N, long n, int r, long d, int t, int p);



void smarsa_MatrixRank (unif01_Gen *gen, sres_Chi2 *res,
                        long N, long n, int r, int s, int L, int k);



void smarsa_Savir2 (unif01_Gen *gen, sres_Chi2 *res,
                    long N, long n, int r, long m, int t);



void smarsa_GCD (unif01_Gen *gen, smarsa_Res2 *res,
                 long N, long n, int r, int s);

 
#endif
 

