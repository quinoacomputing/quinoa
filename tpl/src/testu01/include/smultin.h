 
/* smultin.h for ANSI C */
#ifndef SMULTIN_H
#define SMULTIN_H
 
#include "gdef.h"
#include "fmass.h"
#include "statcoll.h"
#include "gofw.h"
#include "unif01.h"


#define smultin_MAX_DELTA  8



#define smultin_MAXB  10


typedef struct {

   double Maxk;


   double SeuilHash;


   double HashLoad;


   double SeuilEColl;


   double SeuilCOverDense;
   double SeuilCOverNorSup;
   double SeuilCOverNorInf;
   double SeuilCOverSparse;


} smultin_Envir;


extern smultin_Envir smultin_env;


#ifdef USE_LONGLONG
   typedef ulonglong smultin_CellType;
#else
   typedef double smultin_CellType;
#endif



typedef smultin_CellType (*smultin_GenerCellType) (unif01_Gen *, int, int,
                                                   long);



smultin_CellType smultin_GenerCellSerial (unif01_Gen *gen, int r, int t,
                                          long d);



smultin_CellType smultin_GenerCellSerial2 (unif01_Gen *gen, int r, int t,
                                           long d);



smultin_CellType smultin_GenerCellPermut (unif01_Gen *gen, int r, int t,
                                          long junk);



smultin_CellType smultin_GenerCellMax (unif01_Gen *gen, int r, int t,
                                       long junk);



smultin_CellType smultin_GenerCellSerialBits (unif01_Gen *gen, int r, int s,
                                              long L);



typedef struct {

   int NbDelta;
   double ValDelta [smultin_MAX_DELTA];


   smultin_GenerCellType GenerCell;


   int bmax;


} smultin_Param;


smultin_Param * smultin_CreateParam (int NbDelta, double ValDelta[],
                                     smultin_GenerCellType GenerCell,
                                     int bmax);



void smultin_DeleteParam (smultin_Param *par);


extern smultin_Param smultin_ParamDefault;


typedef enum {
   smultin_CollNotInit,           /* Not initialized */
   smultin_CollExact,             /* Exact distribution */
   smultin_CollNormal,            /* Normal approximation */
   smultin_CollPoissonSparse,     /* Poisson approximation: sparse case */
   smultin_CollPoissonDense       /* Poisson approximation: dense case */
} smultin_CollApproxType;



typedef struct {

   lebool Hashing;


   smultin_CellType NbCellsTotal;


   lebool Over;


   smultin_CollApproxType CollApprox;


   double Mu [smultin_MAX_DELTA];
   double Sigma [smultin_MAX_DELTA];


   double EsEmpty;


   long CountSize;
   long Count1Size;
   long *Count;
   long *Count1;
   smultin_CellType *Cell;
   smultin_CellType *Cell1;


   long NbSize;
   long Nb1Size;
   smultin_CellType *Nb;
   smultin_CellType *Nb1;


   smultin_CellType NbCells [smultin_MAXB + 1];
   double EsCells [smultin_MAXB + 1];
   smultin_CellType WbCells [smultin_MAXB + 1];


   double NbCollisions;


   statcoll_Collector *Collector [smultin_MAX_DELTA];


   gofw_TestArray sVal2 [smultin_MAX_DELTA];
   gofw_TestArray pVal2 [smultin_MAX_DELTA];


   double pCollLeft, pCollRight;


   double pColl;


   double pEmpty;


   double pWb [smultin_MAXB + 1];


   int NbDeltaOld;


   double *TabFj[smultin_MAX_DELTA];


   int nLimit;


   lebool flagTab;


} smultin_Res;


smultin_Res * smultin_CreateRes (smultin_Param *par);



void smultin_DeleteRes (smultin_Res *res);


typedef double (*smultin_MNTermeType) (double, double, long);



double smultin_MNTermeColl (double, double, long j);



double smultin_MNTermePowDiv (double Delta, double E, long j);


double smultin_MNTermeKhi2 (double, double E, long j);



double smultin_MNTermeLogLikhood (double, double E, long j);



void smultin_MultinomMuSigma (long n, double k, double theta1,
                              double theta2, smultin_MNTermeType F,
                              double *Mu, double *Sigma);


void smultin_PowDivMomCorChi (double Delta, long n, double k,
                              double *MuC, double *SigmaC);



void smultin_PowDivMom (double Delta, long n, double k,
                        double NbExp, double *Mu, double *Sigma);



fmass_INFO smultin_CreateCollisions (long n, smultin_CellType k);



void smultin_DeleteCollisions (fmass_INFO W);



double smultin_FDistCollisions (fmass_INFO W, long c);



double smultin_FBarCollisions (fmass_INFO W, long c);



double smultin_CollisionsTerm (fmass_INFO W, long c);


void smultin_Multinomial (unif01_Gen *gen, smultin_Param *par,
   smultin_Res *res, long N, long n, int r, long d, int t, lebool Sparse);



void smultin_MultinomialOver (unif01_Gen *gen, smultin_Param *par,
   smultin_Res *res, long N, long n, int r, long d, int t, lebool Sparse);



void smultin_MultinomialBits (unif01_Gen *gen, smultin_Param *par,
   smultin_Res *res, long N, long n, int r, int s, int L, lebool Sparse);



void smultin_MultinomialBitsOver (unif01_Gen *gen, smultin_Param *par,
   smultin_Res *res, long N, long n, int r, int s, int L, lebool Sparse);

 
#endif
 

