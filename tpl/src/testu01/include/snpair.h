 
/* snpair.h  for ANSI C */
#ifndef SNPAIR_H
#define SNPAIR_H
 
#include "gdef.h"
#include "chrono.h"
#include "statcoll.h"
#include "unif01.h"



#define snpair_MAXM 512



#define snpair_MAXREC 12


typedef struct {

   int Seuil1;        /* Recursion threshold for snpair_FindClosePairs */
   int Seuil2;        /* Recursion threshold for snpair_CheckBoundary  */
   int Seuil3;        /* L1 = 1 + (lg (n/Seuil3)) / sqrt(t) */
   int Seuil4;        /* L2 = 1 + (lg (n/Seuil4)) / sqrt(t) */


} snpair_Envir;


extern snpair_Envir snpair_env;



extern long snpair_MaxNumPoints;


typedef double * snpair_PointType;



typedef snpair_PointType * snpair_PointTableType;



typedef struct snpair_Res snpair_Res;

typedef void (*snpair_DistanceType) (snpair_Res *res, snpair_PointType,
                                     snpair_PointType);



typedef void (*snpair_VerifPairsType) (snpair_Res *res, snpair_PointType [],
                                       long, long, int, int);



typedef void (*snpair_MiniProcType) (snpair_Res *res, snpair_PointType [],
                                     long, long, long, long, int, int);



enum snpair_StatType {
   snpair_NP,
   snpair_NPS,
   snpair_NPPR,
   snpair_mNP,
   snpair_mNP1,
   snpair_mNP1S,
   snpair_mNP2,
   snpair_mNP2S,
   snpair_NJumps,
   snpair_BB,
   snpair_BM,
   snpair_StatType_N
};




struct snpair_Res {

   long n;

   lebool CleanFlag;               /* If TRUE, free all memory */
   void *work;


   snpair_PointTableType Points[1 + snpair_MAXREC];


   int NumClose;


   double * CloseDist;


   snpair_DistanceType   Distance;
   snpair_VerifPairsType VerifPairs;
   snpair_MiniProcType   MiniProc;


   statcoll_Collector * Yn;
   statcoll_Collector * Y;
   statcoll_Collector * U;
   statcoll_Collector * V;
   statcoll_Collector * S;
   statcoll_Collector * TheWn;
   statcoll_Collector * TheWni;
   statcoll_Collector * ThepValAD;
   statcoll_Collector * BitMax;


   double sVal [snpair_StatType_N];
   double pVal [snpair_StatType_N];


};



snpair_Res * snpair_CreateRes (void);



void snpair_DeleteRes (snpair_Res *res);


void snpair_QuickSort (snpair_PointType A[], long l, long r, int c);



void snpair_VerifPairs0 (snpair_Res *res, snpair_PointType A[], long r, long s,
                         int np, int c);



void snpair_VerifPairs1 (snpair_Res *res, snpair_PointType A[], long r, long s,
                         int np, int c);



void snpair_MiniProc0 (snpair_Res *res, snpair_PointType A[], long r, long s,
                       long u, long v, int np, int c);



void snpair_MiniProc1 (snpair_Res *res, snpair_PointType A[], long r, long s,
                       long u, long v, int np, int c);



void snpair_CheckBoundary (snpair_Res *res, long r, long s, long u, long v,
                           int nr, int nrb, int np, int c);



void snpair_FindClosePairs (snpair_Res *res, long r, long s, int nr,
                            int np, int c);


void snpair_DistanceCP (snpair_Res *res, snpair_PointType, snpair_PointType);



void snpair_WriteDataCP (unif01_Gen *gen, char *TestName, long N, long n,
                         int r, int t, int p, int m, lebool Torus);



void snpair_WriteResultsCP (unif01_Gen *gen, chrono_Chrono *Timer,
                            snpair_Res *res, long N, long m);



extern lebool snpair_mNP2S_Flag;



void snpair_ClosePairs (unif01_Gen *gen, snpair_Res *res,
                        long N, long n, int r, int t, int p, int m);



void snpair_ClosePairs1 (unif01_Gen *gen, snpair_Res *res,
                         long N, long n, int r, int t, int p, int m);


void snpair_ReTestY (long N, long n, int m, double t0, double t1);


void snpair_DistanceCPBitM (snpair_Res *res, snpair_PointType,
                            snpair_PointType);




void snpair_ClosePairsBitMatch (unif01_Gen *gen, snpair_Res *res,
                                long N, long n, int r, int t);


extern lebool snpair_TimeBB;



void snpair_DistanceBB (snpair_Res *res, snpair_PointType, snpair_PointType);




void snpair_WriteDataBB (unif01_Gen *gen, char *TestName, long N, long n,
                         int r, int t, int p, lebool Tor, int L1, int L2);



void snpair_WriteResultsBB (unif01_Gen *gen, chrono_Chrono *Timer,
                            snpair_Res *res, long N);



void snpair_BickelBreiman (unif01_Gen *gen, snpair_Res *res, long N,
                           long n, int r, int t, int p, lebool Torus);

#endif

