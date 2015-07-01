 
#ifndef SVARIA_H
#define SVARIA_H
 
#include "unif01.h"
#include "sres.h"



extern lebool svaria_Timer;


void svaria_SampleMean (unif01_Gen *gen, sres_Basic *res,
                        long N, long n, int r);



void svaria_SampleCorr (unif01_Gen *gen, sres_Basic *res,
                        long N, long n, int r, int k);



void svaria_SampleProd (unif01_Gen *gen, sres_Basic *res,
                        long N, long n, int r, int t);



void svaria_SumLogs (unif01_Gen *gen, sres_Chi2 *res,
                     long N, long n, int r);



void svaria_WeightDistrib (unif01_Gen *gen, sres_Chi2 *res, long N, long n,
                           int r, long k, double alpha, double beta);



void svaria_CollisionArgMax (unif01_Gen *gen, sres_Chi2 *res,
                             long N, long n, int r, long k, long m);



void svaria_SumCollector (unif01_Gen *gen, sres_Chi2 *res,
                          long N, long n, int r, double g);



void svaria_AppearanceSpacings (unif01_Gen *gen, sres_Basic *res,
                               long N, long Q, long K, int r, int s, int L);

 
#endif
 

