 
/*  utouzin.h  for ANSI C */
#ifndef UTOUZIN_H
#define UTOUZIN_H
 
#include "unif01.h"


unif01_Gen * utouzin_CreateMRG00a (long s1, long s2, long s3, long s4,
                                   long s5);



unif01_Gen * utouzin_CreateMRG00b (long s1, long s2, long s3, long s4,
                                   long s5, long s6);



unif01_Gen * utouzin_CreateMRG00c (long s1, long s2, long s3, long s4,
                                   long s5, long s6, long s7);



unif01_Gen * utouzin_CreateMRG00d (long s1, long s2, long s3, long s4,
                                   long s5, long s6, long s7, long s8);




unif01_Gen * utouzin_CreateMRG00e (long s10, long s11, long s12, 
                                   long s20, long s21, long s22);




unif01_Gen * utouzin_CreateMRG00f (long s10, long s11, long s12, 
                                   long s20, long s21, long s22);



unif01_Gen * utouzin_CreateMRG00g (long s10, long s11, long s12, 
                                   long s20, long s21, long s22,
                                   long s30, long s31, long s32);



unif01_Gen * utouzin_CreateMRG00h (long s10, long s11, long s12, long s13,
                                   long s20, long s21, long s22, long s23);



void utouzin_DeleteGen (unif01_Gen * gen);
 
#endif
 

