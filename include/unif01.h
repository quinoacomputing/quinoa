 
/*  unif01.h  for ANSI C  */
#ifndef UNIF01_H
#define UNIF01_H
 
#include "gdef.h"


typedef struct {
   void *state;
   void *param;
   char *name;
   double (*GetU01) (void *param, void *state);
   unsigned long (*GetBits) (void *param, void *state);
   void (*Write) (void *state);
} unif01_Gen;


#define unif01_MASK32  0xffffffffUL



#define unif01_NORM32  4294967296.0
#define unif01_INV32   2.328306436538696289e-10



extern lebool unif01_WrLongStateFlag;



double unif01_StripD (unif01_Gen *gen, int r);



long unif01_StripL (unif01_Gen *gen, int r, long d);



unsigned long unif01_StripB (unif01_Gen *gen, int r, int s);



void unif01_WriteNameGen (unif01_Gen *gen);



void unif01_WriteState (unif01_Gen *gen);



void unif01_WrLongStateDef (void);



unif01_Gen * unif01_CreateDummyGen (void);



void unif01_DeleteDummyGen (unif01_Gen *gen);



void unif01_DeleteGen (unif01_Gen *gen); 



unif01_Gen * unif01_CreateDoubleGen (unif01_Gen *gen, int s);



unif01_Gen * unif01_CreateDoubleGen2 (unif01_Gen *gen, double h);



unif01_Gen * unif01_CreateLacGen (unif01_Gen *gen, int k, long I[]);



unif01_Gen * unif01_CreateLuxGen (unif01_Gen *gen, int k, int L);



unif01_Gen * unif01_CreateBiasGen (unif01_Gen *gen, double a, double p);



unif01_Gen * unif01_CreateTruncGen (unif01_Gen *gen, int s);



unif01_Gen * unif01_CreateBitBlockGen (unif01_Gen *gen, int r, int s,
                                       int w);



void unif01_DeleteDoubleGen (unif01_Gen *gen);
void unif01_DeleteLacGen    (unif01_Gen *gen);
void unif01_DeleteLuxGen    (unif01_Gen *gen);
void unif01_DeleteBiasGen   (unif01_Gen *gen);
void unif01_DeleteTruncGen  (unif01_Gen *gen);
void unif01_DeleteBitBlockGen (unif01_Gen *gen);


typedef struct {
   unif01_Gen *gen1;
   unif01_Gen *gen2;
} unif01_Comb2_Param;



unif01_Gen * unif01_CreateCombAdd2 (unif01_Gen *gen1, unif01_Gen *gen2,
                                    char *name);



unif01_Gen * unif01_CreateCombAdd3 (unif01_Gen *gen1, unif01_Gen *gen2,
                                    unif01_Gen *gen3, char *name);



unif01_Gen * unif01_CreateCombXor2 (unif01_Gen *gen1, unif01_Gen *gen2,
                                    char *name);



unif01_Gen * unif01_CreateCombXor3 (unif01_Gen *gen1, unif01_Gen *gen2,
                                    unif01_Gen *gen3, char *name);



void unif01_DeleteCombGen (unif01_Gen *gen);


unif01_Gen * unif01_CreateParallelGen (int k, unif01_Gen *gen[], int L);



void unif01_DeleteParallelGen (unif01_Gen *gen);



unif01_Gen *unif01_CreateExternGen01 (char *name, double (*gen01)(void*,void*),
                                      unsigned long (*gen01_bits)(void*,void*));



unif01_Gen *unif01_CreateExternGenBits (char *name,
                                        unsigned int (*genB)(void));



unif01_Gen *unif01_CreateExternGenBitsL (char *name,
                                         unsigned long (*genB)(void));



void unif01_DeleteExternGen01 (unif01_Gen * gen);
void unif01_DeleteExternGenBits (unif01_Gen * gen);
void unif01_DeleteExternGenBitsL (unif01_Gen * gen);


typedef struct {
   unif01_Gen *gen;
   long n;
   double time;
   double mean;
   lebool fU01;
   } unif01_TimerRec;



void unif01_TimerGen (unif01_Gen *gen, unif01_TimerRec *timer, long n,
                      lebool fU01);



void unif01_TimerSumGen (unif01_Gen *gen, unif01_TimerRec *timer, long n,
                         lebool fU01);



void unif01_WriteTimerRec (unif01_TimerRec *timer);



void unif01_TimerGenWr (unif01_Gen *gen, long n, lebool fU01);



void unif01_TimerSumGenWr (unif01_Gen *gen, long n, lebool fU01);

 
#endif
 

