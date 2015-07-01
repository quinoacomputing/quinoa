 
/* fcho.h for ANSI C */
#ifndef FCHO_H
#define FCHO_H
 


typedef struct {
   void *param;
   double (*Choose) (void *param, long, long);
   void (*Write) (void *param, long, long);
   char *name;
} fcho_Cho;



long fcho_ChooseParamL (fcho_Cho *cho, long min, long max, long i, long j);


typedef struct {
   fcho_Cho *Chon;
   fcho_Cho *Chop2;
} fcho_Cho2;



fcho_Cho2 * fcho_CreateCho2 (fcho_Cho *Chon, fcho_Cho *Chop2);



void fcho_DeleteCho2 (fcho_Cho2 *cho);


typedef double (*fcho_FuncType) (double);



double fcho_Linear (double x);



double fcho_LinearInv (double x);



double fcho_2Pow (double x);



fcho_Cho * fcho_CreateSampleSize (double a, double b, double c,
                                  fcho_FuncType F, char *name);



void fcho_DeleteSampleSize (fcho_Cho *cho);



int fcho_Chooses (int r, int s, int resol);

 
#endif
 

