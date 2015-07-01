 
/* num.h for ANSI C */

#ifndef NUM_H
#define NUM_H
 
#include "gdef.h"


#define num_Pi     3.14159265358979323846


#define num_ebase  2.7182818284590452354


#define num_Rac2   1.41421356237309504880


#define num_1Rac2  0.70710678118654752440


#define num_Ln2    0.69314718055994530941


#define num_1Ln2   1.44269504088896340737


#define num_MaxIntDouble   9007199254740992.0


#define num_MaxTwoExp   64


extern double num_TwoExp[];   


#define num_MAXTENNEGPOW   16


extern double num_TENNEGPOW[];


#define num_Log2(x) (num_1Ln2 * log(x)) 



long num_RoundL (double x);



double num_RoundD (double x);



int num_IsNumber (char S[]);



void num_IntToStrBase (long k, long b, char S[]);



void num_Uint2Uchar (unsigned char output[], unsigned int input[], int L);



void num_WriteD (double x, int i, int j, int k);



void num_WriteBits (unsigned long x, int k);



long num_MultModL (long a, long s, long c, long m);



double num_MultModD (double a, double s, double c, double m);



long num_InvEuclid (long m, long z);



unsigned long num_InvExpon (int E, unsigned long z);

 
#endif
 

