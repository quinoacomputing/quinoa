/*************************************************************************\
 *
 * Package:        ProbDist
 * File:           fbar.c
 * Environment:    ANSI C
 *
 * Copyright (c) 2002 Pierre L'Ecuyer, DIRO, Université de Montréal.
 * e-mail: lecuyer@iro.umontreal.ca
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted without a fee for private, research,
 * academic, or other non-commercial purposes.
 * Any use of this software in a commercial environment requires a
 * written licence from the copyright owner.
 *
 * Any changes made to this package must be clearly identified as such.
 *
 * In scientific publications which used this software, a reference to it
 * would be appreciated.
 *
 * Redistributions of source code must retain this copyright notice
 * and the following disclaimer.
 *
 * THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
\*************************************************************************/

#include "fbar.h"
#include "fdist.h"

#include "num.h"
#include "num2.h"
#include "util.h"
#include "gdef.h"

#include <math.h>
#include <float.h>

double fdist_belog(double);
extern const double fdist_XINF;
extern const double fdist_XBIG;
extern const double fdist_XBIGM;

/* EpsArray[j]: Epsilon required for j decimal degits of precision */
static const double EpsArray[] = {
   0.5, 0.5E-1, 0.5E-2, 0.5E-3, 0.5E-4, 0.5E-5, 0.5E-6, 0.5E-7, 0.5E-8,
   0.5E-9, 0.5E-10, 0.5E-11, 0.5E-12, 0.5E-13, 0.5E-14, 0.5E-15, 0.5E-16,
   0.5E-17, 0.5E-18, 0.5E-19, 0.5E-20, 0.5E-21, 0.5E-22, 0.5E-23, 0.5E-24,
   0.5E-25, 0.5E-26, 0.5E-27, 0.5E-28, 0.5E-29, 0.5E-30, 0.5E-31, 0.5E-32,
   0.5E-33, 0.5E-34, 0.5E-35
};

/* Compute IMAX extra terms in the tails of discrete distributions */
static const long IMAX = 20;



/*=========================================================================*/

double fbar_Unif (double x)
{
   if (x <= 0.0)
      return 1.0;
   if (x >= 1.0)
      return 0.0;
   return 1.0 - x;
}


/*=========================================================================*/

double fbar_Expon (double x)
{
   if (x <= 0.0)
      return 1.0;
   if (x >= fdist_XBIGM)
      return 0.0;
   return exp (-x);
}


/*=========================================================================*/

double fbar_Weibull (double c, double x)
{
   double temp;
   util_Assert (c > 0.0, "fbar_Weibull:   c <= 0");
   if (x <= 0.0)
      return 1.0;
   if (x >= DBL_MAX_EXP * FLT_RADIX && c >= 1.0)
      return 0.0;
   temp = c*log(x);
   if (temp >= DBL_MAX_EXP * num_Ln2)
      return 0.0;
   temp = exp(temp);
   return (exp (-temp));
}


/*=========================================================================*/

double fbar_Logistic (double x)
{
   if (x <= -fdist_XBIG) {
      return 1.0;
   }
   if (x >= fdist_XBIG) {
      return exp (-x);
   }
   return 1.0 / (1.0 + exp (x));
}


/*=========================================================================*/

double fbar_Pareto (double c, double x)
{

   util_Assert (c > 0.0, "fbar_Pareto:   c <= 0");
   if (x <= 1.0)
      return 1.0;
   return (pow (x, -c));
}

/**************************************************************************/

double fbar_Normal1 (double x)
/*
 * Returns P[X >= x] = 1 - F(x) where F is the normal distribution by
 * computing the complementary distribution directly; it is thus more
 * precise in the tail.
 */
{
   static const double A[25] = {
      6.10143081923200418E-1,
     -4.34841272712577472E-1,
      1.76351193643605501E-1,
     -6.07107956092494149E-2,
      1.77120689956941145E-2,
     -4.32111938556729382E-3,
      8.54216676887098679E-4,
     -1.27155090609162743E-4,
      1.12481672436711895E-5,
      3.13063885421820973E-7,
     -2.70988068537762022E-7,
      3.07376227014076884E-8,
      2.51562038481762294E-9,
     -1.02892992132031913E-9,
      2.99440521199499394E-11,
      2.60517896872669363E-11,
     -2.63483992417196939E-12,
     -6.43404509890636443E-13,
      1.12457401801663447E-13,
      1.7281533389986098E-14,
     -4.2641016949424E-15,
     -5.4537197788E-16,
      1.5869760776E-16,
      2.08998378E-17,
     -0.5900E-17
   };
   const double kk = 5.30330085889910643300;      /* 3.75 Sqrt(2) */
   double y, t;
   int Neg;

   if (x >= fdist_XBIG) {
      return 0.0;
   }
   if (x <= -fdist_XBIG) {
      return 1.0;
   }

   if (x >= 0.0)
      Neg = 0;
   else {
      Neg = 1;
      x = -x;
   }

   t = (x - kk) / (x + kk);
   y = num2_EvalCheby (A, 24, t);
   y = y * exp (-x * x / 2.0) / 2.0;

   if (Neg == 1)
      return (1.0 - y);
   else
      return (y);
}


/*=========================================================================*/

double fbar_Normal2 (double x)
{
   static const double V[121] = {
        1.2533141373155,      1.137490921203605,      1.037824575853727,
      0.951527192071207,     0.8763644564536924,     0.8105337152790306,
     0.7525711790634081,     0.7012808218544303,     0.6556795424187987,
       0.61495459615093,     0.5784303460476312,     0.5455421356582171,
     0.5158156382179634,     0.4888504415275737,     0.4643069280394423,
     0.4418957328326002,     0.4213692292880546,     0.4025146181296722,
     0.3851482907984348,     0.3691112106902635,     0.3542651113297938,
     0.3404893532870847,     0.3276783146905521,       0.31573921586941,
     0.3045902987101033,     0.2941592970402893,      0.284382146748493,
     0.2752018941576065,     0.2665677689682238,     0.2584343943120386,
     0.2507611114439651,      0.243511400615456,     0.2366523829135607,
      0.230154390478801,     0.2239905946538289,     0.2181366833614714,
     0.2125705804420318,     0.2072722008565011,     0.2022232366330547,
     0.1974069692375194,     0.1928081047153158,     0.1884126285076003,
     0.1842076773079702,     0.1801814257143918,     0.1763229857571027,
     0.1726223176578506,     0.1690701504076941,     0.1656579109468773,
     0.1623776608968675,     0.1592220399363674,     0.1561842150339759,
      0.153257834853479,     0.1504369887362691,     0.1477161697413935,
      0.145090241289131,     0.1425544070104023,     0.1401041834530503,
     0.1377353753382303,     0.1354440530967635,     0.1332265324471292,
     0.1310793558044918,     0.1289992753343376,      0.126983237485437,
     0.1250283688553504,     0.1231319632579323,     0.1212914698765462,
      0.119504482399253,     0.1177687290432979,     0.1160820633859823,
     0.1144424559276431,      0.112847986320103,     0.1112968362007359,
     0.1097872825783083,     0.1083176917221132,     0.1068865135106745,
     0.1054922762005562,     0.1041335815795983,     0.1028091004723001,
     0.1015175685681028,     0.1002577825460485,    0.09902859647173194,
    0.09782891844465691,    0.09665770747608191,    0.09551397057921558,
    0.09439676005522439,    0.09330517095996169,    0.09223833873763035,
    0.09119543700877471,    0.09017567550106469,    0.08917829811230435,
    0.08820258109597616,    0.08724783136042988,    0.08631338487354936,
    0.08539860516539227,    0.08450288192189578,    0.08362562966329139,
    0.08276628650136918,    0.08192431297018954,    0.08109919092525536,
    0.08029042250654048,    0.07949752916111721,    0.07872005072144664,
    0.07795754453568722,    0.07720958464664668,    0.07647576101624852,
    0.07575567879261112,    0.07504895761704659,    0.07435523096847724,
    0.07367414554294564,    0.07300536066605566,    0.07234854773633338,
    0.07170338969763433,    0.07106958053885212,    0.07044682481930167,
    0.06983483721825942,    0.06923334210724434,    0.06864207314371742,
    0.06806077288496332,     0.0674891924209997,    0.06692709102543307,
    0.06637423582325017
};

   int j;
   lebool negatif;
   double t, u, z, h;
   double r, r1, r2, r3, r4, r5, r6, r7, r8;

   if (x >= fdist_XBIG) {
      return 0.0;
   }
   if (x <= -fdist_XBIG) {
      return 1.0;
   }
   if (x < 0.0) {
      negatif = TRUE;
      x = -x;
   } else {
      negatif = FALSE;
   }
   j = (int) (8.0 * x + 0.5);
   if (j > 120)
      j = 120;
   z = 0.125 * j;
   h = x - z;
   r = V[j];
   r1 = r * z - 1.0;
   r2 = 0.5 * (r + z * r1);
   r3 = (r1 + z * r2) / 3.0;
   r4 = 0.25 * (r2 + z * r3);
   r5 = 0.2 * (r3 + z * r4);
   r6 = (r4 + z * r5) / 6.0;
   r7 = (r5 + z * r6) / 7.0;
   r8 = 0.125 * (r6 + z * r7);
   t = r + h * (r1 + h * (r2 + h * (r3 + h * (r4 + h * (r5 + h * (r6 +
                h * (r7 + h * r8)))))));
   u = t * exp (-0.5 * x * x - 0.9189385332046727);
   if (negatif)
      return 1.0 - u;
   else
      return u;
}


/*=========================================================================*/
#ifdef HAVE_ERF

double fbar_Normal3 (double x)
{
   return 0.5 * erfc (x * num_1Rac2);
}

#endif
/*=========================================================================*/

double fbar_BiNormal1 (double x, double y, double rho, int ndig)
{
   return fdist_BiNormal1(-x, -y, rho, ndig);
}


/*=========================================================================*/

double fbar_BiNormal2 (double x, double y, double rho)
{
   return fdist_BiNormal2 (-x, -y, rho);
}


/*=========================================================================*/

double fbar_LogNormal (double mu, double sigma, double x)
{
   util_Assert (sigma > 0.0, "fbar_LogNormal:  sigma  <= 0");
   if (x <= 0.0)
      return 1.0;
   return fbar_Normal1 ((log (x) - mu) / sigma);
}


/*=========================================================================*/

double fbar_JohnsonSB (double alpha, double beta, double a, double b,
   double x)
{
   util_Assert (beta > 0.0, "fbar_JohnsonSB:  beta  <= 0");
   util_Assert (b > a, "fbar_JohnsonSB:  b  <= a");
   if (x <= a)
      return 1.0;
   if (x >= b)
      return 0.0;
   return fbar_Normal1 (alpha + beta * log ((x - a) / (b - x)));
}


/*=========================================================================*/

double fbar_JohnsonSU (double alpha, double beta, double x)
{
   const double XLIM = 1.0e10;
   double r;
   lebool negative = FALSE;
   util_Assert (beta > 0.0, "fbar_JohnsonSU:  beta  <= 0");
   if (x < 0.0) {
      negative = TRUE;
      x = -x;
   }
   /* compute r = x + sqrt (x * x + 1) */
   if (x < XLIM)
      r = x + sqrt (x * x + 1.0);
   else
      r = 2.0 * x;
   if (negative)
      r = 1.0 / r;

   if (r > 0.0)
      return fbar_Normal1 (alpha + beta * log (r));
   else
      return 1.0;
}


/*=========================================================================*/

static double F2AD[103];          /* Tables for the approximation of the */
static double CoAD[103];          /* Anderson-Darling distribution */

static void AndersonDarlingInit (void)
{
   F2AD[0] = 0.0;                F2AD[1] = 1.7315E-10;
   F2AD[2] = 2.80781E-5;         F2AD[3] = 1.40856E-3;
   F2AD[4] = 9.58772E-3;         F2AD[5] = 2.960552E-2;
   F2AD[6] = 6.185146E-2;        F2AD[7] = 1.0357152E-1;
   F2AD[8] = 1.5127241E-1;       F2AD[9] = 2.0190317E-1;
   F2AD[10] = 2.5318023E-1;      F2AD[11] = 3.0354278E-1;
   F2AD[12] = 3.5200015E-1;      F2AD[13] = 3.9797537E-1;
   F2AD[14] = 4.4117692E-1;      F2AD[15] = 4.8150305E-1;
   F2AD[16] = 5.1897375E-1;      F2AD[17] = 5.5368396E-1;
   F2AD[18] = 5.8577199E-1;      F2AD[19] = 6.1539864E-1;
   F2AD[20] = 6.4273362E-1;      F2AD[21] = 6.6794694E-1;
   F2AD[22] = 6.9120359E-1;      F2AD[23] = 7.126605E-1;
   F2AD[24] = 7.3246483E-1;      F2AD[25] = 7.507533E-1;
   F2AD[26] = 7.6765207E-1;      F2AD[27] = 7.8327703E-1;
   F2AD[28] = 7.9773426E-1;      F2AD[29] = 8.1112067E-1;
   F2AD[30] = 8.2352466E-1;      F2AD[31] = 8.3502676E-1;
   F2AD[32] = 8.4570037E-1;      F2AD[33] = 8.5561231E-1;
   F2AD[34] = 8.6482346E-1;      F2AD[35] = 8.7338931E-1;
   F2AD[36] = 8.8136046E-1;      F2AD[37] = 8.8878306E-1;
   F2AD[38] = 8.9569925E-1;      F2AD[39] = 9.0214757E-1;
   F2AD[40] = 9.081653E-1;       F2AD[41] = 9.1378043E-1;
   F2AD[42] = 9.1902284E-1;      F2AD[43] = 9.2392345E-1;
   F2AD[44] = 9.2850516E-1;      F2AD[45] = 9.3279084E-1;
   F2AD[46] = 9.3680149E-1;      F2AD[47] = 9.4055647E-1;
   F2AD[48] = 9.440736E-1;       F2AD[49] = 9.4736933E-1;
   F2AD[50] = 9.5045883E-1;      F2AD[51] = 9.5335611E-1;
   F2AD[52] = 9.5607414E-1;      F2AD[53] = 9.586249E-1;
   F2AD[54] = 9.6101951E-1;      F2AD[55] = 9.6326825E-1;
   F2AD[56] = 9.6538067E-1;      F2AD[57] = 9.6736563E-1;
   F2AD[58] = 9.6923135E-1;      F2AD[59] = 9.7098548E-1;
   F2AD[60] = 9.7263514E-1;      F2AD[61] = 9.7418694E-1;
   F2AD[62] = 9.7564704E-1;      F2AD[63] = 9.7702119E-1;
   F2AD[64] = 9.7831473E-1;      F2AD[65] = 9.7953267E-1;
   F2AD[66] = 9.8067966E-1;      F2AD[67] = 9.8176005E-1;
   F2AD[68] = 9.827779E-1;       F2AD[69] = 9.8373702E-1;
   F2AD[70] = 9.8464096E-1;      F2AD[71] = 9.8549304E-1;
   F2AD[72] = 9.8629637E-1;      F2AD[73] = 9.8705386E-1;
   F2AD[74] = 9.8776824E-1;      F2AD[75] = 9.8844206E-1;
   F2AD[76] = 9.8907773E-1;      F2AD[77] = 9.8967747E-1;
   F2AD[78] = 9.9024341E-1;      F2AD[79] = 9.9077752E-1;
   F2AD[80] = 9.9128164E-1;      F2AD[81] = 9.9175753E-1;
   F2AD[82] = 9.9220682E-1;      F2AD[83] = 9.9263105E-1;
   F2AD[84] = 9.9303165E-1;      F2AD[85] = 9.9340998E-1;
   F2AD[86] = 9.9376733E-1;      F2AD[87] = 9.9410488E-1;
   F2AD[88] = 9.9442377E-1;      F2AD[89] = 9.9472506E-1;
   F2AD[90] = 9.9500974E-1;      F2AD[91] = 9.9527876E-1;
   F2AD[92] = 9.95533E-1;        F2AD[93] = 9.9577329E-1;
   F2AD[94] = 9.9600042E-1;      F2AD[95] = 9.9621513E-1;
   F2AD[96] = 9.964181E-1;       F2AD[97] = 0.99661;
   F2AD[98] = 9.9679145E-1;      F2AD[99] = 9.9696303E-1;
   F2AD[100] = 9.9712528E-1;     F2AD[101] = 9.9727872E-1;
   F2AD[102] = 9.9742384E-1;

   CoAD[0] = 0.0;
   CoAD[1] = 0.0;                 CoAD[2] = 0.0;            
   CoAD[3] = 0.0;                 CoAD[4] = 0.0;            
   CoAD[5] = -1.87E-3;            CoAD[6] = 0.00898;        
   CoAD[7] = 0.0209;              CoAD[8] = 0.03087;        
   CoAD[9] = 0.0377;              CoAD[10] = 0.0414;        
   CoAD[11] = 0.04386;            CoAD[12] = 0.043;         
   CoAD[13] = 0.0419;             CoAD[14] = 0.0403;        
   CoAD[15] = 0.038;              CoAD[16] = 3.54804E-2;    
   CoAD[17] = 0.032;              CoAD[18] = 0.0293;        
   CoAD[19] = 2.61949E-2;         CoAD[20] = 0.0228;        
   CoAD[21] = 0.0192;             CoAD[22] = 1.59865E-2;    
   CoAD[23] = 0.0129;             CoAD[24] = 0.0107;        
   CoAD[25] = 8.2464E-3;          CoAD[26] = 0.00611;       
   CoAD[27] = 0.00363;            CoAD[28] = 1.32272E-3;    
   CoAD[29] = -5.87E-4;           CoAD[30] = -2.75E-3;      
   CoAD[31] = -3.95248E-3;        CoAD[32] = -5.34E-3;      
   CoAD[33] = -6.892E-3;          CoAD[34] = -8.10208E-3;   
   CoAD[35] = -8.93E-3;           CoAD[36] = -9.552E-3;     
   CoAD[37] = -1.04605E-2;        CoAD[38] = -0.0112;       
   CoAD[39] = -1.175E-2;          CoAD[40] = -1.20216E-2;   
   CoAD[41] = -0.0124;            CoAD[42] = -1.253E-2;     
   CoAD[43] = -1.27076E-2;        CoAD[44] = -0.0129;       
   CoAD[45] = -1.267E-2;          CoAD[46] = -1.22015E-2;   
   CoAD[47] = -0.0122;            CoAD[48] = -1.186E-2;     
   CoAD[49] = -1.17218E-2;        CoAD[50] = -0.0114;       
   CoAD[51] = -1.113E-2;          CoAD[52] = -1.08459E-2;   
   CoAD[53] = -0.0104;            CoAD[54] = -9.93E-3;      
   CoAD[55] = -9.5252E-3;         CoAD[56] = -9.24E-3;      
   CoAD[57] = -9.16E-3;           CoAD[58] = -8.8004E-3;    
   CoAD[59] = -8.63E-3;           CoAD[60] = -8.336E-3;     
   CoAD[61] = -8.10512E-3;        CoAD[62] = -7.94E-3;      
   CoAD[63] = -7.71E-3;           CoAD[64] = -7.55064E-3;   
   CoAD[65] = -7.25E-3;           CoAD[66] = -7.11E-3;      
   CoAD[67] = -6.834E-3;          CoAD[68] = -0.0065;       
   CoAD[69] = -6.28E-3;           CoAD[70] = -6.11008E-3;   
   CoAD[71] = -5.86E-3;           CoAD[72] = -5.673E-3;     
   CoAD[73] = -5.35008E-3;        CoAD[74] = -5.11E-3;      
   CoAD[75] = -4.786E-3;          CoAD[76] = -4.59144E-3;   
   CoAD[77] = -4.38E-3;           CoAD[78] = -4.15E-3;      
   CoAD[79] = -4.07696E-3;        CoAD[80] = -3.93E-3;      
   CoAD[81] = -3.83E-3;           CoAD[82] = -3.74656E-3;   
   CoAD[83] = -3.49E-3;           CoAD[84] = -3.33E-3;      
   CoAD[85] = -3.20064E-3;        CoAD[86] = -3.09E-3;      
   CoAD[87] = -2.93E-3;           CoAD[88] = -2.78136E-3;   
   CoAD[89] = -2.72E-3;           CoAD[90] = -2.66E-3;      
   CoAD[91] = -2.56208E-3;        CoAD[92] = -2.43E-3;      
   CoAD[93] = -2.28E-3;           CoAD[94] = -2.13536E-3;   
   CoAD[95] = -2.083E-3;          CoAD[96] = -1.94E-3;      
   CoAD[97] = -1.82E-3;           CoAD[98] = -1.77E-3;      
   CoAD[99] = -1.72E-3;           CoAD[100] = -1.71104E-3;  
   CoAD[101] = -1.741E-3;         CoAD[102] = -0.0016;

}

double fbar_AndersonDarling (long N, double X)
{
   /* This function is not very precise for x < 0.05 */
   const double h = 0.05;         /* the step of the interpolation table */
   static int ADFlag = 0;
   double q;
   double Res, Cor;
   int i;

   if (N == 1) {
      if (X <= 0.38629436111989)
         return 1.0;
      if (X >= fdist_XBIGM)
         return 0.0;
      if (X < 6.0) {
         q = 1.0 - 4.0 * exp(-X - 1.0);
         return 1.0 - sqrt (q);
      } else {
         q = 4.0 * exp(-X - 1.0);
         return 0.5*q*(1.0 + 0.25*q*(1.0 + 0.5*q*(1.0 + 0.125*q*(5.0 + 3.5*q))));
      }
   }

   if (N <= 0) {
      util_Warning (1, "fbar_AndersonDarling:   N < 1");
      return -1.0;
   }

   if (X > 10.0)
      /* Sinclair-Spurr upper tail approximation (3.5) */
      return 1.732 * exp(-X) / sqrt(num_Pi * X);

   if (X > 5.0) {
      /* asymptotic X:  our empirical fit */
      Res = exp (-0.56 - 1.06 * X);
      q = exp (-1.03 - 1.06 * X);         /* Empirical correction in 1/N */
      return Res + q / N;
   }

   if (X <= 0.2)
      return 1.0 - fdist_AndersonDarling (N, X);

   if (ADFlag == 0) {
      AndersonDarlingInit ();
      ADFlag = 1;
   }

   i = 1 + (int) (X / h);
   q = X / h - i;

   /* Newton backwards quadratic interpolation */
   Res = (F2AD[i - 2] - 2.0 * F2AD[i - 1] + F2AD[i]) * q * (q + 1.0) / 2.0
      + (F2AD[i] - F2AD[i - 1]) * q + F2AD[i];

   /* Empirical correction in 1/N */
   Cor = (CoAD[i] * (q + 1.0) - CoAD[i - 1] * q) / N;

   Res = 1.0 - Res - Cor;
   if (Res >= 1.0)
      return 1.0;
   if (Res <= 0.0)
      return 0.0;
   return Res;
}


/*=========================================================================*/

double fbar_ChiSquare1 (long N, double x)
/*
 * Returns an approximation of the complementary Chi square cdf (N degrees
 * of freedom). Similar to p:116 of W.J.Kennedy Jr and J.E.Gentle.
 * Statistical computing, Dekker, New York, 1980. More precise in the
 * tail than simply returning  1 - fdist_ChiSquare.
 */
{
   const double XBIG_CHI = 2000.0;
   const double tiers = 0.33333333333333333;
   const double pt2 = 0.22222222222222222;
   const double moinshuit = -8.3;
   const double gam = 0.8862269254527579825931;
   double H, E, DemiX, Terme, Sommation, Y;
   long i;

   util_Assert (N > 0, "Calling fbar_ChiSquare1 with N < 1");
   if (x <= 0.0)
      return 1.0;
   if (N >= 150) {
      if (x >= N * fdist_XBIG)
         return 0.0;
   } else {
      if (x >= XBIG_CHI)
         return 0.0;
   }

   if (N > 1000) {
      if (x < 2.0)
         return 1.0;
      x = (pow ((x / N), tiers) - (1.0 - pt2 / N)) / sqrt (pt2 / N);
      if (x > 35.0)
         return 0.0;
      if (x <= moinshuit)
         return 1.0;
      return fbar_Normal1 (x);
   }

   DemiX = x / 2.0;

   if (!(N & 1)) {             /* even N */
      Terme = exp (-DemiX);
      Sommation = Terme;
      for (i = 1; i < N / 2; i++) {
	 Terme = Terme * DemiX / i;
	 Sommation += Terme;
      }
      Y = Sommation;

   } else {
      H = 2.0 * fbar_Normal1 (sqrt (x));
      if (N == 1)
	 return H;

      E = exp (-DemiX);
      Terme = sqrt (DemiX) * E / gam;
      for (i = 3; i < N; i += 2) {
	 H += Terme;
	 Terme = Terme * DemiX * 2.0 / i;
      }
      Y = H + Terme;
   }

   if (Y > 1.0)
      return 1.0;
   else 
      return Y;
}


/*=========================================================================*/

double fbar_ChiSquare2 (long n, int d, double x)
{
   util_Assert (n > 0, "fbar_ChiSquare2:   n <= 0");
   if (x <= 0.0)
      return 1.0;
   return fbar_Gamma (n / 2.0, d, x / 2.0);
}


/*=========================================================================*/

double fbar_Gamma (double alpha, int d, double x)
{
   const double aLIM = 1.0E5;
   const double RENORM = 1.0E100;
   const double EPS = EpsArray[d];
   double V[6];
   double v, res, A, B, R, term, dif;
   int i;

   util_Assert (alpha > 0.0, "fbar_Gamma:   a <= 0");
   util_Assert (d > 0, "fbar_Gamma:   d <= 0");
   util_Assert (d < 16, "fbar_Gamma:   d > 15");
   if (x <= 0.0)
      return 1.0;
   if (1.0 == alpha)
      return fbar_Expon (x);

   if (alpha >= 70.0) {
      if (x >= alpha * fdist_XBIG)
         return 0.0;
   } else {
      if (x >= fdist_XBIGM)
         return 0.0;
   }

   if (alpha >= aLIM) {
      double d2 = x + 1.0/3.0 - alpha - 0.02/alpha;
      double S = alpha - 1.0/2.0;
      double z = d2 * sqrt((1 + fdist_belog(S/x))/x);
      return fbar_Normal1 (z);
   }

   if (x <= 1.0 || x < alpha)
      return 1.0 - fdist_Gamma (alpha, d, x);

   v = exp (alpha * log (x) - x - num2_LnGamma (alpha));

   A = 1.0 - alpha;
   B = A + x + 1.0;
   term = 0.0;
   V[0] = 1.0;
   V[1] = x;
   V[2] = x + 1.0;
   V[3] = x * B;
   res = V[2] / V[3];

   do {
      A += 1.0;
      B += 2.0;
      term += 1.0;
      V[4] = B * V[2] - A * term * V[0];
      V[5] = B * V[3] - A * term * V[1];
      if (V[5] != 0.0) {
         R = V[4] / V[5];
         dif = fabs (res - R);
         if (dif <= EPS * R)
            return (v * res);
         res = R;
      }
      for (i = 0; i < 4; i++)
         V[i] = V[i + 2];
      if (fabs (V[4]) >= RENORM) {
         for (i = 0; i < 4; i++)
            V[i] /= RENORM;
      }
   } while (1);

   /* to eliminate a warning from the compiler; never reached */
   return 0.0;
}


/*=========================================================================*/

static double KSPlusbarAsymp (long n, double x)
{
   /* Compute the probability of the KSPlus distribution using 
      an asymptotic formula */
   double t = (6.0*n*x + 1);
   double z = t*t/(18.0*n);
   double v = 1.0 - (2.0*z*z - 4.0*z - 1.0)/(18.0*n);
   if (v <= 0.0)
      return 0.0;
   v = v*exp(-z);
   if (v >= 1.0)
      return 1.0;
   return v;
}


/*-------------------------------------------------------------------------*/

static double KSPlusbarUpper (long n, double x)
{
   /* Compute the probability of the KSPlus distribution in the upper
      tail using Smirnov's stable formula */
   const double EPSILON = 1.0E-10;
   double q;
   double Sum = 0.0;
   double term;
   double t;
   double LogCom;
   double LOGJMAX;
   int j;
   int jmax = (int)(n - n*x);

   /* We must avoid log(0) for j = jmax and q ~ 1.0 */
   if ((1.0 - x - (double)jmax / n) <= 0.0)
      jmax--;

   j = jmax/2;
   LogCom = num2_LnFactorial((int)n) - num2_LnFactorial(j) -
            num2_LnFactorial((int)(n-j));
   LOGJMAX = LogCom;

   while (j > 0) {
      q = (double)j / n + x;
      term = LogCom + (j - 1)*log (q) + (n - j)*num2_log1p (-q);
      t = exp (term);
      Sum += t;
      LogCom += log ((double)j / (n - j + 1));
      if (t <= Sum*EPSILON)
         break;
      j--;
   }

   j = jmax/2;
   LogCom = LOGJMAX + log ((double)(n - j)/(j + 1));
   j++;

   while (j <= jmax) {
      q = (double)j / n + x;
      term = LogCom + (j - 1)*log(q) + (n - j)*num2_log1p(-q);
      t = exp (term);
      Sum += t;
      LogCom += log ((double)(n - j)/(j + 1));
      if (t <= Sum*EPSILON)
         break;
      j++;
   }

   Sum *= x;
   /* add the term j = 0 */
   Sum += exp (n*num2_log1p (-x));
   return Sum;
}


/*=========================================================================*/

double fbar_KSPlus (long N, double x)
{
   const double NxParam = 6.5;    /* frontier: alternating series */
   const long NParam = 4000;      /* frontier: non-alternating series */
   const long NAsymp = 200000;    /* frontier: asymptotic */

   util_Assert (N > 0, "Calling fbar_KSPlus with N < 1");
   if (x <= 0.0)
      return 1.0;
   if ((x >= 1.0) || (N*x*x >= 370.0))
      return 0.0;
   if (N == 1)
      return 1.0 - x;

   if (N * x <= NxParam)
      return 1.0 - fdist_KSPlus (N, x);

   if (N >= NAsymp)
      return KSPlusbarAsymp (N, x);

   if ((N <= NParam) || (N*x*x > 1.0))
      return KSPlusbarUpper(N, x);

/*   return (1.0 - 2.0*x/3.0)*exp(-2.0*N*x*x);  */
   return KSPlusbarAsymp (N, x);
}


/*=========================================================================*/

static double KSSpecial (long n, double x)
{
#define NLIM 20

   if ((n * x * x >= 370.0) || (x >= 1.0))
      return 0.0;
   if (x <= 0.5 / n)
      return 1.0;
   if (n == 1)
      return 2.0 - 2.0 * x;

   if (x <= 1.0 / n) {
      double w;
      double t = 2.0 * x - 1.0 / n;
      if (n <= NLIM) {
         w = num2_Factorial ((int) n);
         return 1.0 - w * pow (t, (double) n);
      }
      w = num2_LnFactorial ((int) n) + n * log (t);
      return 1.0 - exp (w);
   }

   if (x >= 1.0 - 1.0 / n) {
      return 2.0 * pow (1.0 - x, (double) n);
   }
   return -1.0;
}

#undef NLIM
/*-------------------------------------------------------------------------*/

double fbar_KS1 (long n, double x)
{
   double v = KSSpecial(n, x);
   if (v >= 0.0)
      return v;

   if (n <= 400) {
      if (n*x*x < 4.0)
         return 1.0 - fdist_KS1(n, x);
      else 
         return 2.0 * KSPlusbarUpper(n, x);
   }

   if (n*x*x >= 2.2) {
      if (n <= 200000)
         return 2.0 * KSPlusbarUpper(n, x);
      return 2.0*KSPlusbarAsymp (n, x);
   }
 
   return 1.0 - fdist_KS1(n, x);
}


/*=========================================================================*/

double fbar_CramerMises (long N, double x)
{
   return 1.0 - fdist_CramerMises (N, x);
}


double fbar_WatsonG (long N, double x)
{
   return 1.0 - fdist_WatsonG (N, x);
}


/*=========================================================================*/

double fbar_WatsonU (long N, double x)
{
/*
 * Only the asymptotic form has been implemented. In the trivial case
 * N = 1, we simply return 0.5
 */
   const double xSepare = 0.15;
   if (x <= 0.0)
      return 1.0;
   if (x >= fdist_XBIG)
      return 0.0;

   if (N == 1)                    /* N = 1, degenerate case */
      return 0.5;

   if (x > xSepare) {
      /* this series converges rapidly for x > 0.15 */
      const int JMAX = 10;
      int j;
      double signe;
      double v;
      double terme;
      double somme;
      v = exp (-(x * 2.0 * num_Pi * num_Pi));
      signe = 1.0;
      somme = 0.0;
      j = 1;
      do {
         terme = pow (v, (double) j * j);
         somme += signe * terme;
         signe = -signe;
         ++j;
      } while (!(terme < DBL_EPSILON || j > JMAX));
      util_Warning (j > JMAX, "fbar_WatsonU:  sum1 has not converged");
      v = 2.0 * somme;
      if (v <= 0.0)
         return 0.0;
      return v;
   }

   return 1.0 - fdist_WatsonU (N, x);
}


/*=========================================================================*/




/******************************\
 *
 *  DISCRETE DISTRIBUTIONS
 *
\******************************/


/*=========================================================================*/

static const double epsilonScan = 1.0E-7;

static double ScanGlaz (long N, double d, long m)
{
   long j, jmoy;
   double temp;
   double jr, jm1r, Nr = N;
   int signe;
   double q = 1.0 - d;
   double Q4, Q3, Q2, Q1;
   double Bin, BinMoy;

   jmoy = (long) ((N + 1) * d);    /* max term of the Binomial */
   if (jmoy < m - 1)
      jmoy = m - 1;

   /*---------------------------------------------------------*/
   /* Compute Q1: formula (2.5) in Glaz (1989)                */
   /* Compute Q2: formula (A.6) in Berman and Eagleson (1985) */
   /* Compute Q3, Q4 : Theorem (3.2) in Glaz (1989)           */
   /*---------------------------------------------------------*/

   /* compute the probability of term j = jmoy */
   Q1 = 0.0;
   for (j = 1; j <= jmoy; j++) {
      jr = j;
      Q1 += log (Nr - jr + 1.0) - log (jr);
   }
   Q1 += jmoy * log (d) + (Nr - jmoy) * log (q);
   BinMoy = exp (Q1);
   Q1 = BinMoy;
   jm1r = jmoy - m + 1;
   if ((jmoy - m + 1) & 1)
      signe = -1;
   else
      signe = 1;
   Q2 = signe * BinMoy;
   Q3 = signe * BinMoy * (2.0 - jm1r * jm1r + jm1r);
   Q4 = signe * BinMoy * (jm1r + 1.0) * (jm1r + 2.0) * (6.0 + jm1r * jm1r -
      5.0 * jm1r);

   /* compute the probability of terms j > jmoy */
   if ((jmoy - m + 1) & 1)
      signe = -1;
   else
      signe = 1;

   jm1r = jmoy - m + 1;
   Bin = BinMoy;
   for (j = jmoy + 1; j <= N; j++) {
      jr = j;
      jm1r += 1.0;
      signe = -signe;
      Bin = (Bin * (Nr - jr + 1.0) * d) / (jr * q);
      if (Bin < epsilonScan)
         break;
      Q1 += Bin;
      Q2 += signe * Bin;
      Q3 += signe * Bin * (2.0 - jm1r * jm1r + jm1r);
      Q4 += signe * Bin * (jm1r + 1.0) * (jm1r + 2.0) * (6.0 + jm1r * jm1r -
         5.0 * jm1r);
   }

   Q1 = 1.0 - Q1;
   Q3 /= 2.0;
   Q4 /= 12.0;
   if (m == 3) {
      /* Problem with this formula; I do not get the same results as Glaz */
      Q4 = ((Nr * (Nr - 1.0) * d * d * pow (q, Nr - 2.0)) / 8.0
         + Nr * d * 2.0 * pow (1.0 - 2.0 * d, Nr - 1.0))
         - 4.0 * pow (1.0 - 2.0 * d, Nr);
      if (d < 1.0 / 3.0) {
         Q4 += Nr * d * 2.0 * pow (1.0 - 3.0 * d, Nr - 1.0)
               + 4.0 * pow (1.0 - 3.0 * d, Nr);
      }
   }
   /* compute probability: Glaz, equations (3.2) and (3.3) */
   Q3 = Q1 - Q2 - Q3;
   Q4 = Q3 - Q4;
   /* when the approximation is bad, avoid overflow */
   temp = log (Q3) + (Nr - m - 2.0) * log (Q4 / Q3);
   if (temp >= 0.0)
      return 0.0;
   if (temp < (-30.0))
      return 1.0;
   Q4 = exp (temp);
   return 1.0 - Q4;
}

/*----------------------------------------------------------------------*/

static double ScanWNeff (long N, double d, long m)
{
   double q = 1.0 - d;
   double temp;
   double Bin;
   double Sum;
   long j;

   /*--------------------------------------*/
   /* Anderson-Titterington: equation (4)  */
   /*--------------------------------------*/

   /* compute the probability of term j = m */
   Sum = 0.0;
   for (j = 1; j <= m; j++) {
      Sum += log ((double) (N - j + 1)) - log ((double) j);
   }
   Sum += m * log (d) + (N - m) * log (q);
   Bin = exp (Sum);
   temp = (m / d - N - 1.0) * Bin;
   Sum = Bin;

   /* compute the probability of terms j > m */
   for (j = m + 1; j <= N; j++) {
      Bin *= (N - j + 1) * d / (j * q);
      if (Bin < epsilonScan)
         break;
      Sum += Bin;
   }
   Sum = 2.0 * Sum + temp;
   return Sum;
}

/*----------------------------------------------------------------------*/

static double ScanAsympt (long N, double d, long m)
{
   double Kappa;
   double temp;
   double Theta;
   double Sum;

   /*--------------------------------------------------------------*/
   /* Anderson-Titterington: asymptotic formula after equation (4) */
   /*--------------------------------------------------------------*/

   Theta = sqrt (d / (1.0 - d));
   temp = sqrt ((double) N);
   Kappa = m / (d * temp) - temp;
   temp = Theta * Kappa;
   temp = temp * temp / 2.0;
   Sum = 2.0 * fbar_Normal1 (Theta * Kappa) +
      (Kappa * Theta * exp (-temp)) / (d * sqrt (2.0 * num_Pi));
   return Sum;
}

/*----------------------------------------------------------------------*/

double fbar_Scan (long N, double d, long m)
{
   double mu;
   double prob;

   util_Assert (N >= 2, "Calling fbar_Scan with N < 2");
   util_Assert (d > 0.0 && d < 1.0,
      "Calling fbar_Scan with d outside (0,1)");
   if (m > N)
      return 0.0;
   if (m <= 1)
      return 1.0;
   if (m <= 2) {
      if ((N - 1) * d >= 1.0)
         return 1.0;
      return (1.0 - pow (1.0 - (N - 1) * d, (double) N));
   }
   if (d >= 0.5 && m <= (N + 1) / 2.0)
      return 1.0;
   if (d > 0.5)
      return (-1.0);              /* Error */
   /* util_Assert (d <= 0.5, "Calling fbar_Scan with d > 1/2"); */

   mu = N * d;                    /* mean of a binomial */
   if (m <= mu + d)
      return 1.0;
   if (mu <= 10.0)
      return ScanGlaz (N, d, m);
   prob = ScanAsympt (N, d, m);
   if ((d >= 0.3 && N >= 50.0) || (N * d * d >= 250.0 && d < 0.3)) {
      if (prob <= 0.4)
         return prob;
   }
   prob = ScanWNeff (N, d, m);
   if (prob <= 0.4)
      return prob;
   prob = ScanGlaz (N, d, m);
   if (prob > 0.4 && prob <= 1.0)
      return prob;
   return 1.0;
}


/*=========================================================================*/

double fbar_Geometric (double p, long n)
{
   util_Assert (p >= 0.0 && p <= 1.0, "fbar_Geometric:   p not in [0, 1]");
   if (n <= 0)
      return 1.0;
   if (p >= 1.0)                  /* In fact, p == 1 */
      return 0.0;
   if (p <= 0.0)                  /* In fact, p == 0 */
      return 1.0;

   return pow (1.0 - p, (double) n);
}


/*=========================================================================*/

double fbar_Poisson1 (double lam, long s)
{
   const double lamlim = 150.0;
   long i;
   double term, sum;

   util_Assert (lam >= 0.0, "fbar_Poisson1:   lambda < 0");
   if (s <= 0)
      return 1.0;

   /* If lam > lamlim, we use the Chi2 distribution according to the exact
      relation, with 2s + 2 degrees of freedom

      fdist_Poisson (lam, s) = 1 - fdist_ChiSquare (2s + 2, 2*lam)

      which also equals   1 - fdist_Gamma (s + 1, lam) */
   if (lam > lamlim)
      return fdist_Gamma ((double) s, 15, lam);

   if (s <= lam)
      return 1.0 - fdist_Poisson1 (lam, s - 1);

   /* Sum at least IMAX prob. terms from i = s to i = oo */
   sum = term = fmass_PoissonTerm1 (lam, s);
   i = s + 1;
   while (term > fmass_Epsilon || i <= s + IMAX) {
      term *= lam / i;
      sum += term;
      i++;
   }
   return sum;
}


/*=========================================================================*/

double fbar_Poisson2 (fmass_INFO W, long s)
/*
 * fbar_Poisson (lam, s) = 1 - fdist_Poisson (lam, s - 1)
 */
{
   double lam;

   util_Assert (W != NULL, "fbar_Poisson2:   fmass_INFO is NULL pointer");
   lam = W->paramR[0];

   if (s <= 0)
      return 1.0;

   /* For large lam,  we use the Chi2 distribution according to the exact
      relation, with 2s + 2 degrees of freedom

      fdist_Poisson (lam, s) = 1 - fdist_ChiSquare (2s + 2, 2*lam)
      fdist_Poisson (lam, s) = 1 - fdist_Gamma (s + 1, lam)
    */

   if (W->cdf == NULL)
      return fdist_Gamma ((double) s, 15, lam);

   if (s > W->smax)
      return fbar_Poisson1 (lam, s);

   if (s < W->smin)
      return 1.0;

   if (s > W->smed)
      /* We keep the complementary distribution in the upper part of cdf */
      return W->cdf[s - W->smin];
   else
      return 1.0 - W->cdf[s - 1 - W->smin];
}


/*=========================================================================*/

double fbar_Binomial2 (fmass_INFO W, long s)
{
   double p;
   long n;

   util_Assert (W != NULL, "fbar_Binomial2:   fmass_INFO is NULL pointer");
   n = W->paramI[0];
   p = W->paramR[0];
   util_Assert (p >= 0.0 && p <= 1.0, "fbar_Binomial2:   p not in [0, 1]");

   if (0 == n)
      return 1.0;
   if (s < 1)
      return 1.0;
   if (s > n)
      return 0.0;
   if (p == 0.0)
      return 0.0;
   if (p == 1.0)
      return 1.0;

   if (W->cdf != NULL) {
      if (s >= W->smax) {
         /* Add IMAX dominant terms to get a few decimals in the tail */
         const double q = 1.0 - p;
         double z, sum, term;
         long i;
         sum = term = fmass_BinomialTerm3 (n, p, s);
         if (fabs (q) > 0.0) {
            z = p / q;
         } else {
            z = 0.0;
            util_Warning (1, "fbar_Binomial2:   p / q = infinite");
         }
         i = s;
         while (i < n && i < s + IMAX) {
            term = term * z * (n - i) / (i + 1);
            sum += term;
            i++;
         }
         return sum;
         /* return fdist_Beta (s, n - s + 1, 10, p); */
      }

      if (s <= W->smin)
         return 1.0;

      if (s > W->smed)
         /* We keep the complementary distribution in the upper part of cdf */
         return W->cdf[s - W->smin];
      else
         return 1.0 - W->cdf[s - 1 - W->smin];

   } else {
      return 1.0 - fdist_Binomial1 (n, p, s - 1);
   }
}


/*=========================================================================*/

double fbar_NegaBin2 (fmass_INFO W, long s)
{
   double p;
   long n;

   util_Assert (W != NULL, "fbar_NegaBin2:   fmass_INFO is NULL pointer");
   n = W->paramI[0];
   p = W->paramR[0];
   util_Assert (p >= 0.0 && p <= 1.0, "fbar_NegaBin2:   p not in [0, 1]");

   if (s < 1)
      return 1.0;
   if (p >= 1.0)                  /* In fact, p == 1 */
      return 0.0;
   if (p <= 0.0)                  /* In fact, p == 0 */
      return 1.0;

   if (W->cdf == NULL)
      return fdist_Binomial1 (s - 1 + n, p, n - 1);

   if (s >= W->smax)
      return fdist_Binomial1 (s - 1 + n, p, n - 1);
   if (s <= W->smin)
      return 1.0;
   if (s > W->smed)
      /* We keep the complementary distribution in the upper part of cdf */
      return W->cdf[s - W->smin];
   else
      return 1.0 - W->cdf[s - 1 - W->smin];

}


/*=========================================================================*/
