/*************************************************************************\
 *
 * Package:        ProbDist
 * File:           fdist.c
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

#include "fdist.h"
#include "fmass.h"
#include "fbar.h"

#include "num.h"
#include "num2.h"
#include "util.h"
#include "tables.h"

#include <stddef.h>
#include <limits.h>
#include <float.h>
#include <math.h>


double fdist_belog (double);
void fdist_CalcB4 (double, double *, double *, double *, double *);
static double Pomeranz (long n, double x);


const double fdist_XINF  = DBL_MAX; /* x infinity for some distributions */
const double fdist_XBIG  = 100.0;   /* x infinity for some distributions */
const double fdist_XBIGM = 1000.0;  /* x infinity for some distributions */

/* EpsArray[j]: Epsilon required for j decimal degits of precision */
static const double EpsArray[] = {
   0.5, 0.5E-1, 0.5E-2, 0.5E-3, 0.5E-4, 0.5E-5, 0.5E-6, 0.5E-7, 0.5E-8,
   0.5E-9, 0.5E-10, 0.5E-11, 0.5E-12, 0.5E-13, 0.5E-14, 0.5E-15, 0.5E-16,
   0.5E-17, 0.5E-18, 0.5E-19, 0.5E-20, 0.5E-21, 0.5E-22, 0.5E-23, 0.5E-24,
   0.5E-25, 0.5E-26, 0.5E-27, 0.5E-28, 0.5E-29, 0.5E-30, 0.5E-31, 0.5E-32,
   0.5E-33, 0.5E-34, 0.5E-35
};

/* static const double EpsilonLR = 1.0E-15;  */
#define X_EPSILON 1.0e-3          /* For x --> 0 */

#define TWOPI 6.28318530717958647688



/*-------------------------------------------------------------------------*/

double fdist_belog (double x)
/*
 * This is the function   (1 - x*x + 2*x*log(x)) / ((1 - x)*(1 - x))
 */
{
   if (x > 1.0)
      return -fdist_belog(1.0/x);
   if (x < 1.0e-20)
      return 1.0;
   if (x < 0.9)
      return (1.0 - x * x + 2.0 * x * log (x)) / ((1.0 - x) * (1.0 - x));
   if (x == 1.0)
      return 0.0;
   {
      /* For x near 1, use a series expansion to avoid loss of precision. */
      double term;
      const double EPS = 1.0e-12;
      const double Y = 1.0 - x;
      double ypow = 1.0;
      double sum = 0.0;
      int j = 2;
      do {
         ypow *= Y;
         term = ypow / (j * (j + 1));
         sum += term;
         j++;
      } while (fabs (term / sum) > EPS);

      return 2.0 * sum;
   }
}


/*=========================================================================*/

double fdist_Unif (double x)
{
   if (x <= 0.0)
      return 0.0;
   if (x >= 1.0)
      return 1.0;
   return x;
}


/*=========================================================================*/

double fdist_Expon (double x)
{
   if (x <= 0.0)
      return 0.0;
   if (x >= fdist_XBIG)
      return 1.0;
   if (x > X_EPSILON)
      return 1.0 - exp (-x);

   /* Avoid loss of precision for small x */
   return x * (1.0 - x * (0.5 - x*(1.0 / 6.0 - x/24.0)));
}


/*=========================================================================*/

double fdist_Weibull (double c, double x)
{
   double y;
   util_Assert (c > 0.0, "fdist_Weibull:   c <= 0");
   if (x <= 0.0)
      return 0.0;
   if (x >= fdist_XBIG && c >= 1.0)
      return 1.0;

   y = c*log(x);
   if (y >= 5.0)
      return 1.0;
   y = exp(y);

   if (y > X_EPSILON)
      return (1.0 - exp (-y));

   /* Avoid loss of precision for small y */
   return y * (1.0 - y * (0.5 - y*(1.0/6.0 - y/24.0)));
}


/*=========================================================================*/

double fdist_ExtremeValue (double x)
{
   if (x <= -10.0)
      return 0.0;
   if (x >= fdist_XBIG)
      return 1.0;
   return exp (-exp (-x));
}


/*=========================================================================*/

double fdist_Logistic (double x)
{
   if (x <= -fdist_XBIG)
      return exp (x);
   if (x >= fdist_XBIG)
      return 1.0;
   return 1.0 / (1.0 + exp (-x));
}


/*=========================================================================*/

double fdist_Pareto (double c, double x)
{
   double y = c*log(x);
   util_Assert (c > 0.0, "fdist_Pareto:   c <= 0");
   if (x <= 1.0)
      return 0.0;
   if (y >= 50.0)
      return 1.0;
   y = exp(y);
   return (1.0 - 1.0 / y);
}


/*=========================================================================*/

double fdist_Normal1 (double x)
/*
 * Returns P[X <= x] for the normal distribution.
 * As in p:90 of W.J.Kennedy Jr and J.E.Gentle. Statistical computing.
 * Dekker, New York, 1980.
 */
{
   static const double Racinedeux = 1.4142135623730950488;
   static const double racineunsurpi = 0.56418958354775628694;

   static const double p10 = 2.4266795523053175e2;
   static const double p11 = 2.1979261618294152e1;
   static const double p12 = 6.9963834886191355;
   static const double p13 = -3.5609843701815385e-2;

   static const double p20 = 3.004592610201616005e2;
   static const double p21 = 4.519189537118729422e2;
   static const double p22 = 3.393208167343436870e2;
   static const double p23 = 1.529892850469404039e2;
   static const double p24 = 4.316222722205673530e1;
   static const double p25 = 7.211758250883093659e0;
   static const double p26 = 5.641955174789739711e-1;
   static const double p27 = -1.368648573827167067e-7;

   static const double p30 = -2.99610707703542174e-3;
   static const double p31 = -4.94730910623250734e-2;
   static const double p32 = -2.26956593539686930e-1;
   static const double p33 = -2.78661308609647788e-1;
   static const double p34 = -2.23192459734184686e-2;

   static const double q10 = 2.1505887586986120e2;
   static const double q11 = 9.1164905404514901e1;
   static const double q12 = 1.5082797630407787e1;
   static const double q13 = 1.0;

   static const double q20 = 3.004592609569832933e2;
   static const double q21 = 7.909509253278980272e2;
   static const double q22 = 9.313540948506096211e2;
   static const double q23 = 6.389802644656311665e2;
   static const double q24 = 2.775854447439876434e2;
   static const double q25 = 7.700015293522947295e1;
   static const double q26 = 1.278272731962942351e1;
   static const double q27 = 1.0;

   static const double q30 = 1.06209230528467918e-2;
   static const double q31 = 1.91308926107829841e-1;
   static const double q32 = 1.05167510706793207e0;
   static const double q33 = 1.98733201817135256e0;
   static const double q34 = 1.0;

   static const double xasymp = 40.0;
   double Ycarre, unsurY2, Y, R, erf;

   if (x < -xasymp)
      return 0.0;
   if (x > xasymp)
      return 1.0;

   if (x < 0.0)
      return 1.0 - fdist_Normal1 (-x);

   Y = x / Racinedeux;
   Ycarre = x * x / 2.0;
   if (Y < 0.447) {
      R = (p10 + Ycarre * (p11 + Ycarre * (p12 + Ycarre * p13))) /
         (q10 + Ycarre * (q11 + Ycarre * (q12 + Ycarre * q13)));
      erf = Y * R;
   } else {
      if (Y <= 4.0) {
         R = (p20 + Y * (p21 + Y * (p22 + Y * (p23 + Y * (p24 + Y * (p25 +
                           Y * (p26 + Y * p27))))))) / (q20 + Y * (q21 +
               Y * (q22 + Y * (q23 + Y * (q24 + Y * (q25 + Y * (q26 +
                              Y * q27)))))));
         if (-Ycarre < DBL_MIN_EXP * num_Ln2)
            erf = 1.0;
         else
            erf = 1.0 - exp (-Ycarre) * R;
      } else {
         double temp;
         unsurY2 = 1.0 / Ycarre;
         R = (p30 + unsurY2 * (p31 + unsurY2 * (p32 + unsurY2 *
                  (p33 + unsurY2 * p34)))) / (q30 + unsurY2 *
            (q31 + unsurY2 * (q32 + unsurY2 * (q33 + unsurY2 * q34))));
         if (-Ycarre < DBL_MIN_EXP * num_Ln2)
            temp = 0.0;
         else
            temp = exp (-Ycarre);
         erf = 1.0 - (temp / Y) * (racineunsurpi + R / Ycarre);
      }
   }
   return ((1.0 + erf) / 2.0);
}


/**************************************************************************/
/* 
 * The precision of double is 16 decimals; we shall thus use COEFFMAX = 24
 * coefficients. But the approximation is good to 30 decimals of precision
 * with 44 coefficients.
 */
#define COEFFMAX 24

static const double Normal2_A[44] = {
   6.10143081923200417926465815756e-1,
   -4.34841272712577471828182820888e-1,
   1.76351193643605501125840298123e-1,
   -6.0710795609249414860051215825e-2,
   1.7712068995694114486147141191e-2,
   -4.321119385567293818599864968e-3,
   8.54216676887098678819832055e-4,
   -1.27155090609162742628893940e-4,
   1.1248167243671189468847072e-5,
   3.13063885421820972630152e-7,
   -2.70988068537762022009086e-7,
   3.0737622701407688440959e-8,
   2.515620384817622937314e-9,
   -1.028929921320319127590e-9,
   2.9944052119949939363e-11,
   2.6051789687266936290e-11,
   -2.634839924171969386e-12,
   -6.43404509890636443e-13,
   1.12457401801663447e-13,
   1.7281533389986098e-14,
   -4.264101694942375e-15,
   -5.45371977880191e-16,
   1.58697607761671e-16,
   2.0899837844334e-17,
   -5.900526869409e-18,
   -9.41893387554e-19,
   2.14977356470e-19,
   4.6660985008e-20,
   -7.243011862e-21,
   -2.387966824e-21,
   1.91177535e-22,
   1.20482568e-22,
   -6.72377e-25,
   -5.747997e-24,
   -4.28493e-25,
   2.44856e-25,
   4.3793e-26,
   -8.151e-27,
   -3.089e-27,
   9.3e-29,
   1.74e-28,
   1.6e-29,
   -8.0e-30,
   -2.0e-30
};


double fdist_Normal2 (double x)
/*
 * Returns P[X < x] for the normal distribution.
 * As in J. L. Schonfelder, Math. of Computation, Vol. 32,
 * pp 1232--1240, (1978).
 */
{
   double t, r;
   if (x <= -fdist_XBIG)
      return 0.0;
   if (x >= fdist_XBIG)
      return 1.0;

   x = -x / num_Rac2;
   if (x < 0) {
      x = -x;
      t = (x - 3.75) / (x + 3.75);
      r = 1.0 - 0.5 * exp (-x * x) * num2_EvalCheby (Normal2_A, COEFFMAX, t);
   } else {
      t = (x - 3.75) / (x + 3.75);
      r = 0.5 * exp (-x * x) * num2_EvalCheby (Normal2_A, COEFFMAX, t);
   }
   return (r);
}


/*=========================================================================*/
#ifdef HAVE_ERF

double fdist_Normal3 (double x)
{
   return (erfc (-x * num_1Rac2)) / 2.0;
}

#endif


/*=========================================================================*/

double fdist_Normal4 (double x)
{
   static const double V[121] = {
      1.2533141373155, 1.137490921203605, 1.037824575853727,
      0.951527192071207, 0.8763644564536924, 0.8105337152790306,
      0.7525711790634081, 0.7012808218544303, 0.6556795424187987,
      0.61495459615093, 0.5784303460476312, 0.5455421356582171,
      0.5158156382179634, 0.4888504415275737, 0.4643069280394423,
      0.4418957328326002, 0.4213692292880546, 0.4025146181296722,
      0.3851482907984348, 0.3691112106902635, 0.3542651113297938,
      0.3404893532870847, 0.3276783146905521, 0.31573921586941,
      0.3045902987101033, 0.2941592970402893, 0.284382146748493,
      0.2752018941576065, 0.2665677689682238, 0.2584343943120386,
      0.2507611114439651, 0.243511400615456, 0.2366523829135607,
      0.230154390478801, 0.2239905946538289, 0.2181366833614714,
      0.2125705804420318, 0.2072722008565011, 0.2022232366330547,
      0.1974069692375194, 0.1928081047153158, 0.1884126285076003,
      0.1842076773079702, 0.1801814257143918, 0.1763229857571027,
      0.1726223176578506, 0.1690701504076941, 0.1656579109468773,
      0.1623776608968675, 0.1592220399363674, 0.1561842150339759,
      0.153257834853479, 0.1504369887362691, 0.1477161697413935,
      0.145090241289131, 0.1425544070104023, 0.1401041834530503,
      0.1377353753382303, 0.1354440530967635, 0.1332265324471292,
      0.1310793558044918, 0.1289992753343376, 0.126983237485437,
      0.1250283688553504, 0.1231319632579323, 0.1212914698765462,
      0.119504482399253, 0.1177687290432979, 0.1160820633859823,
      0.1144424559276431, 0.112847986320103, 0.1112968362007359,
      0.1097872825783083, 0.1083176917221132, 0.1068865135106745,
      0.1054922762005562, 0.1041335815795983, 0.1028091004723001,
      0.1015175685681028, 0.1002577825460485, 0.09902859647173194,
      0.09782891844465691, 0.09665770747608191, 0.09551397057921558,
      0.09439676005522439, 0.09330517095996169, 0.09223833873763035,
      0.09119543700877471, 0.09017567550106469, 0.08917829811230435,
      0.08820258109597616, 0.08724783136042988, 0.08631338487354936,
      0.08539860516539227, 0.08450288192189578, 0.08362562966329139,
      0.08276628650136918, 0.08192431297018954, 0.08109919092525536,
      0.08029042250654048, 0.07949752916111721, 0.07872005072144664,
      0.07795754453568722, 0.07720958464664668, 0.07647576101624852,
      0.07575567879261112, 0.07504895761704659, 0.07435523096847724,
      0.07367414554294564, 0.07300536066605566, 0.07234854773633338,
      0.07170338969763433, 0.07106958053885212, 0.07044682481930167,
      0.06983483721825942, 0.06923334210724434, 0.06864207314371742,
      0.06806077288496332, 0.0674891924209997, 0.06692709102543307,
      0.06637423582325017
   };

   int j;
   lebool negatif;
   double t, u, z, h;
   double r, r1, r2, r3, r4, r5, r6, r7, r8;

   if (x <= -fdist_XBIG)
      return 0.0;
   if (x >= fdist_XBIG)
      return 1.0;
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
      return u;
   else
      return 1.0 - u;
}


/*=========================================================================*/

static double InitBiNormal (double x, double y, double rho)
{
   /* The special cases of the BiNormal */
   if (fabs (rho) > 1.0) {
      util_Error ("fdist_BiNormal:   |rho| > 1");
      return -1.0;
   }
   if (x == 0.0 && y == 0.0)
      return 0.25 + asin(rho)/TWOPI;
   if (rho == 1.0) {
      x = util_Min(x,y);
      return fdist_Normal2 (x);
   }
   if (rho == 0.0) {
      return fdist_Normal2 (x) * fdist_Normal2 (y);
   }
   if (rho == -1.0) {
      if (y <= -x)
         return 0.0;
      else
         return fdist_Normal2 (x) - fdist_Normal2 (-y);
   }
   if ((x <= -fdist_XBIG) || (y <= -fdist_XBIG))
      return 0.0;
   if (x >= fdist_XBIG)
      return fdist_Normal2 (y);
   if (y >= fdist_XBIG)
      return fdist_Normal2 (x);

   return -2.0;
}


/*=========================================================================*/

double fdist_BiNormal1 (double x, double y, double rho, int ndig)
{
   double a2, ap, b, cn, conex, ex, g2, gh, gk, gw, h2, h4, rr, s1, s2,
      sgn, sn, sp, sqr, t, w2, wh, wk;
   int is = -1;
   int flag = 1;
   const double ah = -x;
   const double ak = -y;
   const double con = num_Pi * num_TENNEGPOW[ndig];
   const double EPSILON = 0.5 * num_TENNEGPOW[ndig];

   util_Assert (ndig <= 15, "fdist_BiNormal1:   ndig > 15");

   b = InitBiNormal (x, y, rho);
   if (b >= 0.0)
      return b;

   gh = fdist_Normal2 (-ah) / 2.0;
   gk = fdist_Normal2 (-ak) / 2.0;

   b = 0;
   rr = (1 - rho) * (1 + rho);
   sqr = sqrt (rr);
   flag = 1;
   if (ah != 0) {
      b = gh;
      if (ah * ak < 0)
         b = b - .5;
      else if (ah * ak == 0) {
         flag = 0;
      }
   } else if (ak == 0) {
      return asin (rho) / TWOPI + .25;
   }
   if (flag)
      b += gk;
   if (ah != 0) {
      flag = 0;
      wh = -ah;
      wk = (ak / ah - rho) / sqr;
      gw = 2 * gh;
      is = -1;
   }

   do {
      if (flag) {
         wh = -ak;
         wk = (ah / ak - rho) / sqr;
         gw = 2 * gk;
         is = 1;
      }
      flag = 1;
      sgn = -1;
      t = 0;
      if (wk != 0) {
         if (fabs (wk) >= 1) {
            if (fabs (wk) == 1) {
               t = wk * gw * (1 - gw) / 2;
               b = b + sgn * t;
               if (is >= 0)
                  break;
               else
                  continue;
            } else {
               sgn = -sgn;
               wh = wh * wk;
               g2 = fdist_Normal2 (wh);
               wk = 1 / wk;
               if (wk < 0)
                  b = b + .5;
               b = b - (gw + g2) / 2 + gw * g2;
            }
         }
         h2 = wh * wh;
         a2 = wk * wk;
         h4 = h2 * .5;
         ex = 0;
         if (h4 < 150.0)
            ex = exp (-h4);
         w2 = h4 * ex;
         ap = 1;
         s2 = ap - ex;
         sp = ap;
         s1 = 0;
         sn = s1;
         conex = fabs (con / wk);
         do {
            cn = ap * s2 / (sn + sp);
            s1 = s1 + cn;
            if (fabs (cn) <= conex)
               break;
            sn = sp;
            sp = sp + 1;
            s2 = s2 - w2;
            w2 = w2 * h4 / sp;
            ap = -ap * a2;
         } while (1);
         t = (atan (wk) - wk * s1) / TWOPI;
         b = b + sgn * t;
      }
      if (is >= 0)
         break;
   } while (ak != 0);

   if (b < EPSILON)
      b = 0;
   if (b > 1)
      b = 1;
   return b;
}


/*=========================================================================*/

double fdist_BiNormal2 (double dh, double dk, double rho)
{
   const double twopi = 2.0 * num_Pi;
   double W[11][3];
   double X[11][3];
   double h, k, hk, bvn, hs, asr, sn, as, a, b, c, d, sp, rs, ep, bs, xs;
   int i, lg, ng, is;

   bvn = InitBiNormal (dh, dk, rho);
   if (bvn >= 0.0)
      return bvn;
/*
   I have made small changes in Genz's Matlab function to make it compatible
   with module fdist. (R. Simard)
*/

/*
//   Copyright (C) 2005, Alan Genz,  All rights reserved.               
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided the following conditions are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The contributor name(s) may not be used to endorse or promote 
//        products derived from this software without specific prior written 
//        permission.
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
//   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
//   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
//   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
//   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
//   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
//   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
//   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
//   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//   function p = bvnl( dh, dk, r )
//
//  A function for computing bivariate normal probabilities.
//  bvnl calculates the probability that x < dh and y < dk. 
//    parameters  
//      dh 1st upper integration limit
//      dk 2nd upper integration limit
//      r   correlation coefficient
//
//   Author
//       Alan Genz
//       Department of Mathematics
//       Washington State University
//       Pullman, Wa 99164-3113
//       Email : alangenz@wsu.edu
//   This function is based on the method described by 
//        Drezner, Z and G.O. Wesolowsky, (1989),
//        On the computation of the bivariate normal inegral,
//        Journal of Statist. Comput. Simul. 35, pp. 101-107,
//    with major modifications for double precision, for |r| close to 1,
//    and for matlab by Alan Genz - last modifications 7/98.
//
//      p = bvnu( -dh, -dk, r );
//      return
//
//   end bvnl
//
//      function p = bvnu( dh, dk, r )
//
//  A function for computing bivariate normal probabilities.
//  bvnu calculates the probability that x > dh and y > dk. 
//    parameters  
//      dh 1st lower integration limit
//      dk 2nd lower integration limit
//      r   correlation coefficient
//
//   Author
//       Alan Genz
//       Department of Mathematics
//       Washington State University
//       Pullman, Wa 99164-3113
//       Email : alangenz@wsu.edu
//
//    This function is based on the method described by 
//        Drezner, Z and G.O. Wesolowsky, (1989),
//        On the computation of the bivariate normal inegral,
//        Journal of Statist. Comput. Simul. 35, pp. 101-107,
//    with major modifications for double precision, for |r| close to 1,
//    and for matlab by Alan Genz - last modifications 7/98.
//        Note: to compute the probability that x < dh and y < dk, use 
//              bvnu( -dh, -dk, r ). 
//
*/
   if (fabs (rho) < 0.3) {
      ng = 0;
      lg = 3;
/*       Gauss Legendre points and weights, n =  6 */
      W[1][0] = 0.1713244923791705;
      W[2][0] = 0.3607615730481384;
      W[3][0] = 0.4679139345726904;

      X[1][0] = 0.9324695142031522;
      X[2][0] = 0.6612093864662647;
      X[3][0] = 0.2386191860831970;

   } else if (fabs (rho) < 0.75) {
      ng = 1;
      lg = 6;
/*       Gauss Legendre points and weights, n = 12 */
      W[1][1] = 0.4717533638651177e-1;
      W[2][1] = 0.1069393259953183;
      W[3][1] = 0.1600783285433464;
      W[4][1] = 0.2031674267230659;
      W[5][1] = 0.2334925365383547;
      W[6][1] = 0.2491470458134029;

      X[1][1] = 0.9815606342467191;
      X[2][1] = 0.9041172563704750;
      X[3][1] = 0.7699026741943050;
      X[4][1] = 0.5873179542866171;
      X[5][1] = 0.3678314989981802;
      X[6][1] = 0.1252334085114692;

   } else {
      ng = 2;
      lg = 10;
/*       Gauss Legendre points and weights, n = 20 */
      W[1][2] = 0.1761400713915212e-1;
      W[2][2] = 0.4060142980038694e-1;
      W[3][2] = 0.6267204833410906e-1;
      W[4][2] = 0.8327674157670475e-1;
      W[5][2] = 0.1019301198172404;
      W[6][2] = 0.1181945319615184;
      W[7][2] = 0.1316886384491766;
      W[8][2] = 0.1420961093183821;
      W[9][2] = 0.1491729864726037;
      W[10][2] = 0.1527533871307259;

      X[1][2] = 0.9931285991850949;
      X[2][2] = 0.9639719272779138;
      X[3][2] = 0.9122344282513259;
      X[4][2] = 0.8391169718222188;
      X[5][2] = 0.7463319064601508;
      X[6][2] = 0.6360536807265150;
      X[7][2] = 0.5108670019508271;
      X[8][2] = 0.3737060887154196;
      X[9][2] = 0.2277858511416451;
      X[10][2] = 0.7652652113349733e-1;
   }

   h = -dh;
   k = -dk;
   hk = h * k;
   bvn = 0;
   if (fabs (rho) < 0.925) {
      hs = (h * h + k * k) / 2.0;
      asr = asin (rho);
      for (i = 1; i <= lg; ++i) {
         sn = sin (asr * (1.0 - X[i][ng]) / 2.0);
         bvn += W[i][ng] * exp ((sn * hk - hs) / (1.0 - sn * sn));
         sn = sin (asr * (1.0 + X[i][ng]) / 2.0);
         bvn += W[i][ng] * exp ((sn * hk - hs) / (1.0 - sn * sn));
      }
      bvn =
         bvn * asr / (4.0 * num_Pi) + fdist_Normal2 (-h) * fdist_Normal2 (-k);

   } else {
      if (rho < 0.0) {
         k = -k;
         hk = -hk;
      }
      if (fabs (rho) < 1.0) {
         as = (1.0 - rho) * (1.0 + rho);
         a = sqrt (as);
         bs = (h - k) * (h - k);
         c = (4.0 - hk) / 8.0;
         d = (12.0 - hk) / 16.0;
         asr = -(bs / as + hk) / 2.0;
         if (asr > -100.0)
            bvn =
               a * exp (asr) * (1.0 - c * (bs - as) * (1.0 -
                  d * bs / 5.0) / 3.0 + c * d * as * as / 5.0);

         if (-hk < 100.0) {
            b = sqrt (bs);
            sp = sqrt (twopi) * fdist_Normal2 (-b / a);
            bvn = bvn - exp (-hk / 2.0) * sp * b * (1.0 - c * bs * (1.0 -
                  d * bs / 5.0) / 3.0);
         }
         a = a / 2.0;
         for (i = 1; i <= lg; ++i) {
            for (is = -1; is <= 1; is += 2) {
               xs = (a * (is * X[i][ng] + 1.0));
               xs = xs * xs;
               rs = sqrt (1.0 - xs);
               asr = -(bs / xs + hk) / 2.0;
               if (asr > -100.0) {
                  sp = (1.0 + c * xs * (1.0 + d * xs));
                  ep = exp (-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs;
                  bvn += a * W[i][ng] * exp (asr) * (ep - sp);
               }
            }
         }
         bvn = -bvn / twopi;
      }
      if (rho > 0.0) {
         if (k > h)
            h = k;
         bvn += fdist_Normal2 (-h);
      }
      if (rho < 0.0) {
         xs = fdist_Normal2 (-h) - fdist_Normal2 (-k);
         if (xs < 0.0)
            xs = 0.0;
         bvn = -bvn + xs;
      }
   }
   if (bvn <= 0.0)
      return 0.0;
   if (bvn >= 1.0)
      return 1.0;
   return bvn;
}


/*=========================================================================*/

double fdist_LogNormal (double mu, double sigma, double x)
{
   util_Assert (sigma > 0.0, "fdist_LogNormal:  sigma  <= 0");
   if (x <= 0.0)
      return 0.0;
   return fdist_Normal2 ((log (x) - mu) / sigma);
}


/*=========================================================================*/

double fdist_JohnsonSB (double alpha, double beta, double a, double b,
   double x)
{
   util_Assert (beta > 0.0, "fdist_JohnsonSB:  beta  <= 0");
   util_Assert (b > a, "fdist_JohnsonSB:  b  <= a");
   if (x <= a)
      return 0.0;
   if (x >= b)
      return 1.0;
   return fdist_Normal2 (alpha + beta * log ((x - a) / (b - x)));
}


/*=========================================================================*/

double fdist_JohnsonSU (double alpha, double beta, double x)
{
   const double XLIM = 1.0e10;
   double r;
   lebool negative = FALSE;
   util_Assert (beta > 0.0, "fdist_JohnsonSU:  beta  <= 0");
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
      return fdist_Normal2 (alpha + beta * log (r));
   else
      return 0.0;
}


/**************************************************************************/

double fdist_ChiSquare1 (long N, double x)
/*
 * Returns an approximation of the Chi square cdf (N degrees of freedom)
 * As in p:116 of W.J.Kennedy Jr and J.E.Gentle. Statistical computing,
 * Dekker, New York, 1980.
 */
{
   const double tiers = 0.33333333333333333;
   const double pt2 = 0.22222222222222222;
   const double moinsdixhuit = -18.8055;
   const double gam = 0.8862269254527579825931;
   double H, H2, E, DemiX, Terme, Sommation, Y = 0;
   long i;

   util_Assert (N > 0, "fdist_ChiSquare1:   k < 1");
   if (x <= 0.0)
      return 0.0;
   if (x >= fdist_XBIG * N)
      return 1.0;

   if (N > 1000) {
      if (x < 2.0)
         return 0.0;
      x = (pow ((x / N), tiers) - (1.0 - pt2 / N)) / sqrt (pt2 / N);
      if (x > 5.0)
         return 1.0;
      if (x < moinsdixhuit)
         return 0.0;
      return fdist_Normal2 (x);

   } else {
      DemiX = x / 2.0;
      if (!(N & 1)) {             /* even N */
         if (-DemiX < DBL_MIN_EXP * num_Ln2)
            Terme = 0.0;
         else
            Terme = exp (-DemiX);
         Sommation = Terme;
         for (i = 1; i < N / 2; i++) {
            Terme = Terme * DemiX / ((double) i);
            Sommation += Terme;
         }
         Y = 1.0 - Sommation;
      } else {
         H2 = -1.0 + 2.0 * fdist_Normal2 (sqrt (x));
         if (N == 1)
            return H2;
         if (-DemiX < DBL_MIN_EXP * num_Ln2)
            E = 0.0;
         else
            E = exp (-DemiX);
         Terme = sqrt (DemiX) * E / gam;
         H = H2;
         for (i = 3; i < N; i += 2) {
            H -= Terme;
            Terme = Terme * DemiX * 2.0 / ((double) i);
         }
         Y = H - Terme;
      }
   }
   if (Y < 0.0)
      return 0.0;
   else
      return Y;
}


/*=========================================================================*/

double fdist_ChiSquare2 (long n, int d, double x)
{
   util_Assert (n > 0, "fdist_ChiSquare2:   n <= 0");
   if (x <= 0.0)
      return 0.0;
   if (x >= fdist_XBIG * n)
      return 1.0;
   return fdist_Gamma (n / 2.0, d, x / 2.0);
}


/*=========================================================================*/

#define Student_n1 20
#define Student_x1 8.01
#define Student_kmax 200
#define Student_eps 0.5E-16

double fdist_Student1 (long n, double x)
{
   double a, u, b, y, z, z2, prec;
   long k;

   util_Assert (n > 0, "fdist_Student1:   n <= 0");
   if (n == 1) {
      if (x < -0.5)
         return atan(-1.0/x) / num_Pi;
      return 0.5 + (atan (x)) / num_Pi;
   }

   if (n == 2) {
      z = 1.0 + x * x / 2.0;
      if (x >= 0.)
         return 0.5 + x / (2.0 * sqrt (z) * num_Rac2);
      else
         return 0.25 / (z * (0.5 - x /(2.0*sqrt(z)*num_Rac2)));
   }

   /* first case: small n and small x */
   if (n <= Student_n1 && x <= Student_x1) {
      b = 1.0 + x * x / n;
      y = x / sqrt ((double) n);
      z = 1.0;
      for (k = n - 2; k >= 2; k -= 2) {
         z = 1.0 + z * (k - 1.0) / (k * b);
      }
      if (n % 2 == 0) {
         u = (1.0 + z * y / sqrt (b)) / 2.0;
         if (u >= 0.)
            return u;
         else
            return 0.;
      } else {
         if (y > -1.0)
            return (0.5 + (atan (y) + z * y / b) / num_Pi);
         else {
            u = (atan (-1.0 / y) + z * y / b) / num_Pi;
            if (u >= 0.)
               return u;
            else
               return 0.;
         }
      }

   /* second case: large n and small x */
   } else if (x < Student_x1) {
      a = n - 0.5;
      b = 48.0 * a * a;
      z2 = a * num2_log1p (x * x / n);
      z = sqrt (z2);
      y = (((((64.0 * z2 + 788.0) * z2 + 9801.0) * z2 + 89775.0) * z2 +
            543375.0) * z2 + 1788885.0) * z / (210.0 * b * b * b);
      y -=
         (((4.0 * z2 + 33.0) * z2 + 240.0) * z2 + 855.0) * z / (10.0 * b * b);
      y += z + (z2 + 3.0) * z / b;
      if (x >= 0.0)
         return fbar_Normal1 (-y);
      else
         return fbar_Normal1 (y);

   /* third case: large x */
   } else {
      /* Compute the Student probability density */
      b = 1.0 + x * x / n;
      /* to avoid overflow with the 2 Gamma functions, use their logarithm.
         However, for large n, there will be some loss of precision */
      y = num2_LnGamma ((n + 1) / 2.0) - num2_LnGamma (n / 2.0);
      y = exp (y);
      y *= pow (b, -(n + 1) / 2.0) / sqrt (num_Pi * n);

      y *= 2.0 * sqrt (n * b);
      z = y / n;
      k = 2;
      z2 = prec = 10.0;
      while (k < Student_kmax && prec > Student_eps) {
         y *= (k - 1) / (k * b);
         z += y / (n + k);
         prec = fabs (z - z2);
         z2 = z;
         k += 2;
      }
      util_Warning (k >= Student_kmax, "fdist_Student1: k >= Student_kmax");
      if (x >= 0.0)
         return 1.0 - z / 2.0;
      else
         return z / 2.0;
   }
}

/*=========================================================================*/

double fdist_Student2 (long n, int d, double x)
{
   util_Assert (n > 0, "fdist_Student2:   n <= 0");
   util_Assert (d > 0, "fdist_Student2:   d <= 0");
   util_Assert (d <= 15, "fdist_Student2:   d > 15");
   if (x <= -fdist_XBIG)
      return 0.0;
   if (x >= fdist_XBIG)
      return 1.0;

   if (x >= 0.0)
      return 0.5 * (1.0 + fdist_Beta (0.5, 0.5 * n, d, x * x / (n + x * x)));
   else
      return 0.5 * (fdist_Beta (0.5 * n, 0.5, d, n / (n + x * x)));
}

/*=========================================================================*/

double fdist_Gamma (double alpha, int d, double x)
{
   const double ALIM = 1.0e5;
   const double EPS = EpsArray[d];

   util_Assert (alpha > 0.0, "fdist_Gamma:   a <= 0");
   util_Assert (d > 0, "fdist_Gamma:   d <= 0");
   util_Assert (d < 16, "fdist_Gamma:   d > 15");
   if (x <= 0.0)
      return 0.0;
   if (1.0 == alpha)
      return fdist_Expon (x);

   if (alpha >= ALIM) {
      double d2 = x + 1.0/3.0 - alpha - 0.02/alpha;
      double S = alpha - 1.0/2.0;
      double z = d2 * sqrt((1 + fdist_belog(S/x))/x);
      return fdist_Normal2 (z);
   }

   if (x <= 1.0 || x < alpha) {
      double v, z, an, term;
      v = exp (alpha * log (x) - x - num2_LnGamma (alpha));
      z = 1.0;
      term = 1.0;
      an = alpha;
      do {
         an += 1.0;
         term *= x / an;
         z += term;
      } while (term >= EPS * z);
      return z * v / alpha;

   } else
      return 1.0 - fbar_Gamma (alpha, d, x);
}


/*=========================================================================*/

static double Isubx_pq_small (double p, double q, double x, int d)
/* 
 * Evaluates fdist_Beta (p, q, d, x) when 0 < p <= 1 and 0 < q <= 2 to a
 * precision of d = -log10 (2 epsilon) decimal digits. Uses a series
 * expansion in powers of x.
 */
{

   int k = 0;
   double s, u, v;
   double epsilon;
   util_Assert (p > 0.0 && p <= 1.0, "Isubx_pq_small:   p not in (0, 1] ");
   util_Assert (q > 0.0 && q <= 2.0, "Isubx_pq_small:   q not in (0, 2] ");

   epsilon = EpsArray[d];
   u = pow (x, p);
   s = u / p;
   do {
      u = (k + 1 - q) * x * u / (k + 1);
      v = u / (k + 1 + p);
      s += v;
      k++;
   } while ((fabs (v) / s) > epsilon);

   v = num2_LnGamma (p + q) - num2_LnGamma (p) - num2_LnGamma (q);
   return s * exp (v);
}

/*-------------------------------------------------------------------------*/

static void forward (double p, double q, double x, double I0, double I1,
   int nmax, double I[])
/* 
 * Given I0 = fdist_Beta (p, q, x) and I1 = fdist_Beta (p, q + 1, x),
 * generates fdist_Beta (p, q + n, x) for n = 0, 1, 2, ..., nmax, and
 * stores the result in I.
 */
{

   int n;

   I[0] = I0;
   if (nmax > 0)
      I[1] = I1;
   for (n = 1; n < nmax; n++)
      I[n + 1] = (1 + (n - 1 + p + q) * (1. - x) / (n + q)) * I[n]
         - (n - 1 + p + q) * (1. - x) * I[n - 1] / (n + q);
}

/*-------------------------------------------------------------------------*/

static void backward (double p, double q, double x, double I0, int d,
   int nmax, double I[])
/*
 * Given I0 = fdist_Beta (p, q, x), generates fdist_Beta (p + n, q, x)
 * for n = 0, 1, 2,..., nmax to d significant digits, using a variant of
 * J.C.P. Miller's backward recurrence algorithm. Stores the result in I. 
 */
{

   int n, nu, m, again, ntab;
   double *Itemp, *Iapprox, *Rr;
   double epsilon, r;

   I[0] = I0;
   if (nmax == 0)
      return;

   epsilon = EpsArray[d];
   nu = 2 * nmax + 5;
   ntab = 64;
   while (ntab <= nu)
      ntab *= 2;

   Rr = (double *) util_Calloc ((size_t) ntab, sizeof (double));
   Iapprox = (double *) util_Calloc ((size_t) ntab, sizeof (double));
   Itemp = (double *) util_Calloc ((size_t) ntab, sizeof (double));

   for (n = 1; n <= nmax; n++)
      Iapprox[n] = 0.0;
   for (n = 0; n <= nmax; n++)
      Itemp[n] = I[n];

   do {
      n = nu;
      r = 0.0;
      do {
         r = (n - 1 + p + q) * x / (n + p + (n - 1 + p + q) * x - (n + p) * r);
         if (n <= nmax)
            Rr[n - 1] = r;
         n--;
      } while (n >= 1);

      for (n = 0; n < nmax; n++)
         Itemp[n + 1] = Rr[n] * Itemp[n];

      again = 0;
      for (n = 1; n <= nmax; n++) {
         if (fabs ((Itemp[n] - Iapprox[n])/Itemp[n]) > epsilon) {
            again++;
            for (m = 1; m <= nmax; m++)
               Iapprox[m] = Itemp[m];
            nu += 5;
            if (ntab <= nu) {
               ntab *= 2;
               Rr = (double *) util_Realloc (Rr, ntab * sizeof (double));
               Iapprox = (double *) util_Realloc (Iapprox, ntab * sizeof (double));
               Itemp = (double *) util_Realloc (Itemp, ntab * sizeof (double));
            }
            break;
         }
      }
   } while (again);

   for (n = 0; n <= nmax; n++)
      I[n] = Itemp[n];
   util_Free (Rr);
   util_Free (Iapprox);
   util_Free (Itemp);
}

/*-------------------------------------------------------------------------*/
static const double RENORM = 1.0e300;

static void Isubx_q_fixed (double p, double q, double x, int d, int nmax,
   double I[])
/* 
 * Generates fdist_Beta (p + n, q, x), 0 < p <= 1, for n = 0, 1, 2,...,
 * nmax to d significant digits, using procedure backward. First reduces
 * q modulo 1 to q0, where 0 < q0 <= 1.
 */
{

   int m, mmax;
   double s, q0, Iq0, Iq1;
   double *Iq;

   util_Assert (p > 0.0 && p <= 1.0, "Isubx_q_fixed:   p not in (0, 1] ");
   m = (int) q;                   /* integer part of q */
   s = q - m;                     /* fractionnal part of q */
   if (s > 0) {
      q0 = s;
      mmax = m;
   } else {
      q0 = s + 1;
      mmax = m - 1;
   }
   Iq0 = RENORM * Isubx_pq_small (p, q0, x, d);
   if (mmax > 0)
      Iq1 = RENORM * Isubx_pq_small (p, q0 + 1.0, x, d);

   Iq = (double *) util_Calloc ((size_t) mmax + 1, sizeof (double));
   forward (p, q0, x, Iq0, Iq1, mmax, Iq);
   backward (p, q, x, Iq[mmax], d, nmax, I);
   for (m = 0; m <= nmax; m++)
      I[m] /= RENORM;
   util_Free (Iq);
}

/*-------------------------------------------------------------------------*/

static void Isubx_p_fixed (double p, double q, double x, int d, int nmax,
   double I[])
/* 
 * Generates fdist_Beta (p, q + n, x), 0 < q <= 1, for n = 0, 1, 2,...,
 * nmax to d significant digits, using procedure forward.
 */
{

   int m, mmax;
   double s, p0, I0, Iq0, I1, Iq1;
   double *Ip;

   util_Assert (q > 0.0 && q <= 1.0, "Isubx_p_fixed:   q not in (0, 1] ");

   m = (int) p;                   /* integer part of p */
   s = p - m;                     /* fractionnal part of p */
   if (s > 0) {
      p0 = s;
      mmax = m;
   } else {
      p0 = s + 1;
      mmax = m - 1;
   }
   I0 = RENORM * Isubx_pq_small (p0, q, x, d);
   I1 = RENORM * Isubx_pq_small (p0, q + 1.0, x, d);

   Ip = (double *) util_Calloc ((size_t) mmax + 1, sizeof (double));
   backward (p0, q, x, I0, d, mmax, Ip);
   Iq0 = Ip[mmax];
   backward (p0, q + 1.0, x, I1, d, mmax, Ip);
   Iq1 = Ip[mmax];
   forward (p, q, x, Iq0, Iq1, nmax, I);
   for (m = 0; m <= nmax; m++)
      I[m] /= RENORM;
   util_Free (Ip);
}

/*-------------------------------------------------------------------------*/

static void Beta_q_fixed (double p, double q, double x, int d, int nmax,
   double I[])
{
   int n;
   util_Assert (p > 0.0 && p <= 1.0, "Beta_q_fixed:   p not in (0, 1]");
   util_Assert (q > 0.0, "Beta_q_fixed:   q <= 0");
   util_Assert (nmax >= 0, "Beta_q_fixed:   nmax < 0");
   if (x == 0.0 || x == 1.0) {
      for (n = 0; n <= nmax; n++)
         I[n] = x;
      return;
   }
   if (x <= 0.5)
      Isubx_q_fixed (p, q, x, d, nmax, I);
   else {
      Isubx_p_fixed (q, p, 1.0 - x, d, nmax, I);
      for (n = 0; n <= nmax; n++)
         I[n] = 1.0 - I[n];
   }
}

/*-------------------------------------------------------------------------*/

static void Beta_p_fixed (double p, double q, double x, int d, int nmax,
   double I[])
{
   int n;
   util_Assert (q > 0.0 && q <= 1.0, "Beta_p_fixed:  q not in (0, 1]");
   util_Assert (p > 0.0, "Beta_p_fixed:   p <= 0");
   util_Assert (nmax >= 0, "Beta_p_fixed:  nmax < 0");
   if (x == 0.0 || x == 1.0) {
      for (n = 0; n <= nmax; n++)
         I[n] = x;
      return;
   }
   if (x <= 0.5)
      Isubx_p_fixed (p, q, x, d, nmax, I);
   else {
      Isubx_q_fixed (q, p, 1.0 - x, d, nmax, I);
      for (n = 0; n <= nmax; n++)
         I[n] = 1.0 - I[n];
   }
}


/*-------------------------------------------------------------------------*/
/*
 * The exact section of fdist_Beta below is very slow for large parameters.
 * It is an old algorithm of Gautschi of 1964. There is an algorithm
 * for fdist_Beta (1994) that is recent and is supposed to be very fast
 * (I MUST write the exact reference for later; I think it may have been in
 * Mathematics of Computations???) 
 */
double fdist_Beta (double p, double q, int d, double x)
/*
 * I[j] will contain either the values of fdist_Beta (p0 + j, q, d, x),
 * where 0 < p0 <= 1, for j = 0, 1, 2, ..., n,  with p = p0 + n; or the
 * values of fdist_Beta (p, q0 + j, d, x), where 0 < q0 <= 1, for j = 0,
 * 1, 2, ..., n, with q = q0 + n.
 */
{

   const double pqmax = 1000.0;
   const double pqlim = 30.0;
   int n, flag;
   double p0, q0, u, temp, yd, gam, h1, h3, y;
   double *I;

   util_Assert (p > 0.0, "fdist_Beta:   p <= 0");
   util_Assert (q > 0.0, "fdist_Beta:   q <= 0");
   util_Assert (d > 0, "fdist_Beta:   d <= 0");
   util_Assert (d < 16, "fdist_Beta:   d > 15");
   if (x <= 0.0)
      return 0.0;
   if (x >= 1.0)
      return 1.0;

   if (util_Max (p, q) <= pqmax) {
      if (p < q) {
         n = (int) p;             /* integer part of p */
         p0 = p - n;              /* fractionnal part of p */
         if (p0 <= 0.0) {         /* p0 == 0 not allowed */
            p0 = 1.0;
            n--;
         }
         I = (double *) util_Calloc ((size_t) n + 1, sizeof (double));
         Beta_q_fixed (p0, q, x, d, n, I);
         u = I[n];
         util_Free (I);
         /* There may be numerical errors far in the tails giving very small
            negative values instead of 0. */
         if (u <= 0.0)
            return 0.0;
         else if (u <= 1.0)
            return u;
         else
            return 1.0;

      } else {
         n = (int) q;             /* integer part of q */
         q0 = q - n;              /* fractionnal part of q */
         if (q0 <= 0.0) {         /* q0 == 0 not allowed */
            q0 = 1.0;
            n--;
         }
         I = (double *) util_Calloc ((size_t) n + 1, sizeof (double));
         Beta_p_fixed (p, q0, x, d, n, I);
         u = I[n];
         util_Free (I);
         /* There may be numerical errors far in the tails giving very small
            negative values instead of 0. */
         if (u <= 0.0)
            return 0.0;
         else if (u <= 1.0)
            return u;
         else
            return 1.0;
      }
   }

   if ((p > pqmax && q < pqlim) || (q > pqmax && p < pqlim)) {
      /* Bol'shev approximation for large max(p, q) and small min(p, q) */
      if (x > 0.5)
         return 1.0 - fdist_Beta (q, p, d, 1.0 - x);

      if (p < q) {
         u = p;
         p = q;
         q = u;
         flag = 0;
      } else {
         flag = 1;
      }
      u = p + 0.5 * q - 0.5;
      if (!flag)
         temp = x / (2.0 - x);
      else
         temp = (1.0 - x) / (1.0 + x);
      yd = 2.0 * u * temp;
      gam =
         (exp (q * log (yd) - yd - num2_LnGamma (q)) * (2.0 * yd * yd - (q -
               1.0) * yd - (q * q - 1.0))) / (24.0 * u * u);
      if (flag) {
         yd = fbar_Gamma (q, d, yd);
         return yd - gam;
      } else {
         yd = fdist_Gamma (q, d, yd);
         return yd + gam;
      }
   }

   /* Normal approximation of Peizer and Pratt */
   h1 = p + q - 1.0;
   y = 1.0 - x;
   h3 = sqrt ((1.0 + y * fdist_belog ((p - 0.5) / (h1 * x))
         + x * fdist_belog ((q - 0.5) / (h1 * y)))
      / ((h1 + 1.0 / 6.0) * x * y))
      * ((h1 + 1.0 / 3.0 + 0.02 * (1.0 / p + 1.0 / q + 1.0 / (p + q)))
      * x - p + 1.0 / 3.0 - 0.02 / p - 0.01 / (p + q));

   return fdist_Normal2 (h3);

}


/*=========================================================================*/
#define EPSILON  1.0e-15          /* Tolerance */
#define EPSBETA  0.5e-10          /* < 0.75 sqrt(DBL_EPSILON) */
#define ALPHALIM 100000.0         /* Limiting alpha for normal approx. */
#define MAXJ  2000                /* Max number of terms in series */
#define INV2PI 0.6366197723675813 /* 2 / PI */
#define LOG4  1.38629436111989062 /* Ln(4) */
#define OneRac2  0.70710678118654752     /* 1/sqrt(2) */
#define SQPI_2  0.88622692545275801   /* Sqrt(Pi) / 2 */
#define LOG_SQPI_2 -0.1207822376352453  /* Ln(Sqrt(Pi) / 2) */



/*------------------------------------------------------------------------*/

static double series1 (double alpha, double x)
/* 
 * Compute the series for F(x).
 * This series is used for alpha < 1 and x close to 0.
 */
{
   int j;
   double sum, term;
   double poc;
   poc = 1.0;
   sum = 1.0 / alpha;
   j = 1;
   do {
      poc *= x * (j - alpha) / j;
      term = poc / (j + alpha);
      sum += term;
      ++j;
   } while ((term > sum * EPSILON) && (j < MAXJ));

   return sum * pow (x, alpha);
}


/*------------------------------------------------------------------------*/

static double series2 (double alpha, double y)
/* 
 * Compute the series for G(y).   y = 0.5 - x.
 * This series is used for alpha < 1 and x close to 1/2.
 */
{
   int j;
   double term, sum;
   double poc;
   const double z = 4.0 * y * y;

   /* Compute the series for G(y) */
   poc = sum = 1.0;
   j = 1;
   do {
      poc *= z * (j - alpha) / j;
      term = poc / (2 * j + 1);
      sum += term;
      ++j;
   } while ((term > sum * EPSILON) && (j < MAXJ));

   return sum * y;
}


/*------------------------------------------------------------------------*/

static double series3 (double alpha, double x)
/* 
 * Compute the series for F(x).
 * This series is used for alpha > 1 and x close to 0.
 */
{
   int j;
   double sum, term;
   const double z = -x / (1.0 - x);

   sum = term = 1.0;
   j = 1;
   do {
      term *= z * (j - alpha) / (j + alpha);
      sum += term;
      ++j;
   } while ((fabs (term) > sum * EPSILON) && (j < MAXJ));

   return sum * x;
}


/*------------------------------------------------------------------------*/

static double series4 (double alpha, double y)
/* 
 * Compute the series for G(y).   y = 0.5 - x.
 * This series is used for alpha > 1 and x close to 1/2.
 */
{
   int j;
   double term, sum;
   const double z = 4.0 * y * y;

   term = sum = 1.0;
   j = 1;
   do {
      term *= z * (j + alpha - 0.5) / (0.5 + j);
      sum += term;
      ++j;
   } while ((term > sum * EPSILON) && (j < MAXJ));

   return sum * y;
}


/*-------------------------------------------------------------------------*/

static double Peizer (double alpha, double x)
/*
 * Normal approximation of Peizer and Pratt
 */
{
   const double y = 1.0 - x;
   double z;
   z = sqrt ((1.0 - y * fdist_belog (2.0 * x) - x * fdist_belog (2.0 * y))
      / ((2.0*alpha - 5.0 / 6.0) * x * y)) * 
      (2.0*x - 1.0) * (alpha - 1.0 / 3.0 + 0.025 / alpha);

   return fdist_Normal2 (z);
}


/*------------------------------------------------------------------------*/

void fdist_CalcB4 (double alpha, double *pB, double *plogB, double *pC,
                   double *plogC)
{
   /* Compute Beta(alpha, alpha) and Beta(alpha, alpha)*4^(alpha-1). */
   double temp;

   if (alpha <= EPSBETA) {
      /* For a -> 0, B(a,a) = (2/a)*(1 - 1.645*a^2 + O(a^3)) */
      *pB = 2.0 / alpha;
      *pC = *pB / (4.0*(1.0 - alpha*LOG4));

   } else if (alpha <= 1.0) {
      *plogB = 2.0 * num2_LnGamma (alpha) - num2_LnGamma (2.0*alpha);
      *plogC = *plogB + (alpha - 1.0)*LOG4;
      *pC = exp(*plogC);
      *pB = exp(*plogB);

   } else if (alpha <= 10.0) {
      *plogC = num2_LnGamma (alpha) - num2_LnGamma (0.5 + alpha) + LOG_SQPI_2;
      *plogB = *plogC - (alpha - 1.0)*LOG4;      

   } else if (alpha <= 200.0) {
      /* Convergent series for Gamma(x + 0.5) / Gamma(x) */
      double term = 1.0;
      double sum = 1.0;
      int i = 1;
      while (term > EPSILON*sum) {
         term *= (i - 1.5)*(i - 1.5) /(i*(alpha + i - 1.5));
         sum += term;
         i++;
      }
      temp = SQPI_2 / sqrt ((alpha - 0.5)*sum);
      *plogC = log(temp);
      *plogB = *plogC - (alpha - 1.0)*LOG4;

   } else {
      /* Asymptotic series for Gamma(a + 0.5) / (Gamma(a) * Sqrt(a)) */
      double u = 1.0 / (8.0*alpha);
      temp = 1.0 + u*(-1.0 + u*(0.5 + u*(2.5 - u*(2.625 + 49.875*u))));
      /* This is 4^(alpha - 1)*B(alpha, alpha) */
      temp = SQPI_2 / (sqrt(alpha) * temp);
      *plogC = log(temp);
      *plogB = *plogC - (alpha - 1.0)*LOG4;
   }
}


/*------------------------------------------------------------------------*/

double fdist_BetaSymmetric (double alpha, double x)
/* 
 * Compute the cumulative probability of the symmetrical beta distribution.
 * Returns a negative value on error, otherwise returns u in [0, 1].
 */
{
   double temp, u, logB, logC;
   int isUpper;                   /* True if x > 0.5 */
   double B = 0.;                 /* Beta(alpha, alpha) */
   double C, x0;

   if (alpha <= 0.0) {
      util_Assert (1, "fdist_BetaSymmetric:   p <= 0\n");
      return -1.0;
   }
   if (x <= 0.0) return 0.0;
   if (x >= 1.0) return 1.0;
   if (x == 0.5) return 0.5;
   if (alpha == 1.0) return x;         /* alpha = 1 is the uniform law */
   if (alpha == 0.5)                   /* alpha = 1/2 is the arcsin law */
      return INV2PI * asin(sqrt(x));

   if (alpha > ALPHALIM)
      return Peizer (alpha, x);

   if (x > 0.5) {
      x = 1.0 - x;
      isUpper = 1;
   } else
      isUpper = 0;

   fdist_CalcB4 (alpha, &B, &logB, &C, &logC);

   if (alpha <= 1.0) {
      /* For x = x0, both series use the same number of terms to get the
         required precision */
      if (x > 0.25) {
         temp = -log (alpha);
         if (alpha >= 1.0e-6)
            x0 = 0.25 + 0.005 * temp;
         else
            x0 = 0.13863 + .01235 * temp;
      } else
         x0 = 0.25;

      if (x <= x0)
         u = (series1 (alpha, x)) / B;
      else
         u = 0.5 - (series2 (alpha, 0.5 - x)) / C;

   } else {                       /* 1 < alpha < ALPHALIM */
      if (alpha < 400.0)
         x0 = 0.5 - 0.45 / sqrt(alpha);
      else
         x0 = 0.5 - 1.0 / sqrt(alpha);
      if (x0 < 0.25)
         x0 = 0.25;

      if (x <= x0) {
         temp = (alpha - 1.0) * log (x * (1.0 - x)) - logB;
         u = series3 (alpha, x) * exp(temp) / alpha;

      } else {      
         const double y = 0.5 - x;
         temp = num2_log1p(-4.0*y*y);
         temp = alpha * temp - logC;
         u = 0.5 - (series4 (alpha, y)) * exp(temp);
      }
   }

   if (isUpper)
      return 1.0 - u;
   else
      return u;
}


/*=========================================================================*/
#define NLIM 20

static double KSSpecial (long n, double x)
{
   /* For nx^2 > 18, fbar_KS(n, x) is smaller than DBL_EPSILON */
   if ((n*x*x >= 18.0) || (x >= 1.0))
      return 1.0;

   if (x <= 0.5 / n)
      return 0.0;

   if (n == 1)
      return 2.0 * x - 1.0;

   if (x <= 1.0 / n) {
      double t = 2.0 * x - 1.0 / n;
      double w;
      if (n <= NLIM) {
         w = num2_Factorial ((int) n);
         return w * pow (t, (double) n); 
      }
      w = num2_LnFactorial ((int) n) + n * log (t);
      return exp (w);
   }

   if (x >= 1.0 - 1.0 / n) {
      return 1.0 - 2.0 * pow (1.0 - x, (double) n);
   }

   return -1.0;
}

#undef NLIM
/*-------------------------------------------------------------------------*/

static double Pelz (long n, double x)
{
   /*
      Approximating the Lower Tail-Areas of the Kolmogorov-Smirnov
         One-Sample Statistic,
      Wolfgang Pelz and I. J. Good,
      Journal of the Royal Statistical Society, Series B.
      Vol. 38, No. 2 (1976), pp. 152-156
    */

   const int JMAX = 20;
   const double EPS = 1.0e-10;
   const double C = 2.506628274631001;  /* sqrt(2*Pi) */
   const double C2 = 1.2533141373155001;  /* sqrt(Pi/2) */
   const double PI2 = num_Pi * num_Pi;
   const double PI4 = PI2 * PI2;
   const double RACN = sqrt((double)n);
   const double z = RACN*x;
   const double z2 = z * z;
   const double z4 = z2 * z2;
   const double z6 = z4 * z2;
   const double w = PI2 / (2.0 * z*z);
   double ti, term, tom;
   double sum;
   int j;

   term = 1;
   j = 0;
   sum = 0;
   while (j <= JMAX && term > EPS * sum) {
      ti = j + 0.5;
      term = exp (-ti * ti * w);
      sum += term;
      j++;
   }
   sum *= C / z;

   term = 1;
   tom = 0;
   j = 0;
   while (j <= JMAX && fabs(term) > EPS * fabs(tom)) {
      ti = j + 0.5;
      term = (PI2 * ti * ti - z2) * exp (-ti * ti * w);
      tom += term;
      j++;
   }
   sum += tom * C2 / (RACN * 3.0 * z4);

   term = 1;
   tom = 0;
   j = 0;
   while (j <= JMAX && fabs(term) > EPS * fabs(tom)) {
      ti = j + 0.5;
      term = 6*z6 + 2*z4 + PI2*(2*z4 - 5*z2)*ti*ti +
             PI4*(1 - 2*z2)*ti*ti*ti*ti;
      term *= exp (-ti * ti * w);
      tom += term;
      j++;
   }
   sum += tom * C2 / (n * 36.0 * z * z6);

   term = 1;
   tom = 0;
   j = 1;
   while (j <= JMAX && term > EPS * tom) {
      ti = j;
      term = PI2 * ti * ti * exp (-ti * ti * w);
      tom += term;
      j++;
   }
   sum -= tom * C2 / (n * 18.0 * z * z2);

   term = 1;
   tom = 0;
   j = 0;
   while (j <= JMAX && fabs(term) > EPS * fabs(tom)) {
      ti = j + 0.5;
      ti = ti * ti;
      term = -30*z6 -90*z6*z2 + PI2*(135*z4 - 96*z6)*ti +
         PI4*(212*z4 - 60*z2)*ti*ti + PI2*PI4*ti*ti*ti*(5 - 30*z2);
      term *= exp (-ti * w);
      tom += term;
      j++;
   }
   sum += tom * C2 / (RACN * n * 3240.0 * z4 * z6);

   term = 1;
   tom = 0;
   j = 1;
   while (j <= JMAX && fabs(term) > EPS * fabs(tom)) {
      ti = j*j;
      term = (3*PI2 * ti * z2 - PI4*ti*ti) * exp (-ti * w);
      tom += term;
      j++;
   }
   sum += tom * C2 / (RACN * n * 108.0 * z6);

   return sum;
}


/*=========================================================================*/

static void mMultiply (double *A, double *B, double *C, int m)
{
   int i, j, k;
   double s;
   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++) {
         s = 0.0;
         for (k = 0; k < m; k++)
            s += A[i * m + k] * B[k * m + j];
         C[i * m + j] = s;
      }
}


/*-------------------------------------------------------------------------*/

static void mPower (double *A, int eA, double *V, int *eV, int m, int n)
{
   double *B;
   int eB, i;
   if (n == 1) {
      for (i = 0; i < m * m; i++)
         V[i] = A[i];
      *eV = eA;
      return;
   }
   mPower (A, eA, V, eV, m, n / 2);
   B = (double *) malloc (m * m * sizeof (double));
   mMultiply (V, V, B, m);
   eB = 2 * (*eV);

   if (n % 2 == 0) {
      for (i = 0; i < m * m; i++)
         V[i] = B[i];
      *eV = eB;
   } else {
      mMultiply (A, B, V, m);
      *eV = eA + eB;
   }

   if (V[(m / 2) * m + (m / 2)] > 1.0e140) {
      for (i = 0; i < m * m; i++)
         V[i] = V[i] * 1.0e-140;
      *eV += 140;
   }
   free (B);
}


/*-------------------------------------------------------------------------*/

double fdist_KS2 (long N0, double x)
{
   int k, m, i, j, g, eH, eQ;
   const int n = N0;
   const double d = x;
   double h, s, *H, *Q;

   /* OMIT NEXT 3 LINES IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL */
#if 1
   s = d * d * n;
   if (s > 7.24 || (s > 3.76 && n > 99))
      return 1 - 2 * exp (-(2.000071 + 0.331 / sqrt ((double) n) +
            1.409 / n) * s);
#endif
   k = (int) (n * d) + 1;
   m = 2 * k - 1;
   h = k - n * d;
   H = (double *) malloc (m * m * sizeof (double));
   Q = (double *) malloc (m * m * sizeof (double));

   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
         if (i - j + 1 < 0)
            H[i * m + j] = 0;
         else
            H[i * m + j] = 1;

   for (i = 0; i < m; i++) {
      H[i * m] -= pow (h, (double) (i + 1));
      H[(m - 1) * m + i] -= pow (h, (double) (m - i));
   }

   H[(m - 1) * m] += (2 * h - 1 > 0 ? pow (2 * h - 1, (double) m) : 0);

   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
         if (i - j + 1 > 0)
            for (g = 1; g <= i - j + 1; g++)
               H[i * m + j] /= g;

   eH = 0;
   mPower (H, eH, Q, &eQ, m, n);
   s = Q[(k - 1) * m + k - 1];
   for (i = 1; i <= n; i++) {
      s = s * i / n;
      if (s < 1.0e-140) {
         s *= 1.0e140;
         eQ -= 140;
      }
   }
   s *= pow (10.0, (double) eQ);
   free (H);
   free (Q);
   return s;
}


/*=========================================================================*/

static double Pomeranz (long n, double x)
{
   const double EPS = 5.0e-13;   /* for floors and ceilings */
   const int ENO = 350;
   const double RENO = ldexp(1.0, ENO);   /* for renormalization of V */
   const double IRENO = 1.0/RENO;
   int coreno;               /* counter: how many renormalizations */
   const double t = n*x;
   double sum, maxsum;
   int i, j, k, s;
   int r1, r2;                   /* Indices i and i-1 for V[i][] */
   int jlow, jup, klow, kup, kup0;
   double w, z;
   double *A;
   double **V;
   double **H;   /* work variables = pow(w, j-k) / Factorial(j-k) */

   A = (double*) util_Calloc ((size_t) (2*n + 3), sizeof (double));
   V = (double **) tables_CreateMatrixD (2, n + 2);
   H = (double **) tables_CreateMatrixD (4, n + 1);

   A[0] = A[1] = 0;
   z = t - floor(t);
   w = ceil(t) - t;
   if (w < z)
      z = w;
   A[2] = z;
   A[3] = 1.0 - A[2];
   for (i = 4; i <= 2*n + 1; i++)
      A[i] = A[i-2] + 1.0;
   A[2*n + 2] = n;

   for (j = 1; j <= n+1; j++)
      V[0][j] = 0;
   for (j = 2; j <= n+1; j++)
      V[1][j] = 0;
   V[1][1] = RENO;
   coreno = 1;

   /* Precompute H[][] = (A[j] - A[j-1]^(j-k) / (j-k)! for speed */
   H[0][0] = 1;
   w = 2.0 * A[2] / n;
   for (j = 1; j <= n; j++)
      H[0][j] = w * H[0][j - 1] / j;

   H[1][0] = 1;
   w = (1.0 - 2.0*A[2])/n;
   for (j = 1; j <= n; j++)
      H[1][j] = w*H[1][j-1] / j;

   H[2][0] = 1;
   w = A[2]/n;
   for (j = 1; j <= n; j++)
      H[2][j] = w*H[2][j-1] / j;

   H[3][0] = 1;
   for (j = 1; j <= n; j++)
      H[3][j] = 0;

   r1 = 0;
   r2 = 1;
   for (i = 2; i <= 2 * n + 2; i++) {
      jlow = 2 + floor (A[i] - t + EPS);
      if (jlow < 1)
         jlow = 1;
      jup = ceil (A[i] + t - EPS);
      if (jup > n + 1)
         jup = n + 1;

      klow = 2 + floor (A[i - 1] - t + EPS);
      if (klow < 1)
         klow = 1;
      kup0 = ceil (A[i - 1] + t - EPS);

      /* Find to which case it corresponds */
      w = (A[i] - A[i-1])/n;
      s = -1;
      for (j = 0; j < 4; j++) {
         if (fabs(w - H[j][1]) <= EPS) {
            s = j;
            break;
         }
      }
      util_Assert (s >= 0, "Pomeranz:   s < 0");

      maxsum = -1;
      r1 = (r1 + 1) & 1;    /* i - 1  */
      r2 = (r2 + 1) & 1;    /* i */

      for (j = jlow; j <= jup; j++) {
         kup = kup0;
         if (kup > j)
            kup = j;
         sum = 0;
         for (k = kup; k >= klow; k--)
            sum += V[r1][k] * H[s][j - k];
         V[r2][j] = sum;
         if (sum > maxsum)
            maxsum = sum;
      }

      if (maxsum < IRENO) {
         /* V is too small: renormalize to avoid underflow of prob */
         for (j = jlow; j <= jup; j++)
            V[r2][j] *= RENO;
         coreno++;    /* keep track of log of RENO */
      }
   }

   z = V[r2][n+1];
   util_Free (A);
   tables_DeleteMatrixD (&H);
   tables_DeleteMatrixD (&V);

   w = num2_LnFactorial(n) - coreno*ENO*num_Ln2 + log(z);
   if (w >= 0.)
      return 1.;
   return exp(w);
}


/*-------------------------------------------------------------------------*/
#define NSEP  400
#define NSEP2 4000
#define ZSEP  4.0
#define ZSEP2 0.2

double fdist_KS1 (long n, double x)
{
   double u = KSSpecial(n, x);
   if (u >= 0.0)
      return u;

   if (n <= NSEP) {
      if (n*x*x < ZSEP)
         return Pomeranz (n, x);
      else 
         return 1. - fbar_KS1(n, x);
   }

   if (n*x*x <= ZSEP2 && n <= NSEP2)
      return Pomeranz (n, x);

   return Pelz (n, x);
}


/*=========================================================================*/

double fdist_KSPlus (long N, double x)
{
   const double NxParam = 6.5;    /* frontier: alternating series */
   const long NParam = 4000;      /* frontier: non-alternating series */
   double q;
   double Sum;
   double term;

   util_Assert (N > 0, "Calling fdist_KSPlus with N < 1");
   if (x <= 0.0)
      return 0.0;
   if ((x >= 1.0) || (N*x*x >= 25.0))
      return 1.0;
   if (N == 1)
      return x;

   /*--------------------------------------------------------------*/
   /* the alternating series is stable and fast for N*x very small */
   /*--------------------------------------------------------------*/

   if (N * x <= NxParam) {
      const double Epsilon = 1.0E-300;
      double LogCom = log ((double) N);
      int Sign = -1;
      long j;
      long jmax = (long) (N * x);
      Sum = 0.0;

      for (j = 1; j <= jmax; j++) {
         double jreal = j;
         double Njreal = N - j;
         q = jreal / N - x;
         /* we must avoid log(0.0) for j = jmax and N*x near an integer */
         if (-q > Epsilon) {
            term = LogCom + jreal * log (-q) + (Njreal - 1.0) * num2_log1p (-q);
            Sum += Sign * exp (term);
         }
         Sign = -Sign;
         LogCom += log (Njreal / (j + 1));
      }
      /* add the term j = 0 */
      Sum += exp ((N - 1) * num2_log1p (x));
      if (Sum >= 0.0)
         return Sum * x;
      else
         return 0.0;
   }

   if (N <= NParam) {
      double Njreal;
      double jreal;
      long j;
      long jmax;
      double LogCom = log ((double) N);
      Sum = 0.0;
      jmax = (long) (N * (1.0 - x));
      if (1.0 - x - (double) jmax/N <= 0.0)
         --jmax;

      for (j = 1; j <= jmax; j++) {
         jreal = j;
         Njreal = N - j;
         q = jreal / N + x;
         term = LogCom + (jreal - 1.0) * log (q) + Njreal * num2_log1p(-q);
         Sum += exp (term);
         LogCom += log (Njreal / (jreal + 1.0));
      }
      Sum *= x;

      /* add the term j = 0; avoid log(0.0) */
      if (1.0 > x)
         Sum += exp (N * num2_log1p(-x));
      Sum = 1.0 - Sum;
      if (Sum >= 0.0)
         return Sum;
      else
         return 0.0;
   }

   /*---------------------------*/
   /* Use an asymptotic formula */
   /*---------------------------*/

   term = 2.0 / 3.0;
   q = x * x * N;
   Sum = 1.0 - exp (-2.0 * q) * (1.0 - term * x * (1.0 - x * (1.0 - term * q)
                    - term / N * (0.2 - 19.0 / 15.0 * q + term * q * q)));
   if (Sum >= 0.0)
      return Sum;
   else
      return 0.0;
}


/*=========================================================================*/

double fdist_KSPlusJumpOne (long N, double a, double x)
{
   const double EpsilonLR = 1.E-15;
   const double Epsilon = 1.0E-290;
   const double NxaParam = 6.5;   /* frontier: alternating series */
   double LogCom;
   double q, p1, q1;
   double Sum = 0.0;
   double term;
   double Njreal;
   double jreal;
   int Sign;
   long j;
   long jmax;

   util_Assert (N >= 1, "Calling fdist_KSPlusJumpOne with N < 1");
   util_Assert (a < 1.0 && a > 0.0,
      "Calling fdist_KSPlusJumpOne with a outside (0, 1)");
   if (x <= 0.0)
      return 0.0;
   if (x + a >= 1.0)
      return 1.0;
   LogCom = log ((double) N);

   /*--------------------------------------------------------------------*/
   /* the alternating series is stable and fast for N*(x + a) very small */
   /*--------------------------------------------------------------------*/
   if (N * (x + a) < NxaParam && a + x < 0.5) {
      jmax = (long) (N * (x + a));
      for (j = 1; j <= jmax; j++) {
         jreal = j;
         Njreal = N - j;
         q = jreal / N - x;
         if ((q < 0.0 && (j & 1)) || ((q > 1.0) && ((N - j - 1) & 1)))
            Sign = -1;
         else
            Sign = 1;

         /* we must avoid log(0.0) */
         q1 = fabs (q);
         p1 = fabs (1.0 - q);
         if (q1 > Epsilon && p1 > Epsilon) {
            term = LogCom + jreal * log (q1) + (Njreal - 1.0) * log (p1);
            Sum += Sign * exp (term);
         }
         LogCom += log (Njreal / (jreal + 1.0));
      }
      /* add the term j = 0 */
      Sum += exp ((N - 1) * num2_log1p(x));
      return Sum * x;
   }

   /*---------------------------------------------*/
   /* For N(x + a) >= NxaParam or (a + x) > 0.5, */
   /* use the non-alternating series.  */
   /*---------------------------------------------*/

   /* EpsilonLR because the distribution has a jump */
   jmax = (long) (N * (1.0 - a - x - EpsilonLR));
   for (j = 1; j <= jmax; j++) {
      jreal = j;
      Njreal = N - jreal;
      q = jreal / N + x;
      if (1.0 - q > Epsilon) {
         term = LogCom + (jreal - 1.0) * log (q) + Njreal * num2_log1p (-q);
         Sum += exp (term);
      }
      LogCom += log (Njreal / (jreal + 1.0));
   }
   Sum *= x;

   /* add the term j = 0 */
   if (1.0 - x > Epsilon)
      Sum += exp (N * num2_log1p (-x));
   return 1.0 - Sum;
}


/*=========================================================================*/
#if 0
static lebool IsJump (fdist_FUNC_JUMPS * H, double xa, double xb,
   double ya, double yb, int NJ)
   /* Find a more precise value for the position of the jump in (xa, xb). */
   /* Return FALSE if there is no jump, TRUE if there is a jump. */
{
   const double eps = DBL_EPSILON;
   const int imax = DBL_MANT_DIG;
   const double epsY = H->epsY;
   double *par = H->par;
   wdist_CFUNC F = H->F;
   int i = 0;
   double x = 1.0, y;

   /* Binary search to refine the x-coordinate of the jump */
   while ((i < imax) && (xb - xa > eps * x)) {
      i++;
      x = (xb + xa) / 2.0;
      y = F (par, x);
      if (y - ya > epsY) {
         yb = y;
         xb = x;
      } else {
         ya = y;
         xa = x;
      }
   }

   if (yb - ya < epsY)
      return FALSE;
   H->xJump[NJ] = (xb + xa) / 2.0;
   H->yLeftJump[NJ] = ya;
   H->yRightJump[NJ] = yb;
   return TRUE;
}

/*-------------------------------------------------------------------------*/

void fdist_FindJumps (fdist_FUNC_JUMPS * H, int Detail)
{
   int i, NJ;
   double yRight;
   double yLeft;
   double x;
   const double epsX = H->epsX;
   const double epsY = H->epsY;
   double *par = H->par;
   wdist_CFUNC F = H->F;

   /* Assume no more than 30 jumps initially */
   NJ = 30;
   H->xJump = (double *) util_Calloc ((size_t) NJ + 1, sizeof (double));
   H->yLeftJump = (double *) util_Calloc ((size_t) NJ + 1, sizeof (double));
   H->yRightJump = (double *) util_Calloc ((size_t) NJ + 1, sizeof (double));

   i = 0;
   if (H->xa > H->xb) {
      x = H->xa;
      H->xa = H->xb;
      H->xb = x;
   }
   x = H->xa;
   yLeft = F (par, x);
   while (x < H->xb) {
      x += epsX;
      yRight = F (par, x);
      if (yRight - yLeft > epsY) {
         /* this should be a jump */
         ++i;
         if (i > NJ) {
            NJ *= 2;
            H->xJump = (double *) util_Realloc (H->xJump,
               (NJ + 1) * sizeof (double));
            H->yLeftJump = (double *) util_Realloc (H->yLeftJump,
               (NJ + 1) * sizeof (double));
            H->yRightJump = (double *) util_Realloc (H->yRightJump,
               (NJ + 1) * sizeof (double));
         }
         if (IsJump (H, x - epsX, x, yLeft, yRight, i) == FALSE)
            i--;
      }
      yLeft = yRight;
   }
   NJ = i;
   H->xJump = (double *) util_Realloc (H->xJump, (NJ + 1) * sizeof (double));
   H->yLeftJump = (double *) util_Realloc (H->yLeftJump,
                  (NJ + 1) * sizeof (double));
   H->yRightJump = (double *) util_Realloc (H->yRightJump,
                  (NJ + 1) * sizeof (double));
   H->NJumps = NJ;

   if (Detail > 0) {
      printf ("\n=========================================================");
      printf ("\nCalling fdist_FindJumps for function  %-32s\n", H->doc);
      printf ("\nInterval = (%g, %g)\n", H->xa, H->xb);
      printf ("epsX = %10.5g\nepsY = %10.5g\n", epsX, epsY);
      printf ("Number of jumps = %4d\n", NJ);
      if (NJ == 0)
         return;
      printf ("Jumps of the function:\n\n");
      printf ("           x                   yLeft              yRight"
         "         yRight - yLeft\n\n");
      for (i = 1; i <= NJ; i++) {
         printf (" %19.15g %19.15g %19.15g %19.15g\n", H->xJump[i],
            H->yLeftJump[i], H->yRightJump[i],
            H->yRightJump[i] - H->yLeftJump[i]);
      }
      printf
         ("\n=========================================================\n\n");
   }
}


/*=========================================================================*/

void fdist_FreeJumps (fdist_FUNC_JUMPS * H)
{
   util_Free (H->xJump);
   util_Free (H->yLeftJump);
   util_Free (H->yRightJump);
}


/*=========================================================================*/

double fdist_KSMinusJumpsMany (fdist_FUNC_JUMPS * H, double dMoins)
{
   int k, j, i, jsup;
   double comb, temp, y;
   double *C, *B;
   const int M = H->par[0];
   const int NJ = H->NJumps;

   util_Assert (M <= 64, "fdist_KSMinusJumpsMany:   sample N too large");
   jsup = M * (1.0 - dMoins - EpsilonLR);
   util_Assert (jsup >= 0, "fdist_KSMinusJumpsMany:  jsup < 0");
   B = (double *) util_Calloc ((size_t) jsup + 1, sizeof (double));
   C = (double *) util_Calloc ((size_t) jsup + 1, sizeof (double));

   j = 0;
   while (j <= jsup) {
      y = dMoins + ((double) j) / M;
      i = 1;
      while (i <= NJ && y > H->yRightJump[i])
         ++i;
      /* I believe that Conover is wrong here, because this gives a */
      /* distribution that is continuous on the left, while probability */
      /* distributions must be continuous on the right. That could be */
      /* why these KS distributions with jumps don't seem to work.  */
      /* I may also have some bugs in these functions.  */

      if (i > NJ || y < H->yLeftJump[i])
         C[j] = 1.0 - y;
      else
         C[j] = 1.0 - H->yRightJump[i];
      ++j;
   }

   B[0] = 1.0;
   for (k = 1; k <= jsup; k++) {
      if (C[k] <= 0.0)
         B[k] = 0.0;
      else {
         temp = 0.0;
         comb = 1.0;
         for (j = 0; j < k; j++) {
            temp += comb * B[j] * pow (C[j], (double) (k - j));
            comb *= ((double) (k - j)) / (j + 1);
         }
         B[k] = 1.0 - temp;
      }
   }
   temp = 0.0;
   comb = 1.0;
   for (j = 0; j <= jsup; j++) {
      temp += comb * B[j] * pow (C[j], (double) (M - j));
      comb *= ((double) (M - j)) / (j + 1);
   }
   util_Warning (temp > 1.0 || temp < 0.0,
      "fdist_KSMinusJumpsMany:   Probabilities outside [0, 1]");
   util_Free (C);
   util_Free (B);
   return temp;
}


/*=========================================================================*/

double fdist_KSPlusJumpsMany (fdist_FUNC_JUMPS * H, double dPlus)
{
   int k, j, i, jsup;
   double comb, temp, y;
   double *F, *E;
   const int M = H->par[0];
   const int NJ = H->NJumps;

   util_Assert (M <= 64, "fdist_KSPlusJumpsMany:   sample N too large");
   jsup = M * (1.0 - dPlus - EpsilonLR);
   util_Assert (jsup >= 0, "fdist_KSPlusJumpsMany:  jsup < 0");
   E = (double *) util_Calloc ((size_t) jsup + 1, sizeof (double));
   F = (double *) util_Calloc ((size_t) jsup + 1, sizeof (double));

   j = 0;
   while (j <= jsup) {
      y = (1.0 - dPlus) - ((double) j) / M;
      i = 1;
      while (i <= NJ && y >= H->yRightJump[i])
         ++i;

      if (i > NJ || y <= H->yLeftJump[i])
         F[j] = y;
      else
         F[j] = H->yLeftJump[i];
      ++j;
   }

   E[0] = 1.0;
   for (k = 1; k <= jsup; k++) {
      if (F[k] <= 0.0)
         E[k] = 0.0;
      else {
         temp = 0.0;
         comb = 1.0;
         for (j = 0; j < k; j++) {
            temp += comb * E[j] * pow (F[j], (double) (k - j));
            comb *= ((double) (k - j)) / (j + 1);
         }
         E[k] = 1.0 - temp;
      }
   }
   temp = 0.0;
   comb = 1.0;
   for (j = 0; j <= jsup; j++) {
      temp += comb * E[j] * pow (F[j], (double) (M - j));
      comb *= ((double) (M - j)) / (j + 1);
   }
   util_Warning (temp > 1.0 || temp < 0.0,
      "fdist_KSPlusJumpsMany:   Probabilities outside [0, 1]");
   util_Free (E);
   util_Free (F);
   return temp;
}
#endif

/*=========================================================================*/

double fdist_CramerMises (long N, double x)
{
   const double Epsilon = DBL_EPSILON;
   const int jmax = 10;
   int j;
   double Cor, Res, arg;
   double termX, termS, termJ;
   static const double A[10] = {
      1.0,
      1.11803398875,
      1.125,
      1.12673477358,
      1.1274116945,
      1.12774323743,
      1.1279296875,
      1.12804477649,
      1.12812074678,
      1.12817350091
   };

   util_Assert (N > 0, "fdist_CramerMises:   N <= 0");

   if (N == 1) {
      if (x <= 1.0 / 12.0)
         return 0.0;
      if (x >= 1.0 / 3.0)
         return 1.0;
      return 2.0 * sqrt (x - 1.0 / 12.0);
   }

   if (x <= 0.002 || x < 1.0 / (12.0*N))
      return 0.0;
   if (x > 3.95 || x >= N/3.0)
      return 1.0;

   termX = 0.0625 / x;            /* 1 / (16x) */
   Res = 0.0;
   j = 0;
   do {
      termJ = 4 * j + 1;
      arg = termJ * termJ * termX;
      termS = A[j] * exp (-arg) * num2_BesselK025 (arg);
      Res += termS;
      ++j;
   } while (!(termS < Epsilon || j > jmax));

   util_Warning (j > jmax, "fdist_CramerMises: iterations have not converged");
   Res /= num_Pi * sqrt (x);

   /* Empirical correction in 1/N */
   if (x < 0.0092)
      Cor = 0.0;
   else if (x < 0.03)
      Cor = -0.0121763 + x * (2.56672 - 132.571 * x);
   else if (x < 0.06)
      Cor = 0.108688 + x * (-7.14677 + 58.0662 * x);
   else if (x < 0.19)
      Cor = -0.0539444 + x * (-2.22024 + x * (25.0407 - 64.9233 * x));
   else if (x < 0.5)
      Cor = -0.251455 + x * (2.46087 + x * (-8.92836 + x * (14.0988 -
               x * (5.5204 + 4.61784 * x))));
   else if (x <= 1.1)
      Cor = 0.0782122 + x * (-0.519924 + x * (1.75148 +
            x * (-2.72035 + x * (1.94487 - 0.524911 * x))));
   else
      Cor = exp (-0.244889 - 4.26506 * x);

   Res += Cor / N;
   /* This empirical correction is not very precise, so ... */
   if (Res <= 1.0)
      return Res;
   else
      return 1.0;
}


/*=========================================================================*/

double fdist_WatsonU (long N, double x)
/*
 * Only the asymptotic form has been implemented. In the trivial case
 * N = 1, we simply return 0.5
 */
{
   const int JMAX = 10;
   const double xSepare = 0.15;
   int j;
   double v;
   double terme;
   double somme;

   if (x <= 0.0)
      return 0.0;
   if (x >= fdist_XBIG)
      return 1.0;

   if (N == 1)                    /* N = 1, degenerate case */
      return 0.5;

   if (x > xSepare)
      return 1.0 - fbar_WatsonU (N, x);

   /* this series converges rapidly for x <= 0.15 */
   v = exp (-(0.125 / x));
   somme = v;
   j = 2;
   do {
      terme = pow (v, (double) (2 * j - 1) * (2 * j - 1));
      somme += terme;
      ++j;
   } while (!(terme < v * DBL_EPSILON || j > JMAX));
   util_Warning (j > JMAX, "fdist_WatsonU:  sum2 has not converged");

   v = 2.0 * somme / sqrt (2.0 * num_Pi * x);
   if (v >= 1.0)
      return 1.0;
   return v;
}


/*=========================================================================*/

static double YWA[143];           /* Tables for a spline approximation */
static double MWA[143];           /* of the WatsonG distribution */
static double CoWA[143];          /* Empirical correction in 1/sqrt(n) */

static void WatsonGInit (void)
/*
 * Initialization procedure for fdist_WatsonG
 */
{
   int j;

   YWA[0] = 1.8121832847E-39;      YWA[1] = 2.0503176304E-32;
   YWA[2] = 4.6139577764E-27;      YWA[3] = 6.5869745929E-23;
   YWA[4] = 1.2765816107E-19;      YWA[5] = 5.6251923105E-17;
   YWA[6] = 8.0747150511E-15;      YWA[7] = 4.8819994144E-13;
   YWA[8] = 1.4996052497E-11;      YWA[9] = 2.6903519441E-10;
   YWA[10] = 3.1322929018E-9;      YWA[11] = 2.5659643046E-8;
   YWA[12] = 1.5749759318E-7;      YWA[13] = 7.6105096466E-7;
   YWA[14] = 3.0113293541E-6;      YWA[15] = 1.0070166837E-5;
   YWA[16] = 2.9199826692E-5;      YWA[17] = 7.4970409372E-5;
   YWA[18] = 1.7340586581E-4;      YWA[19] = 3.6654236297E-4;
   YWA[20] = 7.165864865E-4;       YWA[21] = 1.3087767385E-3;
   YWA[22] = 2.2522044209E-3;      YWA[23] = 3.6781862572E-3;
   YWA[24] = 5.7361958631E-3;      YWA[25] = 8.5877444706E-3;
   YWA[26] = 1.23988738E-2;        YWA[27] = 1.73320516E-2;      
   YWA[28] = 2.35382479E-2;        YWA[29] = 3.11498548E-2;      
   YWA[30] = 4.02749297E-2;        YWA[31] = 5.09930445E-2;      
   YWA[32] = 6.33528333E-2;        YWA[33] = 7.73711747E-2;      
   YWA[34] = 9.30338324E-2;        YWA[35] = 1.10297306E-1;      
   YWA[36] = 1.290916098E-1;       YWA[37] = 1.493236984E-1;     
   YWA[38] = 1.708812741E-1;       YWA[39] = 1.936367476E-1;     
   YWA[40] = 2.174511609E-1;       YWA[41] = 2.42177928E-1;      
   YWA[42] = 2.676662852E-1;       YWA[43] = 2.937643828E-1;     
   YWA[44] = 3.203219784E-1;       YWA[45] = 3.471927188E-1;     
   YWA[46] = 3.742360163E-1;       YWA[47] = 4.013185392E-1;     
   YWA[48] = 4.283153467E-1;       YWA[49] = 4.551107027E-1;     
   YWA[50] = 4.815986082E-1;       YWA[51] = 5.076830902E-1;     
   YWA[52] = 5.332782852E-1;       YWA[53] = 5.583083531E-1;     
   YWA[54] = 5.827072528E-1;       YWA[55] = 6.064184099E-1;     
   YWA[56] = 6.293943006E-1;       YWA[57] = 6.515959739E-1;     
   YWA[58] = 6.729925313E-1;       YWA[59] = 6.935605784E-1;     
   YWA[60] = 7.132836621E-1;       YWA[61] = 7.321517033E-1;     
   YWA[62] = 7.501604333E-1;       YWA[63] = 7.673108406E-1;     
   YWA[64] = 7.836086337E-1;       YWA[65] = 7.99063723E-1;      
   YWA[66] = 8.136897251E-1;       YWA[67] = 8.275034914E-1;     
   YWA[68] = 8.405246632E-1;       YWA[69] = 8.527752531E-1;     
   YWA[70] = 8.642792535E-1;       YWA[71] = 8.750622738E-1;     
   YWA[72] = 8.851512032E-1;       YWA[73] = 8.945739017E-1;     
   YWA[74] = 9.033589176E-1;       YWA[75] = 9.115352296E-1;     
   YWA[76] = 9.19132015E-1;        YWA[77] = 9.261784413E-1;     
   YWA[78] = 9.327034806E-1;       YWA[79] = 9.387357465E-1;     
   YWA[80] = 9.44303351E-1;        YWA[81] = 9.494337813E-1;     
   YWA[82] = 9.541537951E-1;       YWA[83] = 9.584893325E-1;     
   YWA[84] = 9.624654445E-1;       YWA[85] = 9.661062352E-1;     
   YWA[86] = 9.694348183E-1;       YWA[87] = 9.724732859E-1;     
   YWA[88] = 9.752426872E-1;       YWA[89] = 9.777630186E-1;     
   YWA[90] = 9.800532221E-1;       YWA[91] = 9.821311912E-1;     
   YWA[92] = 9.840137844E-1;       YWA[93] = 9.85716844E-1;      
   YWA[94] = 9.872552203E-1;       YWA[95] = 9.886428002E-1;     
   YWA[96] = 9.898925389E-1;       YWA[97] = 9.910164946E-1;     
   YWA[98] = 9.920258656E-1;       YWA[99] = 9.929310287E-1;     
   YWA[100] = 9.937415788E-1;      YWA[101] = 9.944663692E-1;    
   YWA[102] = 9.95113552E-1;       YWA[103] = 9.956906185E-1;    
   YWA[104] = 9.962044387E-1;      YWA[105] = 9.966613009E-1;    
   YWA[106] = 9.970669496E-1;      YWA[107] = 9.974266225E-1;    
   YWA[108] = 9.977450862E-1;      YWA[109] = 9.980266707E-1;
   YWA[110] = 9.982753021E-1;      YWA[111] = 9.984945338E-1;
   YWA[112] = 9.98687576E-1;       YWA[113] = 9.98857324E-1;
   YWA[114] = 9.990063842E-1;      YWA[115] = 9.991370993E-1;
   YWA[116] = 9.992515708E-1;      YWA[117] = 9.99351681E-1;
   YWA[118] = 9.994391129E-1;      YWA[119] = 9.995153688E-1;
   YWA[120] = 9.995817875E-1;      YWA[121] = 9.996395602E-1;
   YWA[122] = 9.996897446E-1;      YWA[123] = 9.997332791E-1;
   YWA[124] = 9.997709943E-1;      YWA[125] = 9.998036243E-1;
   YWA[126] = 9.998318172E-1;      YWA[127] = 9.998561438E-1;
   YWA[128] = 9.998771066E-1;      YWA[129] = 9.998951466E-1;
   YWA[130] = 9.999106508E-1;      YWA[131] = 9.99923958E-1;
   YWA[132] = 9.999353645E-1;      YWA[133] = 9.999451288E-1;
   YWA[134] = 9.999534765E-1;      YWA[135] = 9.999606035E-1;
   YWA[136] = 9.999666805E-1;      YWA[137] = 9.999718553E-1;
   YWA[138] = 9.999762562E-1;      YWA[139] = 9.999799939E-1;
   YWA[140] = 9.999831643E-1;      YWA[141] = 9.999858E-1;
   YWA[142] = 9.999883E-1;

   MWA[0] = 0.0;                MWA[1] = 6.909E-15;    
   MWA[2] = 2.763E-14;          MWA[3] = 1.036E-13;    
   MWA[4] = 3.792E-13;          MWA[5] = 4.773E-12;    
   MWA[6] = 4.59E-10;           MWA[7] = 2.649E-8;     
   MWA[8] = 7.353E-7;           MWA[9] = 1.14E-5;      
   MWA[10] = 1.102E-4;          MWA[11] = 7.276E-4;    
   MWA[12] = 3.538E-3;          MWA[13] = 0.01342;     
   MWA[14] = 0.04157;           MWA[15] = 0.1088;      
   MWA[16] = 0.2474;            MWA[17] = 0.4999;      
   MWA[18] = 0.913;             MWA[19] = 1.53;        
   MWA[20] = 2.381;             MWA[21] = 3.475;       
   MWA[22] = 4.795;             MWA[23] = 6.3;         
   MWA[24] = 7.928;             MWA[25] = 9.602;       
   MWA[26] = 11.24;             MWA[27] = 12.76;       
   MWA[28] = 14.1;              MWA[29] = 15.18;       
   MWA[30] = 15.98;             MWA[31] = 16.47;       
   MWA[32] = 16.64;             MWA[33] = 16.49;       
   MWA[34] = 16.05;             MWA[35] = 15.35;       
   MWA[36] = 14.41;             MWA[37] = 13.28;       
   MWA[38] = 12.0;              MWA[39] = 10.6;        
   MWA[40] = 9.13;              MWA[41] = 7.618;       
   MWA[42] = 6.095;             MWA[43] = 4.588;       
   MWA[44] = 3.122;             MWA[45] = 1.713;       
   MWA[46] = 0.3782;            MWA[47] = -0.8726;     
   MWA[48] = -2.031;            MWA[49] = -3.091;      
   MWA[50] = -4.051;            MWA[51] = -4.91;       
   MWA[52] = -5.668;            MWA[53] = -6.327;      
   MWA[54] = -6.893;            MWA[55] = -7.367;      
   MWA[56] = -7.756;            MWA[57] = -8.064;      
   MWA[58] = -8.297;            MWA[59] = -8.46;       
   MWA[60] = -8.56;             MWA[61] = -8.602;      
   MWA[62] = -8.591;            MWA[63] = -8.533;      
   MWA[64] = -8.433;            MWA[65] = -8.296;      
   MWA[66] = -8.127;            MWA[67] = -7.93;       
   MWA[68] = -7.709;            MWA[69] = -7.469;      
   MWA[70] = -7.212;            MWA[71] = -6.943;      
   MWA[72] = -6.663;            MWA[73] = -6.378;      
   MWA[74] = -6.087;            MWA[75] = -5.795;      
   MWA[76] = -5.503;            MWA[77] = -5.213;      
   MWA[78] = -4.927;            MWA[79] = -4.646;      
   MWA[80] = -4.371;            MWA[81] = -4.103;      
   MWA[82] = -3.843;            MWA[83] = -3.593;      
   MWA[84] = -3.352;            MWA[85] = -3.12;       
   MWA[86] = -2.899;            MWA[87] = -2.689;      
   MWA[88] = -2.489;            MWA[89] = -2.3;        
   MWA[90] = -2.121;            MWA[91] = -1.952;      
   MWA[92] = -1.794;            MWA[93] = -1.645;      
   MWA[94] = -1.506;            MWA[95] = -1.377;      
   MWA[96] = -1.256;            MWA[97] = -1.144;      
   MWA[98] = -1.041;            MWA[99] = -0.9449;     
   MWA[100] = -0.8564;          MWA[101] = -0.775;   
   MWA[102] = -0.7001;          MWA[103] = -0.6315;  
   MWA[104] = -0.5687;          MWA[105] = -0.5113;  
   MWA[106] = -0.459;           MWA[107] = -0.4114;  
   MWA[108] = -0.3681;          MWA[109] = -0.3289;  
   MWA[110] = -0.2934;          MWA[111] = -0.2614;  
   MWA[112] = -0.2325;          MWA[113] = -0.2064;  
   MWA[114] = -0.183;           MWA[115] = -0.1621;  
   MWA[116] = -0.1433;          MWA[117] = -0.1265;  
   MWA[118] = -0.1115;          MWA[119] = -9.813E-2;
   MWA[120] = -8.624E-2;        MWA[121] = -7.569E-2;
   MWA[122] = -6.632E-2;        MWA[123] = -5.803E-2;
   MWA[124] = -5.071E-2;        MWA[125] = -4.424E-2;
   MWA[126] = -3.855E-2;        MWA[127] = -3.353E-2;
   MWA[128] = -2.914E-2;        MWA[129] = -2.528E-2;
   MWA[130] = -0.0219;          MWA[131] = -1.894E-2;
   MWA[132] = -1.637E-2;        MWA[133] = -1.412E-2;
   MWA[134] = -1.217E-2;        MWA[135] = -1.046E-2;
   MWA[136] = -8.988E-3;        MWA[137] = -7.72E-3;
   MWA[138] = -6.567E-3;        MWA[139] = -5.802E-3;
   MWA[140] = -0.0053;          MWA[141] = -4.7E-4;
   MWA[142] = -4.3E-4;

   for (j = 0; j <= 11; j++)
      CoWA[j] = 0.0;

   CoWA[12] = 1.25E-5;            CoWA[13] = 3.87E-5;      
   CoWA[14] = 1.004E-4;           CoWA[15] = 2.703E-4;     
   CoWA[16] = 6.507E-4;           CoWA[17] = 1.3985E-3;    
   CoWA[18] = 2.8353E-3;          CoWA[19] = 5.1911E-3;    
   CoWA[20] = 8.9486E-3;          CoWA[21] = 1.41773E-2;   
   CoWA[22] = 2.16551E-2;         CoWA[23] = 3.1489E-2;    
   CoWA[24] = 4.34123E-2;         CoWA[25] = 5.78719E-2;   
   CoWA[26] = 7.46921E-2;         CoWA[27] = 9.45265E-2;   
   CoWA[28] = 1.165183E-1;        CoWA[29] = 1.406353E-1;  
   CoWA[30] = 1.662849E-1;        CoWA[31] = 1.929895E-1;  
   CoWA[32] = 2.189347E-1;        CoWA[33] = 2.457772E-1;  
   CoWA[34] = 2.704794E-1;        CoWA[35] = 2.947906E-1;  
   CoWA[36] = 3.169854E-1;        CoWA[37] = 3.377435E-1;  
   CoWA[38] = 3.573555E-1;        CoWA[39] = 3.751205E-1;  
   CoWA[40] = 3.906829E-1;        CoWA[41] = 4.039806E-1;  
   CoWA[42] = 4.142483E-1;        CoWA[43] = 4.22779E-1;   
   CoWA[44] = 4.288013E-1;        CoWA[45] = 4.330353E-1;  
   CoWA[46] = 4.34452E-1;         CoWA[47] = 4.338138E-1;  
   CoWA[48] = 4.31504E-1;         CoWA[49] = 4.272541E-1;  
   CoWA[50] = 4.220568E-1;        CoWA[51] = 4.158229E-1;  
   CoWA[52] = 4.083281E-1;        CoWA[53] = 3.981182E-1;  
   CoWA[54] = 3.871678E-1;        CoWA[55] = 3.755527E-1;  
   CoWA[56] = 3.628823E-1;        CoWA[57] = 3.520135E-1;  
   CoWA[58] = 3.400924E-1;        CoWA[59] = 3.280532E-1;  
   CoWA[60] = 3.139477E-1;        CoWA[61] = 2.997087E-1;  
   CoWA[62] = 2.849179E-1;        CoWA[63] = 2.710475E-1;  
   CoWA[64] = 2.576478E-1;        CoWA[65] = 2.449155E-1;  
   CoWA[66] = 2.317447E-1;        CoWA[67] = 2.193161E-1;  
   CoWA[68] = 2.072622E-1;        CoWA[69] = 1.956955E-1;  
   CoWA[70] = 1.846514E-1;        CoWA[71] = 1.734096E-1;  
   CoWA[72] = 1.622678E-1;        CoWA[73] = 1.520447E-1;  
   CoWA[74] = 1.416351E-1;        CoWA[75] = 1.32136E-1;   
   CoWA[76] = 1.231861E-1;        CoWA[77] = 1.150411E-1;  
   CoWA[78] = 1.071536E-1;        CoWA[79] = 9.9465E-2;    
   CoWA[80] = 9.22347E-2;         CoWA[81] = 8.54394E-2;   
   CoWA[82] = 7.87697E-2;         CoWA[83] = 7.23848E-2;   
   CoWA[84] = 6.6587E-2;          CoWA[85] = 6.15849E-2;        
   CoWA[86] = 5.6573E-2;          CoWA[87] = 5.17893E-2;   
   CoWA[88] = 4.70011E-2;         CoWA[89] = 4.2886E-2;
   CoWA[90] = 3.91224E-2;         CoWA[91] = 3.53163E-2;
   CoWA[92] = 3.20884E-2;         CoWA[93] = 2.92264E-2;
   CoWA[94] = 2.66058E-2;         CoWA[95] = 2.37352E-2;
   CoWA[96] = 2.14669E-2;         CoWA[97] = 1.94848E-2;        
   CoWA[98] = 1.75591E-2;         CoWA[99] = 1.58232E-2;        
   CoWA[100] = 1.40302E-2;        CoWA[101] = 1.24349E-2;       
   CoWA[102] = 1.11856E-2;        CoWA[103] = 9.9765E-3;        
   CoWA[104] = 8.9492E-3;         CoWA[105] = 8.0063E-3;        
   CoWA[106] = 7.1509E-3;         CoWA[107] = 6.3196E-3;        
   CoWA[108] = 5.6856E-3;         CoWA[109] = 5.0686E-3;        
   CoWA[110] = 4.5085E-3;         CoWA[111] = 3.9895E-3;        
   CoWA[112] = 3.4804E-3;         CoWA[113] = 3.0447E-3;        
   CoWA[114] = 2.7012E-3;         CoWA[115] = 2.2984E-3;        
   CoWA[116] = 2.0283E-3;         CoWA[117] = 1.7399E-3;        
   CoWA[118] = 1.5032E-3;         CoWA[119] = 1.3267E-3;        
   CoWA[120] = 1.1531E-3;         CoWA[121] = 9.92E-4;          
   CoWA[122] = 9.211E-4;          CoWA[123] = 8.296E-4;         
   CoWA[124] = 6.991E-4;          CoWA[125] = 5.84E-4;          
   CoWA[126] = 5.12E-4;           CoWA[127] = 4.314E-4;         
   CoWA[128] = 3.593E-4;          CoWA[129] = 3.014E-4;         
   CoWA[130] = 2.401E-4;          CoWA[131] = 2.004E-4;         
   CoWA[132] = 1.614E-4;          CoWA[133] = 1.257E-4;         
   CoWA[134] = 1.112E-4;          CoWA[135] = 9.22E-5;          
   CoWA[136] = 8.77E-5;           CoWA[137] = 6.22E-5;          
   CoWA[138] = 4.93E-5;           CoWA[139] = 3.92E-5;          
   CoWA[140] = 3.15E-5;           CoWA[141] = 1.03E-5;
   CoWA[142] = 9.6E-6;
}


/*-------------------------------------------------------------------------*/

double fdist_WatsonG (long n, double X)
/*
 * Approximation of the cumulative distribution function of the
 * fdist_WatsonG statistics by the cubic spline function.
 *   Y[.]  - tabular value of the statistic;
 *   M[.]  - tabular value of the first derivative;
 */
{
   static int WatsonFlag = 0;
   const double MinArg = 0.15;
   const double MaxArg = 1.5;
   const double MinTab = 0.1;
   const double Step = 0.01;
   int i, j;
   double Tj;
   double Ti;
   double R;
   double P;
   double H;
   double Res;

   util_Assert (n > 0, "fdist_WatsonG:   N <= 0");

   if (n == 1)                    /* n = 1, degenerate case */
      return 0.5;

   if (!WatsonFlag) {
      /* Initialization of the interpolation table */
      WatsonGInit ();
      WatsonFlag = 1;
   }

   if (X <= MinArg)
      return 0.0;
   if (X >= 10.0)
      return 1.0;
   if (X > MaxArg) {
      R = exp (19.0 - 20.0 * X);
      Res = 1.0 - R;
      /* Empirical Correction in 1/sqrt(n) */
      R = exp (13.34 - 15.26 * X) / sqrt ((double) n);
      Res += R;
      /* The correction in 1/sqrt(n) is not always precise */
      if (Res >= 1.0)
         return 1.0;
      else
         return Res;
   }

   /* Search of the correct slot in the interpolation table */
   i = (int) ((X - MinTab) / Step) + 1;
   Ti = MinTab + i * Step;
   Tj = Ti - Step;
   /* Approximation within the slot */
   j = i - 1;
   H = X - Tj;
   R = Ti - X;
   P = Step * Step / 6.0;
   Res = ((MWA[j] * R * R * R + MWA[i] * H * H * H) / 6.0) / Step;
   Res += ((YWA[j] - MWA[j] * P) * R + (YWA[i] - MWA[i] * P) * H) / Step;

   /* Empirical correction in 1/sqrt(n) */
   Res += (CoWA[i] * H + CoWA[j] * R) / (Step * sqrt ((double) n));

   if (Res >= 1.0)
      return 1.0;
   return Res;
}


/*=========================================================================*/
#define AD_X0 0.38629436111989062
#define AD_X1 37.816242111357

static double AD_N_1 (double x)
{
   /* The Anderson-Darling distribution for N = 1 */
   double term;
   if (x <= AD_X0)
      return 0.0;
   if (x >= AD_X1)
      return 1.0;
   if (x - AD_X0 >= 1.0e-3)
      term = 1.0 - 4.0 * exp (-x - 1.0);
   else {
      const double q = x - AD_X0;
      term = q*(1.0 - q*(0.5 - q/6.0));
   }
   return sqrt (term);
}

#undef AD_X0
#undef AD_X1
/*=========================================================================*/

double fdist_AndersonDarling (long N, double x)
{
   if (1 == N)
      return AD_N_1 (x);
   util_Assert (N > 0, "fdist_AndersonDarling:   N <= 0");

   if (x <= 0.0)
      return 0.0;
   if (x >= fdist_XBIG)
      return 1.0;

   if (x <= 0.2) {
      /* Sinclair and Spurr lower tail approximation (3.6) */
      double q;
      q = 1.784 + 0.9936*x + 0.03287/x - (2.018 + 0.2029/x)/sqrt (x);
      if (q < -18.0)
         return exp(q);
      q = 1.0 + exp(q);
      return 1.0 - 1.0 / q;
   }
   return 1.0 - fbar_AndersonDarling (N, x);
}


/*=========================================================================*/
/* The following code is part of Marsaglia's file ADinf.c.
   Very little has been changed to adapt it to ProbDist. The file was 
   downloaded from the site of the Journal of Statistical Software
      http://www.jstatsoft.org/v09/i02/
*/

/*--------------------------------------------------------------------------*/
        /* This is file ADinf.c */
/*
A procedure for evaluating the limiting distribution of the
             Anderson-Darling statistic A_n=
-n-(1/n)[ln(x_1(1-x_n)+3ln(x_2(1-x_{n-1})+5ln(x_3(1-x_{n-2})+...
   +(2n-1)ln(x_n(1-x_1))]
    where x_1<x_2<...<x_n is an ordered set of purported uniform 
  [0,1) variates.
The function is ADinf(z)=lim_{n->infty} Pr[A_n<z]. About 15 digit accuracy.
If you don't need that much accuracy, use the quick-and-easy adinf(z).
ADinf uses a two-term recursion for coefficients in series for which 
 initial values
require the complementary normal integral, included as cPhi(z).
Otherwise, use erfc() if your C compiler has one with adequate accuracy.
*/

static double ADf (double z, int j)
{                                 /* called by ADinf(); see article. */
   double t, f, fnew, a, b, c, r;
   int i;
   t = (4 * j + 1) * (4 * j + 1) * 1.23370055013617 / z;
   if (t > 150.)
      return 0.;
   a = 2.22144146907918 * exp (-t) / sqrt (t);
   /* initialization requires cPhi */
   /* if you have erfc(), replace 2*cPhi(sqrt(2*t)) with erfc(sqrt(t)) */
   b = 3.93740248643060 * 2. * fbar_Normal2 (sqrt (2 * t));

   r = z * .125;
   f = a + b * r;
   for (i = 1; i < 200; i++) {
      c = ((i - .5 - t) * b + t * a) / i;
      a = b;
      b = c;
      r *= z / (8 * i + 8);
      if (fabs (r) < 1e-40 || fabs (c) < 1.e-40)
         return f;
      fnew = f + c * r;
      if (f == fnew)
         return f;
      f = fnew;
   }
   return f;
}


static double ADinf (double z)
{
   int j;
   double ad, adnew, r;
   if (z < .01)
      return 0.;   /* avoids exponent limits; ADinf(.01)=.528e-52 */
   r = 1. / z;
   ad = r * ADf (z, 0);
   for (j = 1; j < 100; j++) {
      r *= (.5 - j) / j;
      adnew = ad + (4 * j + 1) * r * ADf (z, j);
      if (ad == adnew) {
         return ad;
      }
      ad = adnew;
   }
   return ad;

}


/*------------------------------------------------------------------------*/
/* The following code is part of Marsaglia's file AnDarl.c.
   Very little has been changed to adapt it to ProbDist. The file was 
   downloaded from the site of the Journal of Statistical Software
      http://www.jstatsoft.org/v09/i02/
--------------------------------------*/

/*
    Anderson-Darling test for uniformity.   Given an ordered set
              x_1<x_2<...<x_n
 of purported uniform [0,1) variates,  compute
          a = -n-(1/n)*[ln(x_1*z_1)+3*ln(x_2*z_2+...+(2*n-1)*ln(x_n*z_n)]
 where z_1=1-x_n, z_2=1-x_(n-1)...z_n=1-x_1, then find
  v=adinf(a) and return  p=v+errfix(v), which should be uniform in [0,1),
  that is, the p-value associated with the observed x_1<x_2<...<x_n.
*/

/* Short, practical version of full ADinf(z), z>0.   */
static double adinf (double z)
{
   if (z < 2.)
      return exp (-1.2337141 / z) / sqrt (z) * (2.00012 + (.247105 -
            (.0649821 - (.0347962 - (.011672 -
                     .00168691 * z) * z) * z) * z) * z);
   /* max |error| < .000002 for z<2, (p=.90816...) */
   return
      exp (-exp (1.0776 - (2.30695 - (.43424 - (.082433 - (.008056 -
                     .0003146 * z) * z) * z) * z) * z));
   /* max |error|<.0000008 for 4<z<infinity */
}

/*------------------------------------------------------------------------*/
/* The function AD(n,z) returns Prob(A_n<z) where
    A_n = -n-(1/n)*[ln(x_1*z_1)+3*ln(x_2*z_2+...+(2*n-1)*ln(x_n*z_n)]
          z_1=1-x_n, z_2=1-x_(n-1)...z_n=1-x_1, and
    x_1<x_2<...<x_n is an ordered set of iid uniform [0,1) variates.
*/

static double AD (int n, double z, int isFastADinf)
{
   double c, v, x;
   /* If isFastADinf is true, use the fast approximation adinf (z),
      if it is false, use the more exact ADinf (z) */
   if (isFastADinf)
      x = adinf (z);
   else
      x = ADinf (z);

   /* now x=adinf(z). Next, get v=errfix(n,x) and return x+v; */
   if (x > .8) {
      v = (-130.2137 + (745.2337 - (1705.091 - (1950.646 - (1116.360 -
                     255.7844 * x) * x) * x) * x) * x) / n;
      return x + v;
   }
   c = .01265 + .1757 / n;
   if (x < c) {
      v = x / c;
      v = sqrt (v) * (1. - v) * (49 * v - 102);
      return x + v * (.0037 / (n * n) + .00078 / n + .00006) / n;
   }
   v = (x - c) / (.8 - c);
   v = -.00022633 + (6.54034 - (14.6538 - (14.458 - (8.259 -
               1.91864 * v) * v) * v) * v) * v;
   return x + v * (.04213 + .01365 / n) / n;
}

/* You must give the ADtest(int n, double *x) routine a sorted array
       x[0]<=x[1]<=..<=x[n-1]
    that you are testing for uniformity.
   It will return the p-value associated
   with the Anderson-Darling test, using
    the above adinf() and errfix( ,   )
         Not well-suited for n<7,
     (accuracy could drop to 3 digits).
*/


/*=========================================================================*/
#if 0

double fdist_AndersonDarling2 (long N, double x)
{
   /* This version uses the more exact limiting distribution ADinf */
   if (1 == N)
      return AD_N_1 (x);
   return AD (N, x, 0);
}

#else

double fdist_AndersonDarling2 (long N, double x)
{
   /* This version uses the approximate limiting distribution adinf */
   if (1 == N)
      return AD_N_1 (x);
   return AD ((int)N, x, 1);
}

#endif
/*=========================================================================*/


/***************************************\
 *
 *      DISCRETE DISTRIBUTIONS
 *
\***************************************/




/*=========================================================================*/

double fdist_Geometric (double p, long s)
{
   util_Assert (p >= 0.0 && p <= 1.0, "fdist_Geometric:   p not in [0, 1]");
   if (s < 0)
      return 0.0;
   if (p >= 1.0)                  /* In fact, p == 1 */
      return 1.0;
   if (p <= 0.0)                  /* In fact, p == 0 */
      return 0.0;
   return 1.0 - pow (1.0 - p, (double) (s + 1));
}


/*=========================================================================*/

double fdist_Poisson1 (double lam, long s)
/*
 * On our machine, computing a value using fdist_Gamma is faster than the 
 * naive computation for lamlim > 150.0, slower for lamlim < 150.0
 */
{
   const double lamlim = 150.0;
   long i;
   double term, sum;

   util_Assert (lam >= 0.0, "fdist_Poisson1:   lambda < 0");
   if (lam == 0.0)
      return 1.0;
   if (s < 0)
      return 0.0;

   /* If lam > lamlim, we use the Chi2 distribution according to the exact
      relation, with 2s + 2 degrees of freedom fdist_Poisson (lam, s) = 1 -
      fdist_ChiSquare (2s + 2, 2*lam) which also equals 1 - fdist_Gamma (s +
      1, lam) */
   if (lam > lamlim)
      return fbar_Gamma (s + 1.0, 15, lam);

   /* Naive computation: sum all prob. from i = 0 to i = s */
   sum = term = exp (-lam);
   for (i = 1; i <= s; i++) {
      term *= lam / i;
      sum += term;
   }
   return sum;
}


/*=========================================================================*/

double fdist_Poisson2 (fmass_INFO W, long s)
{
   double lam;
   util_Assert (W != NULL, "fdist_Poisson2:   fmass_INFO is NULL pointer");
   lam = W->paramR[0];

   if (s < 0)
      return 0.0;
   if (lam == 0.0)
      return 1.0;

   /* For large lam, we use the Chi2 distribution according to the exact
      relation, with 2s + 2 degrees of freedom

      fdist_Poisson (lam, s) = 1 - fdist_ChiSquare (2s + 2, 2*lam)

      which equals also 1 - fdist_Gamma (s + 1, lam) */
   if (W->cdf == NULL)
      return fbar_Gamma (s + 1.0, 15, lam);

   if (s >= W->smax)
      return 1.0;

   if (s < W->smin) {
      /* Sum RMAX dominant terms to get a few decimals in the lower tail. One
         could also call fbar_Gamma (s + 1.0, 15, lam) */
      const long RMAX = 20;
      long i;
      double term = fmass_PoissonTerm1 (lam, s);
      double Sum = term;
      i = s;
      while (i > 0 && i >= s - RMAX) {
         term = term * i / lam;
         i--;
         Sum += term;
      }
      return Sum;
   }

   if (s <= W->smed)
      return W->cdf[s - W->smin];
   else
      /* We keep the complementary distribution in the upper part of cdf */
      return 1.0 - W->cdf[s + 1 - W->smin];
}


/*=========================================================================*/

double fdist_Binomial1 (long n, double p, long s)
{
   const int nlim1 = 10000;
   const double varlim = 100.0;
   double epsilon = fmass_Epsilon;
   double y, z, q = 1.0 - p;
   double sum, term, termmid;
   long i, mid;
   int flag = 0;

   util_Assert (p >= 0.0 && p <= 1.0, "fdist_Binomial1:   p not in [0, 1]");
   util_Assert (n >= 0, "fdist_Binomial1:   n < 0");

   if (0 == n)
      return 1.0;
   if (s < 0)
      return 0.0;
   if (s >= n)
      return 1.0;
   if (p <= 0.0)
      return 1.0;
   if (p >= 1.0)
      return 0.0;                 /* For any s < n */

   if (n < nlim1) {               /* Exact Binomial */
      /* Sum RMAX terms to get a few decimals in the lower tail */
      const long RMAX = 20;
      mid = (long) ((n + 1) * p);
      if (mid > s)
         mid = s;
      sum = term = termmid = fmass_BinomialTerm3 (n, p, mid);

      z = q / p;
      i = mid;
      while (term >= epsilon || i >= mid - RMAX) {
         term *= z * i / (n - i + 1);
         sum += term;
         i--;
         if (i == 0) break;
      }

      z = p / q;
      term = termmid;
      for (i = mid; i < s; i++) {
         term *= z * (n - i) / (i + 1);
         if (term < epsilon)
            break;
         sum += term;
      }
      return sum;

   } else {
      if ((p > 0.5) || ((p == 0.5) && (s > n / 2))) {
         /* use F(p, n, s) = 1 - F(q, n, n-s-1) */
         p = q;
         q = 1.0 - p;
         flag = 1;
         s = n - s - 1;
      }
      if (n * p * q > varlim) {   /* Normal approximation */
         /* Uses the Camp-Paulson approximation based on the F-distribution.
            Its maximum absolute error is smaller than 0.007 / sqrt (npq).
            Ref: W. Molenaar; Approximations to the Poisson, Binomial,....
            QA273.6 M64, p. 93 (1970) */
         term = pow ((s + 1) * q / ((n - s) * p), 1.0 / 3.0);
         y = term * (9 - 1.0 / (s + 1)) - 9 + 1.0 / (n - s);
         z = 3.0 * sqrt (term * term / (s + 1) + 1.0 / (n - s));
         y /= z;
         if (flag) {
            return fbar_Normal1 (y);
         } else {
            return fdist_Normal2 (y);
         }

      } else {                    /* Poisson approximation */
         /* Uses a Bol'shev approximation based on the Poisson distribution.
            Error is O(1/n^4) as n -> infinity. Ref: W. Molenaar;
            Approximations to the Poisson, Binomial,.... QA273.6 M64, p. 107,
            Table 6.2, Formule lambda_9 (1970). */
         y = (2 * n - s) * p / (2.0 - p);
         z = (2.0 * y * y - s * y - (double)s * s - 2 * s) / (6 * (2 * n -
               (double)s) * (2 * n - s));
         z = y / (1 - z);
         if (flag) {
            return fbar_Poisson1 (z, s - 1);
         } else {
            return fdist_Poisson1 (z, s);
         }
      }
   }
}


/*=========================================================================*/

double fdist_Binomial2 (fmass_INFO W, long s)
{
   double p;
   long n;

   util_Assert (W != NULL, "fdist_Binomial2: fmass_INFO is NULL pointer");
   n = W->paramI[0];
   p = W->paramR[0];
   util_Assert (p >= 0.0 && p <= 1.0, "fdist_Binomial2:   p not in [0, 1]");

   if (0 == n)
      return 1.0;
   if (s < 0)
      return 0.0;
   if (s >= n)
      return 1.0;
   if (p == 0.0)
      return 1.0;
   if (p == 1.0)
      return 0.0;

   if (W->cdf != NULL) {
      if (s >= W->smax)
         return 1.0;
      if (s < W->smin) {
         /* Sum RMAX terms to get a few decimals in the lower tail */
         const long RMAX = 20;
         long i;
         double term = fmass_BinomialTerm3 (n, p, s);
         double Sum = term;
         const double z = (1.0 - p) / p;
         i = s;
         while (i > 0 && i >= s - RMAX) {
            term *= z * i / (n - i + 1);
            i--;
            Sum += term;
         }
         return Sum;
      }
      if (s <= W->smed)
         return W->cdf[s - W->smin];
      else
         /* We keep the complementary distribution in the upper part of cdf */
         return 1.0 - W->cdf[s + 1 - W->smin];

   } else {
      return fdist_Binomial1 (n, p, s);
   }
}


/*=========================================================================*/

double fdist_NegaBin1 (long n, double p, long s)
{
   const double epsilon = fmass_Epsilon;
   const long lim1 = 100000;
   double sum, term, termmode;
   long i, mode;

   util_Assert (p >= 0.0 && p <= 1.0, "fdist_NegaBin1:   p not in [0, 1]");
   util_Assert (n > 0, "fdist_NegaBin1:   n < 1");

   if (s < 0)
      return 0.0;
   if (p >= 1.0)                  /* In fact, p == 1 */
      return 1.0;
   if (p <= 0.0)                  /* In fact, p == 0 */
      return 0.0;

   /* Compute the maximum term */
   mode = 1 + (long) ((n * (1.0 - p) - 1.0) / p);
   if (mode > s)
      mode = s;

   if (mode <= lim1) {
      sum = term = termmode = fmass_NegaBinTerm1 (n, p, mode);
      for (i = mode; i > 0; i--) {
         term *= i / ((1.0 - p) * (n + i - 1));
         if (term < epsilon)
            break;
         sum += term;
      }

      term = termmode;
      for (i = mode; i < s; i++) {
         term *= (1.0 - p) * (n + i) / (i + 1);
         if (term < epsilon)
            break;
         sum += term;
      }
      if (sum <= 1.0)
         return sum;
      else
         return 1.0;

   } else {
      return 1.0 - fdist_Binomial1 (s + n, p, n - 1);
   }
}


/*=========================================================================*/

double fdist_NegaBin2 (fmass_INFO W, long s)
{
   double p;
   long n;

   util_Assert (W != NULL, "fdist_NegaBin2: fmass_INFO is NULL pointer");
   n = W->paramI[0];
   p = W->paramR[0];
   util_Assert (p >= 0.0 && p <= 1.0, "fdist_NegaBin2:   p not in [0, 1]");

   if (s < 0)
      return 0.0;
   if (p >= 1.0)                  /* In fact, p == 1 */
      return 1.0;
   if (p <= 0.0)                  /* In fact, p == 0 */
      return 0.0;

   if (W->cdf != NULL) {
      if (s >= W->smax)
         return 1.0;
      if (s < W->smin)
         return fdist_NegaBin1 (n, p, s);
      if (s <= W->smed)
         return W->cdf[s - W->smin];
      else
         /* We keep the complementary distribution in the upper part of cdf */
         return 1.0 - W->cdf[s + 1 - W->smin];

   } else {
      return fdist_NegaBin1 (n, p, s);
   }
}


/*=========================================================================*/

double fdist_Scan (long N, double d, long m)
{
   return 1.0 - fbar_Scan (N, d, m);
}


/*=========================================================================*/
