/*************************************************************************\
 *
 * Package:        MyLib
 * File:           num2.c
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

#include "num2.h"
#include "util.h"
#include "num.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>

#define EPSILON  1.0E-15
#define MAXI 50

/* The factorials n! from n = 0 to n = 170 */
static double Factorials[] = {
   1,
   1,
   2,
   6,
   24,
   120,
   720,
   5040,
   40320,
   362880,
   3628800,
   39916800,
   479001600,
   6227020800.0,
   87178291200.0,
   1307674368000.0,
   20922789888000.0,
   355687428096000.0,
   6402373705728000.0,
   1.21645100408832e+17,
   2.43290200817664e+18,
   5.109094217170944e+19,
   1.124000727777608e+21,
   2.585201673888498e+22,
   6.204484017332394e+23,
   1.551121004333099e+25,
   4.032914611266057e+26,
   1.088886945041835e+28,
   3.048883446117138e+29,
   8.841761993739701e+30,
   2.65252859812191e+32,
   8.222838654177922e+33,
   2.631308369336935e+35,
   8.683317618811886e+36,
   2.952327990396041e+38,
   1.033314796638614e+40,
   3.719933267899012e+41,
   1.376375309122634e+43,
   5.23022617466601e+44,
   2.039788208119744e+46,
   8.159152832478977e+47,
   3.34525266131638e+49,
   1.40500611775288e+51,
   6.041526306337383e+52,
   2.658271574788449e+54,
   1.196222208654802e+56,
   5.502622159812088e+57,
   2.586232415111682e+59,
   1.241391559253607e+61,
   6.082818640342675e+62,
   3.041409320171338e+64,
   1.551118753287382e+66,
   8.065817517094388e+67,
   4.274883284060025e+69,
   2.308436973392414e+71,
   1.269640335365828e+73,
   7.109985878048635e+74,
   4.052691950487722e+76,
   2.350561331282879e+78,
   1.386831185456899e+80,
   8.320987112741392e+81,
   5.075802138772248e+83,
   3.146997326038794e+85,
   1.98260831540444e+87,
   1.268869321858842e+89,
   8.247650592082472e+90,
   5.443449390774431e+92,
   3.647111091818868e+94,
   2.480035542436831e+96,
   1.711224524281413e+98,
   1.197857166996989e+100,
   8.504785885678622e+101,
   6.123445837688608e+103,
   4.470115461512683e+105,
   3.307885441519386e+107,
   2.480914081139539e+109,
   1.88549470166605e+111,
   1.451830920282858e+113,
   1.132428117820629e+115,
   8.946182130782973e+116,
   7.156945704626378e+118,
   5.797126020747366e+120,
   4.75364333701284e+122,
   3.945523969720657e+124,
   3.314240134565352e+126,
   2.817104114380549e+128,
   2.422709538367272e+130,
   2.107757298379527e+132,
   1.854826422573984e+134,
   1.650795516090845e+136,
   1.485715964481761e+138,
   1.352001527678402e+140,
   1.24384140546413e+142,
   1.156772507081641e+144,
   1.087366156656742e+146,
   1.032997848823905e+148,
   9.916779348709491e+149,
   9.619275968248206e+151,
   9.426890448883242e+153,
   9.33262154439441e+155,
   9.33262154439441e+157,
   9.425947759838354e+159,
   9.614466715035121e+161,
   9.902900716486175e+163,
   1.029901674514562e+166,
   1.08139675824029e+168,
   1.146280563734708e+170,
   1.226520203196137e+172,
   1.324641819451828e+174,
   1.443859583202493e+176,
   1.588245541522742e+178,
   1.762952551090244e+180,
   1.974506857221073e+182,
   2.231192748659812e+184,
   2.543559733472186e+186,
   2.925093693493014e+188,
   3.393108684451897e+190,
   3.969937160808719e+192,
   4.684525849754288e+194,
   5.574585761207603e+196,
   6.689502913449124e+198,
   8.09429852527344e+200,
   9.875044200833598e+202,
   1.214630436702532e+205,
   1.50614174151114e+207,
   1.882677176888925e+209,
   2.372173242880046e+211,
   3.012660018457658e+213,
   3.856204823625803e+215,
   4.974504222477285e+217,
   6.466855489220472e+219,
   8.471580690878817e+221,
   1.118248651196004e+224,
   1.487270706090685e+226,
   1.992942746161518e+228,
   2.69047270731805e+230,
   3.659042881952547e+232,
   5.01288874827499e+234,
   6.917786472619486e+236,
   9.615723196941086e+238,
   1.346201247571752e+241,
   1.89814375907617e+243,
   2.695364137888161e+245,
   3.854370717180071e+247,
   5.550293832739301e+249,
   8.047926057471987e+251,
   1.17499720439091e+254,
   1.727245890454638e+256,
   2.556323917872864e+258,
   3.808922637630567e+260,
   5.71338395644585e+262,
   8.627209774233235e+264,
   1.311335885683452e+267,
   2.006343905095681e+269,
   3.089769613847349e+271,
   4.789142901463391e+273,
   7.471062926282891e+275,
   1.172956879426414e+278,
   1.853271869493734e+280,
   2.946702272495037e+282,
   4.714723635992059e+284,
   7.590705053947215e+286,
   1.229694218739449e+289,
   2.004401576545302e+291,
   3.287218585534294e+293,
   5.423910666131586e+295,
   9.003691705778433e+297,
   1.503616514864998e+300,
   2.526075744973197e+302,
   4.269068009004703e+304,
   7.257415615307994e+306
};


#define MLIM 50

/* The natural logarithm of factorials n! from n = 0 to 50 */
static double LnFactorials[MLIM + 1] = {
   0.,
   0.,
   0.6931471805599453,
   1.791759469228055,
   3.178053830347946,
   4.787491742782046,
   6.579251212010101,
   8.525161361065415,
   10.60460290274525,
   12.80182748008147,
   15.10441257307552,
   17.50230784587389,
   19.98721449566188,
   22.55216385312342,
   25.19122118273868,
   27.89927138384088,
   30.67186010608066,
   33.50507345013688,
   36.39544520803305,
   39.33988418719949,
   42.33561646075348,
   45.3801388984769,
   48.47118135183522,
   51.60667556776437,
   54.7847293981123,
   58.00360522298051,
   61.26170176100199,
   64.55753862700632,
   67.88974313718154,
   71.257038967168,
   74.65823634883016,
   78.09222355331529,
   81.55795945611503,
   85.05446701758153,
   88.58082754219767,
   92.13617560368708,
   95.7196945421432,
   99.3306124547874,
   102.9681986145138,
   106.6317602606434,
   110.3206397147574,
   114.0342117814617,
   117.7718813997451,
   121.5330815154386,
   125.3172711493569,
   129.1239336391272,
   132.9525750356163,
   136.8027226373264,
   140.6739236482343,
   144.5657439463449,
   148.477766951773
};


/*=========================================================================*/

double num2_Factorial (int n)
{
   util_Assert (n >= 0, "num2_Factorial:   n < 0");
   if (n <= 170)
      return Factorials[n];
   util_Warning (1, "num2_Factorial:   n > 170:   return inf");
   return 1.0 / 0.0;
}


/*=========================================================================*/

double num2_LnFactorial (int n)
{
   util_Assert (n >= 0, "num2_LnFactorial:   n < 0");
   if (n <= MLIM) {
      return LnFactorials[n];

   } else {
      double x = (double) (n + 1);
      double y = 1.0 / (x * x);
      double z = ((-(5.95238095238E-4 * y) + 7.936500793651E-4) * y -
         2.7777777777778E-3) * y + 8.3333333333333E-2;
      z = ((x - 0.5) * log (x) - x) + 9.1893853320467E-1 + z / x;
      return z;
   }
}


/*=========================================================================*/
#ifndef HAVE_LGAMMA

/* The new standard ISO_C99 includes the lgamma function in math.h;
   otherwise, we shall have to use our own. */

double num2_LnGamma (double x)
{
   const double xlimbig = 1.0 / DBL_EPSILON;
   const double xlim1 = 18.0;
   const double dk2 = 0.91893853320467274178; /* Ln (sqrt (2 Pi)) */
   const double dk1 = 0.9574186990510627;
   const int N = 15;              /* Degree of Chebyshev polynomial */
   double y = 0, z = 0;
   int i, k;

   /* Chebyshev coefficients for lnGamma (x + 3), 0 <= x <= 1 In Yudell Luke:
      The special functions and their approximations, Vol. II, Academic Press,
      p. 301, 1969. There is an error in the additive constant in the formula:
      (Ln (2)). */
   static const double A[] = {
      0.52854303698223459887,
      0.54987644612141411418,
      0.02073980061613665136,
      -0.00056916770421543842,
      0.00002324587210400169,
      -0.00000113060758570393,
      0.00000006065653098948,
      -0.00000000346284357770,
      0.00000000020624998806,
      -0.00000000001266351116,
      0.00000000000079531007,
      -0.00000000000005082077,
      0.00000000000000329187,
      -0.00000000000000021556,
      0.00000000000000001424,
      -0.00000000000000000095
   };

   util_Assert (x > 0.0, "num2_LnGamma:   accepts only x > 0");
   if (x > xlim1) {
      if (x > xlimbig)
         y = 0.0;
      else
         y = 1.0 / (x * x);
      z = ((-(5.95238095238E-4 * y) + 7.936500793651E-4) * y -
         2.7777777777778E-3) * y + 8.3333333333333E-2;
      z = ((x - 0.5) * log (x) - x) + dk2 + z / x;
      return z;

   } else if (x > 4.0) {
      k = (int) x;
      z = x - k;
      y = 1.0;
      for (i = 3; i < k; i++)
         y *= z + i;
      y = log (y);

   } else if (x <= 0.0) {
      return DBL_MAX;

   } else if (x < 3.0) {
      k = (int) x;
      z = x - k;
      y = 1.0;
      for (i = 2; i >= k; i--)
         y *= z + i;
      y = -log (y);

   } else {                       /* 3 <= x <= 4 */
      z = x - 3.0;
      y = 0.0;
   }

   z = num2_EvalCheby (A, N, 2.0 * z - 1.0);
   return z + dk1 + y;
}

#endif
/*=========================================================================*/

#define NLIM 100                  /* pour eviter les debordements */

double num2_Combination (int n, int s)
{
   double Res;
   int i;
   int Diff;
   if (s == 0 || s == n)
      return 1.0;
   if (s < 0) {
      util_Warning (1, "num2_Combination:   s < 0");
      return 0.0;
   }
   if (s > n) {
      util_Warning (1, "num2_Combination:   s > n");
      return 0.0;
   }
   if (s > (n / 2))
      s = n - s;
   if (n <= NLIM) {
      Res = 1.0;
      Diff = n - s;
      for (i = 1; i <= s; i++) {
         Res = (Res * (double) (Diff + i)) / (double) (i);
      }
      return Res;
   } else {
      Res = (num2_LnFactorial (n) - num2_LnFactorial (s))
         - num2_LnFactorial (n - s);
      return exp (Res);
   }
}


/*=========================================================================*/
#ifndef HAVE_LOG1P

double num2_log1p (double x)
{
   /* returns a value equivalent to log (1 + x) accurate also for small x. */
   if (fabs (x) > 0.1) {
      return log (1.0 + x);
   } else {
      double term = x;
      double sum = x;
      int s = 2;
      while (fabs (term) > EPSILON * fabs (sum) && s < MAXI) {
         term *= -x;
         sum += term / s;
         s++;
      }
      return sum;
   }
}

#endif
/*=========================================================================*/

void num2_CalcMatStirling (double ***M, int m, int n)
/* Calcul des elements de la matrice MatStirling [0..m, 0..n]. */
{
   int i, j, k;
   *M = (double **) util_Calloc ((size_t) (m + 1), sizeof (double *));
   for (i = 0; i <= m; i++)
      (*M)[i] = (double *) util_Calloc ((size_t) (n + 1), sizeof (double));

   for (i = 0; i <= m; i++) {
      for (j = 0; j <= n; j++) {
         (*M)[i][j] = 0.0;
      }
   }

   (*M)[0][0] = 1.0;
   for (j = 1; j <= n; j++) {
      (*M)[0][j] = 0.0;
      if (j <= m) {
         k = j - 1;
         (*M)[j][j] = 1.0;
      } else
         k = m;
      for (i = 1; i <= k; i++) {
         (*M)[i][j] = (double) (i) * (*M)[i][j - 1] + (*M)[i - 1][j - 1];
      }
   }
}


/*=========================================================================*/

void num2_FreeMatStirling (double ***M, int m)
{
   int i;
   for (i = 0; i <= m; i++)
      free ((*M)[i]);
   free (*M);
   *M = NULL;
}


/*=========================================================================*/

double num2_VolumeSphere (double pLR, int k)
/* Returns volume of unit sphere in dimension k, norm p */
{
   const double eps = 2.0 * DBL_EPSILON;
   int p = pLR;
   double kLR = (double) k;
   double Vol;
   int s;

   util_Assert (pLR >= 0.0, "num2_VolumeSphere:   p < 0");
   if (fabs (pLR - p) <= eps) {
      switch (p) {
      case 0:
         return num_TwoExp[k];
         break;
      case 1:
         return num_TwoExp[k] / num2_Factorial (k);
         break;
      case 2:
         if ((k % 2) == 0) {
            return pow (num_Pi, kLR / 2.0) / num2_Factorial (k / 2);
         } else {
            s = (k + 1) / 2;
            return pow (num_Pi, (double) (s) - 1.0) * num2_Factorial (s) *
               num_TwoExp[2 * s] / num2_Factorial (2 * s);
         }
         break;
      default:
         break;
      }
   }
   Vol = kLR * (num_Ln2 + num2_LnGamma (1.0 + 1.0 / pLR)) -
      num2_LnGamma (1.0 + kLR / pLR);
   return exp (Vol);
}


/*=========================================================================*/

double num2_EvalCheby (const double A[], int N, double x)
{
   int j;
   double xx;
   double b0, b1, b2;
   util_Warning (fabs (x) > 1.0,
      "Chebychev polynomial evaluated at x outside [-1, 1]");
   xx = 2.0 * x;
   b0 = 0.0;
   b1 = 0.0;
   for (j = N; j >= 0; j--) {
      b2 = b1;
      b1 = b0;
      b0 = (xx * b1 - b2) + A[j];
   }
   return (b0 - b2) / 2.0;
}


/*=========================================================================*/
#define DEGREE 6

double num2_BesselK025 (double x)
{
   double rac;
   double xx;
   double temp;
   double Res;
   double C;
   double B;
   int j;
   static const double c[8] = {
      32177591145.0,
      2099336339520.0,
      16281990144000.0,
      34611957596160.0,
      26640289628160.0,
      7901666082816.0,
      755914244096.0
   };

   static const double b[8] = {
      75293843625.0,
      2891283595200.0,
      18691126272000.0,
      36807140966400.0,
      27348959232000.0,
      7972533043200.0,
      755914244096.0
   };

   if (x < 1.E-300)
      return DBL_MAX;

   /*------------------------------------------------------------------*/
   /* x > 0.6 => approximation asymptotique rationnelle dans Luke: */
   /* Yudell L.Luke "Mathematical functions and their approximations", */
   /* Academic Press Inc. New York, 1975, p.371 */
   /*------------------------------------------------------------------*/
   if (x >= 0.6) {
      B = b[DEGREE];
      C = c[DEGREE];
      for (j = DEGREE; j >= 1; j--) {
         B = B * x + b[j - 1];
         C = C * x + c[j - 1];
      }
      Res = sqrt (num_Pi / (2.0 * x)) * exp (-x) * (C / B);
      return Res;
   }

   /*------------------------------------------------------------------*/
   /* x < 0.6 => la serie de K_{1/4} = Pi/Sqrt(2) [I_{-1/4} - I_{1/4}] */
   /*------------------------------------------------------------------*/
   xx = x * x;
   rac = pow (x / 2.0, 0.25);
   Res = (((xx / 1386.0 + 1.0 / 42.0) * xx + 1.0 / 3.0) * xx + 1.0) /
      (1.225416702465177 * rac);
   temp = (((xx / 3510.0 + 1.0 / 90.0) * xx + 0.2) * xx + 1.0) * rac /
      0.906402477055477;
   Res = num_Pi * (Res - temp) / num_Rac2;
   return Res;
}

#undef DEGREE
/*=========================================================================*/

double num2_Digamma (double x)
{
   static const double C7[] = {
      1.3524999667726346383e4, 4.5285601699547289655e4,
      4.5135168469736662555e4, 1.8529011818582610168e4,
      3.3291525149406935532e3, 2.4068032474357201831e2,
      5.1577892000139084710, 6.2283506918984745826e-3
   };

   static const double D7[] = {
      6.9389111753763444376e-7, 1.9768574263046736421e4,
      4.1255160835353832333e4, 2.9390287119932681918e4,
      9.0819666074855170271e3, 1.2447477785670856039e3,
      6.7429129516378593773e1, 1.0
   };

   static const double C4[] = {
      -2.728175751315296783e-15, -6.481571237661965099e-1,
      -4.486165439180193579, -7.016772277667586642, -2.129404451310105168
   };

   static const double D4[] = {
      7.777885485229616042, 5.461177381032150702e1,
      8.929207004818613702e1, 3.227034937911433614e1, 1.0
   };

   double prodPj = 0.0;
   double prodQj = 0.0;
   double digX = 0.0;

   if (x >= 3.0) {
      double x2 = 1.0 / (x * x);
      int j;
      for (j = 4; j >= 0; j--) {
         prodPj = prodPj * x2 + C4[j];
         prodQj = prodQj * x2 + D4[j];
      }
      digX = log (x) - (0.5 / x) + (prodPj / prodQj);

   } else if (x >= 0.5) {
      const double X0 = 1.46163214496836234126;
      int j;
      for (j = 7; j >= 0; j--) {
         prodPj = x * prodPj + C7[j];
         prodQj = x * prodQj + D7[j];
      }
      digX = (x - X0) * (prodPj / prodQj);

   } else {
      double f = (1.0 - x) - floor (1.0 - x);
      digX = num2_Digamma (1.0 - x) + num_Pi / tan (num_Pi * f);
   }

   return digX;
}

