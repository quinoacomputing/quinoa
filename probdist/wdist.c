/*************************************************************************\
 *
 * Package:        ProbDist
 * File:           wdist.c
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

#include "wdist.h"
#include "fdist.h"



/*=========================================================================*/

double wdist_Normal (double Junk[], double x)
{
   return fdist_Normal2 (x);
}


double wdist_ChiSquare (double W[], double x)
{
   long N = (long) W[0];
   return fdist_ChiSquare2 (N, 12, x);
}


double wdist_Unif (double Junk[], double x)
{
   return fdist_Unif (x);
}

