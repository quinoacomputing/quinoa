/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000-2012, Sandia National Laboratories.                    *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/

#ifndef __HA_OVIS_H
#define __HA_OVIS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#ifdef ZOLTAN_OVIS

#include "ovis.h"

/* Structure for storing OVIS parameters */

struct OVIS_parameters {
  int outputLevel;
  char hello[MAX_PARAM_STRING_LEN];
  char dll[MAX_PARAM_STRING_LEN];
  double minVersion;
};

extern int Zoltan_OVIS_Setup(ZZ *, struct OVIS_parameters *);
extern int Zoltan_OVIS_Set_Param(char *, char *);

#endif /* ZOLTAN_OVIS */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
