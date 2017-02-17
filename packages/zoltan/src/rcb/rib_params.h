/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/


#ifndef __RIB_PARAMS_H
#define __RIB_PARAMS_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "params_const.h"

/*  Parameters structure for RIB method.  Used in  */
/*  Zoltan_RIB_Set_Param and Zoltan_RIB.                      */
static PARAM_VARS RIB_params[] = {
               { "RIB_OVERALLOC", NULL, "DOUBLE", 0 },
               { "CHECK_GEOM", NULL, "INT", 0 },
               { "RIB_OUTPUT_LEVEL", NULL, "INT", 0 },
               { "AVERAGE_CUTS", NULL, "INT", 0 },
               { "KEEP_CUTS", NULL, "INT", 0 },
               { "REDUCE_DIMENSIONS", NULL, "INT", 0 },
               { "DEGENERATE_RATIO", NULL, "DOUBLE", 0 },
               { "FINAL_OUTPUT", NULL,  "INT",    0},
               { NULL, NULL, NULL, 0 } };


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
