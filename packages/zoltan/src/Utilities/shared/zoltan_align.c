/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zoltan_align.h"
#include "zoltan_util.h"

/*****************************************************************************/
/*
 *  Routines for properly aligning data.
 */
/*****************************************************************************/

/* 
 * Plauger alignment algorithm, The Standard C Library.
 * Forces malloc'ed variable size struct alignment.
 * ZOLTAN_ALIGN_VAL is defined in Zoltan/include/zoltan_align.h;
 * values are 0,1,3,7U depending upon machine.
 */

int Zoltan_Align(int a)
{
return((ZOLTAN_ALIGN_VAL + a) & ~ZOLTAN_ALIGN_VAL);
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
