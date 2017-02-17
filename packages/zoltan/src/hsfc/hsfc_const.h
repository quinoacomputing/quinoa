/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef ZOLTAN_HSFC_CONST_H
#define ZOLTAN_HSFC_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* function prototypes */

int  Zoltan_HSFC_Set_Param (char *name, char *val) ;
void Zoltan_HSFC_Print_Structure(ZZ *zz);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
