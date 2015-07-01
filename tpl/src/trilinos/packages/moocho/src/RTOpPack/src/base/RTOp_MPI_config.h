/*
// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/* If the macro RTOp_USE_MPI is defined, then these */
/* declarations will be MPI compatible.  If not then */
/* dummy MPI declarations will be used. */
/* */

#ifndef RTOP_MPI_CONFIG_H
#define RTOP_MPI_CONFIG_H

#include <RTOp_ConfigDefs.hpp>  /* This C++ header also has a proper C mode */

#ifdef HAVE_MPI
#define RTOp_USE_MPI  /* This macro is used in several places so we must keep it */
#endif

#ifdef RTOp_USE_MPI
#include "mpi.h"       /* Use real MPI declarations */
#else
/* #warning causes errors on Atlantis
#warning "Compiling in support for dummy MPI, real MPI will not be available!" */
#include "RTOp_mpi.h"  /* Use dummy MPI declarations */
#endif

/*
#ifdef __cplusplus
extern "C" {
#endif
*/

/** \file RTOp_config.h Platform dependent configuration options for RTOp.
 *
 * These typedefs and macros can be adjusted to the specific requirements of the platform.
 * For example, <tt>long double</tt> could be used instead of \c double if greater
 * precsion is needed.
 *
 * Also included are a few macros for MPI interoperability.  Also included in this default
 * header file is is \c RTOp_mpi.h which contains dummy MPI declarations (which are defined
 * in \c RTOp_mpi.c) for a subset of the MPI functions that are correct for the serial
 * case.
 */
/*@{ */

typedef MPI_Datatype            RTOp_Datatype;   /*< Compatible with MPI_Datatype? */
#define RTOpMPI_VALUE_TYPE      MPI_DOUBLE       /*< (MPI only) Compatible with fortran DOUBLE PRECISION? */
#define RTOpMPI_INDEX_TYPE      MPI_INT          /*< (MPI only) Compatible with fortran INTEGER? */
#define RTOpMPI_CHAR_TYPE       MPI_CHAR         /*< (MPI only) Compatible with fortran CHARACTER? */

/* The maxinum number of characters in a name of a reduction/transformation operator class */
#define RTOp_MAX_REDUCT_TRANS_OP_CLASS_NAME  50

/*@} */

/*
#ifdef __cplusplus
}
#endif
*/

#endif /* RTOP_MPI_CONFIG_H */
