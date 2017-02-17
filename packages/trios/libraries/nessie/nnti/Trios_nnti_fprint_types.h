/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */
/**
 *   @file Trios_nnti_fprint_types.h
 *
 *   @brief Pretty print the different NSSI types.
 *
 */

#ifndef _TRIOS_NNTI_FPRINT_TYPES_H_
#define _TRIOS_NNTI_FPRINT_TYPES_H_

#include <stdio.h>
#include <ostream>
#include "Trios_nnti_xdr.h"


#if defined(__STDC__) || defined(__cplusplus)

extern const char* nnti_err_str(const int rc);

    extern void fprint_NNTI_remote_addr(
                    FILE *fp,
                    const char *name,
                    const char *prefix,
                    const NNTI_remote_addr_t *addr);
    extern void fprint_NNTI_remote_addr(
            std::ostream &out,
                    const char *name,
                    const char *prefix,
                    const NNTI_remote_addr_t *addr);

    extern void fprint_NNTI_peer(
                    FILE *fp,
                    const char *name,
                    const char *prefix,
                    const NNTI_peer_t *addr);
    extern void fprint_NNTI_peer(
            std::ostream &out,
                    const char *name,
                    const char *prefix,
                    const NNTI_peer_t *addr);

    extern void fprint_NNTI_buffer(
                    FILE *fp,
                    const char *name,
                    const char *prefix,
                    const NNTI_buffer_t *addr);
    extern void fprint_NNTI_buffer(
            std::ostream &out,
                    const char *name,
                    const char *prefix,
                    const NNTI_buffer_t *addr);

    extern void fprint_NNTI_status(
                    FILE *fp,
                    const char *name,
                    const char *prefix,
                    const NNTI_status_t *addr);
    extern void fprint_NNTI_status(
            std::ostream &out,
                    const char *name,
                    const char *prefix,
                    const NNTI_status_t *addr);


#else /* K&R C */


#endif /* K&R C */

#endif
