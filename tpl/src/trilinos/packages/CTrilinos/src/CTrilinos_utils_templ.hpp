
/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */


/*! @file CTrilinos_utils_templ.hpp
 * @brief Templated utility functions for CTrilinos. */


#ifndef CTRILINOS_UTILS_TEMPL_HPP
#define CTRILINOS_UTILS_TEMPL_HPP


#include "CTrilinos_config.h"


#include <string>
#include "CTrilinos_enums.h"


namespace CTrilinos {


/* convert struct from specific type to generic CTrilinos_Universal_ID_t
 * but keep the content in tact */
template <typename T>
CTrilinos_Universal_ID_t
abstractType( T id )
{
    CTrilinos_Universal_ID_t newid;

    newid.table = id.table;
    newid.index = id.index;
    newid.is_const = id.is_const;

    return newid;
}

/* convert struct from generic CTrilinos_Universal_ID_t to specific type
 * but keep the content in tact */
template <typename T>
T
concreteType( CTrilinos_Universal_ID_t id )
{
    T newid;

    newid.table = id.table;
    newid.index = id.index;
    newid.is_const = id.is_const;

    return newid;
}


} // namespace CTrilinos


#endif

