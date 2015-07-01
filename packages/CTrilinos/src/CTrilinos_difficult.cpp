
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


#include "CTrilinos_config.h"
#include "CTrilinos_difficult.hpp"


namespace CTrilinos {


#ifdef HAVE_CTRILINOS_IFPACK
Ifpack::EPrecType convert_to_difficult_enum( CT_EPrecType_E_t en )
{
    switch (en) {
    case CT_EPrecType_E_POINT_RELAXATION:
        return Ifpack::POINT_RELAXATION;
    case CT_EPrecType_E_POINT_RELAXATION_STAND_ALONE:
        return Ifpack::POINT_RELAXATION_STAND_ALONE;
    case CT_EPrecType_E_BLOCK_RELAXATION:
        return Ifpack::BLOCK_RELAXATION;
    case CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE:
        return Ifpack::BLOCK_RELAXATION_STAND_ALONE;
    case CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_ILU:
        return Ifpack::BLOCK_RELAXATION_STAND_ALONE_ILU;
#ifdef HAVE_CTRILINOS_AMESOS
    case CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_AMESOS:
        return Ifpack::BLOCK_RELAXATION_STAND_ALONE_AMESOS;
    case CT_EPrecType_E_BLOCK_RELAXATION_AMESOS:
        return Ifpack::BLOCK_RELAXATION_AMESOS;
    case CT_EPrecType_E_AMESOS:
        return Ifpack::AMESOS;
    case CT_EPrecType_E_AMESOS_STAND_ALONE:
        return Ifpack::AMESOS_STAND_ALONE;
#endif /* HAVE_CTRILINOS_AMESOS */
    case CT_EPrecType_E_IC:
        return Ifpack::IC;
    case CT_EPrecType_E_IC_STAND_ALONE:
        return Ifpack::IC_STAND_ALONE;
    case CT_EPrecType_E_ICT:
        return Ifpack::ICT;
    case CT_EPrecType_E_ICT_STAND_ALONE:
        return Ifpack::ICT_STAND_ALONE;
    case CT_EPrecType_E_ILU:
        return Ifpack::ILU;
    case CT_EPrecType_E_ILU_STAND_ALONE:
        return Ifpack::ILU_STAND_ALONE;
    case CT_EPrecType_E_ILUT:
        return Ifpack::ILUT;
    case CT_EPrecType_E_ILUT_STAND_ALONE:
        return Ifpack::ILUT_STAND_ALONE;
#ifdef HAVE_IFPACK_SPARSKIT
    case CT_EPrecType_E_SPARSKIT:
        return Ifpack::SPARSKIT;
#endif /* HAVE_IFPACK_SPARSKIT */
#ifdef HAVE_IFPACK_HIPS
    case CT_EPrecType_E_HIPS:
        return Ifpack::HIPS;
#endif /* HAVE_IFPACK_HIPS */
#ifdef HAVE_HYPRE
    case CT_EPrecType_E_HYPRE:
        return Ifpack::HYPRE;
#endif /* HAVE_HYPRE */
    case CT_EPrecType_E_CHEBYSHEV:
        return Ifpack::CHEBYSHEV;
    default:
        throw CTrilinosMiscException("Likely error in preprocessor macros for Ifpack::EPrecType conversion.\n");
        break;
    }
}
#endif /* HAVE_CTRILINOS_IFPACK */


#ifdef HAVE_CTRILINOS_IFPACK
CT_EPrecType_E_t convert_from_difficult_enum( Ifpack::EPrecType en )
{
    switch (en) {
    case Ifpack::POINT_RELAXATION:
        return CT_EPrecType_E_POINT_RELAXATION;
    case Ifpack::POINT_RELAXATION_STAND_ALONE:
        return CT_EPrecType_E_POINT_RELAXATION_STAND_ALONE;
    case Ifpack::BLOCK_RELAXATION:
        return CT_EPrecType_E_BLOCK_RELAXATION;
    case Ifpack::BLOCK_RELAXATION_STAND_ALONE:
        return CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE;
    case Ifpack::BLOCK_RELAXATION_STAND_ALONE_ILU:
        return CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_ILU;
#ifdef HAVE_CTRILINOS_AMESOS
    case Ifpack::BLOCK_RELAXATION_STAND_ALONE_AMESOS:
        return CT_EPrecType_E_BLOCK_RELAXATION_STAND_ALONE_AMESOS;
    case Ifpack::BLOCK_RELAXATION_AMESOS:
        return CT_EPrecType_E_BLOCK_RELAXATION_AMESOS;
    case Ifpack::AMESOS:
        return CT_EPrecType_E_AMESOS;
    case Ifpack::AMESOS_STAND_ALONE:
        return CT_EPrecType_E_AMESOS_STAND_ALONE;
#endif /* HAVE_CTRILINOS_AMESOS */
    case Ifpack::IC:
        return CT_EPrecType_E_IC;
    case Ifpack::IC_STAND_ALONE:
        return CT_EPrecType_E_IC_STAND_ALONE;
    case Ifpack::ICT:
        return CT_EPrecType_E_ICT;
    case Ifpack::ICT_STAND_ALONE:
        return CT_EPrecType_E_ICT_STAND_ALONE;
    case Ifpack::ILU:
        return CT_EPrecType_E_ILU;
    case Ifpack::ILU_STAND_ALONE:
        return CT_EPrecType_E_ILU_STAND_ALONE;
    case Ifpack::ILUT:
        return CT_EPrecType_E_ILUT;
    case Ifpack::ILUT_STAND_ALONE:
        return CT_EPrecType_E_ILUT_STAND_ALONE;
#ifdef HAVE_IFPACK_SPARSKIT
    case Ifpack::SPARSKIT:
        return CT_EPrecType_E_SPARSKIT;
#endif /* HAVE_IFPACK_SPARSKIT */
#ifdef HAVE_IFPACK_HIPS
    case Ifpack::HIPS:
        return CT_EPrecType_E_HIPS;
#endif /* HAVE_IFPACK_HIPS */
#ifdef HAVE_HYPRE
    case Ifpack::HYPRE:
        return CT_EPrecType_E_HYPRE;
#endif /* HAVE_HYPRE */
    case Ifpack::CHEBYSHEV:
        return CT_EPrecType_E_CHEBYSHEV;
    default:
        throw CTrilinosMiscException("Likely error in preprocessor macros for Ifpack::EPrecType conversion.\n");
        break;
    }
}
#endif /* HAVE_CTRILINOS_IFPACK */


} // namespace CTrilinos

