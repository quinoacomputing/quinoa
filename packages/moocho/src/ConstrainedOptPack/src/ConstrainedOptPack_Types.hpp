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

#ifndef CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
#define CONSTRAINED_OPTIMIZATION_PACK_TYPES_H

#include "NLPInterfacePack_Types.hpp"
#include "NLPInterfacePack_NLP.hpp"

namespace ConstrainedOptPack {

#include "NLPInterfacePack_PublicTypes.ud"

/// Bounds type
enum EBounds { FREE, UPPER, LOWER, EQUALITY };

// concrete classes

class VariableBoundsTester;

// abstract classes

class MatrixSymAddDelUpdateableWithOpFactorized;
class MatrixIdentConcat;
class MeritFuncCalc1D;
class MeritFuncCalc;
class MeritFuncNLP;
class MeritFuncNLE;
class MeritFuncNLF;
class MeritFuncNLPDirecDeriv;
class MeritFuncPenaltyParam;
class MeritFuncPenaltyParams;
class DirectLineSearch_Strategy;

// concrete subclasses

class MeritFuncCalc1DQuadratic;
class MeritFuncCalcNLP;
class MeritFuncNLPL1;
class MeritFuncNLPModL1;
//class MeritFuncCalcNLE;
//class MeritFuncCalcNLF;
//class MatrixHessianSuperBasic;
//class MatrixHessianSuperBasicInitDiagonal;
//class MatrixSymPosDefInvCholFactor;
class MatrixSymPosDefLBFGS;
class MatrixSymAddDelBunchKaufman;
class MatrixSymHessianRelaxNonSing;
class MatrixIdentConcatStd;
class DirectLineSearchArmQuad_Strategy;
class DirectLineSearchArmQuad_StrategySetOptions;
class VarReductOrthogDenseStd_Strategy;

// decomposition classes

class DecompositionSystem;
class DecompositionSystemVarReduct;
class DecompositionSystemVarReductPerm;
class DecompositionSystemVarReductPermStd;
class DecompositionSystemVarReductImp;
class DecompositionSystemCoordinate;
class DecompositionSystemOrthogonal;
class DecompositionSystemTester;
class DecompositionSystemTesterSetOptions;

// Abstract QP solvers

class QPSolverRelaxed;
class QPSolverRelaxedTester;
class QPSolverRelaxedTesterSetOptions;

// Concrete QP solvers

//class QPSchur;
//class QPSolverRelaxedQPSchurRangeSpace;

}	// end namespace ConstrainedOptPack 

#endif // CONSTRAINED_OPTIMIZATION_PACK_TYPES_H
