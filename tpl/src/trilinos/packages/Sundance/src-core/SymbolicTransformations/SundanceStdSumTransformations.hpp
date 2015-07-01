/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_STDSUMTRANSFORMATION_H
#define SUNDANCE_STDSUMTRANSFORMATION_H

#include "SundanceDefs.hpp"
#include "SundanceSumTransformationSequence.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;




/**
 * Apply a standard set of transformations
 */
class StdSumTransformations : public SumTransformationSequence
{
public:
  StdSumTransformations();

  virtual ~StdSumTransformations(){;}
};

/** 
 * Rewrite a sum as a polynomial object, if possible.
 */
class IdentifyPolynomialSum : public SumTransformation
{
public:
  /** */
  IdentifyPolynomialSum() : SumTransformation() {;}

  /** */
  virtual ~IdentifyPolynomialSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Put terms in a standard order
 */
class ReorderSum : public SumTransformation
{
public:
  /** */
  ReorderSum() : SumTransformation() {;}

  /** */
  virtual ~ReorderSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};
    

/** 
 * Simplify sums involving unary minuses
 * \f[
 * x + (-y) \rightarrow x-y. 
 * \f]
 * \f[
 * (-x) + (-y) \rightarrow -(x+y). 
 * \f]
 * \f[
 * (-x) + y \rightarrow y-x. 
 * \f]
 */
class RemoveUnaryMinusFromSum : public SumTransformation
{
public:
  /** */
  RemoveUnaryMinusFromSum() : SumTransformation() {;}

  /** */
  virtual ~RemoveUnaryMinusFromSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};
    
/** 
 * Transform a sum by removing a zero term: 
 * \f[
 * x + 0 \rightarrow x. 
 * \f]
 */
class RemoveZeroFromSum : public SumTransformation
{
public:
  /** */
  RemoveZeroFromSum() : SumTransformation() {;}

  /** */
  virtual ~RemoveZeroFromSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Sum two constant exprs without transformation 
 */
class SumConstants : public SumTransformation
{
public:
  /** */
  SumConstants() : SumTransformation() {;}

  /** */
  virtual ~SumConstants(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Transform a sum by moving any constants to the left:
 * \f[
 * x + a \rightarrow a + x
 * \f]
 **/
class MoveConstantsToLeftOfSum : public SumTransformation
{
public:
  /** */
  MoveConstantsToLeftOfSum() : SumTransformation() {;}

  /** */
  virtual ~MoveConstantsToLeftOfSum(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};


/** 
 * Rearrange a sum whose right operand is a sum including a constant
 * such that constants are grouped on the left:
 * \f[
 * \alpha + s_1 (\beta + s_2 u) \rightarrow (\alpha + s_1 \beta) + s_1 s_2 u
 * \f]
 * \f[
 * u + s_1 (\alpha + s_2 v) \rightarrow s_1 \alpha + (u + s_1 s_2 v)
 * \f]
 **/
class RearrangeRightSumWithConstant : public SumTransformation
{
public:
  /** */
  RearrangeRightSumWithConstant() : SumTransformation() {;}

  /** */
  virtual ~RearrangeRightSumWithConstant(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Rearrange a sum whose left operand is a sum including a constant
 * such that constants are grouped on the left:
 * \f[
 * (\alpha + s_1 u) + s_2 \beta \rightarrow (\alpha + s_2 \beta) + s_1 u 
 * \f]
 * \f[
 * (\alpha + s_1 u) + s_2 v \rightarrow \alpha + (s_1 u + s_2 v)
 * \f]
 **/
class RearrangeLeftSumWithConstant : public SumTransformation
{
public:
  /** */
  RearrangeLeftSumWithConstant() : SumTransformation() {;}

  /** */
  virtual ~RearrangeLeftSumWithConstant(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};


/** 
 * Transform sum of integrals to SumOfIntegral objects
 **/
class SumIntegrals : public SumTransformation
{
public:
  /** */
  SumIntegrals() : SumTransformation() {;}

  /** */
  virtual ~SumIntegrals(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    int sign, RCP<ScalarExpr>& rtn) const ;
};
}


#endif
