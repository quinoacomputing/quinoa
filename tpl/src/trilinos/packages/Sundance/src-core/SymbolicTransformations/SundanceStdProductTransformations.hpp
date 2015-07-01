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

#ifndef SUNDANCE_STDPRODUCTTRANSFORMATION_H
#define SUNDANCE_STDPRODUCTTRANSFORMATION_H

#include "SundanceDefs.hpp"
#include "SundanceProductTransformationSequence.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;




/**
 * Apply a standard set of transformations
 */
class StdProductTransformations : public ProductTransformationSequence
{
public:
  StdProductTransformations();

  virtual ~StdProductTransformations(){;}
};
    
/** 
 * Transform a product by removing a zero term: 
 * \f[
 * x \times 0 \rightarrow 0. 
 * \f]
 * \f[
 * 0 \times x \rightarrow 0. 
 * \f]
 */
class RemoveZeroFromProduct : public ProductTransformation
{
public:
  /** */
  RemoveZeroFromProduct() : ProductTransformation() {;}

  /** */
  virtual ~RemoveZeroFromProduct(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

    
/** 
 * Transform a product by removing multiplication by 1.0: 
 * \f[
 * x \times 1.0 \rightarrow x. 
 * \f]
 * \f[
 * 1.0 \times x \rightarrow x. 
 * \f]
 */
class RemoveOneFromProduct : public ProductTransformation
{
public:
  /** */
  RemoveOneFromProduct() : ProductTransformation() {;}

  /** */
  virtual ~RemoveOneFromProduct(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Transform a product by removing multiplication by -1.0: 
 * \f[
 * x \times (-1.0) \rightarrow -x. 
 * \f]
 * \f[
 * -1.0 \times x \rightarrow -x. 
 * \f]
 */
class RemoveMinusOneFromProduct : public ProductTransformation
{
public:
  /** */
  RemoveMinusOneFromProduct() : ProductTransformation() {;}

  /** */
  virtual ~RemoveMinusOneFromProduct(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Multiply two constant exprs without transformation 
 */
class MultiplyConstants : public ProductTransformation
{
public:
  /** */
  MultiplyConstants() : ProductTransformation() {;}

  /** */
  virtual ~MultiplyConstants(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Transform a product by moving any constants to the left:
 * \f[
 * x \times a \rightarrow a \times x
 * \f]
 **/
class MoveConstantsToLeftOfProduct : public ProductTransformation
{
public:
  /** */
  MoveConstantsToLeftOfProduct() : ProductTransformation() {;}

  /** */
  virtual ~MoveConstantsToLeftOfProduct(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Transform a product by any unary minus operations outside the
 * product
 * \f[
 * (-x) \times y \rightarrow -(x \times y)
 * \f]
 * \f[
 * x \times (-y) \rightarrow -(x \times y)
 * \f]
 * \f[
 * )-x) \times (-y) \rightarrow x \times y
 * \f]
 **/
class MoveUnaryMinusOutsideProduct : public ProductTransformation
{
public:
  /** */
  MoveUnaryMinusOutsideProduct() : ProductTransformation() {;}

  /** */
  virtual ~MoveUnaryMinusOutsideProduct(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/**
 * Transform a product by associating any hungry diff op in the
 * left operand with the right operand:
 * \f[
 * (u D_x) v \rightarrow u D_x u
 * \f]
 */
class AssociateHungryDiffOpWithOperand : public ProductTransformation
{
public:
  /** */
  AssociateHungryDiffOpWithOperand() : ProductTransformation() {;}

  /** */
  virtual ~AssociateHungryDiffOpWithOperand(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/**
 * Kill a diff op acting on a constant
 * \f[
 * D_x \alpha \rightarrow 0
 * \f]
 * \f[
 * D_x (\alpha + u) \rightarrow D_x u
 * \f]
 */
class KillDiffOpOnConstant : public ProductTransformation
{
public:
  /** */
  KillDiffOpOnConstant() : ProductTransformation() {;}

  /** */
  virtual ~KillDiffOpOnConstant(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/**
 * Bring a constant outside a diff op
 * \f[
 * D_x (\alpha u) \rightarrow \alpha D_x u
 * \f]
 */
class BringConstantOutsideDiffOp : public ProductTransformation
{
public:
  /** */
  BringConstantOutsideDiffOp() : ProductTransformation() {;}

  /** */
  virtual ~BringConstantOutsideDiffOp(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};
    
/**
 * Distribute a sum of diff ops over their operand
 * \f[
 * (D_1 + D_2) u \rightarrow D_1 u + D_2 u
 * \f]
 */
class DistributeSumOfDiffOps : public ProductTransformation
{
public:
  /** */
  DistributeSumOfDiffOps() : ProductTransformation() {;}

  /** */
  virtual ~DistributeSumOfDiffOps(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/**
 * Apply a simple diff op
 */
class ApplySimpleDiffOp : public ProductTransformation
{
public:
  /** */
  ApplySimpleDiffOp() : ProductTransformation() {;}

  /** */
  virtual ~ApplySimpleDiffOp(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Rearrange a product whose right operand is 
 * a product including a constant
 * such that constants are grouped on the left:
 * \f[
 * \alpha (\beta u) \rightarrow (\alpha\beta) u
 * \f]
 * \f[
 * u (\alpha v) \rightarrow \alpha (u v)
 * \f]
 **/
class RearrangeRightProductWithConstant : public ProductTransformation
{
public:
  /** */
  RearrangeRightProductWithConstant() : ProductTransformation() {;}

  /** */
  virtual ~RearrangeRightProductWithConstant(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Rearrange a product whose left operand is a product including a constant
 * such that constants are grouped on the left:
 * \f[
 * (\alpha u)\beta \rightarrow \alpha\beta u
 * \f]
 * \f[
 * (\alpha u)v \rightarrow \alpha (u v)
 * \f]
 **/
class RearrangeLeftProductWithConstant : public ProductTransformation
{
public:
  /** */
  RearrangeLeftProductWithConstant() : ProductTransformation() {;}

  /** */
  virtual ~RearrangeLeftProductWithConstant(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};

/** 
 * Rearrange a product of a constant and an integral so that
 * the constant is under the integral sign:
 * \f[
 * \alpha \int u \rightarrow \int \alpha u
 * \f]
 **/
class TakeConstantUnderIntegralSign : public ProductTransformation
{
public:
  /** */
  TakeConstantUnderIntegralSign() : ProductTransformation() {;}

  /** */
  virtual ~TakeConstantUnderIntegralSign(){;}

  /** */
  virtual bool doTransform(const RCP<ScalarExpr>& left, 
    const RCP<ScalarExpr>& right,
    RCP<ScalarExpr>& rtn) const ;
};
}

#endif
