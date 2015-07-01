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

#include "SundanceExpr.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceStdProductTransformations.hpp"

#include "SundanceSumExpr.hpp"
#include "SundanceProductExpr.hpp"
#include "SundanceConstantExpr.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceDiffOp.hpp"


#include "SundanceUnaryMinus.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceDerivOfSymbFunc.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


StdProductTransformations::StdProductTransformations()
  : ProductTransformationSequence()
{
  append(rcp(new RemoveZeroFromProduct()));
  append(rcp(new RemoveOneFromProduct()));
  append(rcp(new RemoveMinusOneFromProduct()));
  append(rcp(new KillDiffOpOnConstant()));
  append(rcp(new BringConstantOutsideDiffOp()));
  append(rcp(new MoveUnaryMinusOutsideProduct()));
  append(rcp(new AssociateHungryDiffOpWithOperand()));
  append(rcp(new DistributeSumOfDiffOps()));
  append(rcp(new ApplySimpleDiffOp()));
  append(rcp(new TakeConstantUnderIntegralSign()));
  append(rcp(new MultiplyConstants()));
  append(rcp(new MoveConstantsToLeftOfProduct()));
  append(rcp(new RearrangeRightProductWithConstant()));
  append(rcp(new RearrangeLeftProductWithConstant()));
}

bool RemoveZeroFromProduct::doTransform(const RCP<ScalarExpr>& left, 
                                        const RCP<ScalarExpr>& right,
                                        RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, 
               "trying RemoveZerofromProduct");
  
  /* Check for the trivial case of multiplication by zero */
  const ConstantExpr* cl = dynamic_cast<const ConstantExpr*>(left.get());
  const ConstantExpr* cr = dynamic_cast<const ConstantExpr*>(right.get());

  if (cl != 0)
    {
      if (cl->value()==0.0 || cl->value()==-0.0)
        {
          if (verb() > 1)
            {
              Out::println("RemoveOneFromProduct::doTransform "
                           "identified multiplication "
                           "by zero. Applying transformation 0*u --> 0");
            }
          rtn = rcp(new ZeroExpr());
          return true;
        }
    }
  if (cr != 0)
    {
      if (cr->value()==0.0 || cr->value()==-0.0)
        {
          if (verb() > 1)
            {
              Out::println("RemoveOneFromProduct::doTransform "
                           "identified multiplication "
                           "by zero. Applying transformation u*0 --> u");
            }
          rtn = rcp(new ZeroExpr());
          return true;
        }
    }
  return false;
}

bool RemoveOneFromProduct::doTransform(const RCP<ScalarExpr>& left, 
                                       const RCP<ScalarExpr>& right,
                                       RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, 
               "trying RemoveOnefromProduct");

  /* Check for the trivial case of multiplication by one */
  const ConstantExpr* cl = dynamic_cast<const ConstantExpr*>(left.get());
  const ConstantExpr* cr = dynamic_cast<const ConstantExpr*>(right.get());

  if (cl != 0)
    {
      if (cl->value()==1.0)
        {
          if (verb() > 1)
            {
              Out::println("RemoveOneFromProduct::doTransform "
                           "identified multiplication "
                           "by one. Applying transformation 1*u --> u");
            }
          rtn = getScalar(Expr::handle(right));
          return true;
        }
    }
  if (cr != 0)
    {
      if (cr->value()==1.0)
        {
          if (verb() > 1)
            {
              Out::println("RemoveOneFromProduct::doTransform "
                           "identified multiplication "
                           "by one. Applying transformation u*1 --> u");
            }
          rtn = getScalar(Expr::handle(left));
          return true;
        }
    }
  return false;
}


bool RemoveMinusOneFromProduct::doTransform(const RCP<ScalarExpr>& left, 
                                            const RCP<ScalarExpr>& right,
                                            RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, 
               "trying RemoveOnefromProduct");

  /* Check for the trivial case of multiplication by minus one */
  const ConstantExpr* cl = dynamic_cast<const ConstantExpr*>(left.get());
  const ConstantExpr* cr = dynamic_cast<const ConstantExpr*>(right.get());

  if (cl != 0)
    {
      if (cl->value()==-1.0)
        {
          if (verb() > 1)
            {
              Out::println("RemoveMinusOneFromProduct::doTransform "
                           "identified multiplication "
                           "by one. Applying transformation -1*u --> -u");
            }
          rtn = getScalar(-Expr::handle(right));
          return true;
        }
    }
  if (cr != 0)
    {
      if (cr->value()==-1.0)
        {
          if (verb() > 1)
            {
              Out::println("RemoveMinusOneFromProduct::doTransform "
                           "identified multiplication "
                           "by one. Applying transformation u*(-1) --> -u");
            }
          rtn = getScalar(-Expr::handle(left));
          return true;
        }
    }
  return false;
}

bool MoveConstantsToLeftOfProduct::doTransform(const RCP<ScalarExpr>& left, 
                                               const RCP<ScalarExpr>& right,
                                               RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, 
               "trying MoveConstantsToLeftOfProduct");

  /* If the left operand is non-constant and
   * the right operand is a constant, 
   * transform u*constant --> constant*u */
  if (!left->isConstant() && right->isConstant())
    {
      if (verb() > 1)
        {
          Out::println("MoveConstantsToLeftOfProduct::doTransform "
                       "identified right operand "
                       "as a constant. Applying transformation u*alpha "
                       "--> alpha*u.");
        }
      rtn = getScalar(Expr::handle(right) * Expr::handle(left));
      return true;
    }
  return false;
}

bool MoveUnaryMinusOutsideProduct::doTransform(const RCP<ScalarExpr>& left, 
                                               const RCP<ScalarExpr>& right,
                                               RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, 
               "trying MoveUnaryMinusOutsideProduct");

  /* If one of the operands is a unary minus, apply it to the whole
   * product. If both are unary minuses, multiply the operands, removing
   * the unary minuses. */
  const UnaryMinus* ul = dynamic_cast<const UnaryMinus*>(left.get());
  const UnaryMinus* ur = dynamic_cast<const UnaryMinus*>(right.get());
  if (ur != 0 && ul != 0)
    {
      if (verb() > 1)
        {
          Out::println("MoveUnaryMinusOutsideProduct::doTransform "
                       "identified both operands "
                       "as unary minuses. Applying transformation (-x)*(-y) "
                       "--> x*y.");
        }
      rtn = getScalar(ul->arg() * ur->arg());
      return true;
    }
  else if (ur != 0)
    {
      if (verb() > 1)
        {
          Out::println("MoveUnaryMinusOutsideProduct::doTransform "
                       "identified right operand "
                       "as a unary minus. Applying transformation x*(-y) "
                       "--> -(x*y).");
        }
      Expr prod = Expr::handle(left) * ur->arg();
      rtn = rcp(new UnaryMinus(getScalar(prod)));
      return true;
    }
  else if (ul != 0)
    {
      if (verb() > 1)
        {
          Out::println("MoveUnaryMinusOutsideProduct::doTransform "
                       "identified left operand "
                       "as a unary minus. Applying transformation (-x)*y "
                       "--> -(x*y).");
        }
      Expr prod = ul->arg() * Expr::handle(right);
      rtn = rcp(new UnaryMinus(getScalar(prod)));
      return true;
    }
  return false;
}

bool MultiplyConstants::doTransform(const RCP<ScalarExpr>& left, 
  const RCP<ScalarExpr>& right,
  RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, 
               "trying MultiplyConstants");

  /* If both operands are constant, just multiply them */
  if (left->isConstant() && right->isConstant())
    {
      if (left->isImmutable() && right->isImmutable())
        {
          const ConstantExpr* cl = dynamic_cast<const ConstantExpr*>(left.get());
          const ConstantExpr* cr = dynamic_cast<const ConstantExpr*>(right.get());
          TEUCHOS_TEST_FOR_EXCEPTION(cl==0 || cr==0, std::logic_error,
                             "MultiplyConstants::doTransform() logic error: "
                             "L and R identified as immutable, but could "
                             "not be cast to ConstantExprs");
          rtn = rcp(new ConstantExpr(cl->value() * cr->value()));
          return true;
        }
      if (verb() > 1)
        {
          Out::println("MultiplyConstants::doTransform "
                       "identified both operands "
                       "as constants. Forming the product without any "
                       "transformation");
        }
      rtn = rcp(new ProductExpr(left, right));
      return true;
    }
  return false;
}

bool AssociateHungryDiffOpWithOperand::doTransform(const RCP<ScalarExpr>& left, 
  const RCP<ScalarExpr>& right,
  RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, 
               "trying AssociateHungryDiffOpWithOperand");

  if (left->isHungryDiffOp())
    {
      /* if the left operand is a product including 
       * a hungry diff op, rotate the
       * tree such that the diff op associates with the right operand */
      const ProductExpr* pLeft 
        = dynamic_cast<const ProductExpr*>(left.get());
      if (pLeft != 0)
        {
          Expr ll = pLeft->left();
          Expr lr = pLeft->right();
          if (verb() > 1)
            {
              Out::println("AssociateHungryDiffOpWithOperand::doTransform "
                           "identified left "
                           "operand as a product with a hungry diff op "
                           "as the last factor. "
                           "Applying (u*D)*v --> u*(D*v).");
            }
          rtn = getScalar(ll*(lr*Expr::handle(right)));
          return true;
        }
    }
  return false;
}

bool KillDiffOpOnConstant::doTransform(const RCP<ScalarExpr>& left, 
                                       const RCP<ScalarExpr>& right,
                                       RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, "trying KillDiffOpOnConstant");

  if (left->isHungryDiffOp())
    {
      /* first check that the operand is not a constant, in case someone
       * differentiated a constant for some unimaginable reason */
      if (right->isConstant())
        {
          rtn = rcp(new ZeroExpr());
          if (verb() > 1)
            {
              Out::println("KillDiffOpOnConstant::doTransform "
                           "identified constant "
                           "as operand of diff op. Applying D*alpha --> 0");
            }
          return true;
        }
      
      /* transform op*(constant + u) --> op*u */
      const SumExpr* sRight = dynamic_cast<const SumExpr*>(right.get());
      if (sRight != 0 && sRight->leftScalar()->isConstant())
        {
          if (verb() > 1)
            {
              Out::println("KillDiffOpOnConstant::doTransform "
                           "identified constant "
                           "term in operand of diff op. "
                           "Applying D*(alpha+u) --> D*u");
            }
          rtn = getScalar(chooseSign(sRight->sign(), Expr::handle(left)) * sRight->right());
          return true;
        }
    }
  return false;
}

bool BringConstantOutsideDiffOp::doTransform(const RCP<ScalarExpr>& left, 
                                             const RCP<ScalarExpr>& right,
                                             RCP<ScalarExpr>& rtn) const
{
  
  SUNDANCE_OUT(this->verb() > 1,  "trying BringConstantOutsideDiffOp");

  if (left->isHungryDiffOp())
    {
      /* transform op*(constant*u) --> constant*op*u */
      const ProductExpr* pRight 
        = dynamic_cast<const ProductExpr*>(right.get());
      if (pRight != 0 && pRight->leftScalar()->isConstant())
        {
          if (verb() > 1)
            {
              Out::println("BringConstantOutsideDiffOp::doTransform "
                           "identified constant "
                           "coefficient in operand of diff op. "
                           "Applying D*(alpha*u) --> alpha*D*u");
            }
          rtn = getScalar(pRight->left() * (Expr::handle(left) * pRight->right()));
          return true;
        }
    }
  return false;
}

bool DistributeSumOfDiffOps::doTransform(const RCP<ScalarExpr>& left, 
                                         const RCP<ScalarExpr>& right,
                                         RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1,"trying DistributeSumOfDiffOps");

  if (left->isHungryDiffOp())
    {
      /* if the left operand is a sum of hungry diff ops, distribute this
       * multiplication over the sum */
      const SumExpr* sLeft = dynamic_cast<const SumExpr*>(left.get());
      if (sLeft != 0)
        {
          Expr ll = sLeft->left();
          Expr lr = sLeft->right();
          if (verb() > 1)
            {
              Out::println("DistributeSumOfDiffOps::doTransform "
                           "identified left "
                           "operand as a sum of hungry diff ops. "
                           "Applying (D1 + D2)*u --> D1*u + D2*u");
            }
          rtn = getScalar(ll*Expr::handle(right) + chooseSign(sLeft->sign(), lr)*Expr::handle(right));
          return true;
        }
    }
  return false;
}

bool ApplySimpleDiffOp::doTransform(const RCP<ScalarExpr>& left, 
                                    const RCP<ScalarExpr>& right,
                                    RCP<ScalarExpr>& rtn) const
{
  SUNDANCE_OUT(this->verb() > 1, "trying ApplySimpleDiffOp");

  if (left->isHungryDiffOp())
    {
      const Derivative* dLeft 
        = dynamic_cast<const Derivative*>(left.get());
      if (dLeft != 0)
        {
          const SymbolicFuncElement* sf 
            = dynamic_cast<const SymbolicFuncElement*>(right.get());
          if (sf != 0 && optimizeFunctionDiffOps())
            {
              rtn = rcp(new DerivOfSymbFunc(dLeft->multiIndex(), right));
            }
          else
            {
              rtn = rcp(new DiffOp(dLeft->multiIndex(), right));
            }
          return true;
        }
    }
  return false;
}

bool RearrangeRightProductWithConstant::doTransform(const RCP<ScalarExpr>& left, 
                                                    const RCP<ScalarExpr>& right,
                                                    RCP<ScalarExpr>& rtn) const
{
  /* Transform several cases in which the right operand is a product
   * involving a constant. */
  const ProductExpr* pRight 
    = dynamic_cast<const ProductExpr*>(right.get());

  /* By this point, we should have already transformed out any cases
   * in which the right operand is a product having a constant as a
   * right operand, because its constants should have been rotated
   * left. Do a paranoia check to be safe */
  TEUCHOS_TEST_FOR_EXCEPTION(pRight != 0 && pRight->rightScalar()->isConstant(),
                     std::logic_error,
                     "unexpected case in "
                     "RearrangeRightProductWithConstant::doTransform: "
                     "the right operand "
                     << pRight->right() 
                     << "of the right operand " << right->toString()
                     << " is a constant. That case should have been "
                     "transformed out by now.");

  if (pRight != 0 && !pRight->isConstant() 
      && pRight->leftScalar()->isConstant())
    {
      /* if left operand is a constant, and the right operand is a 
       * product involving a constant,
       * transform alpha*(beta*u) --> (alpha*beta)*u */
      if (left->isConstant())
        {
          if (verb() > 1)
            {
              Out::println("RearrangeRightProductWithConstant::doTransform: "
                           "identified left operand "
                           "as a constant, and right operand as a product "
                           "involving a constant. Applying transformation "
                           "alpha*(beta*u) --> (alpha*beta)*u");
            }
          rtn = getScalar((Expr::handle(left) * pRight->left()) * pRight->right());
          return true;
        }
      else
        /* if the left operand is non-constant and the right operand
         * is a product involving a constant,
         * transform u * (alpha*v) --> alpha*(u*v) */
        {
          if (verb() > 1)
            {
              Out::println("RearrangeRightProductWithConstant::doTransform: "
                           "identified left operand "
                           "as non-constant, and right operand as a product "
                           "involving a constant. Applying transformation "
                           "u * (alpha*v) --> alpha*(u*v)");
            }
          rtn = getScalar(pRight->left() * (Expr::handle(left) * pRight->right()));
          return true;
        }
    }
  return false;
}


bool RearrangeLeftProductWithConstant::doTransform(const RCP<ScalarExpr>& left, 
                                                   const RCP<ScalarExpr>& right,
                                                   RCP<ScalarExpr>& rtn) const
{
  /* transform cases in which the left operand is a product, exactly
   * one of whose operands is a constant. Because of the preceding 
   * transformation rules,
   * the constant should be on the left */
  const ProductExpr* pLeft 
    = dynamic_cast<const ProductExpr*>(left.get());

  if (pLeft != 0 && !pLeft->isConstant() 
      && pLeft->leftScalar()->isConstant())
    {

      /* Paranoid check to make sure we don't have the case
       * (u*alpha)*right */
      TEUCHOS_TEST_FOR_EXCEPTION(pLeft != 0 && pLeft->rightScalar()->isConstant(),
                         std::logic_error,
                         "RearrangeLeftProductWithConstant::doTransform: "
                         "the right operand "
                         << pLeft->right() 
                         << "of the left operand " << left->toString()
                         << " is a constant. That case should have been "
                         "transformed out by now.");
      /* if the right operand is a constant, 
       * transform (alpha*u)*beta --> (alpha*beta)*u */
      if (right->isConstant())
        {
          if (verb() > 1)
            {
              Out::println("RearrangeLeftProductWithConstant::doTransform: "
                           "identified right operand "
                           "as a constant, and left operand as a product "
                           "involving a constant. Applying transformation "
                           "(alpha*u)*beta --> (alpha*beta)*u");
            }
          rtn =  getScalar((pLeft->left() * Expr::handle(right)) * pLeft->right());
          return true;
        }
      else
        /* if the right operand is non-constant, 
         * transform (alpha*u)*v --> alpha*(u*v) */
        {
          if (verb() > 1)
            {
              Out::println("RearrangeLeftProductWithConstant::doTransform: "
                           "identified right operand "
                           "as non-constant, and left operand as a product "
                           "involving a constant. Applying transformation "
                           "(alpha*u)*v --> alpha*(u*v)");
            }
          rtn = getScalar(pLeft->left() * (pLeft->right() * Expr::handle(right)));
          return true;
        }
    }
  return false;
}


bool TakeConstantUnderIntegralSign::doTransform(const RCP<ScalarExpr>& left, 
                                         const RCP<ScalarExpr>& right,
                                         RCP<ScalarExpr>& rtn) const
{
  const SumOfIntegrals* sLeft 
    = dynamic_cast<const SumOfIntegrals*>(left.get());
  const SumOfIntegrals* sRight 
    = dynamic_cast<const SumOfIntegrals*>(right.get());

  TEUCHOS_TEST_FOR_EXCEPTION(sLeft != 0 && sRight != 0, std::logic_error,
                     "Product of integrals detected: left=" 
                     << left->toString() << " right=" << right->toString());

  if (sLeft != 0 || sRight != 0)
    {
      if (sLeft != 0)
        {
          SumOfIntegrals* l = new SumOfIntegrals(*sLeft);
          const SpatiallyConstantExpr* cRight 
            = dynamic_cast<const SpatiallyConstantExpr*>(right.get());
          TEUCHOS_TEST_FOR_EXCEPTION(cRight == 0, std::logic_error,
                             "Attempting to multiply non-constant expression "
                             << right->toString() << " with an integral");
          l->multiplyByConstant(cRight);
          
          rtn = rcp(l);
          return true;
        }

      if (sRight != 0)
        {
          SumOfIntegrals* r = new SumOfIntegrals(*sRight);
          const SpatiallyConstantExpr* cLeft
            = dynamic_cast<const SpatiallyConstantExpr*>(left.get());
          TEUCHOS_TEST_FOR_EXCEPTION(cLeft == 0, std::logic_error,
                             "Attempting to multiply non-constant expression "
                             << left->toString() << " with an integral");
          r->multiplyByConstant(cLeft);
          
          rtn = rcp(r);
          return true;
        }
    }
  else
    {
      return false;
    }
  return false;
}
