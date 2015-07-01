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
#include "SundanceStdSumTransformations.hpp"

#include "SundanceSumExpr.hpp"
#include "SundanceProductExpr.hpp"
#include "SundanceConstantExpr.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceDiffOp.hpp"


#include "SundanceUnaryMinus.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceFunctionalPolynomial.hpp"
#include "SundanceSumOfIntegrals.hpp"
#include "SundanceSumOfBCs.hpp"
#include "SundanceNullCellFilterStub.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;

static Time& polysumTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("IdentifyPolynomialSum"); 
  return *rtn;
}
static Time& reorderSumTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("ReorderSum"); 
  return *rtn;
}

static Time& removeUnaryMinusTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("RemoveUnaryMinusFromSum"); 
  return *rtn;
}

static Time& removeZeroTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("RemoveZeroFromSum"); 
  return *rtn;
}

static Time& sumConstantsTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("SumConstants"); 
  return *rtn;
}

static Time& moveConstantsTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("MoveConstantsToLeftOfSum"); 
  return *rtn;
}

static Time& rearrangeRightSumWithConstantTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("RearrangeRightSumWithConstant"); 
  return *rtn;
}

static Time& rearrangeLeftSumWithConstantTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("RearrangeLeftSumWithConstant"); 
  return *rtn;
}

static Time& sumIntegralsTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("SumIntegrals"); 
  return *rtn;
}




StdSumTransformations::StdSumTransformations()
  : SumTransformationSequence()
{
  append(rcp(new RemoveZeroFromSum()));
//  append(rcp(new RemoveUnaryMinusFromSum()));
  //  append(rcp(new ReorderSum()));

  append(rcp(new SumConstants()));

    append(rcp(new MoveConstantsToLeftOfSum()));
    append(rcp(new RearrangeRightSumWithConstant()));
    append(rcp(new RearrangeLeftSumWithConstant()));
  append(rcp(new IdentifyPolynomialSum()));
  append(rcp(new SumIntegrals()));
}


bool IdentifyPolynomialSum::doTransform(const RCP<ScalarExpr>& left, 
                                        const RCP<ScalarExpr>& right,
                                        int sign, 
                                        RCP<ScalarExpr>& rtn) const
{
  TimeMonitor timer(polysumTimer());
  if (useOptimizedPolynomials())
    {
      /* see if this sum can be written as a polynomial in 
       * symbolic functions */
      if (FunctionalPolynomial::isConvertibleToPoly(left.get())
          && FunctionalPolynomial::isConvertibleToPoly(right.get()))
        {
          RCP<FunctionalPolynomial> lp = FunctionalPolynomial::toPoly(left);
          RCP<FunctionalPolynomial> rp = FunctionalPolynomial::toPoly(right);
          rtn = lp->addPoly(rp.get(), sign);
          return true;
        }
    }
  /* otherwise, do no transformation */
  return false;
}


bool ReorderSum::doTransform(const RCP<ScalarExpr>& left, 
                             const RCP<ScalarExpr>& right,
                             int sign, RCP<ScalarExpr>& rtn) const
{
  TimeMonitor timer(reorderSumTimer());
  /* first we check to see whether the terms are already in order.
   * The left and right trees are already ordered, so that if the
   * first term on the right is after the last term on the left, the
   * combination of both is already ordered. In that case, nothing more needs
   * to be done. 
   */
  Expr L = Expr::handle(left);
  Expr R = Expr::handle(right);
  Sundance::Map<Expr, int> tree = L.getSumTree();
  Sundance::Map<Expr, int>::const_reverse_iterator iL = tree.rbegin();
  Expr endOfLeft = iL->first;
  Sundance::Map<Expr, int> rightTree = R.getSumTree();
  Sundance::Map<Expr, int>::const_iterator iR = rightTree.begin();
  Expr startOfRight = iR->first;

  if (endOfLeft.lessThan(startOfRight))
    {
      Tabs tab1;
      SUNDANCE_VERB_MEDIUM(tab1 << "Terms are already ordered, doing nothing");
      return false;
    }
  else
    {
      Tabs tab1;

      for (Map<Expr, int>::const_iterator i=rightTree.begin(); i!=rightTree.end(); i++)
        {
          int leftCount = 0;
          if (tree.containsKey(i->first))
            {
              leftCount = tree[i->first];
            }
          int count = leftCount + sign * i->second;
          tree.put(i->first, count);
        }
      SUNDANCE_VERB_MEDIUM(tab1 << "Ordered terms are: " << tree);      


      Expr sum = new ZeroExpr();

      /* add up all terms. If no terms are left after cancellations, the result will 
       * be zero. */
      for (Map<Expr, int>::const_iterator i=tree.begin(); i!=tree.end(); i++)
        {
          const Expr& term = i->first;
          int count = i->second;
          if (count==0) continue;
          if (count==1) sum = sum + term;
          else if (count==-1) sum = sum - term;
          else sum = sum + count*term;
        }
      
      rtn = getScalar(sum);
      return true;
    }
}

#ifdef OLD_CODE

bool ReorderSum::doTransform(const RCP<ScalarExpr>& left, 
                             const RCP<ScalarExpr>& right,
                             int sign, RCP<ScalarExpr>& rtn) const
{
  Tabs tabs;
  SUNDANCE_VERB_LOW(tabs << "trying ReorderSum");

  Tabs tab0;
  SUNDANCE_VERB_MEDIUM(tab0 << "L=" << l->toString());
  SUNDANCE_VERB_MEDIUM(tab0 << "R=" << r->toString());
  const SumExpr* sLeft = dynamic_cast<const SumExpr*>(l.get());
  const SumExpr* sRight = dynamic_cast<const SumExpr*>(r.get());

  if (sLeft==0 && sRight==0)
    {
      Tabs tab1;
      if (Expr::handle(r).lessThan(Expr::handle(l)))
        {
          if (sign>0)
            {
              SUNDANCE_VERB_MEDIUM(tab1 << "Both are non-sums: R < L, so reordering to R+L");
              rtn = getScalar(Expr::handle(r) + Expr::handle(l));
            }
          else
            {
              SUNDANCE_VERB_MEDIUM(tab1 << "Both are non-sums: R < L, reordering to -R+L");
              rtn = getScalar(-Expr::handle(r) + Expr::handle(l));
            }
          return true;
        }
    }


  if (sLeft != 0)
    {
      Tabs tab1;
      Expr ll = sLeft->left();
      Expr lr = sLeft->right();
      int lSign = sLeft->sign();
      if (Expr::handle(r).lessThan(ll))
        {
          if (sign > 0)
            {
              SUNDANCE_VERB_MEDIUM(tab1 << "L is a sum: R < LL, reordering to R+L");
              rtn = getScalar(Expr::handle(r) + Expr::handle(l));
            }
          else
            {
              SUNDANCE_VERB_MEDIUM(tab1 << "L is a sum: R < LL, reordering to -R+L");
              rtn = getScalar(-Expr::handle(r) + Expr::handle(l));
            }
          return true;
        }
      if (Expr::handle(r).lessThan(lr)) 
        {
          if (sign > 0)
            {
              if (lSign > 0)
                {
                  SUNDANCE_VERB_MEDIUM(tab1 << "L is a sum: LL < R < LR, reordering to LL + R + LR");
                  rtn = getScalar(ll + Expr::handle(r) + lr);
                }
              else
                {
                  SUNDANCE_VERB_MEDIUM(tab1 << "L is a sum: LL < R < LR, reordering to LL + R - LR");
                  rtn = getScalar(ll + Expr::handle(r) - lr);
                }
            }
          else
            {
              if (lSign > 0)
                {
                  SUNDANCE_VERB_MEDIUM(tab1 << "L is a sum: LL < R < LR, reordering to LL - R + LR");
                  rtn = getScalar(ll - Expr::handle(r) + lr);
                }
              else
                {
                  SUNDANCE_VERB_MEDIUM(tab1 << "L is a sum: LL < R < LR, "
                                       "reordering to LL - R - LR");
                  rtn = getScalar(ll - Expr::handle(r) - lr);
                }
            }
          return true;
        }
    }
  if (sRight != 0)
    {
      Tabs tab1;
      Expr rl = sRight->left();
      Expr rr = sRight->right();
      int rSign = sRight->sign();

      if (rr.lessThan(Expr::handle(l)))
        {
          if (sign > 0)
            {
              SUNDANCE_VERB_MEDIUM(tab1 << "R is a sum: R < L, reordering to R+L");
              rtn = getScalar(Expr::handle(r) + Expr::handle(l));
            }
          else
            {
              SUNDANCE_VERB_MEDIUM(tab1 << "R is a sum: R < L, reordering to -R+L");
              rtn = getScalar(-Expr::handle(r) + Expr::handle(l));
            }
          return true;
        }

      if (rl.lessThan(Expr::handle(l)))
        {
          if (sign > 0)
            {
              if (rSign > 0)
                {
                  SUNDANCE_VERB_MEDIUM(tab1 << "R is a sum: RL < L < RR, reordering to RL+L+RR");
                  rtn = getScalar(rl + Expr::handle(l) + rr);
                }
              else
                {
                  SUNDANCE_VERB_MEDIUM(tab1 << "R is a sum: RL < L < RR, reordering to RL+L-RR");
                  rtn = getScalar(rl + Expr::handle(l) - rr);
                }
            }
          else
            {
              if (rSign > 0)
                {
                  SUNDANCE_VERB_MEDIUM(tab1 << "R is a sum: RL < L < RR, reordering to -RL+L-RR");
                  rtn = getScalar(-rl + Expr::handle(l) - rr);
                }
              else
                {
                  SUNDANCE_VERB_MEDIUM(tab1 << "R is a sum: RL < L < RR, reordering to -RL+L+RR");
                  rtn = getScalar(-rl + Expr::handle(l) + rr);
                }
            }
          return true;
        }
    }
  
  return false;
}

#endif


bool RemoveZeroFromSum::doTransform(const RCP<ScalarExpr>& left, const RCP<ScalarExpr>& right,
                                    int sign, RCP<ScalarExpr>& rtn) const
{
  TimeMonitor timer(removeZeroTimer());
  SUNDANCE_OUT(this->verb() > 1, 
               "trying RemoveZeroFromSum");

  /* Check for the trivial case of operation with zero */
  
  /* If I'm constant and my value is zero, return other */
  const ConstantExpr* cl = dynamic_cast<const ConstantExpr*>(left.get());
  if (cl != 0 && (cl->value()==0.0 || cl->value()==-0.0))
    {
      if (verb() > 1)
        {
          Out::println("RemoveZeroFromSum identified left operand as zero.");
          Out::println("Applying transformation 0 + x --> x");
        }
      rtn = chooseSign(sign, right);
      return true;
    }

  /* If other is zero, return me */
  const ConstantExpr* cr = dynamic_cast<const ConstantExpr*>(right.get());
  if (cr != 0 && (cr->value()==0.0 || cr->value()==-0.0)) 
    {
      if (verb() > 1)
        {
          Out::println("RemoveZeroFromSum identified right operand as zero.");
          Out::println("Applying transformation x + 0 --> x");
        }
      rtn = left;
      return true;
    }

  /* otherwise, do no transformation */
  return false;
}

bool MoveConstantsToLeftOfSum::doTransform(const RCP<ScalarExpr>& left, const RCP<ScalarExpr>& right,
                                      int sign, RCP<ScalarExpr>& rtn) const
{
  TimeMonitor timer(moveConstantsTimer());

  /* if the right operand is a constant, 
   * transform u +/- alpha --> +/- alpha + u */
  if (right->isConstant())
    {
      if (verb() > 1)
        {
          Out::println("MoveConstantsToLeftOfSum identified right "
                       "operand as constant.");
        }
      rtn = getScalar(Expr::handle(chooseSign(sign, right)) 
        + Expr::handle(left));
      return true;
    }

  return false;
}


bool RemoveUnaryMinusFromSum::doTransform(const RCP<ScalarExpr>& left,
                                          const RCP<ScalarExpr>& right,
                                          int sign, 
                                          RCP<ScalarExpr>& rtn) const
{
  TimeMonitor timer(removeUnaryMinusTimer());
  SUNDANCE_OUT(this->verb() > 1, 
               "trying RemoveUnaryMinusFromSum");

  /* if the right operand is a unary minus, 
   * transform u +/- (-v) --> u -/+ v */
  const UnaryMinus* ul = dynamic_cast<const UnaryMinus*>(left.get());
  const UnaryMinus* ur = dynamic_cast<const UnaryMinus*>(right.get());
  if (ul != 0 && ur != 0)
    {
      if (verb() > 1)
        {
          Out::println("RemoveUnaryMinusFromSum identified both "
                       "operands as unary minuses.");
        }
      rtn = getScalar(-(Expr::handle(chooseSign(sign, getScalar(ur->arg()))) 
                        + ul->arg()));
      return true;
    }
  else if (ul != 0)
    {
      if (verb() > 1)
        {
          Out::println("RemoveUnaryMinusFromSum identified left "
                       "operand as unary minus.");
        }
      if (sign > 0)
        {
          rtn = getScalar(-(ul->arg() - Expr::handle(right)));
        }
      else
        {
          rtn = getScalar(-(ul->arg() + Expr::handle(right)));
        }
      return true;
    }
  else if (ur != 0)
    {
      if (verb() > 1)
        {
          Out::println("RemoveUnaryMinusFromSum identified right "
                       "operand as unary minus.");
        }
      if (sign > 0)
        {
          rtn = getScalar(Expr::handle(left) - ur->arg());
        }
      else
        {
          rtn = getScalar(Expr::handle(left) + ur->arg());
        }
      return true;
    }

  return false;
}

bool SumConstants::doTransform(const RCP<ScalarExpr>& left, const RCP<ScalarExpr>& right,
                               int sign, RCP<ScalarExpr>& rtn) const
{
  
  TimeMonitor timer(sumConstantsTimer());
  SUNDANCE_OUT(this->verb() > 1, 
               "trying SumConstants");

  /* Check to see if both are constant. If so, sum them and return */
  if (left->isConstant() && right->isConstant())
    {
      if (verb() > 1)
        {
          Out::println("SumConstants identified both "
                       "operands as constant. No transformations applied.");
        }
      rtn = rcp(new SumExpr(left, right, sign));
      return true;
    }
  return false;
}

bool RearrangeRightSumWithConstant::doTransform(const RCP<ScalarExpr>& left, 
                                                const RCP<ScalarExpr>& right,
                                                int sign, RCP<ScalarExpr>& rtn) const
{
  TimeMonitor timer(rearrangeRightSumWithConstantTimer());
  const SumExpr* sRight = dynamic_cast<const SumExpr*>(right.get());

  if (sRight != 0)
    {
      /* The case in which the right operand of right is a constant
       * should have been transformed away by now. Do a paranoid 
       * check to make sure this hasn't happened */
      TEUCHOS_TEST_FOR_EXCEPTION(sRight->rightScalar()->isConstant(),
                         std::logic_error,
                         "RearrangeRightSumWithConstant: unexpected case, "
                         "constant expr"
                         << sRight->right() << " found as right operand "
                         "in sum " << right->toString());

      if (sRight->leftScalar()->isConstant())
        {      
          /* If left operand is a constant, transform
           * alpha + s1*(beta + s2*u) --> (alpha + s1*beta) + s1*s2*u */
          if (left->isConstant())
            {
              if (verb() > 1)
                {
                  Out::println("RearrangeRightSumWithConstant::doTransform "
                               "identified right "
                               "operand as sum involving a constant, "
                               "and left operand as a constant. Applying "
                               "transformation alpha + (beta+u) "
                               "--> (alpha + beta) + u.");
                }
              int s1 = sign;
              int s2 = sRight->sign();
              Expr alpha = Expr::handle(left);
              Expr beta = sRight->left();
              Expr u = sRight->right();
              rtn = getScalar((alpha + chooseSign(s1, beta)) + chooseSign(s1*s2,u));
            }
          else  /* if left operand is non-constant, transform
                 * u + s1*(alpha + s2*v) --> s1*alpha + (u + s1*s2*v) */
            {
              if (verb() > 1)
                {
                  Out::println("RearrangeRightSumWithConstant::doTransform "
                               "identified right "
                               "operand as sum involving a constant, "
                               "and left operand as non-constant. Applying "
                               "transformation u + (alpha + v) "
                               "--> alpha + (u + v)");
                }
              int s1 = sign;
              int s2 = sRight->sign();
              Expr u = Expr::handle(left);
              Expr alpha = sRight->left();
              Expr v = sRight->right();
              rtn = getScalar(chooseSign(s1, alpha) + (u + chooseSign(s1*s2, v)));
            }
          return true;
        }
    }
  return false;
}


bool RearrangeLeftSumWithConstant::doTransform(const RCP<ScalarExpr>& left, 
                                                const RCP<ScalarExpr>& right,
                                                int sign, RCP<ScalarExpr>& rtn) const
{
  TimeMonitor timer(rearrangeLeftSumWithConstantTimer());
  const SumExpr* sLeft = dynamic_cast<const SumExpr*>(left.get());

  if (sLeft != 0 && !left->isConstant())
    {
      /* The case in which the right operand of left is a constant
       * should have been transformed away by now. Do a paranoid 
       * check to make sure this hasn't happened */
      TEUCHOS_TEST_FOR_EXCEPTION(sLeft->rightScalar()->isConstant(),
                         std::logic_error,
                         "RearrangeLeftSumWithConstant::doTransform "
                         ": unexpected case, constant expr"
                         << sLeft->right() << " found as right operand "
                         "in sum " << left->toString());
      
      if (sLeft->leftScalar()->isConstant())
        {
          /* If right operand is a constant, transform 
           * (alpha + s1*u) + s2*beta --> (alpha + s2*beta) + s1*u */
          if (right->isConstant())
            {
              if (verb() > 1)
                {
                  Out::println("RearrangeLeftSumWithConstant::doTransform "
                               "identified right "
                               "operand as constant, "
                               "and left operand as sum involving "
                               "a constant. Applying "
                               "transformation (alpha + u) + beta "
                               "--> (alpha + beta) + u.");
                }
              int s2 = sign;
              int s1 = sLeft->sign();
              Expr alpha = sLeft->left();
              Expr beta = Expr::handle(right);
              Expr u = sLeft->right();
              rtn =  getScalar((alpha + chooseSign(s2, beta)) + chooseSign(s1, u));
            }
          else /* if right operand is a non-constant, transform 
                * (alpha + s1*u) + s2*v --> alpha + (s1*u + s2*v) */
            {
              if (verb() > 1)
                {
                  Out::println("RearrangeLeftSumWithConstant::doTransform "
                               "identified right "
                               "operand as non-constant, "
                               "and left operand as sum involving "
                               "a constant. Applying "
                               "transformation (alpha + u) + v "
                               "--> alpha + (u + v).");
                }
              int s2 = sign;
              int s1 = sLeft->sign();
              Expr alpha = sLeft->left();
              Expr u = sLeft->right();
              Expr v = Expr::handle(right);
              rtn =  getScalar(alpha 
                + (chooseSign(s1, u) + chooseSign(s2, v)));
            }
          return true;
        }
    }
  return false;
}

bool SumIntegrals::doTransform(const RCP<ScalarExpr>& left, 
                               const RCP<ScalarExpr>& right,
                               int sign, RCP<ScalarExpr>& rtn) const
{
  TimeMonitor timer(sumIntegralsTimer());
  SUNDANCE_OUT(this->verb() > 1, 
               "trying SumIntegrals");

  const SumOfIntegrals* sLeft 
    = dynamic_cast<const SumOfIntegrals*>(left.get());
  const SumOfIntegrals* sRight 
    = dynamic_cast<const SumOfIntegrals*>(right.get());

  if (sLeft != 0 || sRight != 0)
    {
      /* make sure we don't have a case where one is an essential BC and
       * the other is not */
      bool leftIsBC = (dynamic_cast<const SumOfBCs*>(sLeft) != 0);
      bool rightIsBC = (dynamic_cast<const SumOfBCs*>(sRight) != 0);
      TEUCHOS_TEST_FOR_EXCEPTION((leftIsBC && !rightIsBC)
                         || (!leftIsBC && rightIsBC), std::runtime_error,
                         "Attempting to add EssentialBC and non-EssentialBC "
                         "integrals: L=" << left->toString() << ", R="
                         << right->toString());

      if (sLeft != 0 && sRight != 0)
        {
          SumOfIntegrals* l;
          if (!leftIsBC) l = new SumOfIntegrals(*sLeft);
          else l = new SumOfBCs(*dynamic_cast<const SumOfBCs*>(sLeft));
          l->merge(sRight, sign);
          rtn = rcp(l);
          return true;
        }

      /* at this point, one of the terms is a global equation. BCs should
       * not be involved at this point */
      TEUCHOS_TEST_FOR_EXCEPTION(leftIsBC, std::runtime_error,
                         "Attempting to add a BC " << left->toString()
                         << " and a global expression " << right->toString());

      TEUCHOS_TEST_FOR_EXCEPTION(rightIsBC, std::runtime_error,
                         "Attempting to add a BC " << right->toString()
                         << " and a global expression " << left->toString());

      if (sLeft != 0 && sRight == 0)
        {
          SumOfIntegrals* l = new SumOfIntegrals(*sLeft);
          const SpatiallyConstantExpr* cRight 
            = dynamic_cast<const SpatiallyConstantExpr*>(right.get());

          TEUCHOS_TEST_FOR_EXCEPTION(cRight == 0, std::logic_error,
                             "Attempting to add non-constant expression "
                             << right->toString() << " to an integral");

          Expr r = Integral(l->nullRegion(), Expr::handle(right));
          const SumOfIntegrals* sr 
            = dynamic_cast<const SumOfIntegrals*>(r.ptr().get());
          l->merge(sr, sign);
          rtn = rcp(l);
          return true;
        }
      else
        {
          SumOfIntegrals* r = new SumOfIntegrals(*sRight);
          if (sign < 0) r->changeSign();

          const SpatiallyConstantExpr* cLeft 
            = dynamic_cast<const SpatiallyConstantExpr*>(left.get());

          TEUCHOS_TEST_FOR_EXCEPTION(cLeft == 0, std::logic_error,
                             "Attempting to add non-constant expression "
                             << left->toString() << " to an integral");

          Expr l = Integral(r->nullRegion(), Expr::handle(right));
          const SumOfIntegrals* sl 
            = dynamic_cast<const SumOfIntegrals*>(l.ptr().get());
          r->merge(sl, 1);
          rtn = rcp(r);
          return true;
        }

    }
  else
    {
      return false;
    }
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "this should not happen");
  return false; // -Wall;
}
