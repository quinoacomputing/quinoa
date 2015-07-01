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

#include "SundanceEvalVector.hpp"
#include "SundanceTempStack.hpp"
#include "PlayaExceptions.hpp"
#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;




EvalVector::EvalVector(TempStack* s)
  : s_(s),
    data_(s->popVectorData()),
    str_()
{
  data_->resize(s->vecSize());
}


EvalVector::EvalVector(TempStack* s, const RCP<Array<double> >& data,
                       const std::string& str)
  : s_(s),
    data_(s->popVectorData()),
    str_(str)
{
  //TimeMonitor t(evalVecTimer());
  data_->resize(data->size());
  int n = data_->size();

  if (n > 0)
    {
      double* x = &((*data_)[0]);
      const double* y = &((*data)[0]);
      for (int i=0; i<n; i++)
        {
          x[i] = y[i];
        }
    }
}

EvalVector::~EvalVector()
{
  s_->pushVectorData(data_);
}


void EvalVector::resize(int n)
{
  //TimeMonitor t(evalVecTimer());
  data_->resize(n);
}

RCP<EvalVector> EvalVector::clone() const
{
  return rcp(new EvalVector(s_, data_, str_));
}

void EvalVector::setToConstant(const double& alpha) 
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();
  if (n > 0)
    {
      double* x = &((*data_)[0]);
      for (int i=0; i<n; i++)
        {
          x[i] = alpha;
        }
    }
  if (shadowOps()) str_ = Teuchos::toString(alpha);
}


void EvalVector::add_SV(const double& alpha,
                        const EvalVector* B)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      const double* const Bx = B->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] += alpha*Bx[i];
        }

      addFlops(2*n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = "(" + str_ + "+" 
        + Teuchos::toString(alpha) + "*" + B->str_ + ")";
      else str_ = Teuchos::toString(alpha) + "*" + B->str_;
    }
}

void EvalVector::add_S(const double& alpha)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      
      for (int i=0; i<n; i++)
        {
          x[i] += alpha;
        }
      addFlops(n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = "(" + str_ + "+" 
        + Teuchos::toString(alpha) + ")";
      else str_ = Teuchos::toString(alpha);
    }
}


void EvalVector::add_V(const EvalVector* A)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] += Ax[i];
        }
      addFlops(n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = "(" + str_ + " + " +  A->str_ + ")";
      else str_ = A->str_;
    }
}

void EvalVector::add_SVV(const double& alpha,
                         const EvalVector* B,
                         const EvalVector* C)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      const double* const Bx = B->start();
      const double* const Cx = C->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] += alpha*Bx[i]*Cx[i];
        }
      addFlops(3*n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = "(" + str_ + " + " 
        + Teuchos::toString(alpha) + "*" + B->str() + "*" + C->str() + ")";
      else str_ = Teuchos::toString(alpha) + "*" + B->str() + "*" + C->str();
    }
}

void EvalVector::add_VV(const EvalVector* A,
                        const EvalVector* B)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      const double* const Bx = B->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] += Ax[i]*Bx[i];
        }
      addFlops(2*n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = "(" + str_ + " + " + A->str() 
        + "*" + B->str() + ")";
      else str_ =A->str() + "*" + B->str();
    }
}


void EvalVector::multiply_S_add_SV(const double& alpha, 
                                   const double& beta,
                                   const EvalVector* C)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();
  
  if (n > 0)
    {
      double* const x = start();
      const double* const Cx = C->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= alpha;
          x[i] += beta*Cx[i];
        }
      addFlops(3*n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = "(" + Teuchos::toString(alpha) + "*" + str_ + "+"
        + Teuchos::toString(beta) + "*" + C->str_ + ")";
      else str_ = Teuchos::toString(beta) + "*" + C->str_;
    }
}


void EvalVector::multiply_S_add_S(const double& alpha, 
                                  const double& beta)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();
  
  if (n > 0)
    {
      double* const x = start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= alpha;
          x[i] += beta;
        }
      addFlops(2*n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = "(" + Teuchos::toString(alpha) + "*" + str_
        + " + " + Teuchos::toString(beta) + ")";
      else str_ = Teuchos::toString(beta);
    }
}

void EvalVector::multiply_V(const EvalVector* A)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();
  
  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= Ax[i];
        }
      addFlops(n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = str() + "*" + A->str();
    }
}

void EvalVector::multiply_V_add_VVV(const EvalVector* A,
                                    const EvalVector* B,
                                    const EvalVector* C,
                                    const EvalVector* D)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      const double* const Bx = B->start();
      const double* const Cx = C->start();
      const double* const Dx = D->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= Ax[i];
          x[i] += Bx[i]*Cx[i]*Dx[i];
        }
      addFlops(4*n);
    }

  if (shadowOps())
    {
      if (str_ != "0") str_ = "(" + str() + "*" + A->str() + " + " 
        + B->str() + "*" + C->str() + "*" + D->str() + ")";
      else str_ = B->str() + "*" + C->str() + "*" + D->str();
    }
}

void EvalVector::multiply_V_add_SVV(const EvalVector* A,
                                    const double& beta,
                                    const EvalVector* C,
                                    const EvalVector* D)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      const double* const Cx = C->start();
      const double* const Dx = D->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= Ax[i];
          x[i] += beta*Cx[i]*Dx[i];
        }
      addFlops(4*n);
    }

  if (shadowOps())
    {
      if (beta != 0.0)
        {
          str_ = "(" + str() + "*" + A->str();
        }
      else
        {
          str_ = str() + "*" + A->str();
        }
      if (beta != 0.0)
        {
          str_ += " + ";
          if (beta != 1.0)
            {
              str_ += Teuchos::toString(beta) + "*";
            }
          str_ += C->str() + "*" + D->str();
        }
    }
}

void EvalVector::multiply_V_add_SV(const EvalVector* A,
                                   const double& beta,
                                   const EvalVector* C)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      const double* const Cx = C->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= Ax[i];
          x[i] += beta*Cx[i];
        }
      addFlops(3*n);
    }

  if (shadowOps())
    {
      str_ = "(" + str() + "*" + A->str() + " + " + Teuchos::toString(beta)
        + "*" + C->str() + ")";
    }
}

void EvalVector::multiply_VV(const EvalVector* A,
                             const EvalVector* B)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      const double* const Bx = B->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= Ax[i]*Bx[i];
        }
      addFlops(2*n);
    }

  if (shadowOps())
    {
      str_ = str() + "*" + A->str() + "*" + B->str();
    }
}


void EvalVector::multiply_SV(const double& alpha,
                             const EvalVector* B)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      
      const double* const Bx = B->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= alpha*Bx[i];
        }
      addFlops(2*n);
    }

  if (shadowOps())
    {
      str_ = Teuchos::toString(alpha) + "*" + str_ + "*" + B->str();
    }
}

void EvalVector::multiply_S(const double& alpha)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  if (n > 0)
    {
      double* const x = start();
      
      for (int i=0; i<n; i++)
        {
          x[i] *= alpha;
        }
      addFlops(n);
    }

  if (shadowOps())
    {
      str_ = Teuchos::toString(alpha) + "*" + str_;
    }
}


void EvalVector::setTo_S_add_SVV(const double& alpha,
                                 const double& beta,
                                 const EvalVector* C,
                                 const EvalVector* D)
{
  //TimeMonitor t(evalVecTimer());
  int n = C->data_->size();
  resize(n);

  if (n > 0)
    {
      double* const x = start();
      const double* const Cx = C->start();
      const double* const Dx = D->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] = alpha + beta*Cx[i]*Dx[i];
        }
      addFlops(3*n);
    }

  if (shadowOps())
    {
      str_ = "(" + Teuchos::toString(alpha) + " + " 
        + Teuchos::toString(beta) + "*" 
        + C->str() + "*" + D->str() + ")";
    }
}

void EvalVector::setTo_S_add_VV(const double& alpha,
                                const EvalVector* B,
                                const EvalVector* C)
{
  //TimeMonitor t(evalVecTimer());
  int n = B->data_->size();
  resize(n);

  if (n > 0)
    {
      double* const x = start();
      const double* const Bx = B->start();
      const double* const Cx = C->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] = alpha + Bx[i]*Cx[i];
        }
      addFlops(2*n);
    }

  if (shadowOps())
    {
      str_ = "(" + Teuchos::toString(alpha) + " + " 
        + B->str() + "*" + C->str() + ")";
    }
}

void EvalVector::setTo_S_add_SV(const double& alpha,
                                const double& beta,
                                const EvalVector* C)
{
  //TimeMonitor t(evalVecTimer());
  int n = C->data_->size();
  resize(n);

  if (n > 0)
    {
      double* const x = start();
      const double* const Cx = C->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] = alpha + beta*Cx[i];
        }
      addFlops(2*n);
    }

  if (shadowOps())
    {
      str_ = "(" + Teuchos::toString(alpha) + " + " 
        + Teuchos::toString(beta) + "*" 
        + C->str() + ")";
    }
}



void EvalVector::setTo_S_add_V(const double& alpha,
                               const EvalVector* B)
{
  //TimeMonitor t(evalVecTimer());
  int n = B->data_->size();
  resize(n);

  if (n > 0)
    {
      double* const x = start();
      const double* const Bx = B->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] = alpha + Bx[i];
        }
      addFlops(n);
    }

  if (shadowOps())
    {
      str_ = "(" + Teuchos::toString(alpha) + " + " 
        + B->str() + ")";
    }
}


void EvalVector::setTo_V(const EvalVector* A)
{
  //TimeMonitor t(evalVecTimer());
  int n = A->data_->size();
  resize(n);

  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      
      for (int i=0; i<n; i++)
        {
          x[i] = Ax[i];
        }
    }

  if (shadowOps())
    {
      str_ = A->str();
    }
}



void EvalVector::setTo_SVV(const double& alpha,
                           const EvalVector* B,
                           const EvalVector* C)
{
  //TimeMonitor t(evalVecTimer());
  int n = B->data_->size();
  resize(n);

  if (n > 0)
    {
      double* const x = start();
      const double* const Bx = B->start();
      const double* const Cx = C->start();

      
      for (int i=0; i<n; i++)
        {
          x[i] = alpha*Bx[i]*Cx[i];
        }
      addFlops(2*n);
    }

  if (shadowOps())
    {
      str_ = Teuchos::toString(alpha) + "*" 
        + B->str() + "*" + C->str();
    }
}

void EvalVector::setTo_VV(const EvalVector* A,
                          const EvalVector* B)
{
  //TimeMonitor t(evalVecTimer());
  int n = A->data_->size();
  resize(n);

  if (n > 0)
    {
      double* const x = start();
      const double* const Ax = A->start();
      const double* const Bx = B->start();

      
      for (int i=0; i<n; i++)
        {
          x[i] = Ax[i]*Bx[i];
        }

      addFlops(n);
    }

  if (shadowOps())
    {
      str_ = A->str() + "*" + B->str();
    }
}

void EvalVector::setTo_SV(const double& alpha,
                          const EvalVector* B)
{
  //TimeMonitor t(evalVecTimer());
  int n = B->data_->size();
  resize(n);

  if (n > 0)
    {
      double* const x = start();
      const double* const Bx = B->start();

      
      for (int i=0; i<n; i++)
        {
          x[i] = alpha*Bx[i];
        }
      addFlops(n);
    }

  if (shadowOps())
    {
      str_ = Teuchos::toString(alpha) + "*" 
        + B->str();
    }
}

void EvalVector::applyUnaryOperator(const UnaryFunctor* func, 
                                    Array<RCP<EvalVector> >& opDerivs)
{
  //TimeMonitor t(evalVecTimer());
  int n = data_->size();

  int order = opDerivs.size()-1;

  TEUCHOS_TEST_FOR_EXCEPTION(order < 0 || order > 2, std::runtime_error,
                     "illegal order=" << order << " in "
                     "EvalVector::applyUnaryOperator()");

  opDerivs[0] = s_->popVector();
  opDerivs[0]->resize(n);
  if (order > 0)
    {
      opDerivs[1] = s_->popVector();
      opDerivs[1]->resize(n);
    }
  if (order > 1)
    {
      opDerivs[2] = s_->popVector();
      opDerivs[2]->resize(n);
    }
  
  if (n > 0)
    {
      double* const x = start();
      double* const f = opDerivs[0]->start();
      if (order==0)
        {
          func->eval0(x, n, f);
        }
      else if (order==1)
        {
          double* const df = opDerivs[1]->start();
          func->eval1(x, n, f, df);
        }
      else 
        {
          double* const df = opDerivs[1]->start();
          double* const df2 = opDerivs[2]->start();
          func->eval2(x, n, f, df, df2);
        }

    }
  
  if (shadowOps())
    {
      opDerivs[0]->setString(func->name() + "(" + str() + ")");
      if (order > 0)
        {
          opDerivs[1]->setString(func->name() + "'(" + str() + ")");
        }
      if (order > 1)
        {
          opDerivs[2]->setString(func->name() + "\"(" + str() + ")");
        }
    }
}

void EvalVector::print(std::ostream& os) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(shadowOps() && str_.size()==0, std::runtime_error, "empty eval vector result string!");
  os << str_;

  if (data_->size() > 0)
    {
      os << ", " << *data_;
    }
}








