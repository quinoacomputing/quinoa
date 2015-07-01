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

#ifndef SUNDANCE_EVALVECTOR_H
#define SUNDANCE_EVALVECTOR_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceUnaryFunctor.hpp"
#include "SundanceNoncopyable.hpp"


namespace Sundance
{
using namespace Teuchos;
class TempStack;
/**
 *
 */
class EvalVector : public Noncopyable,
                   public ObjectWithClassVerbosity<EvalVector>
{
  friend class EvalManager;
  friend class TempStack;

private:

  /** */
  EvalVector(TempStack* s);

  /** */
  EvalVector(TempStack* s, const RCP<Array<double> >& data,
    const std::string& str);


public:
  /** 
   * EvalVector has a nontrivial destructor. Upon destruction, 
   * the vector's underlying data object is not destroyed, but rather
   * is put back on the stack of temporary vectors. 
   */
  ~EvalVector();

  /** \name Mathematical operations */
  //@{

  /** */
  void add_SV(const double& alpha, 
    const EvalVector* B) ;

  /**
   * Perform the operation 
   * \f[ 
   * this = this + alpha*B*C
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void add_SVV(const double& alpha,
    const EvalVector* B,
    const EvalVector* C) ;

  /** */
  void add_V(const EvalVector* A) ;

  /** */
  void add_S(const double& alpha);

  /**
   * Perform the operation 
   * \f[ 
   * this = this + A*B
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void add_VV(const EvalVector* A,
    const EvalVector* B) ;


  /**  
   * Perform a scaled addition with another vector,
   * \f[ 
   * this = \alpha this + \beta C
   * \f]
   * The operation is done in-place, overwriting the old values of the
   * vector. 
   */
  void multiply_S_add_SV(const double& alpha, 
    const double& beta,
    const EvalVector* C) ;

  /** Scale and add a constant to this vector. 
   * The operation is done in-place, overwriting the old values of
   * the vector. Each element x[i] is updated as:
   * \f[
   * this = alpha * this + beta
   * \f]
   */
  void multiply_S_add_S(const double& alpha,
    const double& beta) ;

  /**
   * Perform the operation 
   * \f[ 
   * this = this*A + B*C*D
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void multiply_V_add_VVV(const EvalVector* A,
    const EvalVector* B,
    const EvalVector* C,
    const EvalVector* D) ;


  /**
   * Perform the operation 
   * \f[ 
   * this = this*A + beta*C*D
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void multiply_V_add_SVV(const EvalVector* A,
    const double& beta,
    const EvalVector* C,
    const EvalVector* D) ;

  /**
   * Perform the operation 
   * \f[ 
   * this = this*A + beta*C
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void multiply_V_add_SV(const EvalVector* A,
    const double& beta,
    const EvalVector* C) ;

  /**
   * Perform the operation 
   * \f[ 
   * this = this*A*B
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void multiply_VV(const EvalVector* A,
    const EvalVector* B) ;

  /**
   * Perform the operation 
   * \f[ 
   * this = this*alpha*B
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void multiply_SV(const double& alpha,
    const EvalVector* B) ;

  /**
   * Perform the operation 
   * \f[ 
   * this = this*A
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void multiply_V(const EvalVector* A) ;

  /**
   * Perform the operation 
   * \f[ 
   * this = this*alpha
   * \f]
   * which shows up in the chain rule expansion of a second derivative.
   * 
   */
  void multiply_S(const double& alpha) ;

  /**
   *
   */
  void setTo_S_add_SVV(const double& alpha,
    const double& beta,
    const EvalVector* C,
    const EvalVector* D);

  /**
   *
   */
  void setTo_S_add_VV(const double& alpha,
    const EvalVector* B,
    const EvalVector* C);

  /**
   *
   */
  void setTo_S_add_SV(const double& alpha,
    const double& beta,
    const EvalVector* C);

  /** 
   *
   */
  void setTo_S_add_V(const double& alpha,
    const EvalVector* B);


  /**
   *
   */
  void setTo_V(const EvalVector* A);

  /**
   *
   */
  void setTo_VV(const EvalVector* A,
    const EvalVector* B);

  /**
   *
   */
  void setTo_SV(const double& alpha,
    const EvalVector* B);

  /**
   *
   */
  void setTo_SVV(const double& alpha,
    const EvalVector* B,
    const EvalVector* C);

      


  /**
   * Set every element to a constant value
   */
  void setToConstant(const double& alpha) ;

  /** 
   * Apply a unary function
   */
  void applyUnaryOperator(const UnaryFunctor* func, 
    Array<RCP<EvalVector> >& opDerivs);
      
      
  /** */
  RCP<EvalVector> clone() const ;

  /** */
  void resize(int n);

  /** */
  int length() const {return data_->size();}
      
  /** */
  void print(std::ostream& os) const ;

  /** */
  const double * start() const {return &((*data_)[0]);}

  /** */
  double * start() {return &((*data_)[0]);}

  const std::string& str() const {return str_;}

  void setString(const std::string& str) {str_ = str;}

  inline static bool& shadowOps() {static bool rtn = false; return rtn;}

  bool isValid() const {return data_.get() != 0 && s_ != 0;}
  //@}

      

  inline static double& totalFlops() {static double rtn = 0; return rtn;}

private:

  inline static void addFlops(const double& flops) {totalFlops() += flops;}

  mutable TempStack* s_;

  RCP<Array<double> > data_;

  std::string str_;

};



}

namespace std
{
inline ostream& operator<<(std::ostream& os, 
  const Sundance::EvalVector& vec)
{
  vec.print(os);
  return os;
}
}

                  
#endif
