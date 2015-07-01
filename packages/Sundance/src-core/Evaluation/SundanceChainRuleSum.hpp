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

#ifndef SUNDANCE_CHAINRULESUM_H
#define SUNDANCE_CHAINRULESUM_H

#include "SundanceDefs.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvaluator.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance 
{

using namespace Teuchos;


/** */
class IndexPair 
{
public:
  /** */
  IndexPair(int argIndex, int valueIndex)
    : argIndex_(argIndex), valueIndex_(valueIndex) {;}
      
  /** */
  int argIndex() const {return argIndex_;}
      
  /** */
  int valueIndex() const {return valueIndex_;}
      
private:
  int argIndex_;
  int valueIndex_;
};

/** */
class DerivProduct
{
public:
  /** */
  DerivProduct() : coeff_(1.0), constants_(), variables_() {}
  /** */
  DerivProduct(const double& coeff) : coeff_(coeff), constants_(), variables_() {}

  /** */
  void addConstantFactor(const IndexPair& p) {constants_.append(p);}

  /** */
  void addVariableFactor(const IndexPair& p) {variables_.append(p);}

  /** */
  bool isConstant() const {return numVariables()==0;}

  /** */
  int numConstants() const {return constants_.size();}

  /** */
  int numVariables() const {return variables_.size();}

  /** */
  const double& coeff() const {return coeff_;}

  /** */
  const IndexPair& constant(int i) const {return constants_[i];}

  /** */
  const IndexPair& variable(int i) const {return variables_[i];}
        
private:

  double coeff_;

  Array<IndexPair> constants_;

  Array<IndexPair> variables_;
};


/** */
class ChainRuleSum : public ObjectWithClassVerbosity<Evaluator>
{
public:
  /** */
  ChainRuleSum(const MultipleDeriv& md, 
    int resultIndex,
    bool resultIsConstant);

  /** */
  void addTerm(int argDerivIndex, 
    bool argDerivIsConstant,
    const Array<DerivProduct>& sum);


  /** */
  void evalConstant(const EvalManager& mgr,
    const Array<RCP<Array<double> > >& constantArgResults,
    const Array<double>& constantArgDerivs,
    double& constResult) const ;

  /** */
  void evalVar(const EvalManager& mgr,
    const Array<RCP<Array<double> > >& constantArgResults,
    const Array<RCP<Array<RCP<EvalVector> > > > & vArgResults,
    const Array<double>& constantArgDerivs,
    const Array<RCP<EvalVector> >& varArgDerivs,
    RCP<EvalVector>& varResult) const ;

  /** */
  int resultIndex() const {return resultIndex_;}

  /** */
  bool resultIsConstant() const {return resultIsConstant_;}

  /** */
  int numTerms() const {return terms_.size();}

  /** */
  bool argDerivIsConstant(int i) const {return argDerivIsConstant_[i];}

  /** */
  int argDerivIndex(int i) const {return argDerivIndex_[i];}

  /** */
  const Array<DerivProduct>& terms(int i) const {return terms_[i];}

  /** */
  const MultipleDeriv& deriv() const {return md_;}

private:
  MultipleDeriv md_;
  int resultIndex_;
  bool resultIsConstant_;

  Array<int> argDerivIndex_;
  Array<int> argDerivIsConstant_;
  Array<Array<DerivProduct> > terms_;
};

}
               
#endif
