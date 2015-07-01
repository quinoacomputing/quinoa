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

#ifndef SUNDANCE_CELLVECTOREXPR_H
#define SUNDANCE_CELLVECTOREXPR_H

#include "SundanceEvaluatorFactory.hpp"
#include "SundanceCellVectorEvaluator.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Sundance;

/** */
enum CellVectorExprType {CellNormalVector, CellTangentSpace};

/**
 * Expression that returns a component of a cell-determined vector for
 * each cell on 
 * which it is evaluated. This makes sense only for boundary cells. 
 */
class CellVectorExpr
  : public EvaluatableExpr,
    public GenericEvaluatorFactory<CellVectorExpr, CellVectorEvaluator>
{
public:
  /** */
  CellVectorExpr(
    int normalComponentIndex, 
    int dimension,
    const std::string& name
    );

  /** */
  CellVectorExpr(
    int tangentBasisIndex, 
    int tangentComponentIndex,
    int dimension,
    const std::string& name
    );


    
  /** */
  virtual ~CellVectorExpr() {;}

  /** */
  bool isTangent() const {return type_ == CellTangentSpace;}

  /** */
  bool isNormal() const {return !isTangent();}

  /** */
  int componentIndex() const {return componentIndex_;}

  /** */
  int basisMemberIndex() const {return basisMemberIndex_;}

  /** */
  int dimension() const {return dim_;}

  /** */
  virtual XMLObject toXML() const ;

  const std::string& name() const {return name_;}

  /** Write a simple text description suitable 
   * for output to a terminal */
  virtual std::ostream& toText(std::ostream& os, bool paren) const ;
    

  /** */
  virtual Set<MultipleDeriv> 
  internalFindW(int order, const EvalContext& context) const ;

  /** */
  virtual RCP<ExprBase> getRcp() {return rcp(this);}

  /** Ordering operator for use in transforming exprs to standard form */
  virtual bool lessThan(const ScalarExpr* other) const ;
private:
  std::string name_;
  int dim_;
  CellVectorExprType type_;
  int basisMemberIndex_;
  int componentIndex_;
  
};


/** \relates Expr \relates CellVectorExpr */
Expr CellNormalExpr(int dimension, const std::string& name);

/** \relates Expr \relates CellVectorExpr */
Expr CellTangentExpr(int dimension, const std::string& name);

}

#endif
