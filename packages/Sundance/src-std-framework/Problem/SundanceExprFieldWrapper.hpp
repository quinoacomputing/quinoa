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

#ifndef SUNDANCE_EXPRFIELDWRAPPER_H
#define SUNDANCE_EXPRFIELDWRAPPER_H


#include "SundanceDefs.hpp"
#include "PlayaHandleable.hpp"
#include "SundanceFieldBase.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{
using namespace Teuchos;

    
/**
 *
 */
class ExprFieldWrapper : public FieldBase
{
public:
  /** */
  ExprFieldWrapper(const Expr& expr) ;

  /** virtual dtor */
  virtual ~ExprFieldWrapper(){;}

  /** */
  virtual double getData(int cellDim, int cellID, int elem) const ;

  /** */
  virtual bool isDefined(int cellDim, int cellID, int elem) const ;

  /** */
  virtual int numElems() const {return Expr_size_;}

  /** */
  virtual bool isPointData() const {return isPointData_;}

  /* */
  GET_RCP(FieldBase);
  /**
   * Return the cell filter on which this field is defined 
   */
  virtual const CellFilter& domain() const 
    { // here we return only the first element ()
      return discreteSpace_.cellFilters(indices_[0][0]);
    }

public:
  Expr expr_;

  // this field should be unique for all the variables
  const DiscreteFunctionData* df_;

  DiscreteSpace discreteSpace_;

  // ---- this field can be expicitly asked -----
  //RCP<DOFMapBase> map_;

  Array< Array<int> > indices_;

  int Expr_size_;

  bool isPointData_;
};
}



#endif
