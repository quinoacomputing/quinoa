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

#ifndef SUNDANCE_DISCRETEFUNCTION_H
#define SUNDANCE_DISCRETEFUNCTION_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceDiscreteFunctionData.hpp"
#include "SundanceFuncWithBasis.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "PlayaVectorDecl.hpp"

namespace Sundance
{
using namespace Teuchos;
  

/** 
 * DiscreteFunction represents a function that is discretized
 * on a finite-element space.
 */
class DiscreteFunction : public DiscreteFunctionStub,
                         public FuncWithBasis
{
public:
  /** */
  DiscreteFunction(const DiscreteSpace& space, const std::string& name="");

  /** */
  DiscreteFunction(const DiscreteSpace& space, const Vector<double>& vec, 
    const std::string& name="");

  /** */
  DiscreteFunction(const DiscreteSpace& space, const double& constantValue,
    const std::string& name="");
  /** */
  DiscreteFunction(const DiscreteSpace& space, const Array<string>& names);

  /** */
  DiscreteFunction(const DiscreteSpace& space, const Vector<double>& vec, 
    const Array<string>& names);

  /** */
  DiscreteFunction(const DiscreteSpace& space, const double& constantValue,
    const Array<string>& name);

  /** */
  static const DiscreteFunction* discFunc(const Expr& expr);


  /** */
  static DiscreteFunction* discFunc(Expr& expr);

  /** */
  void updateGhosts() const ;

  /** */
  void setVector(const Vector<double>& vec);

  /** */
  const Vector<double>& getVector() const 
    {return data_->getVector();}

  /** */
  const DiscreteSpace& discreteSpace() const 
    {return data_->discreteSpace();}

  /** */
  const Mesh& mesh() const {return discreteSpace().mesh();}

  /** */
  const RCP<DOFMapBase>& map() const {return discreteSpace().map();}


  RCP<GhostView<double> >  ghostView() const 
    {return data_->ghostView();}

  const DiscreteFunctionData* data() const {return data_.get();}


  /** virtual destructor */
  virtual ~DiscreteFunction() {;}

  /* boilerplate */
  GET_RCP(ExprBase);


  /** */
  RCP<const MapStructure> getLocalValues(int cellDim, 
    const Array<int>& cellLID,
    Array<Array<double> >& localValues) const ;


private:
  /** */
  RCP<DiscreteFuncDataStub> getRCP(DiscreteFunctionData* ptr);

  RCP<DiscreteFunctionData> data_;

};


/** \relates DiscreteFunction
 * Replace the vector in oldVals with the vector from newVals.
 */
void updateDiscreteFunction(const Expr& newVals, Expr oldVals);


/** \relates DiscreteFunction
 * Make a copy of the discrete function u0. The copy will have a shallow
 * copy of u0's space, and a deep copy of u0's vector. 
 */
Expr copyDiscreteFunction(const Expr& u0, const string& name = "");


/** \relates DiscreteFunction
 * Add a vector v to the vector underlying the discrete function u.
 */
void addVecToDiscreteFunction(Expr u, const Vector<double>& v);

/** \relates DiscreteFunction
 * Get a shallow copy of the vector underlying a discrete function 
 */
Vector<double> getDiscreteFunctionVector(const Expr& u);


/** \relates DiscreteFunction
 * Set the vector underlying a discrete function 
 */
void setDiscreteFunctionVector(Expr u, const Vector<double>& v);


/** \relates DiscreteFunction
 * Get the mesh underlying a discrete function 
 */
Mesh getDiscreteFunctionMesh(const Expr& u);


/** \relates DiscreteFunction
 * Get the discrete space on which a discrete function is defined 
 */
DiscreteSpace getDiscreteSpace(const Expr& u);


}



#endif
