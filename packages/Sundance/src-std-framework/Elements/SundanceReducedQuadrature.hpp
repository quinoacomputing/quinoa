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

#ifndef SUNDANCE_REDUCED_QUADRATURE_H
#define SUNDANCE_REDUCED_QUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"

namespace Sundance
{
using namespace Teuchos;


/** 
 * Reduced quadrature regards all coefficients as constants, using the midpoint value on
 * each element as representative of the entire element. The product of test and trial functions
 * is then integrated using reference integration.
 */
class ReducedQuadrature : public QuadratureFamilyBase
{
public:
  /** */
  ReducedQuadrature();

  /** */
  virtual ~ReducedQuadrature(){;}


  /** */
  virtual XMLObject toXML() const ;

  /** Describable interface */
  virtual std::string description() const 
    {return "ReducedQuadrature()";}

  /* handleable boilerplate */
  GET_RCP(QuadratureFamilyStub);

  /** */
  virtual int getNumPoints( const CellType &cellType ) const ;


  /** Get the quadrature points and weights for the given cell type */
  virtual void getPoints(const CellType& cellType, 
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** This methos is for the ACI integration */
  virtual void getAdaptedWeights(const CellType& cellType ,
  									 int cellDim,
  	                                 int celLID ,
  	                	             int facetIndex ,
  	                                 const Mesh& mesh ,
  	                                 const ParametrizedCurve& globalCurve ,
  	                                 Array<Point>& quadPoints ,
  	                                 Array<double>& quadWeights ,
  	                                 bool &isCut) const ;

protected:
};
}


#endif
