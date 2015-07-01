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
/*
 * SundanceSurfQuadrature.hpp
 *
 *  Created on: Oct 24, 2011
 *      Author: benk
 */

#ifndef SUNDANCE_SURFQUADRATURE_H
#define SUNDANCE_SURFQUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceQuadratureFamily.hpp"

namespace Sundance
{
using namespace Teuchos;


/**
 * The surface integral class. The surface integral is build from triangles.
 * The maximum number of triangles is 4. <br>
 * IMPORTANT: this quadrature class should only be used for Surface Integrals in 3D with Brick cells
 */
class SurfQuadrature : public QuadratureFamilyBase
{
public:
  /** */
  SurfQuadrature( const QuadratureFamily& quad );

  /** */
  virtual ~SurfQuadrature(){;}


  /** */
  virtual XMLObject toXML() const ;

  /** Describable interface */
  virtual std::string description() const
    {return "SurfQuadrature[order=" + Teuchos::toString(order())
        +  "]";}

  /* handleable boilerplate */
  GET_RCP(QuadratureFamilyStub);

  /** return the maximal number of line segments inside one cell */
  static int getNrMaxTrianglePerCell() { return 4; }

protected:

  /** for surface integral integrals only this method should be used */
  virtual void getQuadRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;


  /** for surface integral integrals only this method should be used */
  virtual void getTriangleRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

private:

  /** The quadrature which will be used for the surface integration */
  const QuadratureFamily& quad_;

};
}


#endif
