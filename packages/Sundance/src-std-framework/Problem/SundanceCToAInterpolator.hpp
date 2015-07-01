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

#ifndef SUNDANCE_CTOAINTERPOLATOR_H
#define SUNDANCE_CTOAINTERPOLATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceAToCPointLocator.hpp"

namespace Sundance
{
  using namespace Teuchos;
  using namespace Playa;

  /**
   * CToAInterpolator interpolates a discrete function at particle positions
   *
   * Note: not tested in parallel.
   */
  class CToAInterpolator
  {
  public:
    /** */
    CToAInterpolator(const AToCPointLocator& locator,
                     const Expr& field);

    /** */
    void interpolate(const Teuchos::Array<double>& positions,
                     Teuchos::Array<double>& results) const ;

   

    /** */
    void interpolate(const std::vector<double>& positions,
      std::vector<double>& results) const 
      {
        Teuchos::Array<double> in(positions.size());
        for (int i=0; i<in.size(); i++) in[i] = positions[i];
        Teuchos::Array<double> out;
        interpolate(in, out);
        for (int i=0; i<out.size(); i++) results[i] = out[i];
      }

    /** */
    void updateField(const Expr& field) ;


  private:

    int dim_;
    int nFacets_;
    int rangeDim_;
    RCP<Array<double> > elemToVecValuesMap_;
    AToCPointLocator locator_;
  };
}


#endif
