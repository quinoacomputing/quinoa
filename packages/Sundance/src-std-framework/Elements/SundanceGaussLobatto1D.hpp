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

#ifndef SUNDANCE_GAUSSLOBATTO1D_H
#define SUNDANCE_GAUSSLOBATTO1D_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  using namespace Teuchos;

    /**
     * Get abscissas and weights for Gauss-Lobatto-Legendre quadrature
     * on a line segment.
     */
    class GaussLobatto1D
    {
    public:
      /** create an n-point rule on the interval [-1, 1] */
      GaussLobatto1D(int n);
      /** create an n-point rule on the interval [a, b] */
      GaussLobatto1D(int n, double a, double b);

      /** return the number of points in the rule */
      int nPoints() const {return nodes_.length();}
      /** get the abscissas */
      const Array<double>& nodes() const {return nodes_;}
      /** get the weights */
      const Array<double>& weights() const {return weights_;}

      static bool unitTest();
    private:
      void computeWeights(int n, double a, double b);

      Array<double> nodes_;
      Array<double> weights_;
    };
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
