// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubatureDirectDef.hpp
    \brief  Definition file for the Intrepid2::CubatureDirect class.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid2 {

template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getCubatureData(ArrayPoint  &                cubPoints,
                                                                    ArrayWeight &                cubWeights,
                                                                    const CubatureTemplate *     cubData) const {

  int numCubPoints = getNumPoints();
  int cellDim      = getDimension();
  // check size of cubPoints and cubWeights
  TEUCHOS_TEST_FOR_EXCEPTION( ( ( (int)cubPoints.size() < numCubPoints*cellDim ) || ( (int)cubWeights.size() < numCubPoints ) ),
                      std::out_of_range,
                      ">>> ERROR (CubatureDirect): Insufficient space allocated for cubature points or weights.");

  for (int pointId = 0; pointId < numCubPoints; pointId++) {
    for (int dim = 0; dim < cellDim; dim++) {
      cubPoints(pointId,dim) = cubData->points_[pointId][dim];
    }
    cubWeights(pointId) = cubData->weights_[pointId];
  }
} // end getCubatureData


template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint  & cubPoints,
                                                                ArrayWeight & cubWeights) const {
  getCubatureData( cubPoints, cubWeights, &(exposeCubatureData()[degree_]) );
} // end getCubature

template<class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getCubature(ArrayPoint& cubPoints,
		                                                ArrayWeight& cubWeights,
                                                                ArrayPoint& cellCoords) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (CubatureDirect): Cubature defined in reference space calling method for physical space cubature.");
}



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getNumPoints() const {
  return exposeCubatureData()[degree_].numPoints_;
} // end getNumPoints



template <class Scalar, class ArrayPoint, class ArrayWeight>
int CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getDimension() const {
  return dimension_;
} // end dimension



template <class Scalar, class ArrayPoint, class ArrayWeight>
void CubatureDirect<Scalar,ArrayPoint,ArrayWeight>::getAccuracy(std::vector<int> & accuracy) const {
  accuracy.assign(1, degree_);
} // end getAccuracy

    
} // end namespace Intrepid2
