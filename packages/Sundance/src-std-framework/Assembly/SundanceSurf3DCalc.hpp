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
 * SundanceSurf3DCalc.hpp
 *
 *  Created on: Oct 25, 2011
 *      Author: benk
 */

#ifndef SUNDANCESURF3DCALC_HPP_
#define SUNDANCESURF3DCALC_HPP_

#include "SundanceParametrizedCurve.hpp"
#include "SundanceOut.hpp"
#include "SundancePoint.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceCellType.hpp"
#include "SundanceMesh.hpp"

namespace Sundance {

/** This class contains the computational methods which are necessary for the 3D surface and cut-cell
 * calculations.*/
class SundanceSurf3DCalc {
public:

	/** empty Ctor */
	SundanceSurf3DCalc() {;}

	/** empty Dtor */
	virtual ~SundanceSurf3DCalc() {;}

	/** method to compute the intersection surface between one curve and a brick cell.
	 * The resulted surface is represented by a triangle surface.
	 * @param maxCellType [IN]
	 * @param maxCellLID [IN]
	 * @param mesh [IN]
	 * @param paramCurve [IN] the surface object
	 * @param intersectPoints [OUT] the array with the intersection points in the reference coordinate
	 * @param brickPoints [OUT] 12 points of the brick cell
	 * @param edgeIndex [OUT] the intersection points are by default on the edges of the brick
	 * @param triangleIndex [OUT] contains the triangles, the size of this vector is 3*nrTriags. */
    static void getCurveQuadPoints(CellType  maxCellType ,
    		                       int maxCellLID ,
    		                       const Mesh& mesh ,
								   const ParametrizedCurve& paramCurve,
								   const Array<Point>& brickPoints,
								   Array<Point>& intersectPoints ,
								   Array<int>& edgeIndex,
								   Array<int>& triangleIndex);


	static const int edegIndex[12][2];

	static const int faceEdges[6][4];

	static const int edgeFaces[12][2];

	static const int edgeNeighboredges[12][12];

};


}
#endif /* SUNDANCESURF3DCALC_HPP_ */
