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
 * SundanceCurveIntagralCalc.hpp
 *
 *  Created on: Sep 2, 2010
 *      Author: benk
 */

#ifndef SUNDANCECURVEINTEGRALCALC_HPP_
#define SUNDANCECURVEINTEGRALCALC_HPP_

#include "SundanceParametrizedCurve.hpp"
#include "SundanceOut.hpp"
#include "SundancePoint.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceCellType.hpp"
#include "SundanceMesh.hpp"

namespace Sundance {

/** Class to compute the intersection/quadrature points of a cell with a curve in 2/3D */
class CurveIntegralCalc {
public:

	/** Empty Ctor*/
	CurveIntegralCalc();

	virtual ~CurveIntegralCalc() {;}

    /**
         [in]  maxCellType <br>
         [in]  const Array<Point>& cellPoints (the physical points of the cell) <br>
         [in]  paramCurve <br>
         [in]  Quadraturefamily , quadrature 1D for curve, or 2D for surface <br>
         [out] Array<Point>& curvePoints <br>
         [out] Array<Point>& curveDerivs  <br>
         [out] Array<Point>& curveNormals <br> */

    static void getCurveQuadPoints(CellType  maxCellType ,
    		                       int maxCellLID ,
    		                       const Mesh& mesh ,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );

    /** this method shows if special methods for polygon can be used */
    static bool usePolygoneCurveQuadrature(
    		const CellType& cellType ,
    		const Mesh& mesh ,
    		const ParametrizedCurve& paramCurve);

    /** test if we are in 3D for the Brick surface quadrature */
    static bool use3DSurfQuadrature(
    		const CellType& cellType ,
    		const Mesh& mesh ){
    	if ( (mesh.cellType(mesh.spatialDim()) == BrickCell) && (cellType == BrickCell) ) { return true; }
    	else {return false;}
    }

private :

    /** the simplest method which considers only the intersection points on the
     * facets, and inside the element we use */
    static void getCurveQuadPoints_line(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );

    /** projects the point from 2D surf to 3D points of the element (needed for surf integrals)*/
    static void get3DRealCoordsOnSurf(const Point &refP ,
    		                          const Array<Point>& cellPoints,
    		                          const Array<int> &triangles ,
    		                          const int nrTriag ,
    		                          const Array<Point> edgeIntersectPoint,
    		                          int &triagIndex ,
    		                          Point &realPoint);


    /** polygons in 1D curves, uses the information from the polygon */
    static void getCurveQuadPoints_polygon(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );

    /** special method for 3D surface quadrature similar to the one in getCurveQuadPoints_line ,
     * but it uses special Quadrature class, up to 4 triangles per brick cell*/
    static void getSurfQuadPoints(
                                   int  maxCellLID ,
                                   CellType  maxCellType ,
                                   const Mesh& mesh ,
                                   const Array<Point>& cellPoints,
								   const ParametrizedCurve& paramCurve,
								   const QuadratureFamily& quad ,
								   Array<Point>& curvePoints ,
								   Array<Point>& curveDerivs ,
								   Array<Point>& curveNormals );
};

}

#endif /* SUNDANCECURVEINTEGRALCALC_HPP_ */
