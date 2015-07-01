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
 * SundanceGaussLobattoQuadrature.hpp
 *
 *  Created on: Jan 20, 2011
 *      Author: benk
 */

#ifndef SUNDANCEGAUSSLOBATTOQUADRATURE_HPP_
#define SUNDANCEGAUSSLOBATTOQUADRATURE_HPP_

#include "SundanceDefs.hpp"
#include "PlayaTabs.hpp"
#include "SundanceQuadratureFamilyBase.hpp"

namespace Sundance {

class GaussLobattoQuadrature : public QuadratureFamilyBase {

public:

	/** In this case the */
	GaussLobattoQuadrature(int order);

	/** */
	virtual ~GaussLobattoQuadrature()
	{;}

	/** */
	virtual XMLObject toXML() const;

	/** Describable interface */
	virtual std::string description() const
	{
		return "GaussLobattoQuadrature[order=" + Teuchos::toString(order()) + "]";
	}

	/* handleable boilerplate */
	GET_RCP(QuadratureFamilyStub);

protected:
	/** compute a rule for the reference line cell */
	virtual void getLineRule(Array<Point>& quadPointsL,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference triangle cell */
	virtual void getTriangleRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference quad cell */
	virtual void getQuadRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/** compute a rule for the reference tet cell */
	virtual void getTetRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const ;

	/** compute a rule for the reference brick cell */
	virtual void getBrickRule(Array<Point>& quadPoints,
			Array<double>& quadWeights) const;

	/**
	 * Compute adapted weights according to curve
	 *
	 * @param cellType
	 * @param cellDim
	 * @param cellLID
	 * @param facetIndex
	 * @param mesh
	 * @param globalCurve
	 * @param quadPoints
	 * @param quadWeights
	 * @param changedWeights
	 */
	virtual void getAdaptedWeights(const CellType& cellType, int cellDim,
			int cellLID, int facetIndex, const Mesh& mesh,
			const ParametrizedCurve& globalCurve,
			Array<Point>& quadPoints, Array<double>& quadWeights,
			bool &weightsChanged) const;

	/** Get the weights for one quad */
	virtual void getAdaptedQuadWeights(int cellLID, const Mesh& mesh,
			const ParametrizedCurve& globalCurve, Array<Point>& quadPoints,
			Array<double>& quadWeights, bool& weightsChanged) const;

	/** Get the weights for one quad knowing that the curve is a 2D polygon */
	virtual void getAdaptedQuadWeights_polygon(int cellLID, const Mesh& mesh,
			const ParametrizedCurve& globalCurve, Array<Point>& quadPoints,
			Array<double>& quadWeights, bool& weightsChanged) const;

	/** Get the weights for one quad knowing that the curve is a 3D surface */
	virtual void getAdaptedQuadWeights_surf(int cellLID, const Mesh& mesh,
			const ParametrizedCurve& globalCurve, Array<Point>& quadPoints,
			Array<double>& quadWeights, bool& weightsChanged) const;

private:

	/** get the triangle quadrature points for the adaptive integration*/
	void getTriangleQuadPoints(Array<Point>& pnt , Array<double>& weight ) const;

	/**
	 * @param i , the i-th point
	 * @param pnt 1D points where we interpolate
	 * @param x the position where we evaluate */
	inline double evalLagrangePol(int i , Array<Point>& pnt , double x) const{
		int m = pnt.size();
		double p = 1.0;
		for (int j = 0 ; (j <= i-1)&&(j<m) ; j++){
		 p = p *( x - pnt[j][0] )/(pnt[i][0] - pnt[j][0]);
		}
		for (int j = i+1 ; j < m ; j++){
		 p = p *( x - pnt[j][0] )/(pnt[i][0] - pnt[j][0]);
		}
	    return p;
	}

	/** this function applies quadrature to the Lagrange functions <br>
	 * all the coordinates are in the reference coordinates */
	inline void makeInterpolantQuad( double px, double py, double ofx, double ofy,
			                         int nr1DPoints , int nr2DPoints ,
			                         Array<Point>& linePoints ,
			                         Array<Point>& quadPoints ,
			                         Array<double>& pointWeights ,
			                         Array<double>& weightsOut ,
			                         double areFakt ) const {
		int xi,yi;
		double xc,yc, intR , valTmp , valBasisX , valBasisY ;

        //SUNDANCE_MSG3(verb_, "px="<<px<<",py="<<py<<",ofx="<<ofx<<",ofy="<<ofy <<
        //		" , nr2DPoints:" << nr2DPoints << " , pointWeights.size():" <<
        //		pointWeights.size() << " , areFakt:" << areFakt);

		// if the area of integration is not significant then just return
		if ( fabs(ofx*ofy) < 1e-14 ) {
			//SUNDANCE_MSG3(verb_, "return");
			return;
		}

		for (int i = 0 ; i < nr2DPoints ; i++){
			// we should add the integration values, and we not set the values //weightsOut[i] = 0.0;
			yi = i % nr1DPoints; xi = i / nr1DPoints; intR = 0.0;
			// for each quadrature point
			for (int q = 0 ; q < pointWeights.size() ; q++){
				xc = px + quadPoints[q][0]*ofx;
				yc = py + quadPoints[q][1]*ofy;
				valBasisX = evalLagrangePol(xi , linePoints , xc);
				valBasisY = evalLagrangePol(yi , linePoints , yc);
				valTmp = pointWeights[q] * ( valBasisX *  valBasisY);
				intR = intR + valTmp;
				//SUNDANCE_MSG3(verb_, "i="<<i<< " , valBasisX=" << valBasisX  << " , valBasisY=" << valBasisY << " , valBasis="
				//		<< valBasisX*valBasisY <<" , val=" << valTmp*::fabs(ofx*ofy/areFakt) );
				//SUNDANCE_MSG3(verb_, "i="<<i<<",xi="<<xi<<",yi="<<yi<<",q="<<q<<",xc="<<xc<<",yc="<<yc<<
				//		",pointWeights[q]=" << pointWeights[q]);
			}
			intR = intR * ::fabs(ofx*ofy/areFakt);
			weightsOut[i] = weightsOut[i] + intR;
		}
	}


	/** this function applies quadrature to the Lagrange functions <br>
	 * all the coordinates are in the reference coordinates */
	inline void makeInterpolantBrick(int nr1DPoints , int nr3DPoints ,
			                         Array<Point>& linePoints ,
			                         Array<Point>& quadPoints ,
			                         Array<double>& pointWeights ,
			                         Array<double>& weightsOut ,
			                         double areFakt ) const {
		int xi,yi,zi;
		double xc,yc,zc, intR , valTmp , valBasisX , valBasisY , valBasisZ;

		for (int i = 0 ; i < nr3DPoints ; i++){
			zi = i/(nr1DPoints*nr1DPoints); yi = (i%(nr1DPoints*nr1DPoints)) / nr1DPoints;
			xi = (i%(nr1DPoints*nr1DPoints)) % nr1DPoints; intR = 0.0;
			// for each quadrature point
			for (int q = 0 ; q < pointWeights.size() ; q++){
				xc = quadPoints[q][0]; yc = quadPoints[q][1]; zc = quadPoints[q][2];
				valBasisX = evalLagrangePol(xi , linePoints , xc);
				valBasisY = evalLagrangePol(yi , linePoints , yc);
				valBasisZ = evalLagrangePol(zi , linePoints , zc);
				valTmp = pointWeights[q] * ( valBasisX * valBasisY * valBasisZ );
				intR = intR + valTmp;
				//SUNDANCE_MSG3(4, "i="<<i<< " , valBasisX=" << valBasisX  << " , valBasisY=" << valBasisY << " , valBasisZ=" << valBasisZ <<
				//		" , valBasis=" << valBasisX*valBasisY <<" , val=" << valTmp );
				//SUNDANCE_MSG3(4, "i="<<i<<",xi="<<xi<<",yi="<<yi<<",q="<<q<<",xc="<<xc<<",yc="<<yc<<",zc="<<zc<<
				//		",pointWeights["<<q<<"]=" << pointWeights[q]);
			}
			intR = intR * areFakt ;
			weightsOut[i] = weightsOut[i] + intR;
		}
	}



	/** nr of points in 1D, the rest should be tensor product */
	int nrPointin1D_;

	/** the verbosity of the object*/
	int verb_;

	static const int quadsEdgesPoints[4][2];

	/** qvasi information for the 3D cut-cell integral*/
	static const int edge3DProjection[12];
};

}

#endif /* SUNDANCEGAUSSLOBATTOQUADRATURE_HPP_ */
