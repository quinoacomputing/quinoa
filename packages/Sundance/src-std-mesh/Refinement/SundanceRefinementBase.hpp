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
 * SundanceRefinementBase.hpp
 *
 *  Created on: Apr 29, 2010
 *      Author: benk
 */

#ifndef SUNDANCEREFINEMENTBASE_HPP_
#define SUNDANCEREFINEMENTBASE_HPP_

#include "SundanceDefs.hpp"
#include "PlayaHandleable.hpp"
#include "SundancePoint.hpp"
#include "SundanceParametrizedCurve.hpp"

namespace Sundance {

using namespace Teuchos;

/** define the predicate */

#define REFINE_MESH_(name, code) \
  class name : public RefinementBase, \
               public Playa::Handleable<RefinementBase> \
  { \
  public:\
    name() : RefinementBase(){;}            \
    virtual ~name(){;} \
    virtual int refine( const int cellLevel , \
						 const Point& cellPos , \
						 const Point& cellDimameter) const code \
    GET_RCP(RefinementBase);\
  }

#define REFINE_MESH(name, code) REFINE_MESH_(name, code);


/** define the predicate also with estimation */
#define REFINE_MESH_ESTIMATE_(name, code , codeEst ) \
  class name : public RefinementBase, \
               public Playa::Handleable<RefinementBase> \
  { \
  public:\
    name() : RefinementBase(){;}            \
    virtual ~name(){;} \
    virtual int refine( const int cellLevel , \
						 const Point& cellPos , \
						 const Point& cellDimameter) const code \
	virtual int estimateRefinementLevel( const Point& cellPos , \
						 const Point& cellDimameter) const codeEst \
    GET_RCP(RefinementBase);\
  }

#define REFINE_MESH_ESTIMATE(name, code, codeEst) REFINE_MESH_ESTIMATE_(name, code, codeEst);


/** Base class for mesh refinement , but also to define computational domain
 * (which must not be the same as the mesh domain) <br>
 * It is important to take a look at the refinement protocol defined:<br>
 *  0 -> no action , 1 -> refine , 2 -> coarse  <br>
 *  <br> <br>
 *  */
class RefinementBase  {

public:

	RefinementBase() {;}

	virtual ~RefinementBase() {;}

	/** Function which has to answer the question per Cell if the
	 * cell needs to be refined or not
	 * @param cellLevel [in] refinement level of the actual cell
	 * @param cellPos [in] the position of the cell (the middle point)
	 * @param cellDimameter [in] resolution of the cell in each direction, this could also be CellDiameter
	 * @returns refinementAction [out] , 0 -> no action , 1 -> refine , 2 -> coarse */
	virtual int refine(const int cellLevel ,
			            const Point& cellPos ,
			            const Point& cellDimameter) const { return 1; }

	/** Function meant to be used for load estimation of a given cell <br>
	 * This function has a "dummy" body so the sub classes must not overwrite this function
	 * @param cellPos [in] the position of the cell (the middle point)
	 * @param cellDimameter [in] resolution of the cell in each direction (this parameter might not be used)
	 * @returns averageLevel [out] returns the estimated (average) refinement level of this cell */
	virtual int estimateRefinementLevel(
			             const Point& cellPos,
			             const Point& cellDimameter ) const { return 1; }

};


class GeometryRefinement : public RefinementBase , public Playa::Handleable<RefinementBase>{

public:

	GeometryRefinement(const ParametrizedCurve& curve , int level)
	: RefinementBase() , curve_(curve) , level_(level) {;}

	virtual ~GeometryRefinement() {;}

	/** see super class */
	virtual int refine(const int cellLevel ,
			            const Point& cellPos ,
			            const Point& cellDimameter) const;

	GET_RCP(RefinementBase);

private:

	/** the curve which will be refined at*/
	const ParametrizedCurve& curve_;

	/** the refinement level at the geometry*/
	const int level_;

};
}

#endif /* SUNDANCEREFINEMENTBASE_HPP_ */
