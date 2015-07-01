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
 * SundanceHNMapBase.hpp
 *
 *  Created on: Mar 18, 2010
 *      Author: benk
 */

#ifndef SUNDANCEHNMAPBASE_HPP_
#define SUNDANCEHNMAPBASE_HPP_

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"

namespace Sundance
{
using namespace Teuchos;

/**
 * The abstract class which extends the functionalities of the DOF map <br>
 * The only additional functionality is that we have a restriction on the DOFs
 * , with the pre-fill transformations these constraints can be build in into the matrix
 *
 */
class HNDoFMapBase
{
public:

	/** Empty Ctor */
	HNDoFMapBase(const Mesh& mesh, int nFuncs, int setupVerb) : mesh_(mesh){;}

	virtual ~HNDoFMapBase() {;}

	/**
	 * @param cellLID [in] the maxCell LID input
	 * @param funcID  [in] the function ID
	 * @param trafoMatrixSize [in/out]
	 * @param doTransform [out]
	 * @param transfMatrix [out] (we assume that the array is already pre-sized )*/
	  virtual void getTrafoMatrixForCell(
		    int cellLID,
		    int funcID,
		    int& trafoMatrixSize,
		    bool& doTransform,
		    Array<double>& transfMatrix ) const=0;

	/** Function to apply transformation for facets
	 * @param cellDim , the facet dimension
	 * @param cellLID , facet LID
	 * @param facetIndex , facet index in the maxCofacet
	 * @param funcID  [in] the function ID
	 * @param trafoMatrixSize [in/out]
	 * @param doTransform [out]
	 * @param transfMatrix [out] (we assume that the array is already pre-sized )*/
	  virtual void getTrafoMatrixForFacet(
			  int cellDim,
			  int cellLID,
			  int facetIndex,
			  int funcID,
			  int& trafoMatrixSize,
			  bool& doTransform,
			  Array<double>& transfMatrix ) const = 0;



	  /** Function used for plotting for hanging node DOFMaps
	   *  Returns for one hanging node (element) the global DoFs which contribute to that hanging local DoF
	   * @param cellDim [in] the dimension
	   * @param cellLID [in] the LID of the cell
	   * @param funcID  [in] the function ID, (to wchich the DOFs belong)
	   * @param dofs    [out] the global DoF s
	   * @param coefs   [out] the coefficient of each global DoF */
	  virtual void getDOFsForHNCell(
			int cellDim,
			int cellLID,
	        int funcID,
	        Array<int>& dofs ,
	        Array<double>& coefs ) const=0;

	  /** Returns the dimension where the DoF map is defined
	   *  For HN we do transformation only for the maxCell type */
	  int getSpacialMeshDim() const { return mesh_.spatialDim();}

protected:

private:

	  const Mesh mesh_;

};

}


#endif /* SUNDANCEHNMAPBASE_HPP_ */
