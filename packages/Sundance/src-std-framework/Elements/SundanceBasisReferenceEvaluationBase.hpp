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

#ifndef SUNDANCE_BASISREFERENCEEVALUATIONBASE_H
#define SUNDANCE_BASISREFERENCEEVALUATIONBASE_H

#include "SundanceDefs.hpp"
#include "SundanceCellType.hpp"
#include "Teuchos_Array.hpp"


namespace Sundance {class Point;}
namespace Sundance {class SpatialDerivSpecifier;}
namespace Sundance {class MultiIndex;}

namespace Sundance 
{
using namespace Teuchos;

using Teuchos::Array;


/** 
 * Abstract interface for evaluation of basis functions and their
 * spatial derivatives on reference cells.
 */
class BasisReferenceEvaluationBase
{
public:

  /** \brief Evaluate the basis functions (or some mixed spatial derivative of
   * the basis functions) for an array of points on the "reference cell" for a
   * given cell type.
   *
   * \param  cellType
   *           [in] The type of cell on which the basis is currently being
   *           evaluated. 
   * \param  pts
   *           [in] Array of points on the reference cell (or master cell)
   *           where the basis functions are to be computed. 
   * \param  deriv
   *           [in] Specification of which differential operator is
   *           to be applied to the basis functions.
   * \param  result
   *           [out] On output, gives a triply nested array which contain
   *           the basis functions (or the requested basis function
   *           derivatives) evaluated at the given points <tt>pts</tt>.  The
   *           size of the outer array <tt>results</tt> is either zero
   *           or spatialDim, depending on whether this is a scalar or
   *           vector basis, respectively. The size of the next
   *           array level is equal to the number of evaluation points. 
   *           Finally, the size of the innermost array level is equal to
   *           the number of DOFs visible from the given cell type.
x   *           Specifically,
   *           \code 
   *           results[k][pointIndex][basisIndex] 
   *           \endcode gives the value
   *           of the spatial derivative of the \f$k\f$-th component of
   *           \f[\frac{\partial^{d_x+d_y+d_z}}{\partial x^{d_x} \partial
   *           y^{d_y} \partial z^{d_z}}\psi_i(x,y,z)\f],
   *           where \f$d_x\f$ =
   *           <tt>deriv[0]</tt>, \f$d_y\f$ = <tt>deriv[1]</tt> (in 2D or 3D)
   *           and \f$d_Z\f$
   *           = <tt>deriv[2]</tt> (in 3D) at the point <tt>pointIndex</tt> 
   *           (where
   *           <tt>0 <= pointIndex < pts.size()</tt>) for the basis function
   *           \f$i\f$ = <tt>basisIndex</tt> (where <tt>0 <= basisIndex <
   *           mapStructure.numBasisChunks()</tt>). 
   */
  virtual void refEval(
    const CellType& cellType,
    const Array<Point>& pts,
    const SpatialDerivSpecifier& deriv,
    Array<Array<Array<double> > >& result,
    int verbosity = 0
    ) const = 0 ;  

  /**
   * Computes the constraints for DoFs which are on hanging elements. <br>
   * The child cell is which constrains the hanging local DoF, and the parent cell is needed to find the
   * global DoFs.
   * @param indexInParent  [in] each (child)cell which has one hanging node, has a parent cell which
   * has one facet, where there are global DoFs
   * @param maxCellDim     [in] the dimension of the maximal cell
   * @param maxNrChild     [in] how many children has one parent cell, this tells us if we have trisection or bisection
   * @param facetDim       [in] the hanging element dimension which is a facet of the child cell
   * @param facetIndex     [in] the hanging element facet index in the child cell
   * @param nodeIndex      [in] one element (e.g edge) might have more than one DoF, specify which DoF on the
   * elemnt do we want to constrain
   * @param localDoFs      [out] the local DoFs in the parent cell which contribute to the hanging DoF
   * @param parentFacetDim [out] the facet dimension where the local DoFs (localDoFs) are
   * @param parentFacetIndex[out] the facet index where the local DoFs (localDoFs) are
   * @param parentFacetNode[out] the facet node where the local DoFs (localDoFs) is (e.g.: one edge might have 2 DoFs in P3)
   * @param coefs          [out] the belonging coefficients to the parents local DoF
   *
   */
  virtual void  getConstrainsForHNDoF(
  				    const int indexInParent,
  				    const int maxCellDim,
  				    const int maxNrChild,
  				    const int facetDim,
  				    const int facetIndex,
  				    const int nodeIndex,
  				    Array<int>& localDoFs,
  			    	Array<int>& parentFacetDim,
  			    	Array<int>& parentFacetIndex,
  			    	Array<int>& parentFacetNode,
  				    Array<double>& coefs
  				    ) {};

};

}


#endif
