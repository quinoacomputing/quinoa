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

#ifndef SUNDANCE_BASISDOFTOPOLOGYBASE_H
#define SUNDANCE_BASISDOFTOPOLOGYBASE_H

#include "SundanceDefs.hpp"
#include "SundanceCellType.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Describable.hpp"

namespace Sundance {


using Teuchos::Array;
using Sundance::CellType;

/** 
 * Abstract interface for specification of the topology of degree-of-freedom
 * (DOF) assignments on reference cells in any dimension. Currently,
 * only an enumerated set of cell types are supported (see <tt>CellType</tt>).
 *
 * A function \f$g(x)\f$ defined on a cell is represented as:

 \f[
 g(x) = \sum_{i=0}^{N-1} \bar{g}_i \psi_i(x)
 \f]

 * where \f$x\in\Re^{D}\f$ is the spatial spatial coordinate for spatical dim
 * \f$D\f$ that lives in the cell's domain, \f$\psi_i(x)\f$ is the \f$i\f$th
 * basis function (of order = <tt>order()</tt>), \f$\bar{g}_i\f$ is the
 * (constant) coefficient for the \f$i\f$th basis function, and \f$N\f$ is the
 * number of basis functions.  Therefore, given the coefficients of the basis
 * function on the cell \f$\bar{g}_i\f$, one can compute the value of the
 * function \f$g(x)\f$ on the cell at any point \f$x\f$ in the domain of the
 * cell given the above summation formula.  This interface refers to the
 * coefficients \f$\bar{g}_i\f$ as degrees of freedom (<b>DOFs</b>).
 *
 * This interface allows the specification basis functions and basis
 * coefficients (i.e. DOFs) to be associated with any of the facets of a cell,
 * including the cell itself without restriction.  See the function
 * <tt>getLocalDOFs()</tt> for how this mapping of DOFs to facets for a single
 * function defined on the cell.
 *
 * It is important to note that each cell type, i.e. the values of the enum
 * <tt>CellType</tt>, has a "agreed upon" geometry for the "reference cell"
 * and this geometry must be known by the client of the basis family and the
 * concrete implementations of the basis family.
 *
 * This is an interoperability interface and is therefore
 * deliberately minimalist. 
 */
class BasisDOFTopologyBase :
    public Teuchos::Describable
{
public:

  /** 
   * \brief Inform caller as to whether a given combination 
   * of cell types is supported. 
   *
   * \param maximalCellType
   *         [in] maximal-dimension cell type to which the cell is connected.
   *         For example, a triangle in 3D could be connected to a prism
   *         or to a tetrahedron. 
   *
   * \param cellType 
   *         [in] type of cell for which we want DOF information.
   *
   */
  virtual bool supportsCellTypePair(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const = 0 ;



  /** \brief Get a description of the DOF numbering and distribution scheme
   * for this basis function on the given cell type.
   *
   * \param  cellType
   *           [in] Specification of the cell topology
   * \param  dofs
   *           [out] Array of dof numbering information, to be filled in
   *           during the call.  On output,
   *           <tt>dofs.size()==dimension(cellType)</tt>.  See description of
   *           <tt>dofs</tt> below for more details.
   *
   * The DOF description is returned in the nested array <tt>dofs</tt>, and is
   * to be interpreted as follows: The outer dimension of the description
   * array <tt>dofs.size()</tt> is <tt>cellDim</tt>, where <tt>cellDim</tt> is
   * the spatial dimension of the cell.  The DOFs attached to facets are
   * stored in array entries <tt>dofs[s]</tt> where <tt>s=0...cellDim-1</tt>.
   * Those associated with the cell body are stored in
   * <tt>dofs[cellDim-1]</tt>. For cell dofs attached to facets, the dof
   * <tt>facetDofIndex</tt> associated with the facet <tt>facetIndex</tt> of
   * facet dimension <tt>facetDim</tt> is given by:

   \code
   dofs[facetDim][facetIndex][faceDofIndex] 
   \endcode
   
   * For dofs attached to the cell body, the local DOF within the entire cell
   * is given by dof is given by

   \code
   dofs[cellDim][0][dofIndex]
   \endcode 

   * More specifically:<ul>
   *
   * <li><tt>dof[facetDim].size()</tt> gives the number of facets of the facet
   * dimension <tt>facetDim</tt>, where <tt>0 <= facetDim <= cellDim</tt>
   *
   * <li><tt>dof[facetDim][facetIndex].size()</tt> gives the number of degrees
   * of freedom (DOFs) on the facet <tt>facetIndex</tt> with facet dimension
   * <tt>facetDim</tt>, where <tt>0 <= facetDim <= cellDim</tt> and <tt>0 <=
   * facetIndex < numFacets(cellType,facetDim)</tt>.
   *
   * </ul>
   *
   * For example, the Lagrange basis functions of order 0 through 3 on 2D
   * triangles would have the following dof arrays:

   \verbatim

   Order 0:

   { {}, {}, {{0}} }
   
   Order 1:

   { { {0}, {1}, {2} }, {}, {} }
    
   Order 2:

   { { {0}, {1}, {2} }, { {3}, {4}, {5} }, {} }
    
   Order 3:

   { { {0}, {1}, {2} }, { {3,4}, {5,6}, {7,8} }, {9} }

   \endverbatim

   * Above, we have used the ordering given in Hughes' textbook.
   */
  virtual void getReferenceDOFs(
    const CellType& maximalCellType,
    const CellType& cellType,
    Array<Array<Array<int> > >& dofs
    ) const = 0 ;



  /** \brief Return the total number of degrees of freedom associated with
   * this basis on a specified cell type. Note: the count returned
   * by this function includes DOFs owned by facets of the specified
   * reference cell. 
   * 
   *
   * \param cellType 
   *         [in] type of cell for which we want DOF information.
   *
   * \param maximalCellType
   *         [in] maximal-dimension cell type to which the cell is connected.
   *         For example, a triangle in 3D could be connected to a prism
   *         or to a tetrahedron. 
   * 
   * \return Number of DOFs associated with the cell type and its facets.    
   */
  virtual int nReferenceDOFsWithFacets(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const ;


  /** \brief Return the total number of degrees of freedom associated with
   * this basis on a specified cell type. Note: the count returned
   * by this function DOES NOT include DOFs owned by facets of the specified
   * reference cell. 
   * 
   *
   * \param cellType 
   *         [in] type of cell for which we want DOF information.
   *
   * \param maximalCellType
   *         [in] maximal-dimension cell type to which the cell is connected.
   *         For example, a triangle in 3D could be connected to a prism
   *         or to a tetrahedron. 
   * 
   * \return Number of DOFs associated with the cell type.    
   */
  virtual int nReferenceDOFsWithoutFacets(
    const CellType& maximalCellType,
    const CellType& cellType
    ) const = 0;



  /** \brief Comparison function allowing use of 
   * OrderedHandle<BasisFamilyBase>
   * in sorted containers. This is needed by the MixedDOFMap ctor when
   * it uses an STL map to group functions having the same bases into chunks. 
   *
   *  Note: this method should normally only be called from within the
   *  comparison operator of OrderedHandle, in which context comparisons
   *  between different derived types have already been resolved by 
   *  comparisons of typeid. Thus, we can require that the lessThan()
   *  function be called only with an argument whose typeid is 
   *  equal to that of <tt>*this.</tt> We recommend that all overriding
   *  implementations check that condition.  
   * 
   * \param other 
   *         [in] Pointer to another basis family object. 
   *         Precondition: <tt>typeid(*this)==typeid(*other).</tt>
   */
  virtual bool lessThan(const BasisDOFTopologyBase* other) const = 0 ;
};


} // namespace Sundance


#endif
