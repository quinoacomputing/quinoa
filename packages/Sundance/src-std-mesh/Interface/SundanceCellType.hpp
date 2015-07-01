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

#ifndef CELLTOPOLOGYCODE_H
#define CELLTOPOLOGYCODE_H


#include "SundanceDefs.hpp"

namespace Sundance 
{
/** \defgroup Sundance_CellType_grp Cell Type Description
 */

/** \brief Enumeration of specific 1D, 2D, and 3D cell types.
 *
 * See the implementation of the nonmember functions <tt>dimension()</tt>,
 * <tt>numFacets()</tt>, and <tt>facetType()</tt> a full description of what
 * these cell types are.
 *
 * ToDo: Consider making CellType an abstract base class so that any type of
 * cell can be implemented?
 *
 * \ingroup Sundance_CellType_grp
 */
enum CellType {
  NullCell        ///< No cell specified
  ,PointCell      ///< 0D vertex cell
  ,LineCell       ///< 1D line, or edge, cell
  ,TriangleCell   ///< 2D triangle
  ,TetCell        ///< 3D tetrahedral cell
  ,QuadCell       ///< 2D quadrilateral cell
  ,BrickCell      ///< 3D "brick" cell
  ,PrismCell      ///< 3D prism cell
};

/** \brief Return a std::string representation of the cell type.
 *
 * \ingroup Sundance_CellType_grp
 */
std::string toString(const CellType& c) ;

/** \brief Return the dimension of the cell type.
 *
 * \ingroup Sundance_CellType_grp
 */
int dimension(const CellType& c) ;

/** \brief Return the number of faces of a given facet dimension for a cell type.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>facetDim <= dimension(c)</tt>
 * </ul>
 *
 * \returns -1 for null cells and point cells
 *
 * \ingroup Sundance_CellType_grp
 */
int numFacets(const CellType& c, int facetDim);

/** \return Return the type of facet of a given facet dimension and a
 * particualar facet.
 *
 * \param  c
 *           [in] Cell type
 * \param  facetDim
 *           [in] The dimension of the facet type requested
 * \param  facetIndex
 *           [in] The particualar index of the facet as defined
 *           by the convension of the cell type (which is actually
 *           partially defined in this function).
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>facetDim <= dimension(c)</tt>
 * <li><tt>0 <= facetIndex < numFacets(c,facetDim)</tt>
 * </ul>
 *
 * \ingroup Sundance_CellType_grp
 */
CellType facetType(const CellType& c, int facetDim, int facetIndex);



/** \brief output stream operator for <tt>CellType</tt>.
 *
 * \ingroup Sundance_CellType_grp
 */
inline std::ostream& operator<<(std::ostream& os, const CellType& c)
{
  os << toString(c);
  return os;
}
  
} // namespace Sundance


#endif


