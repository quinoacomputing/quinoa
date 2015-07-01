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

#ifndef SUNDANCE_ATOCPOINTLOCATOR_H
#define SUNDANCE_ATOCPOINTLOCATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceExpr.hpp"
#include "SundanceCellFilter.hpp"

namespace Sundance
{
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Sundance;
  using namespace Teuchos;

  /**
   * AToCPointLocator finds the cell index for a point within an unstructured
   * mesh.
   *
   * Note: not tested in parallel.
   */
  class AToCPointLocator
  {
  public:
    /** */
    AToCPointLocator(const Mesh& mesh, 
                     const CellFilter& subdomain,
                     const std::vector<int>& nx);



    /** Find the index of a point in an overlaid structured grid. */
    int getGridIndex(const double* x) const ;

    /** Use an overlaid structured grid to estimate the location of the point. */
    int guessCell(const double* x) const 
    {return (*table_)[getGridIndex(x)];}

    /** Find the cell that contains the specified point */
    int findEnclosingCell(int initialGuessLID, const double* x) const ;

    /** Find the cell that contains the specified point, also
     * computing local coordinates within that cell. */
    int findEnclosingCell(int initialGuessLID, const double* x,
                          double* localCoords) const ;

    /** */
    void fillMaximalNeighbors(int cellLID, const int* facetLID) const ;

    /** Test whether a point is within a specified cell */
    bool cellContainsPoint(int cellLID, 
                           const double* x, 
                           const int* facetLID) const ;

    /** Test whether a point is within a specified cell, and if so,
     * compute local coordinates within that cell. */
    bool cellContainsPoint(int cellLID, 
                           const double* x, 
                           const int* facetLID,
                           double* localCoords) const ;

    /** */
    const Mesh& mesh() const {return mesh_;}

    /** */
    const CellFilter& subdomain() const {return subdomain_;}



    /** */
    static Point makePoint(int dim, const double* x) ;
  private:

    /** Find the range of structured grid cells within the bounding box of 
     * a cell. */
    void getGridRange(const Mesh& mesh, int cellDim, int cellLID,
                 Array<int>& lowIndex, Array<int>& highIndex) const;




    int dim_;
    Mesh mesh_;
    int nFacets_;
    std::vector<int> nx_;
    Array<double> low_;
    Array<double> high_;
    Array<double> dx_;
    RCP<Array<int> > table_;
    CellFilter subdomain_;
    mutable Array<RCP<Set<int> > > neighborSet_;
  };
}


#endif
