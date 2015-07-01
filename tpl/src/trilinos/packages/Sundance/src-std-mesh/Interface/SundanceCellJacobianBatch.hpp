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

#ifndef SUNDANCE_CELLJACOBIANBATCH_H
#define SUNDANCE_CELLJACOBIANBATCH_H


#include "SundanceDefs.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
using namespace Teuchos;




/**
 * A CellJacobianBatch is a collection of Jacobian matrices 
 * for many quadrature
 * points distributed over batch of cells. All cells 
 * must have the same dimension
 * and number of quadrature points. Affine cells have constant Jacobians,
 * so in that case the quadrature point index can be ignored.
 * See the ReferenceCellBase documentation
 * for definitions of the coordinate systems used.
 *
 * <H4> Data layout </H4>
 * All Jacobian elements for all points and cells are 
 * packed into a single vector.
 * The length of the vector is 
 * \f$ N_{cell} \times N_{quad} \times D^2, \f$ where
 * \f$D\f$ is the spatial dimension.
 * The indices into the vector cycle in the following order 
 * (from slowest to fastest)
 * <ol>
 * <li> cell number
 * <li> quadrature point number
 * <li> physical coordinate direction
 * <li> reference coordinate direction
 * </ol>
 * Thus, the jacobian values for the \f$q\f$-th quadrature point on the
 * \f$c\f$-th cell are the \f$D^2\f$ entries following 
 * the \f$(q + cN_{quad})D^2\f$-th
 * element.
 *
 */
class CellJacobianBatch 
  : public ObjectWithClassVerbosity<CellJacobianBatch>
{
public:
  /** empty ctor */
  CellJacobianBatch();

  /** get the spatial dimension */
  int spatialDim() const {return spatialDim_;}

  /** get the cell dimension */
  int cellDim() const {return cellDim_;}

  /** get the number of cells in the batch */
  int numCells() const {return numCells_;}

  /** get the number of quad points per cell */
  int numQuadPoints() const {return numQuad_;}

  /** resize the batch */
  void resize(int numCells, int numQuad, int spatialDim, int cellDim);

  /** resize the batch, using one quadrature point per cell 
   * (appropriate for affine elements) */
  void resize(int numCells, int spatialDim, int cellDim);

  /** Get a pointer to the values at the q-th quadrature 
   * point on the c-th cell.
   * @param c the index of the cell in the batch
   * @param q the index of the quadrature point
   */
  double* jVals(int c, int q);

  /** 
   * Get a pointer to the start of the c-th Jacobian in the batch. 
   */
  double* jVals(int c)
    {return &(J_[c*jSize_]);}

  /** Get a constant pointer to start of c-th Jacobian in the batch */
  const double *jVals(int c) const { return &(J_[c*jSize_]); }

  /** */
  double* detJ(int c)
    {return &(detJ_[c]);}

  /** get the vector of determinant values */
  const Array<double>& detJ() const 
    {if (!isFactored_ && cellDim()==spatialDim()) factor(); return detJ_;}
            
  /** 
   * Apply a cell's inverse Jacobian to (possibly) multiple rhs
   * stored in column-major order.
   */
  void applyInvJ(int cell, int q, double* rhs,
    int nRhs, bool trans) const ;
            
  /** 
   * Apply an affine cell's inverse Jacobian to (possibly) multiple rhs
   * stored in column-major order.
   */
  void applyInvJ(int cell, double* rhs,
    int nRhs, bool trans) const 
    {applyInvJ(cell, 0, rhs, nRhs, trans);}
          
  /** 
   * Get the explicit inverse of the Jacobian for the given
   * (cell, quad) combination.
   */
  void getInvJ(int cell, int quad, Array<double>& invJ) const ;

  /** 
   * Get the explicit inverse of the Jacobian for the given
   * affine cell.
   */
  void getInvJ(int cell, Array<double>& invJ) const 
    {getInvJ(cell, 0, invJ);}

          
  /** */
  void print(std::ostream& os) const ;

  static double& totalFlops() {static double rtn = 0; return rtn;}



  static void addFlops(const double& flops) {totalFlops() += flops;}

private:
          
  void factor() const ;

  void computeInverses() const ;

  int spatialDim_;
  int cellDim_;
  int jSize_;
  int numCells_;
  int numQuad_;
  mutable Array<int> iPiv_;
  mutable Array<double> J_;
  mutable Array<double>  detJ_;
  mutable Array<double>  invJ_;
  mutable bool isFactored_;
  mutable bool hasInverses_;
  mutable bool hasDetJ_;
};


inline std::ostream& operator<<(std::ostream& os, 
  const CellJacobianBatch& J)
{
  J.print(os);
  return os;
}
}



#endif
