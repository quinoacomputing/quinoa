/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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

#ifndef PLAYA_LINEAROPERATORDECL_HPP
#define PLAYA_LINEAROPERATORDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaLinearOperatorBaseDecl.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "PlayaRowAccessibleOp.hpp"


namespace Playa
{
using namespace Teuchos;


template <class Scalar>  class LinearSolver;
template <class Scalar>  class VectorSpace;
template <class Scalar>  class Vector;
template <class Scalar>  class VectorType;

/** 
 * User-level linear operator class
 */
template <class Scalar>
class LinearOperator : public Playa::Handle<LinearOperatorBase<Scalar> >
{
public:
  /** \name Constructors, Destructors, and Assignment Operators */
  //@{
  /** Empty constructor*/
  LinearOperator();

  /** Constructor with smart pointer */
  LinearOperator(const RCP<LinearOperatorBase<Scalar> >& smartPtr);
  //@}

  /** Return the domain */
  const VectorSpace<Scalar> domain() const ;

  /** Return the range */
  const VectorSpace<Scalar> range() const ;


  /** 
   * Compute
   * \code
   * out = beta*out + alpha*op*in;
   * \endcode
   **/
  void apply(const Vector<Scalar>& in,
    Vector<Scalar>& out) const ;

  /**  
   * Compute
   * \code
   * out = beta*out + alpha*op^T*in;
   * \endcode
   **/
  void applyTranspose(const Vector<Scalar>& in,
    Vector<Scalar>& out) const ;


  //       /** For the moment this does nothing*/
  LinearOperator<Scalar> form() const {return *this;}
      
      
  /** Get a stopwatch for timing vector operations */
  RCP<Time>& opTimer();

  /**
   * Return a TransposeOperator.
   */
  LinearOperator<Scalar> transpose() const ; 


  /** Return a Loadable Matrix  */
  RCP<LoadableMatrix<Scalar> > matrix();

  /** Get a row of the underlying matrix */     
  void getRow(const int& row, 
    Teuchos::Array<int>& indices, 
    Teuchos::Array<Scalar>& values) const ;
    

  /** \name  Block operations  */
  //@{
      
  /** return number of block rows */
  int numBlockRows() const ;
      

  /** return number of block cols */
  int numBlockCols() const ;
      

  /** get the (i,j)-th block */
  LinearOperator<Scalar> getBlock(const int &i, const int &j) const ;


  /** get a writeable copy of the (i,j)-th block */
  LinearOperator<Scalar> getNonconstBlock(const int &i, const int &j) ;

  /** set the (i,j)-th block 
   *  If the domain and/or the range are not set, then we
   *  are building the operator
   */
  void setBlock(int i, int j, 
    const LinearOperator<Scalar>& sub);

  /** Finalize block assembly */
  void endBlockFill();

  //@}

      

private:

};

}


#endif
