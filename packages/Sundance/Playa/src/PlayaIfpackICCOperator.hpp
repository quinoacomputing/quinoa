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

#ifndef PLAYA_IFPACKICCOPERATOR_HPP
#define PLAYA_IFPACKICCOPERATOR_HPP


#include "PlayaEpetraMatrix.hpp"
#include "PlayaLinearOpWithSpacesDecl.hpp"
#include "Ifpack_ICT.h"


namespace Playa
{
/**
 * This is the operator representation of the inverse upper triangular factor 
 * \f${\tilde R}^{-1}\f$  
 * in the incomplete Cholesky factorization of a matrix \f$ A \approx {\tilde R}{\tilde R}^T\f$. 
 *
 * @author Kimberly Kennedy, with some small modifications by Kevin Long
 */
class IfpackICCOperator : 
  public LinearOpWithSpaces<double>,
  public Printable
{
public:
  /** Construct the operator. During construction, the matrix A will be
   * approximately factored using the options given in the arguments. 
   *
   * \param A The matrix to be factored
   * \param fillLevels Specifies the number of on-processor 
   * nonzeros allowed in each
   * row of the factorization, given as a ratio of the 
   * allowed nnz per row in the factor to nnz per row in A. 
   * \param overlapFill Fill allowed for off-processor elements
   * \param dropTol New elements are dropped if |R_ij| < dropTol*R_jj.
   * \param relaxationValue Fraction of dropped element mass to be added to 
   * diagonal.
   * \param relativeThreshold Fraction of diagonal element to be added to diagonal
   * \param absoluteThreshold Amount to be added to each diagonal element
   */
  IfpackICCOperator(const EpetraMatrix* A,
		    int fillLevels,
		    int overlapFill,
		    double dropTol,
		    double relaxationValue,
		    double relativeThreshold,
		    double absoluteThreshold);

  /** 
   * Apply the operator. 
   */
  virtual void apply(
		     Teuchos::ETransp applyType,
		     const Vector<double>& in,
		     Vector<double> out) const ;
  
  /** \name Diagnostic output */
  //@{
  /** Print the matrix */
  virtual void print(std::ostream& os) const ;
  //@}

  /** */
  std::ostream& describe(
			 std::ostream                         &out
			 ,const Teuchos::EVerbosityLevel      verbLevel
			 ,const std::string                   leadingIndent
			 , const std::string                   indentSpacer
			 ) const 
  {
    out << leadingIndent << indentSpacer << this->description() << std::endl;
    return out;
  }
  /** */
  std::string description() const ;
  //@}



private:
  RCP<Ifpack_ICT> precond_;

};


}

#endif 
