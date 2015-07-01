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


// $Id$ 
// $Source$ 


//   


/*! \file NOX_Playa_Vector.H - Header file for NOX::NOXPlaya::Vector. */

#ifndef NOX_Playa_VECTOR_H
#define NOX_Playa_VECTOR_H

#include <iostream>
#include "NOX_Abstract_Vector.H" // base class
#include "NOX_Common.H" // for #include<vector> 
#include "PlayaVectorDecl.hpp" // definitions for Playa Vectors
#include "PlayaVectorSpaceDecl.hpp"


namespace NOX 
{


namespace NOXPlaya 
{
using Teuchos::RCP;
  
/** Implementation of the NOX::Abstract::Vector interface for 
 * Playa::Vector<double> (using Playa handle to PlayaCore).
 *
 * @author Jill Reese and Kevin Long
 */
class Vector : public NOX::Abstract::Vector 
{

public:	

  //! Copy constructor 
  Vector(const NOX::NOXPlaya::Vector& source, 
    NOX::CopyType type = DeepCopy);

  //! Copy constructor
  Vector(const Playa::Vector<double>& source, 
    NOX::CopyType type = DeepCopy);

  //! Copy constructor when user requests specific precision  
  Vector(const NOX::NOXPlaya::Vector& source, 
    int numdigits, 
    NOX::CopyType type = DeepCopy);

  //! Copy constructor when user requests specific precision 
  Vector(const Playa::Vector<double>& source, 
    int numdigits,
    NOX::CopyType type = DeepCopy);

  /* There's no need to write an explicit dtor for this class, because
   * the Playa vector cleans up its own memory. I've commented out 
   * this method for now, but it should probably just be deleted. 
   * - KL 12 July 2004 */ 
  ////! Destruct Vector.
  //~Vector();

  /** Access to a writable reference to the underlying Playa vector*/
  Playa::Vector<double>& getPlayaVector();

  /** Access to a read-only reference to the underlying Playa vector*/
  const Playa::Vector<double>& getPlayaVector() const;

  /** Access to the vector precision expected by user*/
  int getPrecision() const;

  //@{ \name Initialization methods.

  /** Initialize all elements to a constant */
  NOX::Abstract::Vector& init(double value);

  /** Assignment operator */
  NOX::Abstract::Vector& operator=(const NOX::NOXPlaya::Vector& y);
  /** Assignment operator */
  NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);

  /* Element access should not be implemented in this class; it should
   * be done through the LoadableVector and AccessibleVector 
   * interfaces of Playa.
   * -KL 12 July 2004  */
  // Store value in the i-th element of a vector 
  //void setElement(int i, const double& value);

  // Return the i-th element (const version).
  //const double& getElement(int i) const;
  
  /** Element-wise absolute value */
  NOX::Abstract::Vector& abs(const NOX::NOXPlaya::Vector& y);
  /** Element-wise absolute value */
  NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);
  
  /** Element-wise reciprocal */
  NOX::Abstract::Vector& reciprocal(const NOX::NOXPlaya::Vector& y);
  /** Element-wise reciprocal */
  NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);
  
  //@}
  
  //@{ \name Update methods.
  
  // derived
  NOX::Abstract::Vector& scale(double gamma);
  
  // derived
  NOX::Abstract::Vector& scale(const NOX::NOXPlaya::Vector& a);
  NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);
  
  // derived
  NOX::Abstract::Vector& update(double alpha, const NOX::NOXPlaya::Vector& a, double gamma = 0.0);
  NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0);
  
  // derived
  NOX::Abstract::Vector& update(double alpha, const NOX::NOXPlaya::Vector& a, 
    double beta, const NOX::NOXPlaya::Vector& b,
    double gamma = 0.0);
  NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, 
    double beta, const NOX::Abstract::Vector& b,
    double gamma = 0.0);
  
  //@}
  
  //@{ \name Creating new Vectors. 
  
  // derived
#ifdef TRILINOS_6
  NOX::Abstract::Vector* clone(NOX::CopyType type = NOX::DeepCopy) const;
#else
  RCP<NOX::Abstract::Vector> clone(NOX::CopyType type = NOX::DeepCopy) const ;
#endif
  
  //@}
  
  //@{ \name Norms.
  
  // derived
  double norm(NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm) const;
  
  // derived
  double norm(const NOX::NOXPlaya::Vector& weights) const;
  double norm(const NOX::Abstract::Vector& weights) const;
  
  //@}
  
  //@{ \name Dot products
  
  // derived
  double dot(const NOX::NOXPlaya::Vector& y) const;
  double dot(const NOX::Abstract::Vector& y) const;
  double innerProduct(const NOX::Abstract::Vector& y) const;
  
  //@}
  
  // derived
  //! The length of a Playa vector is known only by the vector space it belongs to
  //! We query the vector space and ask its dimension  
  int length() const;

  //! Prints out the vector to the specified stream. 
  /*! 
    For example, a vector would appear as
    \f[ \left[ \; 0.1 \; 2.34 \; 5 \; \right] \f] 
    It will be all on one line, with a single space between each entry, bracketed on either side.
  */
  std::ostream& leftshift(std::ostream& stream) const;

  // derived
  void print() const;

private:

  /** user-specified precision for printing vectors */
  int precision;

  /** The Playa vector owned by this object. */
  Playa::Vector<double> x;
};

} // namespace Playa
} // namespace NOX


namespace std{
//! Function for printing
ostream& operator<<(std::ostream& stream, const NOX::NOXPlaya::Vector& v);
}

#endif
