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


#include "NOX_Common.H"
#include "NOX_Playa_Vector.hpp"
#include "NOX_Utils.H"
#include "NOX_Random.H" // for Random class

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using namespace Teuchos;
using std::runtime_error;
using std::cout;
using std::ostream;

NOX::NOXPlaya::Vector::Vector(const NOX::NOXPlaya::Vector& source, 
			 NOX::CopyType type)
  :
  precision(3) // 3 digits of accuracy is default
{
 switch (type) 
 {
  case NOX::DeepCopy:
    x = source.x.copy();
    break;

  case NOX::ShapeCopy:
    x = ((source.x).space()).createMember();
    break;

  default:
    std::cerr << "NOX:Playa::Vector - invalid CopyType for copy constructor." << std::endl;
    throw "NOX Playa Error";
  }
}

NOX::NOXPlaya::Vector::Vector(const Playa::Vector<double>& source, 
			 NOX::CopyType type)
  :
  precision(3) // 3 digits of accuracy is default
{
  switch (type) 
 {
    
  case NOX::DeepCopy:
    x = source.copy();
    break;

  case NOX::ShapeCopy:
    x = ((source).space()).createMember();
    break;

  default:
    std::cerr << "NOX:Playa::Vector - invalid CopyType for copy constructor." << std::endl;
    throw "NOX Playa Error";
  }
}


NOX::NOXPlaya::Vector::Vector(const NOX::NOXPlaya::Vector& source, 
                         int numdigits,
			 NOX::CopyType type)
  :
  precision(numdigits)
{
 switch (type) 
 {
    
  case NOX::DeepCopy:
    x = source.x.copy();
    break;

  case NOX::ShapeCopy:
    x = ((source.x).space()).createMember();
    break;

  default:
    std::cerr << "NOX:Playa::Vector - invalid CopyType for copy constructor." << std::endl;
    throw "NOX Playa Error";
  }
}

NOX::NOXPlaya::Vector::Vector(const Playa::Vector<double>& source, 
                         int numdigits,
			 NOX::CopyType type)
  :
  precision(numdigits)
{
  switch (type) 
 {
    
  case NOX::DeepCopy:
    x = source.copy();
    break;

  case NOX::ShapeCopy:
    x = ((source).space()).createMember();
    break;

  default:
    std::cerr << "NOX:Playa::Vector - invalid CopyType for copy constructor." << std::endl;
    throw "NOX Playa Error";
  }
}



Playa::Vector<double>& NOX::NOXPlaya::Vector::getPlayaVector()
{
  return x;
}
 
const Playa::Vector<double>& NOX::NOXPlaya::Vector::getPlayaVector() const
{
  return x;
}

int NOX::NOXPlaya::Vector::getPrecision() const
{
  return precision;
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::operator=(
					   const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const NOX::NOXPlaya::Vector&>(source));
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::operator=(
					   const NOX::NOXPlaya::Vector& source)
{
  // in Playa operator= results in a shallow copy while 
  // acceptCopyOf(source.x) provides the deep copy we want
  x = source.getPlayaVector().copy();
  return *this;
}
  

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::init(double value)
{
  x.setToConstant(value);
  return *this;
}


NOX::Abstract::Vector& NOX::NOXPlaya::Vector::abs(
					     const NOX::Abstract::Vector& base)
{
  return abs(dynamic_cast<const NOX::NOXPlaya::Vector&>(base));
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::abs(
					     const NOX::NOXPlaya::Vector& base)
{
  x.acceptCopyOf(base.x);
  x.selfAbs();
  return *this;
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::reciprocal(
					    const NOX::Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const NOX::NOXPlaya::Vector&>(base));
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::reciprocal(
					    const NOX::NOXPlaya::Vector& base)
{
  x.acceptCopyOf(base.x);
  x.selfReciprocal();
  return *this;
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::scale(double alpha)
{
  x.scale(alpha);
  return *this;
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::update(
					       double alpha, 
					       const NOX::Abstract::Vector& a, 
					       double gamma)
{
  return update( alpha, dynamic_cast<const NOX::NOXPlaya::Vector&>(a), gamma);
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::update(
						 double alpha, 
						 const NOX::NOXPlaya::Vector& a, 
						 double gamma)
{
  x.update(alpha,a.x,gamma);
  return *this;
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::update(
					      double alpha, 
					      const NOX::Abstract::Vector& a, 
					      double beta, 
					      const NOX::Abstract::Vector& b,
					      double gamma)
{
  return update(alpha, dynamic_cast<const NOX::NOXPlaya::Vector&>(a), 
		beta, dynamic_cast<const NOX::NOXPlaya::Vector&>(b), gamma);
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::update(
					       double alpha, 
					       const NOX::NOXPlaya::Vector& a, 
					       double beta, 
					       const NOX::NOXPlaya::Vector& b,
					       double gamma)
{
  x.update(alpha,a.x,beta,b.x,gamma);
  return *this;
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::scale(
					      const NOX::Abstract::Vector& a)
{  
  return scale(dynamic_cast<const NOX::NOXPlaya::Vector&>(a));
}

NOX::Abstract::Vector& NOX::NOXPlaya::Vector::scale(const NOX::NOXPlaya::Vector& a)
{  
  x.selfDotStar(a.x);
  return *this;
}

#ifdef TRILINOS_6
NOX::Abstract::Vector* NOX::NOXPlaya::Vector::clone(NOX::CopyType type) const
{
  return new NOX::NOXPlaya::Vector(*this, type);
}
#else
RCP<NOX::Abstract::Vector> NOX::NOXPlaya::Vector::clone(NOX::CopyType type) const
{
  return rcp(new NOX::NOXPlaya::Vector(*this, type));
}
#endif

double NOX::NOXPlaya::Vector::norm(NOX::Abstract::Vector::NormType type) const
{
  
  if (this->length() == 0)
    return 0.0;

  double value;			// final answer

  switch (type) 
  {
  case MaxNorm:
    value = x.normInf();
    break;
  case OneNorm:
    value = x.norm1();
    break;
  case TwoNorm:
  default:
    value = x.norm2();
   break;
  }

  return value;
}

double NOX::NOXPlaya::Vector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const NOX::NOXPlaya::Vector&>(weights));
}

double NOX::NOXPlaya::Vector::norm(const NOX::NOXPlaya::Vector& weights) const
{
  if (weights.length() != this->length()) 
  {
    std::cerr << "NOX::NOXPlaya::Vector::norm - size mismatch for weights vector" << std::endl;
    throw "NOX::NOXPlaya Error";
  }
  return x.norm2(weights.getPlayaVector());
}

double NOX::NOXPlaya::Vector::dot(const NOX::Abstract::Vector& y) const
{
  return dot(dynamic_cast<const NOX::NOXPlaya::Vector&>(y));
}

double NOX::NOXPlaya::Vector::innerProduct(const NOX::Abstract::Vector& y) const
{
  return dot(dynamic_cast<const NOX::NOXPlaya::Vector&>(y));
}

double NOX::NOXPlaya::Vector::dot(const NOX::NOXPlaya::Vector& y) const
{
  if (y.length() != this->length()) 
  {
    std::cerr << "NOX::NOXPlaya::Vector::dot - size mismatch for y vector" << std::endl;
    throw "NOX::NOXPlaya Error";
  }

  return x.dot(y.x);
}

int NOX::NOXPlaya::Vector::length() const
{
  return (x.space()).dim();
}



ostream& NOX::NOXPlaya::Vector::leftshift(std::ostream& stream) const
{
  x.print(stream);
  return stream;
}

namespace std
{
ostream& operator<<(std::ostream& stream, const NOX::NOXPlaya::Vector& v)
{
  return v.leftshift(stream);
}
}

void NOX::NOXPlaya::Vector::print() const
{
  cout << *this << std::endl;
}
