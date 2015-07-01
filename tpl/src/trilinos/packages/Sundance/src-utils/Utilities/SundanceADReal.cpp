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


#include "SundanceADReal.hpp"
#include "PlayaExceptions.hpp"


using namespace Sundance;
using namespace Teuchos;

ADReal ADReal::operator-() const
{
	ADReal rtn = *this;
	rtn.value_ = -rtn.value_;
  rtn.gradient_ = -gradient_;

  addFlops(1 + gradient_.dim());
	return rtn;
}

ADReal ADReal::operator+(const ADReal& other) const 
{
	ADReal rtn(*this);
	rtn += other;
	return rtn;
}

ADReal ADReal::operator-(const ADReal& other) const 
{
	ADReal rtn(*this);
	rtn -= other;
	return rtn;
}

ADReal ADReal::operator*(const ADReal& other) const 
{
	ADReal rtn(*this);
	rtn *= other;
	return rtn;
}

ADReal ADReal::operator/(const ADReal& other) const 
{
	ADReal rtn(*this);
	rtn /= other;
	return rtn;
}

ADReal& ADReal::operator+=(const ADReal& other) 
{
	value_ += other.value_;
	gradient_ += other.gradient_;

  addFlops(1 + gradient_.dim());

	return *this;
}

ADReal& ADReal::operator-=(const ADReal& other) 
{
	value_ -= other.value_;
	gradient_ -= other.gradient_;

  addFlops(1 + gradient_.dim());

	return *this;
}

ADReal& ADReal::operator*=(const ADReal& other) 
{
	gradient_ = (other.value_*gradient_ + value_*other.gradient_);
	value_ *= other.value_;

  addFlops(1 + 3*gradient_.dim());

	return *this;
}

ADReal& ADReal::operator/=(const ADReal& other) 
{
  TEUCHOS_TEST_FOR_EXCEPTION(other.value_ == 0.0, std::runtime_error,
                     "ADReal::operator/=() division by zero");

	gradient_ = (gradient_/other.value_ 
							 - value_/(other.value_*other.value_) * other.gradient_);
	value_ /= other.value_;

  addFlops(3 + 3*gradient_.dim());

	return *this;
}


// operations with constant reals 

ADReal& ADReal::operator+=(const double& other) 
{
	value_ += other;
  addFlops(1);
	return *this;
}

ADReal& ADReal::operator-=(const double& other) 
{
	value_ -= other;

  addFlops(1);
	return *this;
}

ADReal& ADReal::operator*=(const double& other) 
{
	gradient_ *= other;
	value_ *= other;

  addFlops(1 + gradient_.dim());

	return *this;
}

ADReal& ADReal::operator/=(const double& other) 
{
  TEUCHOS_TEST_FOR_EXCEPTION(other == 0.0, std::runtime_error,
                     "ADReal::operator/=() division by zero");

  addFlops(2 + gradient_.dim());
	gradient_ *= 1.0/other;
	value_ /= other;
	return *this;
}


ADReal ADReal::operator+(const double& other) const 
{
	ADReal rtn(*this);
	rtn += other;
	return rtn;
}

ADReal ADReal::operator-(const double& other) const 
{
	ADReal rtn(*this);
	rtn -= other;
	return rtn;
}

ADReal ADReal::operator*(const double& other) const 
{
	ADReal rtn(*this);
	rtn *= other;
	return rtn;
}

ADReal ADReal::operator/(const double& other) const 
{
	ADReal rtn(*this);
	rtn /= other;
	return rtn;
}


namespace Sundance
{
ADReal operator+(const double& scalar, const ADReal& a)
{
  return a+scalar;
}

ADReal operator-(const double& scalar, const ADReal& a)
{
  return -a+scalar;
}

ADReal operator*(const double& scalar, const ADReal& a)
{
  return a*scalar;
}

ADReal operator/(const double& scalar, const ADReal& a)
{
  ADReal rtn(a);
  rtn.reciprocate();
  return rtn*scalar;
}
}

void ADReal::reciprocate() 
{
	TEUCHOS_TEST_FOR_EXCEPTION(value_==0.0, std::runtime_error,
                     "div by 0 in ADReal::reciprocate()");

  addFlops(1 + gradient_.dim());

	gradient_ /= (value_*value_);
  gradient_ *= -1.0;
	value_ = 1.0/value_;
}


	





