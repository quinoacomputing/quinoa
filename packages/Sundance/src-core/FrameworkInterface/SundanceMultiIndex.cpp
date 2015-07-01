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

#include "SundanceMultiIndex.hpp"
#include "PlayaExceptions.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Sundance;
using namespace Teuchos;

MultiIndex::MultiIndex()
	: m_(maxDim(), 0)
{;}

MultiIndex::MultiIndex(int x, int y, int z)
	: m_(maxDim(), 0)
{
	m_[0] = x;
	m_[1] = y;
	m_[2] = z;
}

MultiIndex MultiIndex::operator+(const MultiIndex& other) const 
{
	MultiIndex rtn;

	for (int i=0; i<maxDim(); i++)
		{
			rtn.m_[i] = m_[i] + other[i];
		}
	return rtn;
}

MultiIndex MultiIndex::operator-(const MultiIndex& other) const 
{
	MultiIndex rtn;

	for (int i=0; i<maxDim(); i++)
		{
			rtn.m_[i] = m_[i] - other[i];
		}
	return rtn;
}

MultiIndex MultiIndex::operator-() const 
{
	MultiIndex rtn;

	for (int i=0; i<maxDim(); i++)
		{
			rtn.m_[i] = -m_[i];
		}
	return rtn;
}

string MultiIndex::toString() const
{
	return "(" + Teuchos::toString(m_[0]) + ","
		+ Teuchos::toString(m_[1]) + ","
		+ Teuchos::toString(m_[2]) + ")";
} 

XMLObject MultiIndex::toXML() const
{
	XMLObject rtn("MultiIndex");
	rtn.addAttribute("indices", toString());
	return rtn;
}

bool MultiIndex::operator==(const MultiIndex& m) const 
{
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] != m[i]) return false;
		}
	return true;
}

bool MultiIndex::operator<(const MultiIndex& m) const 
{
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] > m.m_[i]) return false;
      if (m_[i] < m.m_[i]) return true;
		}
	return false;
}



int MultiIndex::order() const 
{
  int h = 0;
	for (int i=0; i<maxDim(); i++)
		{
			h += m_[i];
		}
	return h;
}

bool MultiIndex::isValid() const 
{
	for (int i=0; i<maxDim(); i++)
		{
      if (m_[i] < 0) return false;
		}
	return true;
}

int MultiIndex::firstOrderDirection() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(order() != 1, std::logic_error,
                     "bad order in MultiIndex::firstOrderDirection() const");
	for (int i=0; i<maxDim(); i++)
		{
			if (m_[i] == 1) return i;
		}
  return -1;
}



string MultiIndex::coordForm() const
{
  std::string rtn;
  
  for (int i=0; i<m_[0]; i++)
		{
      rtn += "x";
		}
  for (int i=0; i<m_[1]; i++)
		{
      rtn += "y";
		}
  for (int i=0; i<m_[2]; i++)
		{
      rtn += "z";
		}
	return rtn;
}


