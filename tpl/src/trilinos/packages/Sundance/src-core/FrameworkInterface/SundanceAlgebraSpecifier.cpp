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

#include "SundanceAlgebraSpecifier.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;
using namespace std;
using namespace Sundance;
using namespace Teuchos;




AlgebraSpecifier::AlgebraSpecifier()
  : EnumTypeField<AlgebraType>(ScalarAT), 
    direction_(-1)
{}

AlgebraSpecifier::AlgebraSpecifier(int direction)
  : EnumTypeField<AlgebraType>(CoordCompAT), 
    direction_(direction)
{}

AlgebraSpecifier::AlgebraSpecifier(
  const AlgebraType& at)
  : EnumTypeField<AlgebraType>(at), direction_(-1)
{
  assertNotType(CoordCompAT);
}

int AlgebraSpecifier::direction() const 
{
  assertType(CoordCompAT);
  return direction_;
}

bool AlgebraSpecifier::operator<(const AlgebraSpecifier& other) const 
{
  if (type() < other.type()) return true;
  if (type() > other.type()) return false;
  
  if (isCoordinateComponent())
    return direction() < other.direction();
  return false;
}

string AlgebraSpecifier::toString() const 
{
  TeuchosOStringStream os;
  os << *this;
  return os.str();
}

namespace std
{

ostream& operator<<(std::ostream& os, 
  const Sundance::AlgebraType& at)
{
  switch(at)
  {
    case ScalarAT:
      os << "Scalar";
      break;
    case VectorAT:
      os << "Vector";
      break;
    case CoordCompAT:
      os << "CoordComp";
      break;
    case NormalAT:
      os << "Normal";
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(1);
  }
  return os;
}

ostream& operator<<(std::ostream& os, const Sundance::AlgebraSpecifier& as)
{
  os << as.type();
  if (as.isCoordinateComponent()) os << "(d=" << as.direction() << ")";
  else os << "()";
  return os;
}

}

namespace Sundance
{

/** \relates AlgebraSpecifier */
AlgebraSpecifier vectorAlgebraSpec()
{
  return AlgebraSpecifier(VectorAT);
}

/** AlgebraSpecifier */
AlgebraSpecifier scalarAlgebraSpec()
{
  return AlgebraSpecifier(ScalarAT);
}

/** AlgebraSpecifier */
AlgebraSpecifier normalAlgebraSpec()
{
  return AlgebraSpecifier(NormalAT);
}

/** AlgebraSpecifier */
AlgebraSpecifier coordAlgebraSpec(int dir)
{
  return AlgebraSpecifier(dir);
}

}
