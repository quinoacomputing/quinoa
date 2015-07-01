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

#ifndef SUNDANCE_ORDEREDTUPLE_H
#define SUNDANCE_ORDEREDTUPLE_H

#include "SundanceDefs.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  using namespace Teuchos;

  /** OrderedPair provides a means of lexigraphic comparison of a pair of
   * objects. The pair {a1, b1} is compared to {a2, b2} by first
   * comparing the most significant entries a1 and a2, and if they are
   * equal, comparing the least significant entries b1 and b2. */
  template<class A, class B>
    class OrderedPair
    {
    public:
      /** */
      OrderedPair(const A& _a, const B& _b)
        : a_(_a), b_(_b) {;}

      /** */
      inline bool operator<(const OrderedPair<A, B>& other) const
        {
          if ( a_ < other.a_ ) 
            {
              return true;
            }
          if ( other.a_ < a_) 
            {
              return false;
            }

          bool rtn = b_ < other.b_;
          return rtn;
        }

      /** */
      const A& first() const {return a_;}

      /** */
      const B& second() const {return b_;}

    private:
      A a_;
      B b_;
    };

  /** Lexigraphically-comparable triple of objects. */
  template<class A, class B, class C>
    class OrderedTriple : public OrderedPair<A, OrderedPair<B, C> >
    {
    public:
      /** */
      OrderedTriple(const A& _a, const B& _b, const C& _c)
        : OrderedPair<A, OrderedPair<B, C> >(_a, OrderedPair<B,C>(_b,_c))
        {;}

      const A& a() const {return this->first();}

      const B& b() const {return this->second().first();}

      const C& c() const {return this->second().second();}
    };

  /** Lexigraphically-comparable quartet of objects. */
  template<class A, class B, class C, class D>
    class OrderedQuartet : public OrderedPair<A, OrderedTriple<B, C, D> >
    {
    public:
      /** */
      OrderedQuartet(const A& _a, const B& _b, const C& _c, const D& _d)
        : OrderedPair<A, OrderedTriple<B, C, D> >(_a, OrderedTriple<B,C,D>(_b,_c,_d))
        {;}

      const A& a() const {return this->first();}
      const B& b() const {return this->second().first();}
      const C& c() const {return this->second().second().first();}
      const D& d() const {return this->second().second().second();}
    };

  /** */
  template <class A, class B>
  inline std::ostream& operator<<(std::ostream& os, const OrderedPair<A,B>& p)
  {
    os << "{" << p.first() << ", " << p.second() << "}";
    return os;
  }

}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
