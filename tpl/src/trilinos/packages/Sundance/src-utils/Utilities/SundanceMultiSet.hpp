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

#ifndef SUNDANCE_MULTISET_H
#define SUNDANCE_MULTISET_H

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Array.hpp"
#include <set>


namespace Sundance
{
using namespace Teuchos;

/** 
 * Extension of STL multiset, adding some nicer syntax 
 * and an iostream insertion operator.
 */
template<class Key>
class MultiSet : public std::multiset<Key>
{
public:
  /** */
  MultiSet() : std::multiset<Key>() {;}

  /** Test whether the specified key is present in the set */
  bool contains(const Key& key) const {return this->find(key) != this->end();}

  /** Put a new entry in the map */
  void put(const Key& key) {this->insert(key);}

  /** Write into an array */
  Array<Key> elements() const ;

  /** Write to stream */
  std::ostream& toStream(std::ostream& os) const ;

  /** Merge with another multiset, returning the merged set */
  MultiSet<Key> merge(const MultiSet<Key>& other) const ;

  /** Take another set and merge into this one, overwriting the original
   * with the merged set */
  void mergeFrom(const MultiSet<Key>& other) ;

  /** Write into a set, i.e., collapsing repeated entries */
  Set<Key> toSet() const ;

  /** Write to a std::string */
  std::string toString() const ;

  /** */
  bool operator==(const MultiSet<int>& other) const 
    {
      return !((*this) < other || other < (*this));
    }
};


template<class Key> inline
Array<Key> MultiSet<Key>::elements() const
{
  Array<Key> rtn;

  typename MultiSet<Key>::const_iterator iter;

  for (iter=this->begin(); iter != this->end(); iter++)
  {
    rtn.append(*iter);
  }
  return rtn;
}

template<class Key> inline
MultiSet<Key> MultiSet<Key>::merge(const MultiSet<Key>& other) const
{
  MultiSet<Key> rtn = *this;

  typename MultiSet<Key>::const_iterator iter;

  for (iter=other.begin(); iter != other.end(); iter++)
  {
    rtn.put(*iter);
  }
  return rtn;
}

template<class Key> inline
void MultiSet<Key>::mergeFrom(const MultiSet<Key>& other) 
{
  typename MultiSet<Key>::const_iterator iter;

  for (iter=other.begin(); iter != other.end(); iter++)
  {
    put(*iter);
  }
}


template<class Key> inline
Set<Key> MultiSet<Key>::toSet() const
{
  Set<int> rtn;
  typename MultiSet<Key>::const_iterator iter;

  for (iter=this->begin(); iter != this->end(); iter++)
  {
    rtn.put(*iter);
  }
  return rtn;
}

  

template<class Key> inline
std::ostream& MultiSet<Key>::toStream(std::ostream& os) const
{
  typename MultiSet<Key>::const_iterator iter;

  int k = 0;
  os << "{";
  for (iter=this->begin(); iter != this->end(); iter++, k++)
  {
    os << *iter;
    if (k<((int) this->size()-1)) os << ", ";
  }
  os << "}";

  return os;
}

template<class Key> inline
string MultiSet<Key>::toString() const
{
  std::ostringstream os;
  os << *this;
  return os.str();
}


/** \relates MultiSet Create a multiset */
template<class Key> inline
MultiSet<Key> makeMultiSet(const Key& k)
{
  MultiSet<Key> rtn;
  rtn.put(k);
  return rtn;
}

/** \relates MultiSet Create a multiset */
template<class Key> inline
MultiSet<Key> makeMultiSet(const Key& k1, const Key& k2)
{
  MultiSet<Key> rtn = makeMultiSet<Key>(k1);
  rtn.put(k2);
  return rtn;
}

/** \relates MultiSet Create a multiset */
template<class Key> inline
MultiSet<Key> makeMultiSet(const Key& k1, const Key& k2, const Key& k3)
{
  MultiSet<Key> rtn = makeMultiSet<Key>(k1, k2);
  rtn.put(k3);
  return rtn;
}

/** \relates MultiSet Create a multiset */
template<class Key> inline
MultiSet<Key> makeMultiSet(const Key& k1, const Key& k2, 
  const Key& k3, const Key& k4)
{
  MultiSet<Key> rtn = makeMultiSet<Key>(k1, k2, k3);
  rtn.put(k4);
  return rtn;
}

/** \relates MultiSet Create a multiset */
template<class Key> inline
MultiSet<Key> makeMultiSet(const Key& k1, const Key& k2, 
  const Key& k3, const Key& k4,
  const Key& k5)
{
  MultiSet<Key> rtn = makeMultiSet<Key>(k1, k2, k3, k4);
  rtn.put(k5);
  return rtn;
}

/** \relates MultiSet Create a multiset */
template<class Key> inline
MultiSet<Key> makeMultiSet(const Key& k1, const Key& k2, 
  const Key& k3, const Key& k4,
  const Key& k5, const Key& k6)
{
  MultiSet<Key> rtn = makeMultiSet<Key>(k1, k2, k3, k4, k5);
  rtn.put(k6);
  return rtn;
}

/** \relates MultiSet Create a multiset */
template<class Key> inline
MultiSet<Key> makeMultiSet(const Array<Key>& k)
{
  MultiSet<Key> rtn;
  for (int i=0; i<k.size(); i++) rtn.put(k[i]);
  return rtn;
}

}

namespace std
{
/** \relates Sundance::MultiSet 
 * Write to a stream
 */
template<class Key> inline
ostream& operator<<(std::ostream& os, const Sundance::MultiSet<Key>& m)
{return m.toStream(os);}
}


#endif

