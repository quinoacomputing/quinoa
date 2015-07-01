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

#ifndef SUNDANCE_MAP_H
#define SUNDANCE_MAP_H

#include "SundanceDefs.hpp"
#include <map>

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  using namespace Teuchos;

  /**
   * Extension of STL map, adding some nicer put/get/contains syntax 
   * and an iostream insertion operator.
   */
  template<class Key, class Value, class Compare = std::less<Key> >
    class Map : public std::map<Key, Value, Compare>
    {
    public:
      /** */
      Map() : std::map<Key,Value,Compare>() {;}

      /** Test whether the specified key is present in the map */
      inline bool containsKey(const Key& key) const {return this->find(key) != this->end();}

      /** Put a new (key, value) entry in the map */
      inline void put(const Key& key, const Value& value)
        {this->operator[](key) = value;}

      /** Look up value and return a read-only reference */
      inline const Value& get(const Key& key) const
        {return (*(this->find)(key)).second;}

      /** Look up value and return a modifiable reference */
      inline Value& get(const Key& key) 
        {return (*(this->find)(key)).second;}
    };

}

namespace std
{
   /** \relates Sundance::Map 
    * Write to a stream
    */
  template<class Key, class Value, class Compare> inline
    std::ostream& operator<<(std::ostream& os, const Sundance::Map<Key, Value, Compare>& m)
    {
      typename Sundance::Map<Key, Value, Compare>::const_iterator iter;

      os << "Map[";
      int k = 0 ;
      for (iter=m.begin(); iter != m.end(); iter++, k++)
        {
          os << "{" << (*iter).first << ", " << (*iter).second << "}";
          if (k < (int) m.size()-1) os << ", ";
        }
      os << "]";
      return os;
    }
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
