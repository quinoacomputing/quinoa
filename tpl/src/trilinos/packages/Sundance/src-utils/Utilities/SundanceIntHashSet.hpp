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

#ifndef SUNDANCE_INTHASHSET_H
#define SUNDANCE_INTHASHSET_H

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "Teuchos_Array.hpp"
#include <list>


namespace Sundance
{
  using namespace Teuchos;

  /** 
   *
   */
  class IntHashSet
  {
  public:
    /** */
    IntHashSet();

    /** */
    void setCapacity(int capacity);

    /** */
    inline void put(int x) 
    {
      
      std::list<int>& d = data_[hashFunc(x)];
      for (std::list<int>::const_iterator i=d.begin(); i != d.end(); i++)
        {
          if (x == *i) return;
        }
      d.push_back(x);
      size_++;
    }

    /** */
    bool contains(int x) const ;

    /** */
    int size() const {return size_;}

    /** */
    void fillArray(int* a) const ;

  private:

    inline int hashFunc( const int key ) const { return (seed() ^ key)%capacity_; }

    static unsigned int seed() {static int rtn = (2654435761U); return rtn;}
    unsigned int capacity_;
    Array<std::list<int> > data_;
    int size_;
  };

}


#endif
