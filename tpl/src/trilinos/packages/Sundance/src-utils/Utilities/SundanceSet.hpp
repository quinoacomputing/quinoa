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

#ifndef SUNDANCE_SET_H
#define SUNDANCE_SET_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include <set>
#include <algorithm>



namespace Sundance
{
using namespace Teuchos;

/** 
 * Extension of STL set, adding some nicer syntax 
 * and an iostream insertion operator.
 */
template<class Key, class Compare = std::less<Key> >
class Set 
{
public:
    
  typedef typename std::set<Key, Compare>::iterator iterator;
  typedef typename std::set<Key, Compare>::const_iterator const_iterator;
  typedef typename std::set<Key, Compare>::reverse_iterator reverse_iterator;
  typedef typename std::set<Key, Compare>::const_reverse_iterator const_reverse_iterator;

  typedef typename std::set<Key, Compare>::size_type size_type;
  typedef typename std::set<Key, Compare>::pointer pointer;
  typedef typename std::set<Key, Compare>::const_pointer const_pointer;
  typedef typename std::set<Key, Compare>::const_reference const_reference;
  typedef typename std::set<Key, Compare>::reference reference;
  typedef typename std::set<Key, Compare>::difference_type difference_type;
    
  /** */
  Set() : set_() {;}

  /** */
  Set(const std::set<Key>& s) : set_(s) {;}

  /** */
  iterator begin() {return set_.begin();}

  /** */
  const_iterator begin() const {return set_.begin();}

  /** */
  iterator end() {return set_.end();}

  /** */
  const_iterator end() const {return set_.end();}

  /** */
  reverse_iterator rbegin() {return set_.rbegin();}

  /** */
  const_reverse_iterator rbegin() const {return set_.rbegin();}

  /** */
  reverse_iterator rend() {return set_.rend();}

  /** */
  const_reverse_iterator rend() const {return set_.rend();}

  /** */
  std::pair<iterator, bool> insert(const Key& x) {return set_.insert(x);}

  /** */
  iterator insert(iterator pos, const Key& x) {return set_.insert(pos,x);}
  
  /** */
  template <typename InputIterator>
  void insert(InputIterator first, InputIterator last) {set_.insert(first,last);}

  /** */
  void erase(iterator pos) {set_.erase(pos);}

  /** */
  void erase(const Key& x) {set_.erase(x);}

  /** */
  void erase(iterator first, iterator last) {set_.erase(first, last);}

  /** */
  void clear() {set_.clear();}

  /** */
  iterator find(const Key& x) {return set_.find(x);}

  /** */
  const_iterator find(const Key& x) const {return set_.find(x);}

  /** */
  iterator lower_bound(const Key& x) {return set_.lower_bound(x);}

  /** */
  const_iterator lower_bound(const Key& x) const {return set_.lower_bound(x);}

  /** */
  iterator upper_bound(const Key& x) {return set_.upper_bound(x);}

  /** */
  const_iterator upper_bound(const Key& x) const {return set_.upper_bound(x);}

  /** */
  std::pair<iterator,iterator> equal_range(const Key& x)
    {return set_.equal_range(x);}

  /** */
  std::pair<const_iterator,const_iterator> equal_range(const Key& x) const 
    {return set_.equal_range(x);}

  /** */
  int size() const {return set_.size();}

  /** */
  int max_size() const {return set_.max_size();}

  /** */
  bool empty() const {return set_.empty();}
    
  /** */
  const std::set<Key, Compare>& set() const {return set_;}
    
  /** */
  std::set<Key, Compare>& set() {return set_;}

  /** Test whether the specified key is present in the set */
  bool contains(const Key& key) const {return this->find(key) != this->end();}

  /** Put a new entry in the set */
  void put(const Key& key) {insert(key);}

  /** Write into an array */
  Array<Key> elements() const ;

  /** */
  void elements(Array<Key>& keys) const ;

  /** */
  void merge(const Set<Key, Compare>& other);

  /** */
  Set<Key, Compare> intersection(const Set<Key, Compare>& other) const ;

  /** */
  Set<Key, Compare> setUnion(const Set<Key, Compare>& other) const ;

  /** */
  Set<Key, Compare> setDifference(const Set<Key, Compare>& other) const ;

  /** */
  std::ostream& toStream(std::ostream& os) const ;

  /** */
  std::string toString() const ;

private:
  std::set<Key, Compare> set_;
};


template<class Key, class Compare> inline
Array<Key> Set<Key, Compare>::elements() const
{
  Array<Key> rtn;

  typename Set<Key, Compare>::const_iterator iter;

  for (iter=this->begin(); iter != this->end(); iter++)
  {
    rtn.append(*iter);
  }
  return rtn;
}


template<class Key, class Compare> inline
void Set<Key, Compare>::elements(Array<Key>& rtn) const
{
  rtn.resize(0);
  typename Set<Key, Compare>::const_iterator iter;

  for (iter=this->begin(); iter != this->end(); iter++)
  {
    rtn.append(*iter);
  }
}

template<class Key, class Compare> inline
void Set<Key, Compare>::merge(const Set<Key, Compare>& other)
{
  typename Set<Key, Compare>::const_iterator iter;

  for (iter=other.begin(); iter != other.end(); iter++)
  {
    put(*iter);
  }
}

template<class Key, class Compare> inline
Set<Key, Compare> Set<Key, Compare>::intersection(const Set<Key, Compare>& other) const
{
  Set<Key, Compare> rtn;

  set_intersection(this->begin(), this->end(),
    other.begin(), other.end(), 
    std::insert_iterator<Set<Key, Compare> >(rtn, rtn.begin())); 
  return rtn;
}

template<class Key, class Compare> inline
Set<Key, Compare> Set<Key, Compare>::setUnion(const Set<Key, Compare>& other) const
{
  Set<Key, Compare> rtn;

  set_union(this->begin(), this->end(),
    other.begin(), other.end(), 
    std::insert_iterator<Set<Key, Compare> >(rtn, rtn.begin())); 
  return rtn;
}

template<class Key, class Compare> inline
Set<Key, Compare> Set<Key, Compare>::setDifference(const Set<Key, Compare>& other) const
{
  Set<Key, Compare> rtn;

  set_difference(this->begin(), this->end(),
    other.begin(), other.end(), 
    std::insert_iterator<Set<Key, Compare> >(rtn, rtn.begin())); 
  return rtn;
}

template<class Key, class Compare> inline
std::ostream& Set<Key, Compare>::toStream(std::ostream& os) const
{
  typename Set<Key, Compare>::const_iterator iter;

  int k = 0;
  os << "{";
  for (iter=this->begin(); iter != this->end(); iter++, k++)
  {
    os << *iter;
    if (k<(this->size()-1)) os << ", ";
  }
  os << "}";

  return os;
}

template<class Key, class Compare> inline
string Set<Key, Compare>::toString() const
{
  std::ostringstream os;
  os << *this;
  return os.str();
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Array<Key>& k)
{
  Set<Key> rtn;
  for (int i=0; i<k.size(); i++) rtn.put(k[i]);
  return rtn;
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Key& k)
{
  Set<Key> rtn;
  rtn.put(k);
  return rtn;
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Key& k1, const Key& k2)
{
  Set<Key> rtn = makeSet<Key>(k1);
  rtn.put(k2);
  return rtn;
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Key& k1, const Key& k2, const Key& k3)
{
  Set<Key> rtn = makeSet<Key>(k1, k2);
  rtn.put(k3);
  return rtn;
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Key& k1, const Key& k2, const Key& k3, const Key& k4)
{
  Set<Key> rtn = makeSet<Key>(k1, k2, k3);
  rtn.put(k4);
  return rtn;
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Key& k1, const Key& k2, const Key& k3, const Key& k4,
  const Key& k5)
{
  Set<Key> rtn = makeSet<Key>(k1, k2, k3, k4);
  rtn.put(k5);
  return rtn;
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Key& k1, const Key& k2, const Key& k3, const Key& k4,
  const Key& k5, const Key& k6)
{
  Set<Key> rtn = makeSet<Key>(k1, k2, k3, k4, k5);
  rtn.put(k6);
  return rtn;
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Key& k1, const Key& k2, const Key& k3, const Key& k4,
  const Key& k5, const Key& k6, const Key& k7)
{
  Set<Key> rtn = makeSet<Key>(k1, k2, k3, k4, k5, k6);
  rtn.put(k7);
  return rtn;
}

/** \relates Set Creates a set */
template<class Key> inline
Set<Key> makeSet(const Key& k1, const Key& k2, const Key& k3, const Key& k4,
  const Key& k5, const Key& k6, const Key& k7, const Key& k8)
{
  Set<Key> rtn = makeSet<Key>(k1, k2, k3, k4, k5, k6, k7);
  rtn.put(k8);
  return rtn;
}


/** \relates Set */
template <typename Key, typename Compare> bool operator==(
  const Set<Key, Compare>& L,
  const Set<Key, Compare>& R)
{
  return L.set() == R.set();
}

/** \relates Set */
template <typename Key, typename Compare> bool operator!=(
  const Set<Key, Compare>& L,
  const Set<Key, Compare>& R)
{
  return L.set() != R.set();
}

/** \relates Set */
template <typename Key, typename Compare> bool operator<=(
  const Set<Key, Compare>& L,
  const Set<Key, Compare>& R)
{
  return L.set() <= R.set();
}

/** \relates Set */
template <typename Key, typename Compare> bool operator<(
  const Set<Key, Compare>& L,
  const Set<Key, Compare>& R)
{
  return L.set() < R.set();
}


/** \relates Set */
template <typename Key, typename Compare> bool operator>(
  const Set<Key, Compare>& L,
  const Set<Key, Compare>& R)
{
  return L.set() > R.set();
}

/** \relates Set */
template <typename Key, typename Compare> bool operator>=(
  const Set<Key, Compare>& L,
  const Set<Key, Compare>& R)
{
  return L.set() >= R.set();
}



  
}

namespace std
{
/** \relates Sundance::Set */
template<class Key, class Compare> inline
ostream& operator<<(std::ostream& os, const Sundance::Set<Key, Compare>& m)
{return m.toStream(os);}
}


#endif
