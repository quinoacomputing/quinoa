/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/

#ifndef util_PairIter_hpp
#define util_PairIter_hpp

#include <utility>
#include <iterator>

namespace phdmesh {

/**
 *  \brief  Pair of begin and end iterators wrapped to provide
 *          a container-like view of the span.
 *
 *  \author H. Carter Edwards  <hcedwar@sandia.gov>
 *  \date   November 2007
 */
template< class IterType ,
          class IterCategory =
             typename std::iterator_traits< IterType >::iterator_category >
class PairIter {};

//----------------------------------------------------------------------
// Specialized for random access iterators, others TBD.

/**
 *  \brief  Pair of begin and end iterators wrapped to provide
 *          a container-like view of the span.
 */
template< class IterType >
class PairIter< IterType , std::random_access_iterator_tag > 
  : public std::pair< IterType , IterType >
{
private:
  typedef std::pair< IterType , IterType > Pair ;
  typedef PairIter< IterType , std::random_access_iterator_tag > Self ;
  typedef std::iterator_traits< IterType > Traits ;
public:

  //--------------------------------

  typedef          IterType                iterator ;
  typedef typename Traits::value_type      value_type ;
  typedef typename Traits::pointer         pointer ;
  typedef typename Traits::reference       reference ;
  typedef typename Traits::difference_type difference_type ;
  typedef          size_t                  size_type ;

  //--------------------------------

  ~PairIter() {}

  PairIter() : Pair() { Pair::second = Pair::first ; }

  PairIter( const Self & rhs ) : Pair( rhs ) {}

  PairIter( const Pair & rhs ) : Pair( rhs ) {}

  Self & operator = ( const Self & rhs )
    { Pair::first = rhs.first ; Pair::second = rhs.second ; return *this ; }

  Self & operator = ( const Pair & rhs )
    { Pair::first = rhs.first ; Pair::second = rhs.second ; return *this ; }

  //--------------------------------

  bool operator == ( const Self & rhs ) const
    { return Pair::first == rhs.first && Pair::second == rhs.second ; }

  bool operator != ( const Self & rhs ) const
    { return Pair::first != rhs.first || Pair::second != rhs.second ; }

  bool operator == ( const Pair & rhs ) const
    { return Pair::first == rhs.first && Pair::second == rhs.second ; }

  bool operator != ( const Pair & rhs ) const
    { return Pair::first != rhs.first || Pair::second != rhs.second ; }

  //--------------------------------

  Self & operator ++ () { ++ Pair::first ; return *this ; }

  Self operator ++ (int) { Self tmp(*this); ++ Pair::first ; return tmp ; }

  reference operator * ()  const { return * Pair::first ; }
  pointer   operator -> () const { return & * Pair::first ; }

  //--------------------------------
  // Container-like functionality for random access iterators.

  reference front() const { return * Pair::first ; }
  reference back()  const { return  Pair::second[-1] ; }

  iterator begin() const { return  Pair::first ; }
  iterator end()   const { return  Pair::second ; }

  template<class Iterator>
  PairIter( Iterator i , Iterator e ) : Pair(i,e) {}

  template<class Container>
  explicit
  PairIter( const Container & c ) : Pair( c.begin() , c.end() ) {}

  template<class Container>
  explicit
  PairIter( Container & c ) : Pair( c.begin() , c.end() ) {}

  bool empty () const { return ! ( Pair::first < Pair::second ) ; }

  operator bool () const { return Pair::first < Pair::second ; }

  reference operator [] ( size_t n ) const { return Pair::first[n] ; }

  size_t size() const
    {
      const difference_type d = std::distance( Pair::first , Pair::second );
      return d < 0 ? 0 : (size_t) d ;
    }
};

} // namespace phdmesh

#endif /* util_PairIter_hpp */

