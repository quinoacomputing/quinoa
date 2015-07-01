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
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 */

#ifndef util_ParallelIndex_hpp
#define util_ParallelIndex_hpp

#include <utility>
#include <vector>

#include <util/Basics.hpp>
#include <util/Parallel.hpp>

namespace phdmesh {

/** Parallel cross-reference index for a collection of
 *  'key_type' keys.  Each processor constructs a
 *  ParallelIndex with its local collection of keys.
 *  The ParallelIndex may be queried for
 *  (1) which other processors that submitted the same keys or
 *  (2) which processors submitted an arbitrary set of keys.
 */
class ParallelIndex {
public:

  typedef uint64_type key_type ;

  typedef std::pair< key_type , key_type > KeyProc ;

  struct LessKeyProc {
    bool operator()( const KeyProc & lhs , const KeyProc & rhs ) const
      { return lhs < rhs ; }

    bool operator()( const KeyProc & lhs , const key_type rhs ) const
      { return lhs.first < rhs ; }
  };

  ~ParallelIndex();

  ParallelIndex( ParallelMachine , const std::vector<key_type> & );

  /** Query which other processors submitted the
   *  same keys that the local processor submitted.
   */
  void query( std::vector<KeyProc> & ) const ;

  /** Query which processors submitted the given keys.
   *  The local processor is in the output if it submitted a queried key.
   */
  void query( const std::vector<key_type> & , std::vector<KeyProc> & ) const ;

private:
  ParallelIndex();
  ParallelIndex( const ParallelIndex & );
  ParallelIndex & operator = ( const ParallelIndex & );

  ParallelMachine      m_comm ;
  std::vector<KeyProc> m_key_proc ;
};

}

//----------------------------------------------------------------------

#endif

