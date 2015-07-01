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
 * @file
 * @author H. Carter Edwards
 * @date   January 2007
 */

#include <iostream>
#include <map>
#include <set>
#include <list>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <util/TPI.hpp>
#include <util/Basics.hpp>
#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <util/OctTreeOps.hpp>

namespace phdmesh {
namespace {

#if defined( HAVE_MPI )

void all_gather( ParallelMachine comm ,
                 const float * local , float * global , unsigned count )
{
  MPI_Allgather( const_cast<float*>( local ) , count , MPI_FLOAT ,
                 global , count , MPI_FLOAT , comm );
}

#else

void all_gather( ParallelMachine ,
                 const float * local , float * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

#endif

}
}


namespace phdmesh {

typedef std::map< OctTreeKey ,
                  std::pair< std::list< IdentProcBox > ,
                             std::list< IdentProcBox > > > SearchTree ;

namespace {

void search_tree_statistics( ParallelMachine  arg_comm ,
                             const SearchTree & s ,
                             unsigned * const data )
{
  const unsigned huge = std::numeric_limits<unsigned>::max();
  unsigned avg[2] = { 0 , 0 };
  unsigned max[2] = { 0 , 0 };
  unsigned min[2] ;
  min[0] = min[1] = huge ;

  SearchTree::const_iterator i ;

  for ( i = s.begin() ; i != s.end() ; ++i ) {
    const SearchTree::value_type  & inode = *i ;
    const unsigned d_size = inode.second.first.size();
    const unsigned r_size = inode.second.second.size();

    avg[0] += d_size ;
    if ( d_size < min[0] ) { min[0] = d_size ; }
    if ( max[0] < d_size ) { max[0] = d_size ; }

    avg[1] += r_size ;
    if ( r_size < min[1] ) { min[1] = r_size ; }
    if ( max[1] < r_size ) { max[1] = r_size ; }
  }

  // Average for this processor

  const unsigned ncells = s.size() ;

  if ( ncells ) {
    avg[0] = ( avg[0] + ncells - 1 ) / ncells ;
    avg[1] = ( avg[1] + ncells - 1 ) / ncells ;
  }

  if ( min[0] == huge ) { min[0] = 0 ; }
  if ( min[1] == huge ) { min[1] = 0 ; }

  all_reduce( arg_comm , Min<2>( min ) , Max<2>( max ) , Sum<2>( avg ) );

  const unsigned p_size = parallel_machine_size( arg_comm );

  // Average among all processors:

  avg[0] = ( avg[0] + p_size - 1 ) / p_size ;
  avg[1] = ( avg[1] + p_size - 1 ) / p_size ;

  data[0] = min[0] ;
  data[1] = max[0] ;
  data[2] = avg[0] ;

  data[3] = min[1] ;
  data[4] = max[1] ;
  data[5] = avg[1] ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool hsfc_box_covering( const float * const global_box ,
                        const float * const small_box ,
                        OctTreeKey * const covering ,
                        unsigned   &       number )
{
  enum { Dimension    = 3 };
  enum { Combinations = 8 };

  const double min = std::numeric_limits<float>::epsilon();
  const double max = 1.0 - min ;

  // Determine the unit-box bounds and bisection depth for the box

  double ubox_low[ Dimension ] ;
  double ubox_up[  Dimension ] ;

  bool valid = true ;

  // Determine unit box and is maximum length
  // The box is bounded by [eps,1-eps].

  double unit_size = 0.0 ;

  for ( unsigned i = 0 ; i < Dimension ; ++i ) {

    const float global_low = global_box[i] ;
    const float global_up  = global_box[i+Dimension] ;
    const float small_low  = small_box[i] ;
    const float small_up   = small_box[i+Dimension] ;

    if ( small_up < global_low ) {
      // Entirely less than 'min'
      ubox_low[i] = ubox_up[i] = min ;
      valid = false ;
    }
    else if ( global_up < small_low ) {
      // Entirely greater than 'max'
      ubox_low[i] = ubox_up[i] = max ;
      valid = false ;
    }
    else {
      const double scale = 1.0 / ( global_up - global_low );

      double unit_low = ( small_low - global_low ) * scale ;
      double unit_up  = ( small_up  - global_low ) * scale ;

      if ( unit_low < min ) {
        unit_low = min ;
        valid = false ;
      }

      if ( max < unit_up ) {
        unit_up = max ;
        valid = false ;
      }

      if ( unit_up < unit_low ) {
        // A negative volume, consider it a point at the lower
        unit_up = unit_low ;
        valid = false ;
      }
      else {
        const double tmp_size = unit_up - unit_low ;
        if ( unit_size < tmp_size ) { unit_size = tmp_size ; }
      }

      ubox_low[i] = unit_low ;
      ubox_up[i]  = unit_up ;
    }
  }

  // Depth is determined by smallest cell depth
  // that could contain the small_box

  unsigned depth = OctTreeKey::MaxDepth ;

  if ( 0 < unit_size ) {
    const double size_inv = 1.0 / unit_size ;
    depth = 1 ;
    while ( depth < OctTreeKey::MaxDepth && ( 1 << depth ) < size_inv ) {
      ++depth ;
    }
    --depth ;
  }

  // Determine the oct-tree nodes for each key

  const unsigned shift    = OctTreeKey::BitsPerWord - depth ;
  const unsigned num_cell = 1 << depth ;

  // At most two cells in each axis at this depth

  unsigned coord_low[ Dimension ];
  unsigned coord_up[  Dimension ];

  for ( unsigned i = 0 ; i < Dimension ; ++i ) {
    const unsigned low = (unsigned)( ubox_low[i] * num_cell );
    const unsigned up  = (unsigned)( ubox_up[i]  * num_cell );

    if ( low + 1 < up ) {
      std::string msg("phdmesh::hsfc_box_covering FAILED : depth determination logic error");
      throw std::logic_error( msg );
    }

    coord_low[i] = low << shift ;
    coord_up[i]  = up  << shift ;
  }

  unsigned n = 0 ;

  for ( unsigned i = 0 ; i < Combinations ; ++i ) {

    const bool duplicate =
      ( ( i & 01 ) && coord_up[0] == coord_low[0] ) ||
      ( ( i & 02 ) && coord_up[1] == coord_low[1] ) ||
      ( ( i & 04 ) && coord_up[2] == coord_low[2] ) ;

    if ( ! duplicate ) {
      unsigned coord[3] ;

      coord[0] = ( i & 01 ) ? coord_up[0] : coord_low[0] ;
      coord[1] = ( i & 02 ) ? coord_up[1] : coord_low[1] ;
      coord[2] = ( i & 04 ) ? coord_up[2] : coord_low[2] ;

      covering[n] = hsfc3d( depth , coord );

      ++n ;
    }
  }

  number = n ;

  return valid ;
}

//----------------------------------------------------------------------

namespace {

bool intersect( const IdentProcBox & a , const IdentProcBox & b )
{
  // If both 'part' values are nonzero and equal then do not intersect

  const bool outside =
    ( a.part && b.part && a.part == b.part ) ||
    ( b.box[3] < a.box[0] ) || ( a.box[3] < b.box[0] ) ||
    ( b.box[4] < a.box[1] ) || ( a.box[4] < b.box[1] ) ||
    ( b.box[5] < a.box[2] ) || ( a.box[5] < b.box[2] ) ;
  return ! outside ;
}

//----------------------------------------------------------------------

template< class S >
struct SetInsertBuffer {
  enum { N = 128 };
  S      & m_set ;
  TPI::ThreadPool m_pool ;
  unsigned m_lock ;
  unsigned m_iter ;
  typename S::value_type m_buffer[ N ] ;

  void overflow();
  void operator()( const typename S::value_type & v );

  SetInsertBuffer( S & s , TPI::ThreadPool p , unsigned l )
    : m_set(s), m_pool(p), m_lock(l), m_iter(0) {}

  ~SetInsertBuffer() { overflow(); }
};

template<class S>
void SetInsertBuffer<S>::overflow()
{
  TPI::LockGuard update_lock( m_pool , m_lock );
  while ( m_iter ) {
    m_set.insert( m_buffer[ --m_iter ] );
  }
}

template<class S>
void SetInsertBuffer<S>::operator()( const typename S::value_type & v )
{
  m_buffer[ m_iter ] = v ;
  if ( N == ++m_iter ) { overflow(); }
}

void proximity_search_symmetric(
  const SearchTree::const_iterator i_beg ,
  const SearchTree::const_iterator i_end ,
  SetInsertBuffer< std::set< std::pair<IdentProc,IdentProc> > > & arg_out )
{
  SearchTree::const_iterator j ;

  std::list<IdentProcBox>::const_iterator id , ir ;

  const SearchTree::value_type  & inode = *i_beg ;
  const std::list<IdentProcBox> & domain_outer = inode.second.first ;

  const std::list<IdentProcBox>::const_iterator
    beg_dom_out = domain_outer.begin(),
    end_dom_out = domain_outer.end();

  // Outer cell vs. itself

  for ( id = beg_dom_out ; id != end_dom_out ; ++id ) {
    const IdentProcBox & d = *id ;
    for ( ir = id ; ++ir != end_dom_out ; ) {
      const IdentProcBox & r = *ir ;
      if ( intersect( d , r ) ) {
        const IdentProc & dip = d ;
        const IdentProc & rip = r ;
        if ( dip < rip ) {
          std::pair<IdentProc,IdentProc> tmp( dip , rip );
          arg_out( tmp );
        }
        else {
          std::pair<IdentProc,IdentProc> tmp( rip , dip );
          arg_out( tmp );
        }
      }
    }
  }

  // Outer cell searching inner cells.
  // Outer cell always precedes inner cells
  // Iterate forward until the cell is not contained.

  const OctTreeKey & outer_key = inode.first ;

  for ( j = i_beg ; ++j != i_end && outer_key.intersect( (*j).first ) ; ) {

    const SearchTree::value_type & jnode = *j ;

    const std::list<IdentProcBox> & domain_inner = jnode.second.first ;

    const std::list<IdentProcBox>::const_iterator
      beg_dom_inn = domain_inner.begin(),
      end_dom_inn = domain_inner.end();

    // Check domain_outer vs. domain_inner,
    // skip if the same box.

    for ( id = beg_dom_out ; id != end_dom_out ; ++id ) {
      const IdentProcBox & d = *id ;
      for ( ir = beg_dom_inn ; ir != end_dom_inn ; ++ir ) {
        const IdentProcBox & r = *ir ;
        if ( d.ident != r.ident || d.proc != r.proc ) {
          if ( intersect( d , r ) ) {
            const IdentProc & dip = d ;
            const IdentProc & rip = r ;
            if ( dip < rip ) {
              std::pair<IdentProc,IdentProc> tmp( dip , rip );
              arg_out( tmp );
            }
            else {
              std::pair<IdentProc,IdentProc> tmp( rip , dip );
              arg_out( tmp );
            }
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void proximity_search_asymmetric(
  const SearchTree::const_iterator i_beg ,
  const SearchTree::const_iterator i_end ,
  SetInsertBuffer< std::set< std::pair<IdentProc,IdentProc> > > & arg_out )
{
  SearchTree::const_iterator j ;

  std::list<IdentProcBox>::const_iterator id , ir ;

  const SearchTree::value_type & inode = *i_beg ;

  const std::list<IdentProcBox> & domain_outer = inode.second.first ;
  const std::list<IdentProcBox> & range_outer  = inode.second.second ;

  const std::list<IdentProcBox>::const_iterator
    beg_dom_out = domain_outer.begin(),
    end_dom_out = domain_outer.end(),
    beg_ran_out = range_outer.begin(),
    end_ran_out = range_outer.end();

  // domain_outer vs. range_outer

  for ( id = beg_dom_out ; id != end_dom_out ; ++id ) {
    const IdentProcBox & d = *id ;
    for ( ir = beg_ran_out ; ir != end_ran_out ; ++ir ) {
      const IdentProcBox & r = *ir ;
      if ( intersect( d , r ) ) {
        std::pair<IdentProc,IdentProc> tmp( d , r );
        arg_out( tmp );
      }
    }
  }

  // Outer cell searching inner cells.
  // Outer cell always precedes inner cells
  // Iterate forward until the cell is not contained.

  const OctTreeKey & outer_key = inode.first ;

  for ( j = i_beg ; ++j != i_end && outer_key.intersect( (*j).first ) ; ) {

    const SearchTree::value_type & jnode = *j ;

    const std::list<IdentProcBox> & domain_inner = jnode.second.first ;
    const std::list<IdentProcBox> & range_inner  = jnode.second.second ;

    const std::list<IdentProcBox>::const_iterator
      beg_dom_inn = domain_inner.begin(),
      end_dom_inn = domain_inner.end(),
      beg_ran_inn = range_inner.begin(),
      end_ran_inn = range_inner.end();

    // Check domain_outer vs. range_inner

    for ( id = beg_dom_out ; id != end_dom_out ; ++id ) {
      const IdentProcBox & d = *id ;
      for ( ir = beg_ran_inn ; ir != end_ran_inn ; ++ir ) {
        const IdentProcBox & r = *ir ;
        if ( intersect( d , r ) ) {
          std::pair<IdentProc,IdentProc> tmp( d , r );
          arg_out( tmp );
        }
      }
    }

    // Check domain_inner vs. range_outer if non-symmetric

    for ( id = beg_dom_inn ; id != end_dom_inn ; ++id ) {
      const IdentProcBox & d = *id ;
      for ( ir = beg_ran_out ; ir != end_ran_out ; ++ir ) {
        const IdentProcBox & r = *ir ;
        if ( intersect( d , r ) ) {
          std::pair<IdentProc,IdentProc> tmp( d , r );
          arg_out( tmp );
        }
      }
    }
  }
}

//----------------------------------------------------------------------

typedef void (*proximity_search_work_routine)(
  const SearchTree::const_iterator i_beg ,
  const SearchTree::const_iterator i_end ,
  SetInsertBuffer< std::set< std::pair<IdentProc,IdentProc> > > & );

class ProximitySearch {
public:
  enum { NLOCKS = 2 };

  const proximity_search_work_routine m_search ;

  std::set< std::pair<IdentProc,IdentProc> > & m_relation ;
  SearchTree::const_iterator m_tree_iter ;
  SearchTree::const_iterator m_tree_end ;
  
  ~ProximitySearch() {}

  ProximitySearch(
    bool symmetric ,
    const SearchTree & search_tree ,
    std::set< std::pair<IdentProc,IdentProc> > & relation );

  void iterate_tree(TPI::ThreadPool);

private:
  ProximitySearch();
  ProximitySearch( const ProximitySearch & );
  ProximitySearch & operator = ( const ProximitySearch & );
};


void ProximitySearch::iterate_tree(TPI::ThreadPool pool)
{
  enum { N_WORK = 32 };

  try {
    SetInsertBuffer< std::set< std::pair<IdentProc,IdentProc> > >
      tmp( m_relation , pool , 1 );

    const SearchTree::const_iterator i_tree_end = m_tree_end ;

    unsigned n_work = N_WORK ;

    while ( n_work ) {

      n_work = 0 ;

      SearchTree::const_iterator i_beg , i_end ;

      { // Get work:
        TPI::LockGuard get_work_lock( pool , 0 );

        i_end = i_beg = m_tree_iter ;

        while ( n_work < N_WORK && i_tree_end != i_end ) {
          ++i_end ; ++n_work ;
        }
    
        m_tree_iter = i_end ;

      } // get_work_lock is released

      // Perform work:

      for ( ; i_beg != i_end ; ++i_beg ) {
        (*m_search)( i_beg , i_tree_end , tmp );
      }
    }
  }
  catch ( const std::exception & x ) {
    std::cerr << x.what() << std::endl ;
    std::cerr.flush();
  }
  catch ( ... ) {
    std::cerr << "ProximitySearch::iterate_tree FAILED" << std::endl ;
    std::cerr.flush();
  }
}

ProximitySearch::ProximitySearch(
  bool symmetric ,
  const SearchTree & search_tree ,
  std::set< std::pair<IdentProc,IdentProc> > & relation )
: m_search( symmetric ? & proximity_search_symmetric
                      : & proximity_search_asymmetric ) ,
  m_relation( relation ),
  m_tree_iter( search_tree.begin() ),
  m_tree_end(  search_tree.end() )
{
  if ( m_tree_iter != m_tree_end ) {
    TPI::Set_lock_size( NLOCKS );
    TPI::Run( *this , & ProximitySearch::iterate_tree );

    if ( m_tree_iter != m_tree_end ) {
      std::string msg("phdmesh::proximity_search FAILED to complete" );
      throw std::runtime_error(msg);
    }
  }
}

//----------------------------------------------------------------------
// Reset the accumulated node weights to only include
// those nodes in the range [ k_first , k_last ]

void accumulate_weights(
  OctTreeKey k_node ,
  OctTreeKey k_first ,
  const unsigned ord_end ,
  const unsigned depth ,
        float * const weights )
{
  const unsigned ord_node_2 = 2 * oct_tree_offset( depth , k_node );

  if ( k_node.depth() < depth ) {

    double w = 0 ;

    const unsigned d1 = k_node.depth() + 1 ;

    unsigned i = k_first.index( d1 );

    if ( i ) {
      k_node.set_index( d1 , i );

      const unsigned ord = oct_tree_offset( depth , k_node );
      const unsigned ord_2 = ord * 2 ;

      accumulate_weights( k_node , k_first , ord_end , depth , weights );

      // Counts of this node and all of its descending nodes
      w += weights[ord_2] + weights[ ord_2 + 1 ] ;

      k_first = OctTreeKey(); // Done with the lower bound
    }

    for ( ++i ; i <= 8 ; ++i ) {

      k_node.set_index( d1 , i );

      const unsigned ord = oct_tree_offset( depth , k_node );
      const unsigned ord_2 = ord * 2 ;

      if ( ord < ord_end ) {
        accumulate_weights( k_node, k_first , ord_end , depth , weights );

        // Counts of this node and all of its descending nodes
        w += weights[ord_2] + weights[ ord_2 + 1 ] ;
      }
    }

    // Descending node weight

    weights[ ord_node_2 + 1 ] = (float) w ; 
  }
}

//----------------------------------------------------------------------

void oct_key_split(
  const OctTreeKey & key ,
  const unsigned     upper_ord ,
        OctTreeKey & key_upper )
{
  // Split key at key.depth() + 1

  unsigned d = key.depth();

  key_upper = key ;

  if ( upper_ord == 1 ) { // key_upper gets it all
    while ( d && 1 == key_upper.index(d) ) {
      key_upper.clear_index(d);
      --d ;
    }
  }
  else if ( 8 < upper_ord ) { // key_upper get none of it, Increment key_upper

    unsigned i = 0 ;
    while ( d && 8 == ( i = key_upper.index(d) ) ) {
      key_upper.clear_index(d);
      --d ;
    }
    if ( d ) { key_upper.set_index( d , i + 1 ); }
  }
  else {
    key_upper.set_index( d + 1 , upper_ord );
  }
}

//----------------------------------------------------------------------

void partition( 
  const OctTreeKey & k_first ,
  const unsigned     i_end ,
  const OctTreeKey & key ,
  const unsigned     depth ,
  const float      * weights ,
  const double tolerance ,
  const double target_ratio ,
  double w_lower ,
  double w_upper ,
  OctTreeKey & k_upper )
{
  const unsigned ord_node = oct_tree_offset( depth , key );
  const float * const w_node = weights + ord_node * 2 ;

  const unsigned d1 = key.depth() + 1 ;

  // Add weights from nested nodes and their descendents
  // Try to achieve the ratio.

  const unsigned i_first = k_first.index( d1 );

  unsigned i = ( i_first ) ? i_first : 1 ;
  unsigned j = 8 ;
  {
    OctTreeKey k_upp = key ;
    k_upp.set_index( d1 , j );
    while ( i_end <= oct_tree_offset( depth , k_upp ) ) {
      k_upp.set_index( d1 , --j );
    }
  }

  w_lower += w_node[0] ;
  w_upper += w_node[0] ;

  // At the maximum depth?

  if ( key.depth() == depth ) {
    // Assume weight from unrepresented nested nodes is
    // evenly distributed among the nodes in the span [i,j]

    const unsigned n = 1 + j - i ;

    const double val = ((double) w_node[1]) / ((double) n);

    // val = val_lower + val_upper
    // ( w_lower + val_lower ) / ( w_upper + val_upper ) == target_ratio

    const double val_lower =
      ( target_ratio * ( w_upper + val ) - w_lower ) /
      ( target_ratio + 1 ) ;

    if ( 0 < val_lower ) {
      // How much of the range does the lower portion get?
      // Roundoff instead of merely truncating:
      i += (unsigned)( 0.5 + ( n * val_lower ) / val );

      // Can only get up to the maximum
      if ( j < i ) { i = j ; }
    }
    oct_key_split( key , i , k_upper );
  }
  else {

    while ( i != j ) {
      OctTreeKey ki = key ; ki.set_index( d1 , i );
      OctTreeKey kj = key ; kj.set_index( d1 , j );

      const float * const vi = weights + 2 * oct_tree_offset( depth , ki );
      const float * const vj = weights + 2 * oct_tree_offset( depth , kj );

      const double vali = vi[0] + vi[1] ;
      const double valj = vj[0] + vj[1] ;

      if ( 0 < vali && 0 < valj ) {

        // Choose between ( w_lower += vali ) vs. ( w_upper += valj )
        // Knowing that the skipped value will be revisited.

        if ( ( w_lower + vali ) < target_ratio * ( w_upper + valj ) ) {
          // Add to 'w_lower' and will still need more later
          w_lower += vali ;
          ++i ;
        }
        else {
           // Add to 'w_upper' and will still need more later
          w_upper += valj ;
          --j ;
        }
      }
      else {
        if ( vali <= 0.0 ) { ++i ; }
        if ( valj <= 0.0 ) { --j ; }
      }
    }

    // If 'i' has not incremented then 'k_first' is still in force
    OctTreeKey nested_k_first ;
    if ( i_first == i ) { nested_k_first = k_first ; }

    // Split node nested[i] ?
    OctTreeKey ki = key ; ki.set_index( d1 , i );

    const float * const vi = weights + 2 * oct_tree_offset( depth , ki );
    const double vali = vi[0] + vi[1] ;

    double diff = 0.0 ;

    if ( vali <= 0.0 ) {
      diff = 0.0 ; // Nothing can be done.  Give 'i' to the upper range
    }
    else if ( w_lower < w_upper * target_ratio ) {
      // Try adding to w_lower
      diff = ((double) (w_lower + vali)) / ((double)w_upper) - target_ratio ;
      ++i ;
    }
    else {
      // Try adding to w_upper
      diff = ((double)w_lower) / ((double)(w_upper + vali)) - target_ratio ;
    }

    if ( - tolerance < diff && diff < tolerance ) {
      oct_key_split( key , i , k_upper );
    }
    else {
      partition( nested_k_first , i_end , ki ,
                 depth , weights ,
                 tolerance , target_ratio ,
                 w_lower , w_upper , k_upper );
    }
  }
}

//----------------------------------------------------------------------

unsigned processor( const OctTreeKey * const cuts_b ,
                    const OctTreeKey * const cuts_e ,
                    const OctTreeKey & key )
{
  const OctTreeKey * const cuts_p = std::upper_bound( cuts_b , cuts_e , key );

  if ( cuts_p == cuts_b ) {
    std::string msg("phdmesh::processor FAILED: Bad cut-key array");
    throw std::runtime_error(msg);
  }

  return ( cuts_p - cuts_b ) - 1 ;
}

//----------------------------------------------------------------------

void pack( CommAll & comm_all ,
           const OctTreeKey * const cuts_b ,
           const SearchTree & send_tree ,
                 SearchTree * const recv_tree )
{
  const unsigned p_rank = comm_all.parallel_rank();
  const unsigned p_size = comm_all.parallel_size();
  const OctTreeKey * const cuts_e = cuts_b + p_size ;

  SearchTree::const_iterator i ;

  for ( i = send_tree.begin() ; i != send_tree.end() ; ++i ) {
    const OctTreeKey & key = (*i).first ;

    unsigned p = processor( cuts_b , cuts_e , key );

    do {
      if ( p != p_rank ) {
        CommBuffer & buf = comm_all.send_buffer(p);

        const std::list< IdentProcBox > & domain = (*i).second.first ;
        const std::list< IdentProcBox > & range  = (*i).second.second ;

        std::list< IdentProcBox >::const_iterator j ;

        const unsigned dsize = domain.size();
        const unsigned rsize = range.size();

        buf.pack<unsigned>( key.value() , OctTreeKey::NWord );
        buf.pack<unsigned>( dsize );
        buf.pack<unsigned>( rsize );

        for ( j = domain.begin() ; j != domain.end() ; ++j ) {
          const IdentProcBox & box = *j ;
          buf.pack<unsigned>( box.ident );
          buf.pack<unsigned>( box.proc );
          buf.pack<unsigned>( box.part );
          buf.pack<float>( box.box , 6 );
        }

        for ( j = range.begin() ; j != range.end() ; ++j ) {
          const IdentProcBox & box = *j ;
          buf.pack<unsigned>( box.ident );
          buf.pack<unsigned>( box.proc );
          buf.pack<unsigned>( box.part );
          buf.pack<float>( box.box , 6 );
        }
      }
      else if ( recv_tree ) {
        // Copy contents of the send node
        (*recv_tree)[ key ] = (*i).second ;
      }

      // If the cut keys are at a finer grantularity than
      // this key then this key may overlap more than one
      // processor's span.  Check for overlap with the
      // beginning key of the next processor.

      ++p ;

    } while( p < p_size && key.intersect( cuts_b[p] ) );
  }
}

void unpack( CommAll & comm_all , SearchTree & tree )
{
  unsigned domain_size ;
  unsigned range_size ;
  unsigned value[ OctTreeKey::NWord ];
  OctTreeKey key ;
  IdentProcBox box ;

  const unsigned p_size = comm_all.parallel_size();

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = comm_all.recv_buffer(p);

    while ( buf.remaining() ) {
      buf.unpack<unsigned>( value , OctTreeKey::NWord );
      buf.unpack<unsigned>( domain_size );
      buf.unpack<unsigned>( range_size );

      // Insert key, get domain and range

      key.set_value( value );

      SearchTree::mapped_type & node = tree[ key ];

      std::list< IdentProcBox > & domain = node.first ;
      std::list< IdentProcBox > & range  = node.second ;

      for ( unsigned j = 0 ; j < domain_size ; ++j ) {
        buf.unpack<unsigned>( box.ident );
        buf.unpack<unsigned>( box.proc );
        buf.unpack<unsigned>( box.part );
        buf.unpack<float>( box.box , 6 );

        domain.push_back( box );
      }

      for ( unsigned j = 0 ; j < range_size ; ++j ) {
        buf.unpack<unsigned>( box.ident );
        buf.unpack<unsigned>( box.proc );
        buf.unpack<unsigned>( box.part );
        buf.unpack<float>( box.box , 6 );

        range.push_back( box );
      }
    }
  }
}

//----------------------------------------------------------------------

bool communicate( 
  ParallelMachine arg_comm ,
  const OctTreeKey * const arg_cuts ,
  const SearchTree & send_tree ,
        SearchTree & recv_tree ,
  const bool local_flag ) 
{
  const unsigned p_size = parallel_machine_size( arg_comm );

  // Communicate search_tree members

  CommAll comm_all( arg_comm );

  // Sizing pass for pack
  pack( comm_all , arg_cuts , send_tree , NULL );

  // If more than 25% then is dense
  const bool global_flag =
    comm_all.allocate_buffers( p_size / 4 , false , local_flag );

  // Actual packing pass, copy local entries too
  pack( comm_all , arg_cuts , send_tree , & recv_tree );

  comm_all.communicate();

  unpack( comm_all , recv_tree );

  return global_flag ;
}

//----------------------------------------------------------------------

void communicate(
  ParallelMachine arg_comm ,
  const std::set< std::pair<IdentProc,IdentProc> > & send_relation ,
        std::set< std::pair<IdentProc,IdentProc> > & recv_relation )
{
  typedef std::pair<IdentProc,IdentProc> ValueType ;

  CommAll comm_all( arg_comm );

  const unsigned p_rank = comm_all.parallel_rank();
  const unsigned p_size = comm_all.parallel_size();

  std::set< ValueType >::const_iterator i ;

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) {
    const ValueType & val = *i ;
    if ( val.first.proc == p_rank || val.second.proc == p_rank ) {
      recv_relation.insert( val );
    }
    if ( val.first.proc != p_rank ) {
      CommBuffer & buf = comm_all.send_buffer( val.first.proc );
      buf.skip<unsigned>( 4 );
    }
    if ( val.second.proc != p_rank && val.second.proc != val.first.proc ) {
      CommBuffer & buf = comm_all.send_buffer( val.second.proc );
      buf.skip<unsigned>( 4 );
    }
  }

  // If more than 25% messages then is dense

  comm_all.allocate_buffers( p_size / 4 , false );

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) {
    const ValueType & val = *i ;
    if ( val.first.proc != p_rank ) {
      CommBuffer & buf = comm_all.send_buffer( val.first.proc );
      buf.pack<unsigned>( val.first.ident );
      buf.pack<unsigned>( val.first.proc );
      buf.pack<unsigned>( val.second.ident );
      buf.pack<unsigned>( val.second.proc );
    }
    if ( val.second.proc != p_rank && val.second.proc != val.first.proc ) {
      CommBuffer & buf = comm_all.send_buffer( val.second.proc );
      buf.pack<unsigned>( val.first.ident );
      buf.pack<unsigned>( val.first.proc );
      buf.pack<unsigned>( val.second.ident );
      buf.pack<unsigned>( val.second.proc );
    }
  }

  comm_all.communicate();

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = comm_all.recv_buffer( p );
    while ( buf.remaining() ) {
      ValueType val ;
      buf.unpack<unsigned>( val.first.ident );
      buf.unpack<unsigned>( val.first.proc );
      buf.unpack<unsigned>( val.second.ident );
      buf.unpack<unsigned>( val.second.proc );
      recv_relation.insert( val );
    }
  }
}

}

//----------------------------------------------------------------------

void oct_tree_partition_private(
  const unsigned p_first ,
  const unsigned p_end ,
  const unsigned depth ,
  const double   tolerance ,
  float * const weights ,
  const unsigned cuts_length ,
  OctTreeKey * const cuts )
{
  // split tree between [ p_first , p_end )
  const unsigned p_size  = p_end - p_first ;
  const unsigned p_upper = ( p_end + p_first ) / 2 ;

  const double target_fraction =
    ((double) ( p_upper - p_first ) ) / ((double) p_size );

  const double target_ratio = target_fraction / ( 1.0 - target_fraction );

  // Determine k_lower and k_upper such that
  //
  // Weight[ k_first , k_lower ] / Weight [ k_upper , k_last ] == target_ratio
  //
  // Within a tollerance

  const OctTreeKey k_first = cuts[ p_first ];

  const unsigned i_end   =
    p_end < cuts_length ? oct_tree_offset( depth , cuts[ p_end ] )
                        : oct_tree_size( depth );

  // Walk the tree [ k_first , k_last ] and accumulate weight

  accumulate_weights( OctTreeKey() , k_first , i_end , depth , weights );

  OctTreeKey k_root ;
  OctTreeKey & k_upper = cuts[ p_upper ] ;

  unsigned w_lower = 0 ;
  unsigned w_upper = 0 ;

  partition( k_first, i_end, k_root ,
             depth, weights,
             tolerance, target_ratio,
             w_lower, w_upper, k_upper );

  const bool nested_lower_split = p_first + 1 < p_upper ;
  const bool nested_upper_split = p_upper + 1 < p_end ;

  // If splitting both lower and upper, and a thread is available
  // then one of the next two calls could be a parallel thread
  // with a local copy of the shared 'weights' array.

  if ( nested_lower_split ) {
    oct_tree_partition_private( p_first, p_upper, depth,
                                tolerance, weights, cuts_length, cuts );
  }

  if ( nested_upper_split ) {
    oct_tree_partition_private( p_upper, p_end, depth,
                                tolerance, weights, cuts_length, cuts );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void box_global_bounds(
  ParallelMachine            arg_comm ,
  const unsigned             arg_domain_boxes_number ,
  const IdentProcBox * const arg_domain_boxes ,
  const unsigned             arg_range_boxes_number ,
  const IdentProcBox * const arg_range_boxes ,
        float        * const arg_global_box )
{
  enum { Dim = 3 };

  const bool symmetric = arg_range_boxes == arg_domain_boxes ;

  Copy<Dim>( arg_global_box ,         std::numeric_limits<float>::max() );
  Copy<Dim>( arg_global_box + Dim , - std::numeric_limits<float>::max() );

  //------------------------------------
  // Trivial loop threading possible:

  for ( unsigned i = 0 ; i < arg_domain_boxes_number ; ++i ) {
    const float * const box = arg_domain_boxes[i].box ;
    Min<Dim>( arg_global_box ,       box );
    Max<Dim>( arg_global_box + Dim , box + Dim );
  }

  if ( ! symmetric ) {
    for ( unsigned i = 0 ; i < arg_range_boxes_number ; ++i ) {
      const float * const box = arg_range_boxes[i].box ;
      Min<Dim>( arg_global_box ,       box );
      Max<Dim>( arg_global_box + Dim , box + Dim );
    }
  }

  //------------------------------------

  all_reduce( arg_comm , Min<Dim>( arg_global_box ) ,
                         Max<Dim>( arg_global_box + Dim ) );

  // Scale up and down by epsilon

  const double eps = std::numeric_limits<float>::epsilon();

  for ( unsigned i = 0 ; i < Dim ; ++i ) {
    float upper = arg_global_box[i+Dim] ;
    float lower = arg_global_box[i] ;

    double delta = eps * ( upper - lower );

    while ( upper <= arg_global_box[i+Dim] ||
                     arg_global_box[i] <= lower ) {
      upper = (float)( upper + delta );
      lower = (float)( lower - delta );
      delta *= 2 ;
    }

    arg_global_box[i+Dim] = (float) upper ;
    arg_global_box[i]     = (float) lower ;
  }
}

//----------------------------------------------------------------------
// Partition a search tree among processors

void oct_tree_partition( 
  ParallelMachine             arg_comm ,
  const SearchTree          & arg_tree ,
  const double                arg_tolerance ,
  std::vector< OctTreeKey > & arg_cuts )
{
  enum { tree_depth  = 4 };
  enum { tree_size   = OctTreeSize< tree_depth >::value };
  enum { tree_size_2 = tree_size * 2 };

  const unsigned p_size = parallel_machine_size( arg_comm );
  const OctTreeKey k_null ;

  arg_cuts.assign( p_size , k_null );

  float local_count[  tree_size_2 ];
  float global_count[ tree_size_2 ];

  for ( unsigned i = 0 ; i < tree_size_2 ; ++i ) {
    local_count[i] = 0.0 ;
  }

  for ( SearchTree::const_iterator i =  arg_tree.begin() ;
                                   i != arg_tree.end() ; ++i ) {

    const OctTreeKey & key = (*i).first ;

    const std::list< IdentProcBox > & domain = (*i).second.first ;
    const std::list< IdentProcBox > & range  = (*i).second.second ;

    const unsigned depth   = key.depth();
    const unsigned ordinal = oct_tree_offset( tree_depth , key );
    const unsigned num_d   = domain.size();
    const unsigned num_r   = range.size();
    const unsigned number  = num_d + num_r ;

    if ( depth <= 4 ) { // Values for this node:
      local_count[ 2 * ordinal ] += number ;
    }
    else { // Values for a deeper node
      local_count[ 2 * ordinal + 1 ] += number ;
    }
  }

  all_reduce_sum( arg_comm , local_count , global_count , tree_size_2 );

  oct_tree_partition_private( 0, p_size, tree_depth,
                              arg_tolerance, global_count,
                              p_size, & arg_cuts[0]);
}

//----------------------------------------------------------------------

bool oct_tree_proximity_search(
  ParallelMachine            arg_comm ,
  const float        * const arg_global_box ,
  const unsigned             arg_domain_boxes_number ,
  const IdentProcBox * const arg_domain_boxes ,
  const unsigned             arg_range_boxes_number ,
  const IdentProcBox * const arg_range_boxes ,
  const OctTreeKey   * const arg_cut_keys ,
  std::vector< std::pair<IdentProc,IdentProc> > & arg_relation ,
  unsigned * const arg_search_tree_stats )
{
  enum { Dim = 3 };

  const bool symmetric = arg_range_boxes == arg_domain_boxes ||
                         arg_range_boxes == NULL ;
  const unsigned p_size = parallel_machine_size( arg_comm );
  const unsigned p_rank = parallel_machine_rank( arg_comm );

  //----------------------------------------------------------------------
  // Search tree defined by oct-tree covering for boxes

  bool local_violations = false ;
  bool global_violations = false ;

  SearchTree search_tree ;

  {
    OctTreeKey covering[8] ;
    unsigned   number = 0 ;

    for ( unsigned i = 0 ; i < arg_domain_boxes_number ; ++i ) {

      IdentProcBox tmp( arg_domain_boxes[i] );

      tmp.proc = p_rank ;

      const bool valid =
        hsfc_box_covering( arg_global_box, tmp.box, covering, number );

      if ( ! valid ) { local_violations = true ; }

      for ( unsigned k = 0 ; k < number ; ++k ) {
        const OctTreeKey key = covering[k] ;
        search_tree[key].first.push_back(tmp);
      }
    }

    if ( ! symmetric ) {
      for ( unsigned i = 0 ; i < arg_range_boxes_number ; ++i ) {

        IdentProcBox tmp( arg_range_boxes[i] );

        tmp.proc = p_rank ;

        const bool valid =
          hsfc_box_covering( arg_global_box, tmp.box, covering, number );

        if ( ! valid ) { local_violations = true ; }

        for ( unsigned k = 0 ; k < number ; ++k ) {
          const OctTreeKey key = covering[k] ;
          search_tree[key].second.push_back(tmp);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // Use a set to provide a unique and sorted result.

  std::set< std::pair<IdentProc,IdentProc> > tmp_relation ;

  if ( p_size == 1 ) {

    global_violations = local_violations ;

    if ( arg_search_tree_stats ) {
      search_tree_statistics( arg_comm , search_tree ,
                              arg_search_tree_stats );
    }

    ProximitySearch( symmetric, search_tree, tmp_relation);
  }
  else {
    // Communicate search_tree members

    SearchTree local_tree ;

    std::set< std::pair<IdentProc,IdentProc> > local_relation ;

    if ( arg_cut_keys ) {
      global_violations =
        communicate( arg_comm , arg_cut_keys , search_tree , local_tree ,
                     local_violations );
    }
    else {
      const double tolerance = 0.001 ;

      std::vector< OctTreeKey > cuts ;

      oct_tree_partition( arg_comm , search_tree , tolerance , cuts );
     
      global_violations =
        communicate( arg_comm , & cuts[0] , search_tree , local_tree ,
                     local_violations );
    }

    // Local proximity search with received members

    if ( arg_search_tree_stats ) {
      search_tree_statistics( arg_comm , local_tree ,
                              arg_search_tree_stats );
    }

    ProximitySearch( symmetric, local_tree, local_relation);

    // Communicate relations back to domain and range processors

    communicate( arg_comm , local_relation , tmp_relation );
  }

  arg_relation.clear();
  arg_relation.reserve( tmp_relation.size() );

  std::set< std::pair<IdentProc,IdentProc> >::iterator ir ;
  for ( ir = tmp_relation.begin() ; ir != tmp_relation.end() ; ++ir ) {
    arg_relation.push_back( *ir );
  }

  return global_violations ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Course and fast algorithm.
// Generate a dense octree at 'tree_depth',
// sum contributions from all processors,
// and then split the tree among processors.

void oct_tree_partition_course(
  ParallelMachine          comm ,
  const unsigned           length ,
  const OctTreeKey * const keys ,
  const float      * const weights ,
        OctTreeKey * const cut_keys ,
  const unsigned           override_size )
{
  static const char method[] = "phdmesh::oct_tree_partition" ;

  enum { tree_depth  = 4 };
  enum { tree_size   = OctTreeSize< tree_depth >::value };
  enum { tree_size_2 = tree_size * 2 };

  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  float local_weight[ tree_size_2 ];
  float global_weight[ tree_size_2 ];

  for ( unsigned i = 0 ; i < tree_size_2 ; ++i ) {
    local_weight[i] = 0 ;
    global_weight[i] = 0 ;
  }

  for ( unsigned i = 0 ; i < length ; ++i ) {

    const OctTreeKey & k = keys[i] ;
    const float        w = weights[i] ;

    if ( 0 < w ) {
      const unsigned depth   = k.depth();
      const unsigned ordinal = oct_tree_offset( tree_depth , k );

      if ( depth <= tree_depth ) { // Values for this node:
        local_weight[ 2 * ordinal ] += w ;
      }
      else { // Values for a deeper node
        local_weight[ 2 * ordinal + 1 ] += w ;
      }
    }
    else {
      // ERROR: Non-positive weight
      // The root node's descending weight is unused at this moment
      local_weight[1] += 1 ;
    }
  }

  all_reduce_sum( comm , local_weight , global_weight , tree_size_2 );

  // Error handling

  if ( 0 < global_weight[1] ) {
    std::ostringstream msg ;
    msg << method << " : FAILED Due to non-positive weight, P" << p_rank ;
    msg << " had " << local_weight[1] ;
    throw std::range_error( msg.str() );
  }

  const double tolerance = 0.001 ;

  const unsigned ncuts = override_size ? override_size : p_size ;

  oct_tree_partition_private( 0, ncuts, tree_depth,
                              tolerance, global_weight, ncuts, cut_keys );
}

//----------------------------------------------------------------------
// Fine and slow algorithm.
// Generate a dense octree at 'tree_depth',
// sum number of contributions from all processors,
// split the tree among processors,
// communicate the contributions to the owning processors,
// split the contributions among processors via parallel prefix.

void oct_tree_partition_fine(
  ParallelMachine          comm ,
  const unsigned           length ,
  const OctTreeKey * const keys ,
  const float      * const weights ,
        OctTreeKey * const cut_keys )
{
  static const char method[] = "phdmesh::oct_tree_partition" ;

  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  {
    enum { tree_depth  = 4 };
    enum { tree_size   = OctTreeSize< tree_depth >::value };
    enum { tree_size_2 = tree_size * 2 };

    float local_count[ tree_size_2 ];
    float global_count[ tree_size_2 ];

    for ( unsigned i = 0 ; i < tree_size_2 ; ++i ) {
      local_count[i] = 0 ;
      global_count[i] = 0 ;
    }

    for ( unsigned i = 0 ; i < length ; ++i ) {

      const OctTreeKey & k = keys[i] ;
      const float        w = (float)( weights ? weights[i] : 1.0 );

      if ( 0 < w ) {
        const unsigned depth   = k.depth();
        const unsigned ordinal = oct_tree_offset( tree_depth , k );

        if ( depth <= tree_depth ) { // Values for this node:
          local_count[ 2 * ordinal ] += 1 ;
        }
        else { // Values for a deeper node
          local_count[ 2 * ordinal + 1 ] += 1 ;
        }
      }
      else {
        // ERROR: Non-positive weight
        // The root node's descending weight is unused at this moment
        local_count[1] += 1 ;
      }
    }

    all_reduce_sum( comm , local_count , global_count , tree_size_2 );

    // Error handling

    if ( 0 < global_count[1] ) {
      std::ostringstream msg ;
      msg << method << " : FAILED Due to non-positive weight, P" << p_rank ;
      msg << " had " << local_count[1] ;
      throw std::range_error( msg.str() );
    }

    const double tolerance = 0.001 ;

    oct_tree_partition_private( 0, p_size, tree_depth,
                                tolerance, global_count, p_size, cut_keys );

    // close scope to reclaim memory from the local and global trees
  }

  //----------------------------------------------------------------------
  // Cut keys are now the working partitioning,
  // communicate keys and weights to the owning entity.
  // Prefix the set and chop it up among processors accordingly.

  std::map< OctTreeKey , float > entities ;

  {
    CommAll comm_all( comm );

    for ( unsigned i = 0 ; i < length ; ++i ) {
      const OctTreeKey & k = keys[i] ;
      const float        w = (float)( weights ? weights[i] : 1.0 );
      const unsigned p = processor( cut_keys , cut_keys + p_size , k );
      CommBuffer & buf = comm_all.send_buffer(p);
      buf.pack<OctTreeKey::value_type>( k.value() , OctTreeKey::NWord );
      buf.pack<float>( w );
    }

    comm_all.allocate_buffers( p_size / 4 , false );

    for ( unsigned i = 0 ; i < length ; ++i ) {
      const OctTreeKey & k = keys[i] ;
      const float        w = (float)( weights ? weights[i] : 1.0 );
      const unsigned p = processor( cut_keys , cut_keys + p_size , k );
      CommBuffer & buf = comm_all.send_buffer(p);
      buf.pack<OctTreeKey::value_type>( k.value() , OctTreeKey::NWord );
      buf.pack<float>( w );
    }

    comm_all.communicate();

    OctTreeKey key ;
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = comm_all.recv_buffer( p );
      while ( buf.remaining() ) {
        OctTreeKey::value_type val[ OctTreeKey::NWord ];
        float w ;
        buf.unpack<OctTreeKey::value_type>( val , OctTreeKey::NWord );
        buf.unpack<float>( w );
        key.set_value( val );
        entities[ key ] = w ;
      }
    }
  }

  // Iterate and sum the local weights,
  // then scatter the results.

  const float fzero = 0.0 ;

  std::vector< float > global_weights( p_size , fzero );

  // Now prefix the global weights to determine the total global weight
  // and which span of cuts fall into this processor's responsibility.

  {
    float local_weight = 0.0 ;

    const std::map< OctTreeKey , float >::const_iterator e = entities.end();
          std::map< OctTreeKey , float >::const_iterator i = entities.begin();

    for ( ; i != e ; ++i ) {
      local_weight += (*i).second ;
    }

    all_gather( comm , & local_weight , & global_weights[0] , 1 );
  }

  for ( unsigned p = 1 ; p < p_size ; ++p ) {
    global_weights[p] += global_weights[p-1] ;
  }

  {
    // This processor will take care of some range of cuts,
    // most likely just its own.

    const double total_weight = global_weights[p_size - 1] ;
    const double cut_weight   = total_weight / ((double) p_size);
    const double first_weight = p_rank ? global_weights[ p_rank - 1 ] : 0.0 ;

    const unsigned p_first =
      1 + (unsigned)( first_weight / cut_weight );

    const unsigned p_end   =
      ( 1 + p_rank == p_size ) ? p_size :
      1 + (unsigned)( global_weights[p_rank] / cut_weight );

    const unsigned uzero = 0 ;
    std::vector<unsigned> local_tmp(  p_size * 2 , uzero );
    std::vector<unsigned> global_tmp( p_size * 2 , uzero );

    const std::map< OctTreeKey , float >::const_iterator e = entities.end();
    const std::map< OctTreeKey , float >::const_iterator b = entities.begin();
          std::map< OctTreeKey , float >::const_iterator i = b ;

    float w = 0.0 ;

    if ( p_rank == 0 ) { cut_keys[0] = OctTreeKey(); }

    for ( unsigned p = p_first ; p < p_end && e != i ; ++p ) {
      // Where to cut for 'p'
      const double wcut = p * cut_weight - first_weight ;

      // Iterate up to 'wcut'
      while ( e != i && w + (*i).second < wcut ) {
        w += (*i).second ;
        ++i ;
      }
      // Which is closer to 'wcut'?  To add or not
      if ( e != i && w + (*i).second - wcut < wcut - w ) {
        ++i ;
      }

      // Cut between 'i' and '--i'

      OctTreeKey tmp_key ;

      if ( b == i ) {
        tmp_key = cut_keys[ p_rank ];
      }
      else if ( e == i ) {
        tmp_key = cut_keys[ p_rank + 1 ];
      }
      else {
        std::map< OctTreeKey , float >::const_iterator j = i ; --j ;
        const OctTreeKey ki = (*i).first ;
        const OctTreeKey kj = (*j).first ;
        unsigned d = 1 ;
        while ( d < OctTreeKey::MaxDepth && ki.index(d) == kj.index(d) ) {
          tmp_key.set_index( d , ki.index(d) );
          ++d ;
        }
        tmp_key.set_index( d , ( ki.index(d) + kj.index(d) + 1 ) / 2 );
      }
      { // Work-around for shoddy pathscale compiler:
              unsigned * const dst = & local_tmp[ p * 2 ];
        const unsigned * const src = tmp_key.value();
        Copy<2>( dst , src );
        // Copy<2>( & local_tmp[ p * 2 ] , tmp_key.value() );
      }
    }

    all_reduce_bor( comm , & local_tmp[0] , & global_tmp[0] , 2 * p_size );

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      cut_keys[p].set_value( & global_tmp[2*p] );
    }
  }
}



}


