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
 * @date   January 2007
 */

#ifndef pddgeom_OctTreeOps_hpp
#define pddgeom_OctTreeOps_hpp

#include <utility>
#include <vector>
#include <util/Parallel.hpp>
#include <util/OctTree.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** A parallel-unique identified cartesion box.
 *    lower bound is box[0..2]
 *    upper bound is box[3..5]
 *
 *  For the oct-tree proximity search two boxes that have the
 *  same nonzero 'part' value are not checked for intesection.
 */
struct IdentProcBox : public IdentProc {
  unsigned part ;
  float box[6] ;

  IdentProcBox() {}

  IdentProcBox( const IdentProcBox & rhs )
    : IdentProc( rhs ), part( rhs.part ) { Copy<6>( box , rhs.box ); }

  IdentProcBox & operator = ( const IdentProcBox & rhs )
    {
      IdentProc::operator =( rhs );
      part = rhs.part ;
      Copy<6>( box , rhs.box );
      return *this ;
    }
};

//----------------------------------------------------------------------
/** Global bounds for a set of boxes.
 *  The lower bound is the minimum of all boxes, decreased by epsilon.
 *  The upper bound is the maximum of all boxes, increased by epsilon.
 *  Thus all input boxes are fully contained within the global box.
 */
void box_global_bounds(
  ParallelMachine            arg_comm ,
  const unsigned             arg_domain_boxes_number ,
  const IdentProcBox * const arg_domain_boxes ,
  const unsigned             arg_range_boxes_number ,
  const IdentProcBox * const arg_range_boxes ,
        float        * const arg_global_box );

//----------------------------------------------------------------------
/** Search for intersection of domain boxes with range boxes
 *  within a given global bounding box.
 *  Output vector of matches with a domain or range box on
 *  the local processor.
 *
 *  If 'arg_cuts' is given it will be used for the parallel
 *  search.  If 'arg_cuts == NULL' then a balanced internal
 *  partitioning will be generated.
 *
 *  The search_tree_stats are for the local search:
 *    [0] = minimum search tree domain cell size
 *    [1] = maximum search tree domain cell size
 *    [2] = average search tree domain cell size
 *    [3] = minimum search tree range cell size
 *    [4] = maximum search tree range cell size
 *    [5] = average search tree range cell size
 *  These statistics require an extra communication to gather.
 *
 *  Returns 'true' if all small boxes on all processors 
 *  had non-negative volumes and were fully contained within
 *  the global box.
 */
bool oct_tree_proximity_search(
  ParallelMachine            arg_comm ,
  const float        * const arg_global_box ,
  const unsigned             arg_domain_boxes_number ,
  const IdentProcBox * const arg_domain_boxes ,
  const unsigned             arg_range_boxes_number ,
  const IdentProcBox * const arg_range_boxes ,
  const OctTreeKey   * const arg_cuts ,
  std::vector< std::pair<IdentProc,IdentProc> > & arg_relation ,
  unsigned * const arg_search_tree_stats = NULL );

//----------------------------------------------------------------------
/** Generate an oct-tree covering of a small box within a global box.
 *  The cartesian space is mapped to an oct-tree via Hilbert space
 *  filling curve.  The covering consists of 1..8 oct-tree cells,
 *  the 'covering' array must be dimensioned to at least eight.
 *  Returns true for a "good" small box: it a non-negative volume
 *  and is fully contained within the global box.
 */
bool hsfc_box_covering( const float * const global_box ,
                        const float * const small_box ,
                        OctTreeKey  * const covering ,
                        unsigned &          number );

//----------------------------------------------------------------------
/** Given an array of oct-tree keys and weights generate
 *  a course partitioning of the oct-tree key space.
 *  The algorithm is fast with a single 'all_reduce' operation.
 */
void oct_tree_partition_course(
  ParallelMachine          comm ,
  const unsigned           length ,
  const OctTreeKey * const keys ,
  const float      * const weights ,
        OctTreeKey * const cut_keys ,
  const unsigned           override_ncuts = 0 );

//----------------------------------------------------------------------
/** Given an array of oct-tree keys and weights generate
 *  a fine partitioning of the oct-tree key space.
 *  The algorithm is likely to produce a better balance than
 *  the previous course algorithm.  However, it performs
 *  more communcations and computations.
 */
void oct_tree_partition_fine(
  ParallelMachine          comm ,
  const unsigned           length ,
  const OctTreeKey * const keys ,
  const float      * const weights ,
        OctTreeKey * const cut_keys );

//----------------------------------------------------------------------
/**  A recursive kernel used within the oct_tree_partitioning algorithms.
 *   Exposed to support unit testing.
 */

void oct_tree_partition_private(
  const unsigned p_first ,
  const unsigned p_end ,
  const unsigned depth ,
  const double   tolerance ,
  float * const weights ,
  const unsigned cuts_length ,
  OctTreeKey * const cuts );

//----------------------------------------------------------------------

}

#endif

