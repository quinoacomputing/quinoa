/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2008) Sandia Corporation                     */
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
 * @date   June 2008
 */

#ifndef phdmesh_element_Declarations_hpp
#define phdmesh_element_Declarations_hpp

#include <mesh/Types.hpp>
#include <element/CellTopology.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** Attach an cell topology to a Part.
 *  There is at most one cell topology allowed.
 */
void set_cell_topology( Part & , const CellTopology * singleton );

/** Attach an element local topology to a Part.
 *  There is at most one element topology allowed.
 */
template< class Traits >
void set_cell_topology( Part & p )
{ return set_cell_topology( p , cell_topology<Traits>() ); }

const CellTopology * get_cell_topology( const Part & );
const CellTopology * get_cell_topology( const Kernel & );
const CellTopology * get_cell_topology( const Entity & );

//----------------------------------------------------------------------
/** Declare an element member of a part with a cell topology
 *  and nodes conformal to that topology.
 */
Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const unsigned elem_id ,
                          const unsigned node_id[] );

Entity & declare_element( BulkData & mesh ,
                          Part & part ,
                          const unsigned elem_id ,
                          Entity * node[] );

//----------------------------------------------------------------------
/* The element must have a topology. */

Entity & declare_element_side( BulkData & mesh ,
                               const unsigned global_side_id ,
                               Entity & elem , const unsigned local_side_id );

//----------------------------------------------------------------------

}

#endif

