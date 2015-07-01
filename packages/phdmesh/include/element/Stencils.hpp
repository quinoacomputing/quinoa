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
 * @date   May 2008
 */

#ifndef phdmesh_element_Stencils_hpp
#define phdmesh_element_Stencils_hpp

#include <util/Basics.hpp>
#include <mesh/Types.hpp>

namespace phdmesh {
namespace { // To prevent multiple copies for the linker

enum { ElementStencils_OK =
         StaticAssert< phdmesh::Node == 0 &&
                       phdmesh::Edge == 1 &&
                       phdmesh::Face == 2 &&
                       phdmesh::Element == 3 >::OK };

//----------------------------------------------------------------------

template< class TopologyTraits >
int element_node_stencil( EntityType , EntityType , unsigned , unsigned );

template<>
int element_node_stencil<void>( EntityType from_type ,
                                EntityType to_type ,
                                unsigned   identifier ,
                                unsigned   kind )
{
  int ordinal = -1 ;

  if ( Element == from_type && Node == to_type && 0 == kind ) {
    ordinal = (int) identifier ;
  }

  return ordinal ;
}

template< class TopologyTraits >
int element_node_stencil( EntityType from_type ,
                          EntityType to_type ,
                          unsigned   identifier ,
                          unsigned   kind )
{
  enum { number_node = TopologyTraits::node_count };

  int ordinal = -1 ;

  if ( Element == from_type &&
       Node    == to_type &&
       0       == kind &&
       identifier < number_node ) {
    ordinal = (int) identifier ;
  }

  return ordinal ;
}

//----------------------------------------------------------------------

}
}

#endif

