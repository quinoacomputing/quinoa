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

#ifndef phdmesh_element_Dimension_hpp
#define phdmesh_element_Dimension_hpp

#include <util/Array.hpp>
#include <mesh/Field.hpp>
#include <mesh/MetaData.hpp>
#include <element/Stencils.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

class ElementNode : public ArrayDimTag {
public:
  const char * name() const ;
  static const ElementNode & tag();
private:
  ElementNode() {}
  ElementNode( const ElementNode & );
  ElementNode & operator = ( const ElementNode & );
};

typedef Field<double*,ElementNode> ElementNodePointerField ;

template< class NodeField >
inline
ElementNodePointerField &
declare_element_node_pointer_field(
  MetaData & md , const std::string & s ,
  NodeField & node_field )
{
  const unsigned num_states = node_field.number_of_states();

  ElementNodePointerField & f =
    md.template declare_field< ElementNodePointerField >( s, num_states );

  for ( unsigned i = 0 ; i < num_states ; ++i ) {
    FieldState state = (FieldState) i;
    md.declare_field_relation(
      f[ state ] , & element_node_stencil<void> , node_field[ state ] );
  }
  
  return f ;
}

//----------------------------------------------------------------------

struct QuadratureTag : public ArrayDimTag {
  const char * name() const ;
  static const QuadratureTag & tag();
private:
  QuadratureTag() {}
  QuadratureTag( const QuadratureTag & );
  QuadratureTag & operator = ( const QuadratureTag & );
};

//----------------------------------------------------------------------

struct BasisTag : public ArrayDimTag {
  const char * name() const ;
  static const BasisTag & tag();
private:
  BasisTag() {}
  BasisTag( const BasisTag & );
  BasisTag & operator = ( const BasisTag & );
};

//----------------------------------------------------------------------

}

#endif

