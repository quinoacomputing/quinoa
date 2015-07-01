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
 * @author H. Carter Edwards
 */

#include <math.h>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <util/TPI.h>
#include <util/ParallelComm.hpp>

#include <mesh/MetaData.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/Comm.hpp>

#include <element/Stencils.hpp>
#include <element/Declarations.hpp>
#include <element/Hexahedron_Topologies.hpp>

#include "Gears.hpp"

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {

//----------------------------------------------------------------------

GearFields::GearFields( MetaData & S )
: gear_coord(
    S.declare_field<CylindricalField>( std::string("gear_coordinates") ) ),
  model_coord(
    S.declare_field<CartesianField>( std::string("model_coordinates") ) ),
  current_coord(
    S.declare_field<CartesianField>( std::string("coordinates") , 2 ) ),
  displacement(
    S.declare_field<CartesianField>( std::string("displacement") , 2 ) ),
  element_attr(
    S.declare_field<AttributeField>( std::string("element_attribute") ) ),
  test_value(
    S.declare_field<CartesianField>( std::string("test_value") ) ),
  elem_node_test_value(
    S.declare_field<ElementNodePointerField>( std::string("test_pointer") ) )
{
  const Part & universe = S.universal_part();

  S.put_field( gear_coord    , Node , universe , SpatialDimension );
  S.put_field( model_coord   , Node , universe , SpatialDimension );
  S.put_field( current_coord , Node , universe , SpatialDimension );
  S.put_field( displacement  , Node , universe , SpatialDimension );
  S.put_field( test_value  ,   Node , universe , SpatialDimension );
}

//----------------------------------------------------------------------

namespace {

unsigned long
identifier( unsigned nthick ,  // Number of entities through the thickness
            unsigned nradius , // Number of entities through the radius
            unsigned iz ,      // Thickness index
            unsigned ir ,      // Radial index
            unsigned ia )      // Angle index
{
  return iz + nthick * ( ir + nradius * ia );
}

}

Gear::Gear( MetaData & S ,
            const std::string & name ,
            const GearFields & gear_fields ,
            const double center[] ,
            const double rad_min ,
            const double rad_max ,
            const unsigned rad_num ,
            const double z_min ,
            const double z_max ,
            const unsigned z_num ,
            const unsigned angle_num ,
            const int      turn_dir )
  : m_mesh_meta_data( S ),
    m_mesh( NULL ),
    m_gear( S.declare_part(std::string("Gear_").append(name)) ),
    m_surf( S.declare_part(std::string("Surf_").append(name)) ),
    m_gear_coord(    gear_fields.gear_coord ),
    m_model_coord(   gear_fields.model_coord ),
    m_current_coord( gear_fields.current_coord ),
    m_displacement(  gear_fields.displacement ),
    m_test_value(    gear_fields.test_value ),
    m_elem_node_test_value( gear_fields.elem_node_test_value )
{
  typedef Hexahedron<> Hex ;
  enum { SpatialDimension = GearFields::SpatialDimension };

  set_cell_topology< Hex >( m_gear );

  // A test value for each element:

  S.put_field( gear_fields.test_value , Element , m_gear , SpatialDimension );

  // An element->node->test_value pointer field:

  S.put_field( gear_fields.elem_node_test_value, Element, m_gear,
               Hex::node_count );

  m_mesh_meta_data.declare_field_relation(
    gear_fields.elem_node_test_value ,
    & element_node_stencil< Hex > ,
    gear_fields.test_value );

  // Meshing parameters for this gear:

  const double TWO_PI = 2.0 * acos( (double) -1.0 );

  m_center[0] = center[0] ;
  m_center[1] = center[1] ;
  m_center[2] = center[2] ;

  m_z_min     = z_min ;
  m_z_max     = z_max ;
  m_z_inc     = (z_max - z_min) / ((double) ( z_num - 1 ));

  m_rad_min   = rad_min ;
  m_rad_max   = rad_max ;
  m_rad_inc   = (rad_max - rad_min) / ((double) ( rad_num - 1 ));

  m_ang_inc   = TWO_PI / (double) angle_num ;

  m_rad_num   = rad_num ;
  m_z_num     = z_num ;
  m_angle_num = angle_num ;
  m_turn_dir  = turn_dir ;
}

//----------------------------------------------------------------------

Entity * Gear::create_node(
  const std::vector<Part*> & parts ,
  const unsigned long node_id_base ,
  const unsigned iz ,
  const unsigned ir ,
  const unsigned ia ) const
{
  const double angle     = m_ang_inc * ia ;
  const double cos_angle = cos( angle );
  const double sin_angle = sin( angle );

  const double radius = m_rad_min + m_rad_inc * ir ;
  const double x = m_center[0] + radius * cos_angle ;
  const double y = m_center[1] + radius * sin_angle ;
  const double z = m_center[2] + m_z_min + m_z_inc * iz ;

  // Create the node and set the model_coordinates

  unsigned long id_gear = identifier( m_z_num, m_rad_num, iz, ir, ia );
  unsigned long id = node_id_base + id_gear ;

  Entity & node = m_mesh->declare_entity( entity_key(Node,id) , parts );

  double * const gear_data    = field_data( m_gear_coord , node );
  double * const model_data   = field_data( m_model_coord , node );
  double * const current_data = field_data( m_current_coord , node );
  double * const disp_data    = field_data( m_displacement , node );

  gear_data[0] = radius ;
  gear_data[1] = angle ;
  gear_data[2] = z - m_center[2] ;

  model_data[0] = x ;
  model_data[1] = y ;
  model_data[2] = z ;

  current_data[0] = x ;
  current_data[1] = y ;
  current_data[2] = z ;

  disp_data[0] = 0.0 ;
  disp_data[1] = 0.0 ;
  disp_data[2] = 0.0 ;

  return & node ;
}

//----------------------------------------------------------------------

void Gear::mesh( BulkData & M )
{
  static const char method[] = "phdmesh::Gear::mesh" ;

  m_mesh = & M ;

  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();

  entity_id_type counts[ EntityTypeEnd ];
  entity_id_type max_id[ EntityTypeEnd ];

  comm_mesh_stats( M , counts , max_id );

  const unsigned long node_id_base = max_id[ Node ] + 1 ;
  const unsigned long face_id_base = max_id[ Face ] + 1 ;
  const unsigned long elem_id_base = max_id[ Element ] + 1 ;

  const unsigned long elem_id_gear_max =
    m_angle_num * ( m_rad_num - 1 ) * ( m_z_num - 1 );

  std::vector<Part*> elem_parts ;
  std::vector<Part*> face_parts ;

  {
    Part * const p_gear = & m_gear ;
    Part * const p_surf = & m_surf ;
    Part * const p_owns = & m_mesh_meta_data.locally_owned_part();

    elem_parts.push_back( p_gear );
    elem_parts.push_back( p_owns );

    face_parts.push_back( p_gear );
    face_parts.push_back( p_surf );
    face_parts.push_back( p_owns );
  }

  for ( unsigned ia = 0 ; ia < m_angle_num ; ++ia ) {
    for ( unsigned ir = 0 ; ir < m_rad_num - 1 ; ++ir ) {
      for ( unsigned iz = 0 ; iz < m_z_num - 1 ; ++iz ) {

        unsigned long elem_id_gear =
          identifier( m_z_num-1 , m_rad_num-1 , iz , ir , ia );

        if ( ( ( elem_id_gear * p_size ) / elem_id_gear_max ) == p_rank ) {

          unsigned long elem_id = elem_id_base + elem_id_gear ;

          // Create the node and set the model_coordinates

          const unsigned ia_1 = ( ia + 1 ) % m_angle_num ;
          const unsigned ir_1 = ir + 1 ;
          const unsigned iz_1 = iz + 1 ;

          Entity * node[8] ;

          node[0] = create_node( elem_parts, node_id_base, iz  , ir  , ia_1 );
          node[1] = create_node( elem_parts, node_id_base, iz  , ir  , ia   );
          node[2] = create_node( elem_parts, node_id_base, iz_1, ir  , ia   );
          node[3] = create_node( elem_parts, node_id_base, iz_1, ir  , ia_1 );
          node[4] = create_node( elem_parts, node_id_base, iz  , ir_1, ia_1 );
          node[5] = create_node( elem_parts, node_id_base, iz  , ir_1, ia   );
          node[6] = create_node( elem_parts, node_id_base, iz_1, ir_1, ia   );
          node[7] = create_node( elem_parts, node_id_base, iz_1, ir_1, ia_1 );

#if 0 /* VERIFY_CENTROID */

          // Centroid of the element for verification

          const double TWO_PI = 2.0 * acos( (double) -1.0 );
          const double angle = m_ang_inc * ( 0.5 + (double) ia );
          const double z = m_center[2] + m_z_min + m_z_inc * (0.5 + (double)iz);

          double c[3] = { 0 , 0 , 0 };

          for ( unsigned j = 0 ; j < 8 ; ++j ) {
            double * const coord_data = field_data( m_model_coord , *node[j] );
            c[0] += coord_data[0] ;
            c[1] += coord_data[1] ;
            c[2] += coord_data[2] ;
          }
          c[0] /= 8 ; c[1] /= 8 ; c[2] /= 8 ;
          c[0] -= m_center[0] ;
          c[1] -= m_center[1] ;

          double val_a = atan2( c[1] , c[0] );
          if ( val_a < 0 ) { val_a += TWO_PI ; }
          const double err_a = angle - val_a ;
          const double err_z = z - c[2] ;

          const double eps = 100 * std::numeric_limits<double>::epsilon();

          if ( err_z < - eps || eps < err_z ||
               err_a < - eps || eps < err_a ) {
            std::string msg ;
            msg.append("problem setup element centroid error" );
            throw std::logic_error( msg );
          }
#endif

          Entity & elem =
            M.declare_entity( entity_key(Element,elem_id), elem_parts );

          for ( unsigned j = 0 ; j < 8 ; ++j ) {
            M.declare_relation( elem , * node[j] , j );
          }
        }
      }
    }
  }

  // Array of faces on the surface

  {
    const unsigned ir = m_rad_num - 1 ;

    for ( unsigned ia = 0 ; ia < m_angle_num ; ++ia ) {
      for ( unsigned iz = 0 ; iz < m_z_num - 1 ; ++iz ) {

        unsigned long elem_id_gear =
          identifier( m_z_num-1 , m_rad_num-1 , iz , ir-1 , ia );

        if ( ( ( elem_id_gear * p_size ) / elem_id_gear_max ) == p_rank ) {

          unsigned long elem_id = elem_id_base + elem_id_gear ;

          unsigned long face_id =
            face_id_base + identifier( m_z_num-1 , 1 , iz , 0 , ia );

          unsigned face_ord = 5 ;
          Entity * node[4] ;

          const unsigned ia_1 = ( ia + 1 ) % m_angle_num ;
          const unsigned iz_1 = iz + 1 ;

          node[0] = create_node( face_parts, node_id_base, iz  , ir  , ia_1 );
          node[1] = create_node( face_parts, node_id_base, iz  , ir  , ia   );
          node[2] = create_node( face_parts, node_id_base, iz_1, ir  , ia   );
          node[3] = create_node( face_parts, node_id_base, iz_1, ir  , ia_1 );

          Entity & face =
            M.declare_entity( entity_key(Face,face_id), face_parts );

          for ( unsigned j = 0 ; j < 4 ; ++j ) {
            M.declare_relation( face , * node[j] , j );
          }

          Entity & elem = * M.get_entity(entity_key(Element,elem_id),method);

          M.declare_relation( elem , face , face_ord );
        }
      }
    }
  }
}

//----------------------------------------------------------------------
// Iterate nodes and turn them by the angle

void Gear::turn( double turn_angle ) const
{
  const unsigned Length = 3 ;

  const KernelSet & ks = m_mesh->kernels( Node );
  const KernelSet::iterator ek = ks.end();
        KernelSet::iterator ik = ks.begin();
  for ( ; ik != ek ; ++ik ) {
    Kernel & k = *ik ;
    if ( k.has_superset( m_gear ) ) {
      const unsigned n = k.size();
      double * const kernel_gear_data    = field_data( m_gear_coord , k );
      double * const kernel_model_data   = field_data( m_model_coord , k );
      double * const kernel_current_data = field_data( m_current_coord , k );
      double * const kernel_disp_data    = field_data( m_displacement , k );

      for ( unsigned i = 0 ; i < n ; ++i ) {
        double * const gear_data    = kernel_gear_data    + i * Length ;
        double * const model_data   = kernel_model_data   + i * Length ;
        double * const current_data = kernel_current_data + i * Length ;
        double * const disp_data    = kernel_disp_data    + i * Length ;
    
        const double radius = gear_data[0] ;
        const double angle  = gear_data[1] + turn_angle * m_turn_dir ;

        current_data[0] = m_center[0] + radius * cos( angle );
        current_data[1] = m_center[1] + radius * sin( angle );
        current_data[2] = m_center[2] + gear_data[2] ;

        disp_data[0] = current_data[0] - model_data[0] ;
        disp_data[1] = current_data[1] - model_data[1] ;
        disp_data[2] = current_data[2] - model_data[2] ;
      }
    }
  }
}

//----------------------------------------------------------------------

}

