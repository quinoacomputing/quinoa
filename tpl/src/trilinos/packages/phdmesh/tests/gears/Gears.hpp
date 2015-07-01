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

#ifndef gears_Gears_hpp
#define gears_Gears_hpp

#include <vector>

#include <util/Parallel.hpp>
#include <mesh/Types.hpp>
#include <mesh/FieldTraits.hpp>

#include <mesh_io/ExoII.hpp>

#include <element/Dimensions.hpp>

namespace phdmesh {

struct GearFields {

  enum { SpatialDimension = 3 };

  typedef Field<double,Cartesian>            CartesianField ;
  typedef Field<double,Cylindrical>          CylindricalField ;
  typedef exodus::FileSchema::AttributeField AttributeField ;

  CylindricalField & gear_coord ;
  CartesianField   & model_coord ;
  CartesianField   & current_coord ;
  CartesianField   & displacement ;
  AttributeField   & element_attr ;
  CartesianField   & test_value ;
  ElementNodePointerField & elem_node_test_value ;

  GearFields( MetaData & S );

private:
  GearFields();
  GearFields( const GearFields & );
  GearFields & operator = ( const GearFields & );
};

class Gear {
public:
  Gear( MetaData & S ,
        const std::string & name ,
        const GearFields & gear_fields ,
        const double   center[] ,
        const double   rad_min ,
        const double   rad_max ,
        const unsigned rad_num ,
        const double   z_min ,
        const double   z_max ,
        const unsigned z_num ,
        const unsigned angle_num ,
        const int      turn_direction );

  void mesh( BulkData & );
  void turn( double turn_angle ) const ;

  MetaData & m_mesh_meta_data ;
  BulkData * m_mesh ;
  Part & m_gear ;
  Part & m_surf ;
  const GearFields::CylindricalField  & m_gear_coord ;
  const GearFields::CartesianField    & m_model_coord ;
  const GearFields::CartesianField    & m_current_coord ;
  const GearFields::CartesianField    & m_displacement ;
  const GearFields::CartesianField    & m_test_value ;
  const ElementNodePointerField       & m_elem_node_test_value ;

private:

  double m_center[3] ;
  double m_z_min ;
  double m_z_max ;
  double m_z_inc ;
  double m_rad_min ;
  double m_rad_max ;
  double m_rad_inc ;
  double m_ang_inc ;
  unsigned m_rad_num ;
  unsigned m_z_num ;
  unsigned m_angle_num ;
  int      m_turn_dir ;

  Entity * create_node( const std::vector<Part*> & ,
                        const unsigned long node_id_base ,
                        const unsigned iz ,
                        const unsigned ir ,
                        const unsigned ia ) const ;
};

}

#endif

