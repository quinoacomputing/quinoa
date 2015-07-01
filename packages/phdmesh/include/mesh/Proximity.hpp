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

#ifndef phdmesh_Proximity_hpp
#define phdmesh_ProximityTypes_hpp

//----------------------------------------------------------------------

#include <mesh/FieldTraits.hpp>
#include <mesh/Types.hpp>

namespace phdmesh {

class ProximitySearch ;

/** Global geometric proximity search for same-mesh and same-type entities.
 *  Only the owned entities are included in the search.
 */
void proximity_search( BulkData & M ,
                       const ProximitySearch & ,
                       const unsigned entity_type ,
                       std::vector< std::pair<IdentProc,IdentProc> > & );

class ProximitySearch {
public:
  typedef Field<double,Cartesian> CoordinateField ;

  /** Generate a cartesion box enclosing the entity.
   *  box[0-2] = min , box[3-5] = max
   */
  virtual void box( const Entity & , float * const ) const ;

  /** For entities in the partset define a search part id.
   *  Entities will not be compared for proximity if their
   *  part ids are not zero and not equal.  This allows 
   *  relationed entities to be filtered out of the search.
   */
  virtual int part_id( const Kernel & ) const ;

  virtual ~ProximitySearch() ;

  /** Default behavior is to define a bounding box based
   *  upon relationed nodes.  If a part has a CSet member
   *  of the ProximitySearch object then the 'part_id'
   *  is that part's ordinal + 1.
   */
  ProximitySearch( const CoordinateField & node_coord ,
                   const float             box_expansion = 0 );

  const CoordinateField & m_node_coord ;
  const float             m_box_expansion ;

private:
  ProximitySearch();
  ProximitySearch( const ProximitySearch & );
  ProximitySearch & operator = ( const ProximitySearch & );
};

} // namespace phdmesh

#endif

