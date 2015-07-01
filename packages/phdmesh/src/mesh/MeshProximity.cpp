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

#include <cctype>
#include <stdexcept>
#include <iostream>
#include <sstream>

#include <util/OctTreeOps.hpp>

#include <mesh/Types.hpp>
#include <mesh/MetaData.hpp>
#include <mesh/BulkData.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/Proximity.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

#if 0

namespace {

struct ProximityBoxes {
  const ProximitySearch   & m_prox ;
  const BulkData              & m_mesh ;
  const MetaData            & m_mesh_meta_data ;
  KernelSet::const_iterator m_iter ;
  KernelSet::const_iterator m_iter_end ;
  std::vector<IdentProcBox> m_boxes ;

  ProximityBoxes( 
    BulkData & M ,
    const ProximitySearch & prox ,
    const unsigned entity_type ,
    std::vector< std::pair<IdentProc,IdentProc> > & proximity );

  void fill_boxes();
};

void ProximityBoxes::fill_boxes()
{
  Part & owns_part = m_mesh_meta_data.locally_owned_part();
  const unsigned p_rank = m_mesh_meta_data.parallel_rank();

  for(;;) {
    Kernel * kernel = NULL ;
    int part_id ;
    std::vector<IdentProcBox>::iterator idbox ;

    { // Get work, requires a lock
      taskpool::log get_work_lock( 0 );
      for ( ; m_iter != m_iter_end && NULL == kernel ; ++m_iter ) {
        if ( m_iter->has_superset( owns_part ) &&
             0 <= ( part_id = prox.part_id(*m_iter) ) ) {
          kernel = & *m_iter ;
          idbox  =  m_box_iter ;
          m_box_iter += kernel->size();
        }
      }
    }

    if ( NULL != kernel ) {
      const Kernel::iterator j_end = kernel->end();
            Kernel::iterator j     = kernel->begin();

      for ( ; j != j_end ; ++j , ++idbox ) {
        const Entity & entity = **j ; ++j ;

        idbox->ident = entity.identifier();
        idbox->proc  = p_rank ;
        idbox->part  = (unsigned) part_id ;

        m_prox.box( entity , idbox->box );
      }
    }
  }
}


ProximityBoxes::ProximityBoxes(
  BulkData & M ,
  const ProximitySearch & prox ,
  const unsigned entity_type ,
  std::vector< std::pair<IdentProc,IdentProc> > & proximity )
  : m_prox( prox ),
    m_mesh( M ),
    m_mesh_meta_data( M.mesh_meta_data() )
{
  const MetaData & S = M.mesh_meta_data();
  const unsigned p_rank = S.parallel_rank();

  const KernelSet & kernels = m_mesh.kernels( entity_type );

  // Iterate surfaces and generate bounding boxes

  {
    Part & owns_part = m_mesh_meta_data.locally_owned_part();

    const KernelSet::const_iterator i_end = kernels.end();
          KernelSet::const_iterator i     = kernels.begin();

    unsigned count = 0 ;

    for ( i = kernels.begin() ; i != i_end ; ) {
      const Kernel & kernel = *i ; ++i ;
      if ( kernel.has_superset( owns_part ) && 0 <= prox.part_id(kernel) ) {
        count += kernel.size();
      }
    }

    m_boxes.assign( count );
  }

  {
    m_iter_end = kernels.end();
    m_iter     = kernels.begin();

    taskpool::run( *this , ProximityBoxes::fill_boxes , 1 );
  }

  float global_box[6] ;

  // Generate a global square box that will contain all of the
  // local boxes.
  {
    box_global_bounds( S.parallel() ,
                       boxes.size() , & boxes[0] , 0 , NULL ,
                       global_box );

    double global_center[3] , global_size[3] ;

    global_center[0] = ( global_box[0] + global_box[3] ) * 0.5 ;
    global_center[1] = ( global_box[1] + global_box[4] ) * 0.5 ;
    global_center[2] = ( global_box[2] + global_box[5] ) * 0.5 ;

    global_size[0] = global_box[3] - global_box[0] ;
    global_size[1] = global_box[4] - global_box[1] ;
    global_size[2] = global_box[5] - global_box[2] ;

    double max_size = global_size[0] ;
    if ( max_size < global_size[1] ) { max_size = global_size[1] ; }
    if ( max_size < global_size[2] ) { max_size = global_size[2] ; }
    max_size *= 0.5 ;

    global_box[0] = (float) ( global_center[0] - max_size );
    global_box[1] = (float) ( global_center[1] - max_size );
    global_box[2] = (float) ( global_center[2] - max_size );
    global_box[3] = (float) ( global_center[0] + max_size );
    global_box[4] = (float) ( global_center[1] + max_size );
    global_box[5] = (float) ( global_center[2] + max_size );
  }

  unsigned search_stats[6] ;
  unsigned search_tasks ;

  oct_tree_proximity_search( S.parallel() ,
                             global_box ,
                             boxes.size() , & boxes[0] , 0 , NULL ,
                             NULL , proximity , search_stats ,
                             & search_tasks );

  return search_tasks ;
}

#endif

//----------------------------------------------------------------------

void proximity_search(
  BulkData & M ,
  const ProximitySearch & prox ,
  const unsigned entity_type ,
  std::vector< std::pair<IdentProc,IdentProc> > & proximity )
{
  const MetaData & S = M.mesh_meta_data();
  const unsigned p_rank = M.parallel_rank();

  // Iterate surfaces and generate bounding boxes

  std::vector<IdentProcBox> boxes ;

  {
    Part & owns_part = S.locally_owned_part();

    const KernelSet & kernels = M.kernels( entity_type );

    const KernelSet::const_iterator i_end = kernels.end();
          KernelSet::const_iterator i     = kernels.begin();

    unsigned count = 0 ;

    for ( i = kernels.begin() ; i != i_end ; ) {
      const Kernel & kernel = *i ; ++i ;
      if ( kernel.has_superset( owns_part ) && 0 <= prox.part_id(kernel) ) {
        count += kernel.size();
      }
    }

    {
      const IdentProcBox empty ;
      boxes.assign( count , empty );
    }

    count = 0 ;
    for ( i = kernels.begin() ; i != i_end ; ) {

      const Kernel & kernel = *i ; ++i ;
      const int part_id = prox.part_id( kernel );

      if ( kernel.has_superset( owns_part ) && 0 <= part_id ) {

        const Kernel::iterator j_end = kernel.end();
              Kernel::iterator j     = kernel.begin();

        while ( j != j_end ) {
          IdentProcBox & idbox  = boxes[ count ] ; ++count ;
          const Entity & entity = **j ; ++j ;

          idbox.ident = entity.identifier();
          idbox.proc  = p_rank ;
          idbox.part  = (unsigned) part_id ;

          prox.box( entity , idbox.box );
        }
      }
    }
  }

  float global_box[6] ;

  // Generate a global square box that will contain all of the
  // local boxes.  The box is padded outward by a proportionately
  // small value.
  {
    box_global_bounds( M.parallel() ,
                       boxes.size() , & boxes[0] , 0 , NULL ,
                       global_box );

    double global_size[3] ;

    global_size[0] = global_box[3] - global_box[0] ;
    global_size[1] = global_box[4] - global_box[1] ;
    global_size[2] = global_box[5] - global_box[2] ;

    double max_size = global_size[0] ;
    if ( max_size < global_size[1] ) { max_size = global_size[1] ; }
    if ( max_size < global_size[2] ) { max_size = global_size[2] ; }

    global_box[3] = (float) ( global_box[0] + max_size );
    global_box[4] = (float) ( global_box[1] + max_size );
    global_box[5] = (float) ( global_box[2] + max_size );
  }

  unsigned search_stats[6] ;

  oct_tree_proximity_search( M.parallel() ,
                             global_box ,
                             boxes.size() , & boxes[0] , 0 , NULL ,
                             NULL , proximity , search_stats );
}

//----------------------------------------------------------------------

ProximitySearch::~ProximitySearch() {}

ProximitySearch::ProximitySearch(
  const ProximitySearch::CoordinateField & node_coord ,
  const float box_expansion )
  : m_node_coord( node_coord ),
    m_box_expansion(
      std::numeric_limits<float>::epsilon() * 2 < box_expansion ?
      box_expansion : std::numeric_limits<float>::epsilon() * 2 )
{}

void ProximitySearch::box(
  const Entity & entity , float * const entity_box ) const
{
  if ( entity.entity_type() != Node ) {
    const float f_max = std::numeric_limits<float>::max();

    entity_box[0] = entity_box[1] = entity_box[2] =  f_max ;
    entity_box[3] = entity_box[4] = entity_box[5] = - f_max ;

    PairIterRelation jnode = entity.relations( Node );

    for ( ; jnode ; ++jnode ) {
      const Entity & node = * jnode->entity();

      const double * const coord = field_data( m_node_coord , node );

      if ( coord[0] < entity_box[0] ) { entity_box[0] = (float) coord[0] ; }
      if ( coord[1] < entity_box[1] ) { entity_box[1] = (float) coord[1] ; }
      if ( coord[2] < entity_box[2] ) { entity_box[2] = (float) coord[2] ; }

      if ( entity_box[3] < coord[0] ) { entity_box[3] = (float) coord[0] ; }
      if ( entity_box[4] < coord[1] ) { entity_box[4] = (float) coord[1] ; }
      if ( entity_box[5] < coord[2] ) { entity_box[5] = (float) coord[2] ; }
    }

    float box_inc[3] ;

    box_inc[0] = entity_box[3] - entity_box[0] ;
    box_inc[1] = entity_box[4] - entity_box[1] ;
    box_inc[2] = entity_box[5] - entity_box[2] ;

    float max_inc = box_inc[0] ;
    if ( max_inc < box_inc[1] ) { max_inc = box_inc[1] ; }
    if ( max_inc < box_inc[2] ) { max_inc = box_inc[2] ; }

    max_inc *= m_box_expansion ;

    entity_box[0] -= max_inc ;
    entity_box[1] -= max_inc ;
    entity_box[2] -= max_inc ;

    entity_box[3] += max_inc ;
    entity_box[4] += max_inc ;
    entity_box[5] += max_inc ;
  }
  else {

    const double * const coord = field_data( m_node_coord , entity );

    entity_box[0] = entity_box[3] = (float) coord[0] ;
    entity_box[1] = entity_box[4] = (float) coord[1] ;
    entity_box[2] = entity_box[5] = (float) coord[2] ;
  }
}

int ProximitySearch::part_id( const Kernel & kernel ) const
{
  PartSet part_set ;

  kernel.supersets( part_set );

  Part * part = NULL ;

  PartSet::iterator i ;
  for ( i = part_set.begin() ; part == NULL && i != part_set.end() ; ++i ) {
    const ProximitySearch * const tmp = (*i)->attribute<ProximitySearch>();
    if ( this == tmp ) { part = *i ; }
  }

  return part != NULL ? part->mesh_meta_data_ordinal() + 1 : 0 ;
}

} // namespace phdmesh

