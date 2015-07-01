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

#ifndef phdmesh_Comm_hpp
#define phdmesh_Comm_hpp

//----------------------------------------------------------------------
/** @file
 *  @brief BulkData entity relation across processor boundaries.
 *
 *  A parallel entity relation matches a domain mesh entity residing on
 *  a domain processor with a range entity residing on a range processor.
 *
 *  { ( ( domain_entity , domain_proc ) , ( range_entity , range_proc ) ) }
 *
 *  This relation is partitioned among processors as follows:
 *
 *    on the domain_proc : { ( domain_entity , range_proc )[i] }
 *    on the range_proc  : { ( range_entity , domain_proc )[i] }
 *
 *  Members of the relation should be ordered by
 *  domain entity identifer and then by range processor.
 *  Thus on the domain processor members are fully sorted
 *  and can be searched accordingly.
 *  However, on the range processor the ordering is conformal
 *  to the domain processor array.  Thus, for an asymmetric
 *  parallel relation the array on the range processor is
 *  unlikely to be sorted.
 */
//----------------------------------------------------------------------

#include <vector>

#include <util/Parallel.hpp>
#include <util/OctTree.hpp>
#include <mesh/Types.hpp>
#include <mesh/EntityComm.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

//----------------------------------------------------------------------
/** Sort and unique an EntityProc array.
 */
void sort_unique( std::vector<EntityProc> & );

/** Sanity check locally for non-null, same-mesh, off-processor,
 *  and proper ordering.  Return an error string if a problem.
 */
bool verify( const std::vector<EntityProc> & , std::string & );

/** Find the first entry corresponding to the given entity.
 *  The array must be properly sorted.
 */
std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & , Entity & );

std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & , unsigned );

std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & , const EntityProc & );

/** Find the first entry corresponding to the given entity.
 *  The array must be properly sorted.
 */
std::vector<EntityProc>::iterator
lower_bound( std::vector<EntityProc> & , Entity & );

std::vector<EntityProc>::iterator
lower_bound( std::vector<EntityProc> & , const EntityProc & );

//----------------------------------------------------------------------
/** Sanity check on existing or potential parallel relation information.
 *  If the result is invalid then outputs a string with an explanation.
 *  Symmetric version of verification.
 */
bool comm_verify( ParallelMachine ,
                  const std::vector<EntityProc> & ,
                  std::string & );

/** Sanity check on existing or potential parallel relation information.
 *  If the result is invalid then outputs a string with an explanation.
 *  Asymmetric version of verification.
 */
bool comm_verify( ParallelMachine ,
                  const std::vector<EntityProc> & ,
                  const std::vector<EntityProc> & ,
                  std::string & );

//----------------------------------------------------------------------
/** Global counts and maximum identifiers for a mesh's entities.
 *  Array sizes:
 *    counts[ EntityTypeEnd ]
 *    max_id[ EntityTypeEnd ]
 */
bool comm_mesh_stats( BulkData & ,
                      entity_id_type * const counts ,
                      entity_id_type * const max_id ,
                      bool local_flag = false );

//----------------------------------------------------------------------
/** Send local mesh entities to 'recv_mesh' according to 'send'.
 *  Received mesh entities are identified in 'recv'.
 */
void comm_copy(
  BulkData & send_mesh ,
  BulkData & recv_mesh ,
  const std::vector<EntityProc> & send ,
        std::vector<EntityProc> & recv );

//----------------------------------------------------------------------
/** Communicate mesh entities from the send mesh to the receive mesh.
 *  The mesh manager is used to map parts between different meshes
 *  and incorporate received mesh entities into the receive mesh.
 */
bool communicate_entities(
  const EntityComm & manager ,
  BulkData & send_mesh ,
  BulkData & recv_mesh ,
  const std::vector<EntityProc> & send ,
        std::vector<EntityProc> & recv ,
  bool local_flag = false );

//----------------------------------------------------------------------
/** Verify that the shared entity values are bit-wise identical */

bool comm_verify_shared_entity_values(
  const BulkData & , EntityType , const FieldBase & );

//----------------------------------------------------------------------
/** Discover the sharing of all mesh entities by searching for
 *  duplicate identifiers on different processors.
 */
void comm_mesh_discover_sharing( BulkData & );

void comm_mesh_add_sharing( BulkData & , const std::vector<EntityProc> & );

/** Scrub shared entities of any that are 
 *  not owned and not used by an owned entity.
 */
bool comm_mesh_scrub_sharing( BulkData & M );

//----------------------------------------------------------------------
/** Verify parallel consistency of mesh entities' identifiers,
 *  owner processor, and mesh parts.
 *  Return false on all processors if an error is detected
 *  and log the error to 'std::cerr'.
 */
bool comm_mesh_verify_parallel_consistency( BulkData & M );

//----------------------------------------------------------------------
/** Generate all aura entities attached to shared nodes.  */
void comm_mesh_regenerate_aura( BulkData & );

//----------------------------------------------------------------------
/** Rebalance the mesh using the HSFC algorithm.
 *  The shared node coordinate field values must be consistent.
 *  BulkData entity sharing and local auras are updated.
 *  Return the HSFC cut keys used for the rebalance.
 */
void comm_mesh_rebalance( BulkData & ,
                          const Field<double,Cartesian> & node_coord_field ,
                          const Field<float> * const elem_weight_field ,
                          std::vector<OctTreeKey> & cut_keys );

//----------------------------------------------------------------------

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

