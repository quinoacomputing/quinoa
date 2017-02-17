/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_search_util/stk_mesh/PrintEntityProc.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>

namespace stk_classic {
namespace search_util {

typedef stk_classic::search::ident::IdentProc<stk_classic::mesh::EntityKey, unsigned> IdentProc;

typedef std::vector<std::pair<IdentProc, IdentProc> > IdentProcRelation;

// Used to output the results of a coarse or direct search to
// verify which entity contains another entity.
void print_entity_map(stk_classic::diag::Writer &writer,
                      const std::vector<std::pair<stk_classic::mesh::Entity*, stk_classic::mesh::Entity*> >& entity_map,
                      const std::string & relation)
{
  if (writer.shouldPrint()) {
    size_t size = entity_map.size();
    for (size_t i=0; i < size; i++) {
      stk_classic::mesh::EntityKey key1 = entity_map[i].first->key();
      stk_classic::mesh::EntityKey key2 = entity_map[i].second->key();
      const stk_classic::mesh::MetaData& meta1 = stk_classic::mesh::MetaData::get(*(entity_map[i].first));
      const stk_classic::mesh::MetaData& meta2 = stk_classic::mesh::MetaData::get(*(entity_map[i].second));

      writer << "[" << i << "] "
             << meta1.entity_rank_name(stk_classic::mesh::entity_rank(key1)) << " "
             << stk_classic::mesh::entity_id(key1) << relation
             << meta2.entity_rank_name(stk_classic::mesh::entity_rank(key2)) << " "
             << stk_classic::mesh::entity_id(key2) << "\n";
    }
  }
}

/**
 * Used to output a sharing or ghosting vector in human readable
 * form.
 * The "action" argument will typically be "Share " or "Ghost "
 * The "to_from" argument will typically be
 * - for sharing " with "
 * - for ghosting " from " or " to "
 *
 * Decodes the entity key and prints as entity type, entity id
 *
 * Example output: "Share NODE 37 with processor 12"
 */
void print_entity_proc_map(stk_classic::diag::Writer &writer,
                           const std::vector<stk_classic::mesh::EntityProc>& entity_proc,
                           const std::string &action,
                           const std::string &to_from)
{
  if (writer.shouldPrint()) {
    size_t size = entity_proc.size();
    for (size_t i=0; i < size; i++) {
      stk_classic::mesh::EntityKey key = entity_proc[i].first->key();
      const stk_classic::mesh::MetaData& meta = stk_classic::mesh::MetaData::get( *(entity_proc[i].first) );

      writer << "[" << i << "] "
             << action
             << meta.entity_rank_name(stk_classic::mesh::entity_rank(key)) << " "
             << stk_classic::mesh::entity_id(key) << " " << to_from << " processor "
             << entity_proc[i].second << "\n";
    }
  }
}

void print_entity_proc_map(stk_classic::diag::Writer &writer,
                           const std::vector<stk_classic::mesh::Entity*>& entity_proc,
                           const std::string &action,
                           const std::string &to_from)
{
  if (writer.shouldPrint()) {
    size_t size = entity_proc.size();
    for (size_t i=0; i < size; i++) {
      stk_classic::mesh::EntityKey key = entity_proc[i]->key();
      const stk_classic::mesh::MetaData& meta = stk_classic::mesh::MetaData::get( entity_proc[i]->bucket() );

      writer << "[" << i << "] "
             << action
             << meta.entity_rank_name(stk_classic::mesh::entity_rank(key)) << " "
             << stk_classic::mesh::entity_id(key) << " " << to_from << " processor "
             << entity_proc[i]->owner_rank() << "\n";
    }
  }
}


void print_entity_proc_map( stk_classic::diag::Writer & writer ,
                            const stk_classic::mesh::BulkData & mesh )
{
  const stk_classic::mesh::MetaData & meta = stk_classic::mesh::MetaData::get(mesh);
  const std::vector<stk_classic::mesh::Entity*> & comm = mesh.entity_comm();
  const std::vector<stk_classic::mesh::Ghosting*> & ghost = mesh.ghostings();

  size_t counter = 0 ;

  for ( size_t ig = 0 ; ig < ghost.size() ; ++ig ) {

    const stk_classic::mesh::Ghosting & g = * ghost[ig] ;

    writer << "P" << mesh.parallel_rank()
           << " " << g.name() << " Communication:" << std::endl ;

    for ( std::vector<stk_classic::mesh::Entity*>::const_iterator
          i = comm.begin() ; i != comm.end() ; ++i ) {

      const stk_classic::mesh::Entity & entity = **i ;

      std::vector<unsigned> procs ;

      stk_classic::mesh::comm_procs( g , entity , procs );

      if ( ! procs.empty() ) {
        writer << "[" << counter << "] "
               << meta.entity_rank_name( entity.entity_rank() )
               << "[" << entity.identifier() << " " ;
        if ( entity.owner_rank() != mesh.parallel_rank() ) {
          writer << "not_" ;
        }
        writer << "owned ] {" ;
        for ( size_t j = 0 ; j < procs.size() ; ++j ) {
          writer << " " << procs[j] ;
        }
        writer << " }" << std::endl ;
      }
    }
  }
}


/**
 * Used to output the results of a relation vector in human
 * readable form.  This function cannot be used as the default
 * output of an IdentProcRelation since it is using knowledge that
 * what is really being stored in the IdentProc is stk_classic::mesh
 * entity keys.
 */
void print_stk_mesh_relation_map(
  stk_classic::diag::Writer &writer,
  const std::vector<std::string> &entity_names,
  IdentProcRelation relation)
{
  if (writer.shouldPrint()) {
    size_t size = relation.size();
    writer << "relation  [size " << size << "]\n";
    for (size_t i=0; i < size; i++) {
      IdentProc domain = relation[i].first;
      IdentProc range  = relation[i].second;

//       stk_classic::mesh::EntityKey domain_entity_key;
//       stk_classic::mesh::EntityKey range_entity_key;
//       domain_entity_key.value(domain.ident);
//       range_entity_key.value(range.ident);

      stk_classic::mesh::EntityKey domain_entity_key(domain.ident);
      stk_classic::mesh::EntityKey range_entity_key(range.ident);

      writer << "[" << i << "] ("
             << entity_names[stk_classic::mesh::entity_rank(domain_entity_key)] << " "
             << stk_classic::mesh::entity_id(domain_entity_key)
             << ", proc " << domain.proc
             << "    ->    "
             << entity_names[stk_classic::mesh::entity_rank(range_entity_key)] << " "
             << stk_classic::mesh::entity_id(range_entity_key)
             << ", proc " << range.proc
             << ")\n";
    }
  }
}
}
}
