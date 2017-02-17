/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/UseCase_Rebal_4.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/CreateAdjacentEntities.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

#include <stk_rebalance_utils/RebalanceUtils.hpp>

//----------------------------------------------------------------------

using namespace stk_classic::mesh::fixtures;

typedef stk_classic::mesh::Field<double> ScalarField ;

namespace stk_classic {
namespace rebalance {
namespace use_cases {

class GreedySideset : public Partition {
  public :
  struct MeshInfo {
    std::vector<mesh::Entity *>      mesh_entities;
    const VectorField * nodal_coord_ref ;
    const ScalarField * elem_weight_ref;
    std::vector<unsigned>            dest_proc_ids ;

    /** Default Constructor. */
    MeshInfo():
      nodal_coord_ref(NULL),
      elem_weight_ref(NULL) {}

    /** Destructor. */
    ~MeshInfo() {}
  };
  explicit GreedySideset(ParallelMachine pm,
                         const stk_classic::mesh::PartVector & surfaces,
                         mesh::BulkData   & bulk_data);
  virtual ~GreedySideset();
  virtual void reset_dest_proc_data();
  virtual void set_mesh_info ( const std::vector<mesh::Entity *> &mesh_entities,
                               const VectorField   * nodal_coord_ref,
                               const ScalarField   * elem_weight_ref=NULL);
  virtual void determine_new_partition(bool &RebalancingNeeded);
  virtual unsigned num_elems() const;
  virtual int get_new_partition(stk_classic::mesh::EntityProcVec &new_partition);
  virtual bool partition_dependents_needed()const;
  bool find_mesh_entity(const mesh::Entity * entity, unsigned & moid) const;
  unsigned destination_proc(const unsigned moid) const;
  void set_destination_proc(const unsigned moid, const unsigned proc );
  MeshInfo  mesh_information_;
  unsigned  total_number_entities_;
  const stk_classic::mesh::PartVector & surfaces_;
  mesh::BulkData   & bulk_data_;
};

GreedySideset::GreedySideset(ParallelMachine pm,
                             const stk_classic::mesh::PartVector & surfaces,
                             mesh::BulkData   & bulk_data) :
  stk_classic::rebalance::Partition(pm),
  mesh_information_(),
  surfaces_(surfaces),
  bulk_data_(bulk_data) {}
GreedySideset::~GreedySideset() {}
void GreedySideset::reset_dest_proc_data() {
  const int proc = parallel_machine_rank(comm_);
  const unsigned size = mesh_information_.mesh_entities.size();
  mesh_information_.dest_proc_ids.assign(size, proc);
}
void GreedySideset::set_mesh_info ( const std::vector<mesh::Entity *> &mesh_entities,
                                    const VectorField   * nodal_coord_ref,
                                    const ScalarField   * elem_weight_ref){
  MeshInfo mesh_info;

  /* Keep track of the total number of elements. */
  total_number_entities_ = mesh_entities.size();

  mesh_info.mesh_entities = mesh_entities;
  mesh_info.nodal_coord_ref = nodal_coord_ref;
  mesh_info.elem_weight_ref = elem_weight_ref;

  /** Default destination for an entity is the processor
      that already owns the entity, which is this processor.
      The length of the dest_proc_ids vector is the same
      length as the mesh_entities vector.
  */
  mesh_info.dest_proc_ids.assign(mesh_entities.size(), stk_classic::parallel_machine_rank(comm_));

  mesh_information_ = mesh_info;
}

unsigned GreedySideset::num_elems() const {return total_number_entities_ ;}
int GreedySideset::get_new_partition(stk_classic::mesh::EntityProcVec &new_partition){
std::vector<mesh::Entity*>::iterator i=mesh_information_.mesh_entities.begin();
std::vector<unsigned>     ::iterator j=mesh_information_.dest_proc_ids.begin();
  for (;i != mesh_information_.mesh_entities.end(),
        j != mesh_information_.dest_proc_ids.end();
        ++i,++j) {
    mesh::Entity * mesh_entity = *i;
    unsigned proc = *j;
    mesh::EntityProc et(mesh_entity, proc);
    new_partition.push_back(et);
  }
  return 0;
}
bool GreedySideset::partition_dependents_needed()const{return true;}

bool GreedySideset::find_mesh_entity(const mesh::Entity * entity, unsigned & moid) const
{
  unsigned len = mesh_information_.mesh_entities.size();
  for(moid = 0; moid < len; ++moid)
  {
    if(mesh_information_.mesh_entities[moid] == entity) return true;
  }
  return false;
}
unsigned GreedySideset::destination_proc(const unsigned moid) const
{
  return mesh_information_.dest_proc_ids[ moid ];
}
void GreedySideset::set_destination_proc(const unsigned moid,
                                         const unsigned proc )
{
  mesh_information_.dest_proc_ids[ moid ] = proc;
}


void GreedySideset::determine_new_partition(bool &RebalancingNeeded) {

  reset_dest_proc_data();

  stk_classic::mesh::fem::FEMMetaData & fem_meta = stk_classic::mesh::fem::FEMMetaData::get(bulk_data_);
  const stk_classic::mesh::EntityRank side_rank = fem_meta.side_rank();
  const stk_classic::mesh::EntityRank elem_rank = fem_meta.element_rank();

  // Select active ghosted side faces.
  stk_classic::mesh::Selector selector(!fem_meta.locally_owned_part() &
                                stk_classic::mesh::selectIntersection(surfaces_));

  mesh::EntityVector sides;
  mesh::get_selected_entities(selector, bulk_data_.buckets(side_rank), sides);

  const unsigned p_rank = bulk_data_.parallel_rank();
  size_t local_changes = 0;
  const unsigned nSide = sides.size();
  for(unsigned iSide = 0; iSide < nSide; ++iSide)
  {
    const mesh::Entity & side = *sides[iSide];
    const unsigned sideProc = side.owner_rank();
    ThrowRequireMsg(sideProc!=p_rank,
     "When iterating Non-locally owned sides, found a locally owned side.");

    stk_classic::mesh::PairIterRelation iElem = side.relations(elem_rank);
    for ( ; iElem.first != iElem.second; ++iElem.first ) {
      const mesh::Entity & elem = *iElem.first->entity();
      unsigned moid;
      const bool mesh_entity_found = find_mesh_entity(&elem, moid);
      if (mesh_entity_found) {
        const unsigned elemProc = elem.owner_rank();
        ThrowRequireMsg(elemProc==p_rank,
          "When iterating locally owned elements, found a non-locally owned element.");
        const unsigned destProc = destination_proc(moid);
        ThrowRequireMsg(destProc==p_rank || destProc==sideProc,
         " Sanity check failed: "
         "It's possible that an element is connected to "
         "two sides that are owned by different procs.  We don't "
         "yet handle that situation here but we can, at least, "
         "detect it. ");
        if(elemProc != sideProc)
        {
          ++local_changes;
          set_destination_proc(moid, sideProc);
        }
      }
    }
  }
  size_t global_changes = 0;
  stk_classic::all_reduce_sum (comm_, &local_changes, &global_changes, 1);
  RebalancingNeeded = global_changes > 0;
}

enum { nx = 3, ny = 3 };

bool test_greedy_sideset ( stk_classic::ParallelMachine comm )
{
  unsigned spatial_dimension = 2;
  std::vector<std::string> rank_names = stk_classic::mesh::fem::entity_rank_names(spatial_dimension);
  stk_classic::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(spatial_dimension, rank_names);
  stk_classic::mesh::MetaData & meta_data = stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk_classic::mesh::BulkData bulk_data( meta_data , comm , 100 );
  const stk_classic::mesh::EntityRank element_rank    = fem_meta.element_rank();
  const stk_classic::mesh::EntityRank node_rank       = fem_meta.node_rank();

  stk_classic::mesh::fem::CellTopology quad_top(shards::getCellTopologyData<shards::Quadrilateral<4> >());
  stk_classic::mesh::fem::CellTopology line_top(shards::getCellTopologyData<shards::Line<2> >());
  stk_classic::mesh::Part & quad_part( fem_meta.declare_part("quad", quad_top ) );
  stk_classic::mesh::Part & side_part( fem_meta.declare_part("line", line_top ) );
  VectorField & coord_field( fem_meta.declare_field< VectorField >( "coordinates" ) );
  ScalarField & weight_field( fem_meta.declare_field< ScalarField >( "element_weights" ) );

  stk_classic::mesh::put_field( coord_field , node_rank , fem_meta.universal_part() );
  stk_classic::mesh::put_field(weight_field , element_rank , fem_meta.universal_part() );

  fem_meta.commit();
  const unsigned p_rank = bulk_data.parallel_rank();
  bulk_data.modification_begin();

  if ( !p_rank ) {

    std::vector<std::vector<stk_classic::mesh::Entity*> > quads(nx);
    for ( unsigned ix = 0 ; ix < nx ; ++ix ) quads[ix].resize(ny);

    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk_classic::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk_classic::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk_classic::mesh::Entity &q = stk_classic::mesh::fem::declare_element( bulk_data , quad_part , elem , nodes );
        quads[ix][iy] = &q;
      }
    }

    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk_classic::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk_classic::mesh::Entity * e = bulk_data.get_entity( element_rank, elem );
        double * const e_weight = stk_classic::mesh::field_data( weight_field , *e );
        *e_weight = 1.0;
      }
    }
    for ( unsigned iy = 0 ; iy <= ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix <= nx ; ++ix ) {
        stk_classic::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk_classic::mesh::Entity * n = bulk_data.get_entity( node_rank, nid );
        double * const coord = stk_classic::mesh::field_data( coord_field , *n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }
  }

  bulk_data.modification_end();

  // create some sides and faces to rebalance.
  stk_classic::mesh::PartVector add_parts;
  stk_classic::mesh::create_adjacent_entities(bulk_data, add_parts);

  bulk_data.modification_begin();

  const stk_classic::mesh::PartVector surfaces(1, &side_part);
  {
    const stk_classic::mesh::PartVector empty_remove_parts;
    stk_classic::mesh::fem::FEMMetaData & fmeta = stk_classic::mesh::fem::FEMMetaData::get(bulk_data);
    const stk_classic::mesh::EntityRank side_rank = fmeta.side_rank();
    stk_classic::mesh::Selector selector2( fmeta.locally_owned_part());
    mesh::EntityVector sides;
    mesh::get_selected_entities(selector2, bulk_data.buckets(side_rank), sides);

    const unsigned nSide = sides.size();
    for(unsigned iSide = 0; iSide < nSide; ++iSide)
    {
      mesh::Entity & side = *sides[iSide];
      if (side.identifier()==7) {
        bulk_data.change_entity_parts(side, surfaces, empty_remove_parts);
      }
    }
  }
  bulk_data.modification_end();

  // Zoltan partition is specialized form a virtual base class, stk_classic::rebalance::Partition.
  // Other specializations are possible.
  Teuchos::ParameterList emptyList;
  stk_classic::rebalance::Zoltan zoltan_partition(comm, spatial_dimension, emptyList);
  stk_classic::mesh::Selector selector3(fem_meta.locally_owned_part());
  stk_classic::rebalance::rebalance(bulk_data, selector3, &coord_field, NULL, zoltan_partition);
  {
    const int  print_stats = 1;
    int        nentity     = 0;
    double     entity_wgt  = 0;
    int        ncuts       = 0;
    double     cut_wgt     = 0;
    int        nboundary   = 0;
    int        nadj        = 0;
    const int ierr = zoltan_partition.evaluate (print_stats, &nentity, &entity_wgt, &ncuts, &cut_wgt, &nboundary, &nadj);
    std::cout <<" Information returned from the Zoltan evaluate function:"<<std::endl;
    std::cout <<" Error Code:             :"<<ierr      <<std::endl;
    std::cout <<" Number of entities:     :"<<nentity   <<std::endl;
    std::cout <<" Number of cuts:         :"<<ncuts     <<std::endl;
    std::cout <<" Cut Weight:             :"<<cut_wgt   <<std::endl;
    std::cout <<" Number on Boundary:     :"<<nboundary <<std::endl;
    std::cout <<" Number Adjancent:       :"<<nadj      <<std::endl;
    {
      stk_classic::mesh::fem::FEMMetaData & fmeta = stk_classic::mesh::fem::FEMMetaData::get(bulk_data);
      const stk_classic::mesh::EntityRank side_rank = fmeta.side_rank();
      const stk_classic::mesh::EntityRank elem_rank = fmeta.element_rank();
      const mesh::Entity *s = bulk_data.get_entity(side_rank,7);
      if (s) {
        const mesh::Entity & side = *s;
        if (p_rank == side.owner_rank()) {
          stk_classic::mesh::PairIterRelation iElem = side.relations(elem_rank);
          for ( ; iElem.first != iElem.second; ++iElem.first ) {
            const mesh::Entity & elem = *iElem.first->entity();
            const unsigned elemProc = elem.owner_rank();
            if (elemProc!=p_rank) {
              std::cout <<p_rank<<" Good: Found element of of side 7 not owned."
                        <<" Element "<<elemProc
                        <<" on processor "<<p_rank
                        <<" owned by "<<elemProc<<std::endl;
            }
          }
        }
      }
    }
  }

  const double imbalance_threshold = 1.5;
  bool do_rebal = imbalance_threshold < stk_classic::rebalance::check_balance(bulk_data, NULL, element_rank);

  if( !p_rank )
    std::cout << std::endl
     << "Use Case 4: imbalance_threshold after rebalance 1 = " << imbalance_threshold <<", "<<do_rebal << std::endl;

  {
    stk_classic::rebalance::use_cases::GreedySideset greedy_sideset(comm, surfaces, bulk_data);
    stk_classic::mesh::Selector selector4(fem_meta.locally_owned_part());
    stk_classic::rebalance::rebalance(bulk_data, selector4, &coord_field, NULL, greedy_sideset);
  }

  do_rebal = imbalance_threshold < stk_classic::rebalance::check_balance(bulk_data, NULL, element_rank);

  if( !p_rank )
    std::cout << std::endl
     << "Use Case 4: imbalance_threshold after rebalance 2 = " << imbalance_threshold <<", "<<do_rebal << std::endl;
  {
    stk_classic::mesh::fem::FEMMetaData & fmeta = stk_classic::mesh::fem::FEMMetaData::get(bulk_data);
    const stk_classic::mesh::EntityRank side_rank = fmeta.side_rank();
    const stk_classic::mesh::EntityRank elem_rank = fmeta.element_rank();
    mesh::Entity *s = bulk_data.get_entity(side_rank,7);
    if (s) {
      mesh::Entity & side = *s;
      if (p_rank == side.owner_rank()) {
        stk_classic::mesh::PairIterRelation iElem = side.relations(elem_rank);
        for ( ; iElem.first != iElem.second; ++iElem.first ) {
          const mesh::Entity & elem = *iElem.first->entity();
          const unsigned elemProc = elem.owner_rank();
          if (elemProc!=p_rank) {
            std::cerr <<p_rank<<" Error: Found element of of side 7 not owned:"<<elemProc<<std::endl;
          }
          ThrowRequireMsg(elemProc==p_rank,
           "Use case 4 error check failed. Found element of of side 7 not owned.");
        }
      }
    }
  }
  // Check that we satisfy our threshhold
  bool result = !do_rebal ;

  // And verify that all dependent entities are on the same proc as their parent element
  {
    stk_classic::mesh::EntityVector entities;
    stk_classic::mesh::Selector selector5 = fem_meta.locally_owned_part();

    get_selected_entities(selector5, bulk_data.buckets(node_rank), entities);
    result &= verify_dependent_ownership(element_rank, entities);
    get_selected_entities(selector5, bulk_data.buckets(fem_meta.side_rank()), entities);
    result &= verify_dependent_ownership(element_rank, entities);
  }



  return result;
}

} //namespace use_cases
} //namespace rebalance
} //namespace stk_classic


