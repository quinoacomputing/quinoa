/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/BoundaryAnalysis.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_mesh/fixtures/GridFixture.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <Shards_BasicTopologies.hpp>

#include <iomanip>
#include <algorithm>

static const size_t NODE_RANK = stk_classic::mesh::fem::FEMMetaData::NODE_RANK;

class UnitTestStkMeshSkinning {
public:
  UnitTestStkMeshSkinning(stk_classic::ParallelMachine pm) : m_comm(pm),  m_num_procs(0), m_rank(0)
  {
    m_num_procs = stk_classic::parallel_machine_size( m_comm );
    m_rank = stk_classic::parallel_machine_rank( m_comm );
  }

  void test_skinning();

  stk_classic::ParallelMachine m_comm;
  int m_num_procs;
  int m_rank;
};

namespace {

STKUNIT_UNIT_TEST( UnitTestStkMeshSkinning , testUnit )
{
  UnitTestStkMeshSkinning unit(MPI_COMM_WORLD);
  unit.test_skinning();
}

STKUNIT_UNIT_TEST( UnitTestStkMeshSkinning , testSingleShell )
{
  const int spatial_dimension = 3;
  stk_classic::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(spatial_dimension, stk_classic::mesh::fem::entity_rank_names(spatial_dimension));
  stk_classic::mesh::MetaData & meta_data = stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk_classic::mesh::BulkData bulk_data( meta_data, MPI_COMM_WORLD );
  const stk_classic::mesh::EntityRank element_rank = fem_meta.element_rank();

  const unsigned p_rank = bulk_data.parallel_rank();

  stk_classic::mesh::Part & skin_part = fem_meta.declare_part("skin_part");

  stk_classic::mesh::fem::CellTopology shell_top(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >());
  stk_classic::mesh::Part & shell_part = fem_meta.declare_part("shell_part", shell_top);

  fem_meta.commit();

  bulk_data.modification_begin();

  if ( p_rank == 0 ) {

    stk_classic::mesh::EntityId elem_node[4] ;

    // Query nodes from this simple grid fixture via the (i,j,k) indices.
    elem_node[0] = 1;
    elem_node[1] = 2;
    elem_node[2] = 3;
    elem_node[3] = 4;

    stk_classic::mesh::EntityId elem_id = 1;

    stk_classic::mesh::fem::declare_element( bulk_data, shell_part, elem_id, elem_node);

  }
  bulk_data.modification_end();


  stk_classic::mesh::skin_mesh( bulk_data, element_rank, &skin_part);

  {
    const unsigned mesh_rank = element_rank;
    const stk_classic::mesh::MetaData & meta = stk_classic::mesh::MetaData::get(bulk_data);
    stk_classic::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;
    const std::vector<stk_classic::mesh::Bucket*>& buckets = bulk_data.buckets( mesh_rank -1);
    int num_skin_entities = stk_classic::mesh::count_selected_entities( select_skin, buckets);


    stk_classic::all_reduce(MPI_COMM_WORLD, stk_classic::ReduceSum<1>(&num_skin_entities));

    // Verify that the correct 6 sides are present.

    STKUNIT_ASSERT_EQUAL( num_skin_entities, 2 );
  }
}

} //end namespace

void UnitTestStkMeshSkinning::test_skinning()
{
  // This test will only work for np=1
  if (m_num_procs > 1) {
    return;
  }

  stk_classic::mesh::fixtures::GridFixture grid_mesh(MPI_COMM_WORLD);

  stk_classic::mesh::BulkData& bulk_data = grid_mesh.bulk_data();
  stk_classic::mesh::fem::FEMMetaData& fem_meta = grid_mesh.fem_meta();
  const stk_classic::mesh::EntityRank element_rank = fem_meta.element_rank();

  // Create a part for the skin and the shells
  stk_classic::mesh::Part & skin_part = fem_meta.declare_part("skin_part");
  stk_classic::mesh::fem::CellTopology line_top(shards::getCellTopologyData<shards::ShellLine<2> >());
  stk_classic::mesh::Part & shell_part = fem_meta.declare_part("shell_part", line_top);
  fem_meta.commit();

  // Begin modification cycle
  grid_mesh.bulk_data().modification_begin();

  // Generate the plain grid
  grid_mesh.generate_grid();

  // Add the shells
  std::vector<unsigned> count;
  stk_classic::mesh::Selector locally_owned(fem_meta.locally_owned_part());
  stk_classic::mesh::count_entities(locally_owned, bulk_data, count);
  const unsigned num_shell_1_faces = 4;
  const unsigned num_shell_2_faces = 2;
  const unsigned num_shell_faces = num_shell_1_faces + num_shell_2_faces;
  const unsigned num_entities = count[NODE_RANK] +
                                count[element_rank];

  stk_classic::mesh::PartVector shell_parts;
  shell_parts.push_back(&shell_part);

  std::vector<stk_classic::mesh::Entity*> shell_faces;
  for (unsigned i = 1; i <= num_shell_faces; ++i) {
    stk_classic::mesh::Entity& new_shell = bulk_data.declare_entity(element_rank,
                                                            num_entities + i,
                                                            shell_parts);
    shell_faces.push_back(&new_shell);
  }

  // Set up relationships for shells

  // declare shell relationships for first shell
  unsigned node_list_1[5] = {21, 26, 31, 36, 41};
  for (unsigned i = 0; i < num_shell_1_faces; ++i) {
    stk_classic::mesh::Entity& shell = *(shell_faces[i]);
    stk_classic::mesh::Entity& node1 = *(bulk_data.get_entity(NODE_RANK, node_list_1[i]));
    stk_classic::mesh::Entity& node2 = *(bulk_data.get_entity(NODE_RANK, node_list_1[i+1]));
    bulk_data.declare_relation(shell, node1, 0);
    bulk_data.declare_relation(shell, node2, 1);
  }

  // declare shell relationships for second shell
  unsigned node_list_2[3] = {31, 36, 41};
  for (unsigned i = 0; i < num_shell_2_faces; ++i) {
    stk_classic::mesh::Entity& shell = *(shell_faces[i + num_shell_1_faces]);
    stk_classic::mesh::Entity& node1 = *(bulk_data.get_entity(NODE_RANK, node_list_2[i]));
    stk_classic::mesh::Entity& node2 = *(bulk_data.get_entity(NODE_RANK, node_list_2[i+1]));
    bulk_data.declare_relation(shell, node1, 0);
    bulk_data.declare_relation(shell, node2, 1);
  }

  grid_mesh.bulk_data().modification_end();

  // skin the boundary
  stk_classic::mesh::skin_mesh(bulk_data, element_rank, &skin_part);

  // Grab the skin entities
  stk_classic::mesh::Selector skin_selector(skin_part);
  const std::vector<stk_classic::mesh::Bucket*>& edge_buckets = bulk_data.buckets(fem_meta.edge_rank());
  std::vector<stk_classic::mesh::Entity*> skin_entities;
  stk_classic::mesh::get_selected_entities(skin_selector, edge_buckets, skin_entities);

  unsigned num_expected_skin_entites = 16;
  STKUNIT_EXPECT_EQUAL(num_expected_skin_entites, skin_entities.size());

  // Map the element id to the number of skins associated with that element

  std::map<unsigned, unsigned> results;
  std::map<unsigned, unsigned> expected_results;

  expected_results[1] = 2;
  expected_results[2] = 1;
  expected_results[3] = 1;
  expected_results[4] = 1;
  expected_results[5] = 1;
  expected_results[9] = 1;
  expected_results[13] = 2;
  expected_results[14] = 1;
  expected_results[15] = 1;
  expected_results[16] = 1;
  expected_results[42] = 1;
  expected_results[43] = 1;
  expected_results[44] = 1;
  expected_results[45] = 1;
  expected_results[46] = 1;
  expected_results[47] = 1;

  // Vector of of vector of entities (shells) that are expected to share a skin

  std::set<std::set<unsigned> > sharing;
  std::set<std::set<unsigned> > expected_sharing;

  std::set<unsigned> temp;
  temp.insert(44);
  temp.insert(46);
  expected_sharing.insert(temp);

  temp.clear();
  temp.insert(45);
  temp.insert(47);
  expected_sharing.insert(temp);

  // map skin-id to ids of elements it is attached to; we will use this to
  // compute sharing
  for (std::vector<stk_classic::mesh::Entity*>::const_iterator
       itr = skin_entities.begin(); itr != skin_entities.end(); ++itr) {
    stk_classic::mesh::PairIterRelation upward_relation_itr =
      (*itr)->relations(element_rank);
    bool has_multiple = upward_relation_itr.size() > 1;
    std::set<unsigned> sharing_elements;
    for ( ; !upward_relation_itr.empty() ; ++upward_relation_itr ) {
      unsigned elem_id = upward_relation_itr->entity()->identifier();
      if (results.find(elem_id) != results.end()) {
        ++results[elem_id];
      }
      else {
        results[elem_id] = 1;
      }

      if (has_multiple) {
        sharing_elements.insert(elem_id);
      }
    }
    if (has_multiple) {
      sharing.insert(sharing_elements);
    }
  }

  STKUNIT_EXPECT_TRUE(results == expected_results);
  STKUNIT_EXPECT_TRUE(sharing == expected_sharing);
}

