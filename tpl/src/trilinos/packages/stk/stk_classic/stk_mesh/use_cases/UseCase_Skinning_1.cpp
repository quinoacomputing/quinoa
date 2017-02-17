/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <use_cases/UseCase_Skinning.hpp>

#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/BoundaryAnalysis.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

namespace {

static const size_t NODE_RANK = stk_classic::mesh::fem::FEMMetaData::NODE_RANK;

unsigned count_skin_entities( stk_classic::mesh::BulkData & mesh, stk_classic::mesh::Part & skin_part, stk_classic::mesh::EntityRank skin_rank )
{
  const stk_classic::mesh::MetaData & meta = stk_classic::mesh::MetaData::get(mesh);

  stk_classic::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;

  const std::vector<stk_classic::mesh::Bucket*>& buckets = mesh.buckets( skin_rank );

  return count_selected_entities( select_skin, buckets);
}

} // empty namespace

bool skinning_use_case_1(stk_classic::ParallelMachine pm)
{
  bool passed = true;
  {
    //setup the mesh
    stk_classic::mesh::fixtures::HexFixture fixture(pm,3,3,3);

    stk_classic::mesh::fem::FEMMetaData & fem_meta = fixture.m_fem_meta;
    stk_classic::mesh::BulkData & mesh = fixture.m_bulk_data;
    const stk_classic::mesh::EntityRank element_rank = fem_meta.element_rank();
    const stk_classic::mesh::EntityRank side_rank    = fem_meta.side_rank();

    stk_classic::mesh::Part & skin_part = fem_meta.declare_part("skin_part");
    fem_meta.commit();

    fixture.generate_mesh();

    skin_mesh(mesh, element_rank, &skin_part);

    std::vector< stk_classic::mesh::EntityId > elements_to_separate;

    //separate out the middle element
    elements_to_separate.push_back(fixture.elem_id(1,1,1));

    separate_and_skin_mesh(
        fem_meta,
        mesh,
        skin_part,
        elements_to_separate,
        element_rank
        );

    // pointer to middle_element after mesh modification.
    stk_classic::mesh::Entity * middle_element = mesh.get_entity(element_rank, fixture.elem_id(1,1,1));

    unsigned num_skin_entities = count_skin_entities(mesh, skin_part, side_rank);

    stk_classic::all_reduce(pm, stk_classic::ReduceSum<1>(&num_skin_entities));

    //there should be 66 faces in the skin part
    //54 on the outside
    //6 on the inside attached to the entire mesh
    //6 on the inside attected to the element that was detached
    bool correct_skin = ( num_skin_entities == 66 );
    bool correct_relations = true;
    bool correct_comm = true;

    //all nodes connected to the single element that has been broken off
    //should have relations.size() == 4 and comm.size() == 0
    if (middle_element != NULL && middle_element->owner_rank() == mesh.parallel_rank()) {

      stk_classic::mesh::PairIterRelation relations = middle_element->relations(NODE_RANK);

      for (; relations.first != relations.second; ++relations.first) {
        stk_classic::mesh::Entity * current_node = (relations.first->entity());
        //each node should be attached to only 1 element and 3 faces
        correct_relations &= ( current_node->relations().size() == 4 );
        //the entire closure of the element should exist on a single process
        correct_comm      &= ( current_node->comm().size() == 0 );
      }
    }
    passed &= (correct_skin && correct_relations && correct_comm);
  }

  //seperate the entire middle layer of the mesh
  {
    //setup the mesh
    stk_classic::mesh::fixtures::HexFixture fixture(pm,3,3,3);

    stk_classic::mesh::fem::FEMMetaData & fem_meta = fixture.m_fem_meta;
    stk_classic::mesh::BulkData & mesh = fixture.m_bulk_data;
    const stk_classic::mesh::EntityRank element_rank = fem_meta.element_rank();
    const stk_classic::mesh::EntityRank side_rank    = fem_meta.side_rank();

    stk_classic::mesh::Part & skin_part = fem_meta.declare_part("skin_part");
    fem_meta.commit();

    fixture.generate_mesh();

    skin_mesh(mesh, element_rank, &skin_part);

    std::vector< stk_classic::mesh::EntityId > elements_to_separate;

    //separate out the middle level
    elements_to_separate.push_back(fixture.elem_id(1,0,0));
    elements_to_separate.push_back(fixture.elem_id(1,0,1));
    elements_to_separate.push_back(fixture.elem_id(1,0,2));
    elements_to_separate.push_back(fixture.elem_id(1,1,0));
    elements_to_separate.push_back(fixture.elem_id(1,1,1));
    elements_to_separate.push_back(fixture.elem_id(1,1,2));
    elements_to_separate.push_back(fixture.elem_id(1,2,0));
    elements_to_separate.push_back(fixture.elem_id(1,2,1));
    elements_to_separate.push_back(fixture.elem_id(1,2,2));

    separate_and_skin_mesh(
        fem_meta,
        mesh,
        skin_part,
        elements_to_separate,
        element_rank
        );

    // pointer to middle_element after mesh modification.
    unsigned num_skin_entities = count_skin_entities(mesh, skin_part, side_rank);

    stk_classic::all_reduce(pm, stk_classic::ReduceSum<1>(&num_skin_entities));

    //there should be 90 faces in the skin part
    //30 attached to each level of the mesh
    bool correct_skin = ( num_skin_entities == 90 );

    passed &= correct_skin;
  }

  return passed;
}
