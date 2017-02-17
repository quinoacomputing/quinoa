/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/fixtures/BoxFixture.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_mesh/fixtures/QuadFixture.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/FieldData.hpp>

using stk_classic::mesh::Part;

// UnitTestBulkData_new is the beginnings of a refactoring of the bulk
// data unit test.  It relies on a customized BoxFixture to rapidly
// create a mesh for testing.

namespace {

void new_insert_transitive_closure( std::set<stk_classic::mesh::EntityProc,stk_classic::mesh::EntityLess> &  ,
					 const stk_classic::mesh::EntityProc & entry );
void new_comm_sync_send_recv(
   stk_classic::mesh::BulkData & mesh ,
   std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess > & new_send ,
   std::set< stk_classic::mesh::Entity * , stk_classic::mesh::EntityLess > & new_recv );

void new_comm_recv_to_send(
  stk_classic::mesh::BulkData & mesh ,
  const std::set< stk_classic::mesh::Entity * , stk_classic::mesh::EntityLess > & new_recv ,
        std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess > & new_send );

/**
 * The customized box fixture used in this file for testing. This fixture
 * is similar to the BoxFixture it inherits from, with the only difference
 * being the extra parts that this fixture declares for testing purposes.
 */
class TestBoxFixture : public stk_classic::mesh::fixtures::BoxFixture
{
 public:
  TestBoxFixture(stk_classic::ParallelMachine pm = MPI_COMM_WORLD,
                 unsigned block_size = 1000) :
    BoxFixture(pm, block_size),
    m_test_part ( m_fem_meta.declare_part ( "Test Part" ) ),
    m_cell_part ( m_fem_meta.declare_part ( "Cell list" , 3 /*max rank*/ ) ),
    m_part_A_0 ( m_fem_meta.declare_part ( "Part A 0", 0 ) ),
    m_part_A_1 ( m_fem_meta.declare_part ( "Part A 1", 1 ) ),
    m_part_A_2 ( m_fem_meta.declare_part ( "Part A 2", 2 ) ),
    m_part_A_3 ( m_fem_meta.declare_part ( "Part A 3", 3 ) ),
    m_part_A_superset ( m_fem_meta.declare_part ( "Part A superset" ) ),
    m_part_B_0 ( m_fem_meta.declare_part ( "Part B 0", 0 ) ),
    m_part_B_1 ( m_fem_meta.declare_part ( "Part B 1", 1 ) ),
    m_part_B_2 ( m_fem_meta.declare_part ( "Part B 2", 2 ) ),
    m_part_B_3 ( m_fem_meta.declare_part ( "Part B 3", 3 ) ),
    m_part_B_superset ( m_fem_meta.declare_part ( "Part B superset" ) )
  {
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_0 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_1 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_2 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_3 );

    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_0 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_1 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_2 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_3 );

    // None of the tests currently need to make any addtional changes
    // to MetaData; if this changes, the line below will have to be
    // removed.
    m_fem_meta.commit();
  }

  Part & get_test_part () { return m_test_part; }
  Part & get_cell_part () { return m_cell_part; }

  Part & get_part_a_0 () { return m_part_A_0; }
  Part & get_part_a_1 () { return m_part_A_1; }
  Part & get_part_a_2 () { return m_part_A_2; }
  Part & get_part_a_3 () { return m_part_A_3; }

  Part & get_part_a_superset () { return m_part_A_superset; }

  Part & get_part_b_0 () { return m_part_B_0; }
  Part & get_part_b_1 () { return m_part_B_1; }
  Part & get_part_b_2 () { return m_part_B_2; }
  Part & get_part_b_3 () { return m_part_B_3; }

  Part & get_part_b_superset () { return m_part_B_superset; }

 private:
  Part     & m_test_part;   // A simple part
  Part     & m_cell_part;   // A part to put cells in

  Part     & m_part_A_0;
  Part     & m_part_A_1;
  Part     & m_part_A_2;
  Part     & m_part_A_3;

  Part     & m_part_A_superset;

  Part     & m_part_B_0;
  Part     & m_part_B_1;
  Part     & m_part_B_2;
  Part     & m_part_B_3;

  Part     & m_part_B_superset;
};

}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyAssertOwnerDeletedEntity )
{
  TestBoxFixture fixture;

  stk_classic::mesh::BulkData         &bulk = fixture.bulk_data();
  stk_classic::mesh::Part             &new_part = fixture.get_test_part ();
  stk_classic::mesh::PartVector        add_part;
  add_part.push_back ( &new_part );

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  STKUNIT_ASSERT(bulk.modification_end());

  // Find a cell owned by this process
  stk_classic::mesh::Entity  *cell_to_delete = NULL;
  stk_classic::mesh::Entity  *cell_to_delete_copy = NULL;
  std::vector<stk_classic::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(3).begin();
  while ( cur_bucket != bulk.buckets(3).end() )
  {
    stk_classic::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() == fixture.comm_rank() )
      {
        cell_to_delete = &*cur_entity;
        break;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  STKUNIT_ASSERT ( cell_to_delete != NULL );
  cell_to_delete_copy = cell_to_delete;
  bulk.modification_begin();
  bulk.destroy_entity ( cell_to_delete );
  // Destroying an already destroyed entity returns false
  STKUNIT_ASSERT( false == bulk.destroy_entity( cell_to_delete_copy ) );
  bulk.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyDetectsBadKey )
{
  TestBoxFixture fixture;

  stk_classic::mesh::BulkData         &bulk = fixture.bulk_data();
  stk_classic::mesh::Part             &new_part = fixture.get_test_part ();
  stk_classic::mesh::PartVector        add_part, empty_vector;
  add_part.push_back ( &new_part );

  stk_classic::mesh::EntityKey bad_key1 ( 45 , 1 );  // Bad entity rank
  stk_classic::mesh::EntityKey bad_key2 ( 1 , 0 );   // Bad id

  STKUNIT_ASSERT_THROW ( bulk.declare_entity(bad_key1.rank(),
                                             bad_key1.id(),
                                             empty_vector),
                         std::logic_error );
  STKUNIT_ASSERT_THROW ( bulk.declare_entity(bad_key2.rank(),
                                             bad_key2.id(),
                                             empty_vector),
                         std::logic_error );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyDetectsNonOwnerChange )
{
  // Set up a mesh where there are shared nodes. Take one of the nodes, and
  // have the non-owning processes try to make a change to that node; this
  // should cause an exception.

  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;
  unsigned p_size = stk_classic::parallel_machine_size(pm);
  unsigned p_rank = stk_classic::parallel_machine_rank(pm);

  stk_classic::mesh::fixtures::QuadFixture fixture(pm, 1 /*nx*/, p_size /*ny*/);
  fixture.m_fem_meta.commit();
  fixture.generate_mesh();
  stk_classic::mesh::BulkData & bulk = fixture.m_bulk_data;

  stk_classic::mesh::PartVector empty_vector;

  stk_classic::mesh::Entity* shared_node = fixture.node(1 /*x*/, 1 /*y*/);
  // Assert that this node is shared
  if ( p_size > 1 && shared_node && (p_rank == 0 || p_rank == 1) ) {
    STKUNIT_ASSERT_GE(shared_node->sharing().size(), 1u);
  }

  bulk.modification_begin();

  // Non-owners of shared_node will attempt to make a change to it; this should
  // cause an exception
  if (shared_node && p_rank != shared_node->owner_rank()) {
    STKUNIT_ASSERT_THROW(bulk.change_entity_parts(*shared_node,
                                                  empty_vector,  //add parts
                                                  empty_vector), //rem parts
                         std::logic_error);
  }

  bulk.modification_end();
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyGetEntityGuards )
{
  TestBoxFixture fixture;

  stk_classic::mesh::BulkData      &bulk = fixture.bulk_data();
  STKUNIT_ASSERT_THROW ( bulk.get_entity ( 1 , 0 ) , std::logic_error );
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyExplicitAddInducedPart )
{
  TestBoxFixture fixture;
  stk_classic::mesh::BulkData     &bulk = fixture.bulk_data ();
  stk_classic::mesh::PartVector    empty_vector;
  stk_classic::mesh::PartVector    cell_part_vector;

  bulk.modification_begin();

  stk_classic::mesh::Entity &new_cell = bulk.declare_entity ( 3 , fixture.comm_rank()+1 , empty_vector );
  stk_classic::mesh::Entity &new_node = bulk.declare_entity ( 0 , fixture.comm_rank()+1 , empty_vector );

  bulk.declare_relation ( new_cell , new_node , 1 );

  cell_part_vector.push_back ( &fixture.get_cell_part () );
  bulk.change_entity_parts ( new_cell , cell_part_vector );
#ifdef SIERRA_MIGRATION
  bulk.change_entity_parts ( new_node , cell_part_vector );
#else
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_node , cell_part_vector ) , std::runtime_error );
#endif
}

/************************
 * This unit test is not possible currently because of the lack of
 * separation between internal part modification routines and public
 * part modification routines.
STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCannotRemoveFromSpecialParts )
{
  stk_classic::mesh::fixtures::BoxFixture fixture;
  stk_classic::mesh::BulkData          &bulk = fixture.bulk_data();
  stk_classic::mesh::PartVector         test_parts;
  stk_classic::mesh::PartVector         out_parts;
  stk_classic::mesh::PartVector         empty_vector;

  stk_classic::mesh::Entity &new_cell = bulk.declare_entity ( 3 , fixture.comm_rank()+1 , empty_vector );
  test_parts.push_back ( &fixture.fem_meta().universal_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
  test_parts.clear();
  test_parts.push_back ( &fixture.fem_meta().locally_owned_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
  test_parts.clear();
  test_parts.push_back ( &fixture.fem_meta().globally_shared_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
}
 */


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyDefaultPartAddition )
{
  TestBoxFixture fixture;
  stk_classic::mesh::BulkData            &bulk = fixture.bulk_data ();

  bulk.modification_begin();
  stk_classic::mesh::Entity &new_cell = fixture.get_new_entity ( 3 , 1 );
  bulk.modification_end();

  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.fem_meta().universal_part() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.fem_meta().locally_owned_part() ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyChangePartsSerial )
{
  TestBoxFixture fixture;
  stk_classic::mesh::BulkData            &bulk = fixture.bulk_data ();
  stk_classic::mesh::PartVector           create_parts , remove_parts , add_parts, empty_parts;

  create_parts.push_back ( &fixture.get_test_part() );
  create_parts.push_back ( &fixture.get_part_a_3() );
  remove_parts.push_back ( &fixture.get_part_a_3() );
  add_parts.push_back ( &fixture.get_part_b_superset() );
  add_parts.push_back ( &fixture.get_cell_part() );

  bulk.modification_begin();
  stk_classic::mesh::Entity &new_cell = fixture.get_new_entity ( 3 , 1 );
  bulk.change_entity_parts ( new_cell , create_parts , empty_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_test_part() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_a_3() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_a_superset() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_part_b_superset() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_cell_part() ) );

  bulk.modification_begin();
  bulk.change_entity_parts ( new_cell , add_parts , remove_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_test_part() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_part_a_3() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_a_superset() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_b_superset() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_cell_part() ) );

  bulk.modification_begin();
  bulk.change_entity_parts ( new_cell , empty_parts , add_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_test_part() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_part_a_3() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.get_part_a_superset() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_part_b_superset() ) );
  STKUNIT_ASSERT ( !new_cell.bucket().member ( fixture.get_cell_part() ) );

  //Verify still a member of default parts
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.fem_meta().universal_part() ) );
  STKUNIT_ASSERT ( new_cell.bucket().member ( fixture.fem_meta().locally_owned_part() ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyParallelAddParts )
{
  TestBoxFixture fixture;
  stk_classic::mesh::BulkData             &bulk = fixture.bulk_data ();
  stk_classic::mesh::PartVector            add_part;

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  add_part.push_back ( &fixture.get_part_a_0() );

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  STKUNIT_ASSERT(bulk.modification_end());

  bulk.modification_begin();

  for ( std::vector<stk_classic::mesh::Entity*>::const_iterator
        cur_entity =  bulk.entity_comm().begin();
        cur_entity != bulk.entity_comm().end() ; ++cur_entity ) {
    stk_classic::mesh::Entity & entity = **cur_entity ;
    if ( entity.entity_rank() == 0 ) {
      if ( entity.owner_rank() == fixture.comm_rank() ) {
        bulk.change_entity_parts ( entity, add_part, stk_classic::mesh::PartVector() );
      }
    }
  }

  bulk.modification_end();

  for ( std::vector<stk_classic::mesh::Entity*>::const_iterator
        cur_entity =  bulk.entity_comm().begin();
        cur_entity != bulk.entity_comm().end() ; ++cur_entity ) {
    stk_classic::mesh::Entity & entity = **cur_entity ;
    if ( entity.entity_rank() == 0 ) {
      STKUNIT_ASSERT ( entity.bucket().member ( fixture.get_part_a_0 () ) );
    }
  }
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyInducedMembership )
{
  TestBoxFixture fixture;
  stk_classic::mesh::BulkData             &bulk = fixture.bulk_data ();
  stk_classic::mesh::PartVector            create_node_parts , create_cell_parts , empty_parts;

  create_node_parts.push_back ( &fixture.get_part_a_0() );
  create_cell_parts.push_back ( &fixture.get_cell_part() );

  bulk.modification_begin();

  stk_classic::mesh::Entity &node = fixture.get_new_entity ( 0 , 1 );
  stk_classic::mesh::Entity &cell = fixture.get_new_entity ( 3 , 1 );

  bulk.modification_begin();

  bulk.change_entity_parts ( node , create_node_parts , stk_classic::mesh::PartVector () );
  bulk.change_entity_parts ( cell , create_cell_parts , stk_classic::mesh::PartVector () );
  // Add node to cell part
  stk_classic::mesh::RelationIdentifier cell_node_rel_id = 0;
  bulk.declare_relation ( cell , node , cell_node_rel_id );
  bulk.modification_end();

  STKUNIT_ASSERT ( node.bucket().member ( fixture.get_cell_part() ) );

  bulk.modification_begin();
  bulk.destroy_relation ( cell , node, cell_node_rel_id );
  bulk.modification_end();

  STKUNIT_ASSERT ( !node.bucket().member ( fixture.get_cell_part() ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCanRemoveFromSetWithDifferentRankSubset )
{
  TestBoxFixture fixture;
  stk_classic::mesh::BulkData           &bulk = fixture.bulk_data ();
  stk_classic::mesh::PartVector          add_parts , remove_parts, empty_parts;

  add_parts.push_back ( &fixture.get_part_b_3() );
  add_parts.push_back ( &fixture.get_part_a_superset() );

  remove_parts.push_back ( &fixture.get_part_a_superset() );

  bulk.modification_begin();

  stk_classic::mesh::Entity  &e = bulk.declare_entity ( 3 , fixture.comm_rank()+1 , add_parts );
  bulk.modification_end();

  bulk.modification_begin();
  bulk.change_entity_parts ( e , empty_parts , remove_parts );
  bulk.modification_end();

  STKUNIT_ASSERT ( e.bucket().member ( fixture.get_part_b_3() ) );
  STKUNIT_ASSERT ( !e.bucket().member ( fixture.get_part_a_superset() ) );
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCommonGhostingName )
{

  TestBoxFixture fixture;
  stk_classic::mesh::BulkData          &bulk = fixture.bulk_data ();

  bulk.modification_begin();

  if ( fixture.comm_size() == 1 ) return;

  if ( fixture.comm_rank() == 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.create_ghosting ( "Name 1" ) , std::runtime_error );
  }
  else
  {
    STKUNIT_ASSERT_THROW ( bulk.create_ghosting ( "Name 2" ) , std::runtime_error );
  }
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyTrivialDestroyAllGhostings )
{
  TestBoxFixture fixture;

  if ( fixture.comm_size() == 1 ) return;

  stk_classic::mesh::BulkData  &bulk = fixture.bulk_data();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  STKUNIT_ASSERT(bulk.modification_end());

  bulk.modification_begin();

  stk_classic::mesh::Ghosting &ghosting = bulk.create_ghosting ( "Ghost 1" );

  // Find a cell owned by this process
  std::vector<stk_classic::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(3).begin();
  unsigned send_rank = 0;

  std::vector<stk_classic::mesh::EntityProc>  to_send;
  std::vector<stk_classic::mesh::Entity *>    empty_vector;
  while ( cur_bucket != bulk.buckets(3).end() )
  {
    stk_classic::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() == fixture.comm_rank() )
      {
        if ( send_rank == fixture.comm_size() ) send_rank = 0;
        if ( send_rank != fixture.comm_rank() )
          to_send.push_back ( std::make_pair ( &*cur_entity , send_rank ) );
        send_rank++;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }
  bulk.change_ghosting ( ghosting , to_send , empty_vector );
  bulk.modification_end();


  {
    std::vector<stk_classic::mesh::EntityProc> send_list ;
    std::vector<stk_classic::mesh::Entity*>    recv_list ;
    ghosting.send_list( send_list );
    ghosting.receive_list( recv_list );

    STKUNIT_ASSERT ( ! send_list.empty()  );
    STKUNIT_ASSERT ( ! recv_list.empty() );
  }

  // Usage of operator << in Ghosting.cpp
  std::ostringstream oss;
  oss << ghosting;

  bulk.modification_begin();
  bulk.destroy_all_ghosting ();
  bulk.modification_end();

  {
    std::vector<stk_classic::mesh::EntityProc> send_list ;
    std::vector<stk_classic::mesh::Entity*>    recv_list ;
    ghosting.send_list( send_list );
    ghosting.receive_list( recv_list );

    STKUNIT_ASSERT ( send_list.empty() );
    STKUNIT_ASSERT ( recv_list.empty() );
  }
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyChangeGhostingGuards )
{
  TestBoxFixture fixture1, fixture2;
  stk_classic::mesh::BulkData & bulk1 = fixture1.bulk_data ();
  stk_classic::mesh::BulkData & bulk2 = fixture2.bulk_data ();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box1[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };
  int local_box2[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk1.modification_begin();
  fixture1.generate_boxes( root_box, local_box1 );
  STKUNIT_ASSERT(bulk1.modification_end());

  bulk2.modification_begin();
  fixture2.generate_boxes( root_box, local_box2 );
  STKUNIT_ASSERT(bulk2.modification_end());

  bulk1.modification_begin();
  bulk2.modification_begin();

  std::vector<stk_classic::mesh::EntityProc>  to_send;
  std::vector<stk_classic::mesh::Entity *>    empty_vector;
  std::vector<stk_classic::mesh::Bucket *>::const_iterator cur_bucket = bulk1.buckets(3).begin();
  unsigned send_rank = 0;
  while ( cur_bucket != bulk1.buckets(3).end() )
  {
    stk_classic::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() == fixture1.comm_rank() )
      {
        if ( send_rank == fixture1.comm_size() ) send_rank = 0;
        if ( send_rank != fixture1.comm_rank() )
          to_send.push_back ( std::make_pair ( &*cur_entity , send_rank ) );
        ++send_rank;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  stk_classic::mesh::Ghosting &ghosting = bulk1.create_ghosting ( "Ghost 1" );
  STKUNIT_ASSERT_THROW ( bulk2.change_ghosting ( ghosting , to_send , empty_vector ) , std::runtime_error );
  STKUNIT_ASSERT_THROW ( bulk1.change_ghosting ( bulk1.shared_aura() , to_send , empty_vector ) , std::runtime_error );

  ghosting.receive_list(empty_vector);
  ghosting.send_list(to_send);

  bulk1.modification_end();
  bulk2.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyOtherGhostingGuards )
{
  TestBoxFixture fixture;
  stk_classic::mesh::BulkData          &bulk = fixture.bulk_data ();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  STKUNIT_ASSERT(bulk.modification_end());

  bulk.modification_begin();

  std::vector<stk_classic::mesh::EntityProc>  to_send_unowned;
  std::vector<stk_classic::mesh::EntityProc>  empty_send;
  std::vector<stk_classic::mesh::Entity *>    to_remove_not_ghosted;
  std::vector<stk_classic::mesh::Entity *>    empty_remove;
  std::vector<stk_classic::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(3).begin();
  unsigned send_rank = 0;
  while ( cur_bucket != bulk.buckets(3).end() )
  {
    stk_classic::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() != fixture.comm_rank() )
      {
        if ( send_rank == fixture.comm_size() ) send_rank = 0;
        if ( send_rank != fixture.comm_rank() )
          to_send_unowned.push_back ( std::make_pair ( &*cur_entity , send_rank ) );
        ++send_rank;
      }
      else
      {
        to_remove_not_ghosted.push_back ( &*cur_entity );
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  stk_classic::mesh::Ghosting &ghosting = bulk.create_ghosting ( "Ghost 1" );
  if ( to_send_unowned.size() > 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.change_ghosting ( ghosting , to_send_unowned , empty_remove ) , std::runtime_error );
  }
  else
  {
    bulk.change_ghosting ( ghosting , to_send_unowned , empty_remove );
  }

  if ( to_remove_not_ghosted.size() > 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.change_ghosting ( ghosting , empty_send , to_remove_not_ghosted ) , std::runtime_error );
  }
  else
  {
    bulk.change_ghosting ( ghosting , empty_send , to_remove_not_ghosted );
  }
  bulk.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyPartsOnCreate )
{
   TestBoxFixture fixture;
   stk_classic::mesh::BulkData           & bulk = fixture.bulk_data ();
   stk_classic::mesh::Part               & part_a = fixture.get_part_a_0 ();
   stk_classic::mesh::Part               & part_b = fixture.get_part_b_0 ();

   stk_classic::mesh::PartVector           create_vector;
   create_vector.push_back ( &part_a );

   bulk.modification_begin();

   stk_classic::mesh::Entity &node = bulk.declare_entity ( 0 , fixture.comm_rank()+1 ,create_vector );
   bulk.modification_end();

   STKUNIT_ASSERT ( node.bucket().member ( part_a ) );

   bulk.modification_begin();
   create_vector.push_back ( &part_b );
   stk_classic::mesh::Entity &node2 = bulk.declare_entity ( 0 , fixture.comm_size() + fixture.comm_rank() + 1 , create_vector );
   bulk.modification_end();

   STKUNIT_ASSERT ( node2.bucket().member ( part_a ) );
   STKUNIT_ASSERT ( node2.bucket().member ( part_b ) );
}

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyBoxGhosting )
{
  const unsigned p_size = stk_classic::parallel_machine_size( MPI_COMM_WORLD );
  if ( 8 < p_size ) { return ; }

  stk_classic::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, 2, 2, 2 );
  fixture.m_fem_meta.commit();
  fixture.generate_mesh();

  for ( size_t iz = 0 ; iz < 3 ; ++iz ) {
  for ( size_t iy = 0 ; iy < 3 ; ++iy ) {
  for ( size_t ix = 0 ; ix < 3 ; ++ix ) {
    stk_classic::mesh::Entity * const node = fixture.node(ix,iy,iz);
    STKUNIT_ASSERT( NULL != node );

    STKUNIT_ASSERT( fixture.node_id(ix,iy,iz) == node->identifier() );
    stk_classic::mesh::fixtures::HexFixture::Scalar * const node_coord =
      stk_classic::mesh::field_data( fixture.m_coord_field , *node );
    STKUNIT_ASSERT( node_coord != NULL );
  }
  }
  }

  for ( size_t iz = 0 ; iz < 2 ; ++iz ) {
  for ( size_t iy = 0 ; iy < 2 ; ++iy ) {
  for ( size_t ix = 0 ; ix < 2 ; ++ix ) {
    stk_classic::mesh::Entity * const elem = fixture.elem(ix,iy,iz);
    STKUNIT_ASSERT( NULL != elem );

    stk_classic::mesh::PairIterRelation elem_nodes = elem->relations();
    STKUNIT_ASSERT_EQUAL( 8u , elem_nodes.size() );
    stk_classic::mesh::fixtures::HexFixture::Scalar ** const elem_node_coord =
      stk_classic::mesh::field_data( fixture.m_coord_gather_field , *elem );
    for ( size_t j = 0 ; j < elem_nodes.size() ; ++j ) {
      STKUNIT_ASSERT_EQUAL( j , elem_nodes[j].identifier() );
      stk_classic::mesh::fixtures::HexFixture::Scalar * const node_coord =
        stk_classic::mesh::field_data( fixture.m_coord_field , *elem_nodes[j].entity() );
      STKUNIT_ASSERT( node_coord == elem_node_coord[ elem_nodes[j].identifier() ] );
    }
    if ( 8u == elem_nodes.size() ) {
      STKUNIT_ASSERT( elem_nodes[0].entity() == fixture.node(ix,iy,iz));
      STKUNIT_ASSERT( elem_nodes[1].entity() == fixture.node(ix+1,iy,iz));
      STKUNIT_ASSERT( elem_nodes[2].entity() == fixture.node(ix+1,iy,iz+1));
      STKUNIT_ASSERT( elem_nodes[3].entity() == fixture.node(ix,iy,iz+1));
      STKUNIT_ASSERT( elem_nodes[4].entity() == fixture.node(ix,iy+1,iz));
      STKUNIT_ASSERT( elem_nodes[5].entity() == fixture.node(ix+1,iy+1,iz));
      STKUNIT_ASSERT( elem_nodes[6].entity() == fixture.node(ix+1,iy+1,iz+1));
      STKUNIT_ASSERT( elem_nodes[7].entity() == fixture.node(ix,iy+1,iz+1));
    }

  }
  }
  }
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , testEntityComm )
{
  //Test on unpack_field_values in EntityComm.cpp
  //code based on ../base/BulkDataGhosting.cpp
  //Create a simple mesh. Add nodes one element and some parts.

  const int spatial_dimension = 3;

  stk_classic::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(spatial_dimension, stk_classic::mesh::fem::entity_rank_names ( spatial_dimension ));

  stk_classic::mesh::fem::CellTopology tet_top(shards::getCellTopologyData<shards::Tetrahedron<4> >());
  stk_classic::mesh::Part & part_a = fem_meta.declare_part( "block_a", tet_top );
  stk_classic::mesh::Part & part_b = fem_meta.declare_part( "block_b", tet_top );

  stk_classic::mesh::fem::CellTopology node_top(shards::getCellTopologyData<shards::Node>());
  stk_classic::mesh::Part & part_a_0 = fem_meta.declare_part( "block_a_0", node_top );

  typedef stk_classic::mesh::Field<double>  ScalarFieldType;

  ScalarFieldType & volume =
     fem_meta.declare_field < ScalarFieldType > ( "volume" , 4 );
  ScalarFieldType & temperature =
     fem_meta.declare_field < ScalarFieldType > ( "temperature" , 4 );
  stk_classic::mesh::Part  & universal     = fem_meta.universal_part ();
  put_field ( volume , 3 , universal );
  put_field ( temperature , 3 , universal );

  fem_meta.commit();

  stk_classic::mesh::PartVector    create_vector;
  stk_classic::mesh::PartVector    empty_vector;
  create_vector.push_back ( &part_a );
  create_vector.push_back ( &part_b );

  stk_classic::mesh::MetaData & meta = stk_classic::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  stk_classic::mesh::BulkData bulk ( meta , MPI_COMM_WORLD , 100 );

  bulk.modification_begin();

  stk_classic::mesh::Ghosting &ghosts = bulk.create_ghosting ( "Ghost 1" );

  unsigned size2 = stk_classic::parallel_machine_size( MPI_COMM_WORLD );
  unsigned rank_count2 = stk_classic::parallel_machine_rank( MPI_COMM_WORLD );
  int new_id2 = size2 + rank_count2;

  stk_classic::mesh::Entity &elem2 = bulk.declare_entity ( 3 , new_id2+1 ,create_vector );
  STKUNIT_ASSERT_EQUAL( elem2.bucket().member ( part_a ), true );

  unsigned size = stk_classic::parallel_machine_size( MPI_COMM_WORLD );
  unsigned rank_count = stk_classic::parallel_machine_rank( MPI_COMM_WORLD );

  int id_base = 0;
  for ( id_base = 0 ; id_base < 99 ; ++id_base )
  {
    int new_id = size * id_base + rank_count;
    stk_classic::mesh::Entity &new_node = bulk.declare_entity( 0 , new_id+1 , empty_vector );
    STKUNIT_ASSERT_EQUAL( new_node.bucket().member ( part_a_0 ), false );
  }

  //Create a bucket of nodes for sending

  std::vector<stk_classic::mesh::EntityProc>  add_send;

  const std::vector<stk_classic::mesh::Bucket*> & buckets = bulk.buckets( 0 );

  std::vector<stk_classic::mesh::Bucket*>::const_iterator cur_bucket;

  cur_bucket = buckets.begin();

  unsigned send_rank = 0;
  while ( cur_bucket != buckets.end() )
  {
    stk_classic::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( cur_entity->owner_rank() == rank_count )
      {
        if ( send_rank == size ) send_rank = 0;
        if ( send_rank != rank_count )
          add_send.push_back ( std::make_pair ( &*cur_entity , send_rank ) );
        ++send_rank;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess > new_send ;
  std::set< stk_classic::mesh::Entity * ,   stk_classic::mesh::EntityLess > new_recv ;

  //  Keep the closure of the remaining received ghosts.
  //  Working from highest-to-lowest key (rank entity type)
  //  results in insertion of the transitive closure.
  //  Insertion will not invalidate the associative container's iterator.

  for ( std::set< stk_classic::mesh::Entity * , stk_classic::mesh::EntityLess >::iterator
        i = new_recv.end() ; i != new_recv.begin() ; ) {
    --i ;

    const unsigned erank = (*i)->entity_rank();

    for ( stk_classic::mesh::PairIterRelation
          irel = (*i)->relations(); ! irel.empty() ; ++irel ) {
      if ( irel->entity_rank() < erank &&
           in_receive_ghost( ghosts , * irel->entity() ) ) {
        new_recv.insert( irel->entity() );
      }
    }
  }

  //  Initialize the new_send from the new_recv
  new_comm_recv_to_send( bulk , new_recv , new_send );

  //------------------------------------
  // Add the specified entities and their closure to the send ghosting

  for ( std::vector< stk_classic::mesh::EntityProc >::const_iterator
        i = add_send.begin() ; i != add_send.end() ; ++i ) {
        new_insert_transitive_closure( new_send , *i );
  }

  // Synchronize the send and receive list.
  // If the send list contains a not-owned entity
  // inform the owner and receiver to ad that entity
  // to their ghost send and receive lists.

  new_comm_sync_send_recv( bulk , new_send , new_recv );

  //------------------------------------
  // Push newly ghosted entities to the receivers and update the comm list.
  // Unpacking must proceed in entity-rank order so that higher ranking
  // entities that have relations to lower ranking entities will have
  // the lower ranking entities unpacked first.  The higher and lower
  // ranking entities may be owned by different processes,
  // as such unpacking must be performed in rank order.

  //Start of CommAll section:
  {
    stk_classic::CommAll comm( MPI_COMM_WORLD );

    for ( std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess >::iterator
          j = new_send.begin(); j != new_send.end() ; ++j ) {
          stk_classic::mesh::Entity & entity = * j->first ;
      if ( ! in_ghost( ghosts , entity , j->second ) ) {
        // Not already being sent , must send it.
        stk_classic::CommBuffer & buf = comm.send_buffer( j->second );
        buf.pack<unsigned>( entity.entity_rank() );
        stk_classic::mesh::pack_entity_info(  buf , entity );
        stk_classic::mesh::pack_field_values( buf , entity );
      }
    }

    comm.allocate_buffers( size / 4 );

    for ( std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess >::iterator
          j = new_send.begin(); j != new_send.end() ; ++j ) {
          stk_classic::mesh::Entity & entity = * j->first ;
      if ( ! in_ghost( ghosts , entity , j->second ) ) {
        // Not already being sent , must send it.
        stk_classic::CommBuffer & buf = comm.send_buffer( j->second );
        buf.pack<unsigned>( entity.entity_rank() );
        stk_classic::mesh::pack_entity_info(  buf , entity );
        stk_classic::mesh::pack_field_values( buf , entity );

      }
    }

    comm.communicate();

    std::ostringstream error_msg ;

    for ( unsigned rank = 0 ; rank < rank_count ; ++rank ) {

      for ( unsigned p = 0 ; p < size ; ++p ) {

        stk_classic::CommBuffer & buf = comm.recv_buffer(p);

        while ( buf.remaining() ) {

          // Only unpack if of the current entity rank.
          // If not the current entity rank, break the iteration
          // until a subsequent entity rank iteration.
          {
            unsigned this_rank = ~0u ;
            buf.peek<unsigned>( this_rank );
            if ( this_rank != rank ) break ;

            buf.unpack<unsigned>( this_rank );
          }

          // FIXME for Carol; the code below did not work with -np 4
          //STKUNIT_ASSERT_EQUAL( stk_classic::mesh::unpack_field_values( buf , elem2 , error_msg ), false);
	  //std::cout << "Error message for unpack_field_values = " << error_msg.str() << std::endl ;

        }
      }

    }
  }//end of CommAll section

  bulk.modification_end ();
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , testUninitializedMetaData )
{
  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;

  stk_classic::mesh::MetaData meta; // Construct, but do not initialize
  stk_classic::mesh::BulkData bulk(meta, pm);

  meta.set_entity_rank_names(stk_classic::mesh::fem::entity_rank_names(2 /*spatial-dim*/));

  meta.commit();

  bulk.modification_begin();

  STKUNIT_ASSERT_THROW( bulk.declare_entity(0, /*rank*/
                                            1, /*id*/
                                            stk_classic::mesh::PartVector() ),
                        std::logic_error);
}

namespace {

void new_insert_transitive_closure( std::set<stk_classic::mesh::EntityProc,stk_classic::mesh::EntityLess> & new_send ,
                                const stk_classic::mesh::EntityProc & entry )
{
  // Do not insert if I can determine that this entity is already
  // owned or shared by the receiving processor.

  if ( entry.second != entry.first->owner_rank() &&
       ! in_shared( * entry.first , entry.second ) ) {

    std::pair< std::set<stk_classic::mesh::EntityProc,stk_classic::mesh::EntityLess>::iterator , bool >
      result = new_send.insert( entry );

    if ( result.second ) {
      // A new insertion, must also insert the closure

      const unsigned etype = entry.first->entity_rank();
      stk_classic::mesh::PairIterRelation irel  = entry.first->relations();

      for ( ; ! irel.empty() ; ++irel ) {
        if ( irel->entity_rank() < etype ) {
          stk_classic::mesh::EntityProc tmp( irel->entity() , entry.second );
          new_insert_transitive_closure( new_send , tmp );
        }
      }
    }
  }
}


// Synchronize the send list to the receive list.

void new_comm_sync_send_recv(
  stk_classic::mesh::BulkData & mesh ,
  std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess > & new_send ,
  std::set< stk_classic::mesh::Entity * , stk_classic::mesh::EntityLess > & new_recv )
{
  const unsigned parallel_rank = mesh.parallel_rank();
  const unsigned parallel_size = mesh.parallel_size();

  stk_classic::CommAll all( mesh.parallel() );

  // Communication sizing:

  for ( std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess >::iterator
        i = new_send.begin() ; i != new_send.end() ; ++i ) {
    const unsigned owner = i->first->owner_rank();
    all.send_buffer( i->second ).skip<stk_classic::mesh::EntityKey>(2);
    if ( owner != parallel_rank ) {
      all.send_buffer( owner ).skip<stk_classic::mesh::EntityKey>(2);
    }
  }

  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  // Communication packing (with message content comments):
  for ( std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess >::iterator
        i = new_send.begin() ; i != new_send.end() ; ) {
    const unsigned owner = i->first->owner_rank();

    // Inform receiver of ghosting, the receiver does not own
    // and does not share this entity.
    // The ghost either already exists or is a to-be-done new ghost.
    // This status will be resolved on the final communication pass
    // when new ghosts are packed and sent.

    const stk_classic::mesh::EntityKey &entity_key = i->first->key();
    const uint64_t &proc = i->second;

    all.send_buffer( i->second ).pack(entity_key).pack(proc);

    if ( owner != parallel_rank ) {
      // I am not the owner of this entity.
      // Inform the owner of this ghosting need.
      all.send_buffer( owner ).pack(entity_key).pack(proc);

      // Erase it from my processor's ghosting responsibility:
      // The iterator passed to the erase method will be invalidated.
      std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess >::iterator jrem = i ; ++i ;
      new_send.erase( jrem );
    }
    else {
      ++i ;
    }
  }

  all.communicate();

  // Communication unpacking:
  for ( unsigned p = 0 ; p < parallel_size ; ++p ) {
    stk_classic::CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {

      stk_classic::mesh::EntityKey entity_key;
      uint64_t proc(0);

      buf.unpack(entity_key).unpack(proc);

      stk_classic::mesh::Entity * const e = mesh.get_entity( entity_key );

      if ( parallel_rank != proc ) {
        //  Receiving a ghosting need for an entity I own.
        //  Add it to my send list.
        STKUNIT_ASSERT( e != NULL );
        stk_classic::mesh::EntityProc tmp( e , proc );
        new_send.insert( tmp );
      }
      else if ( e != NULL ) {
        //  I am the receiver for this ghost.
        //  If I already have it add it to the receive list,
        //  otherwise don't worry about it - I will receive
        //  it in the final new-ghosting communication.
        new_recv.insert( e );
      }
    }
  }
}

void new_comm_recv_to_send(
  stk_classic::mesh::BulkData & mesh ,
  const std::set< stk_classic::mesh::Entity * , stk_classic::mesh::EntityLess > & new_recv ,
        std::set< stk_classic::mesh::EntityProc , stk_classic::mesh::EntityLess > & new_send )
{
  const unsigned parallel_size = mesh.parallel_size();

  stk_classic::CommAll all( mesh.parallel() );

  for ( std::set< stk_classic::mesh::Entity * , stk_classic::mesh::EntityLess >::const_iterator
        i = new_recv.begin() ; i != new_recv.end() ; ++i ) {
    const unsigned owner = (*i)->owner_rank();
    all.send_buffer( owner ).skip<stk_classic::mesh::EntityKey>(1);
  }

  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  for ( std::set< stk_classic::mesh::Entity * , stk_classic::mesh::EntityLess >::const_iterator
        i = new_recv.begin() ; i != new_recv.end() ; ++i ) {
    const unsigned owner = (*i)->owner_rank();
    const stk_classic::mesh::EntityKey key = (*i)->key();
    all.send_buffer( owner ).pack<stk_classic::mesh::EntityKey>( & key , 1 );
  }

  all.communicate();

  for ( unsigned p = 0 ; p < parallel_size ; ++p ) {
    stk_classic::CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {
      stk_classic::mesh::EntityKey key ;
      buf.unpack<stk_classic::mesh::EntityKey>( & key , 1 );
      stk_classic::mesh::EntityProc tmp( mesh.get_entity( entity_rank(key), entity_id(key) ) , p );
      new_send.insert( tmp );
    }
  }
}

}
