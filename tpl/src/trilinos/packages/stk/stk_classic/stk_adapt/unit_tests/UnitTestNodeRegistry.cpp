/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_percept/Util.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/ConstantFunction.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/fixtures/Fixture.hpp>

#include <stk_adapt/NodeRegistry.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/SerializeNodeRegistry.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <stk_util/parallel/Parallel.hpp>
#include <math.h>

namespace stk_classic {
namespace adapt {
namespace unit_tests {

#define EXTRA_PRINT 0

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(nodeRegistry, createAddNodes_serial_and_1st_parallel)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  percept::PerceptMesh eMesh(3u);
  eMesh.new_mesh(percept::GMeshSpec("3x3x12|bbox:0,0,0,1,1,1"));  // create a 3x3x12 hex mesh in the unit cube
  //int scalarDimension = 0; // a scalar
  //int vectorDimension = 3;
  eMesh.commit();
  eMesh.print_info();

  unsigned p_size = eMesh.get_bulk_data()->parallel_size();
  unsigned p_rank = eMesh.get_bulk_data()->parallel_rank();

  std::cout << "TEST::nodeRegistry::createAddNodes_serial: p_size = "<< p_size << " rank= " << p_rank << std::endl;
  if (p_size >= 1)  // FIXME
    return;

  if (p_size >= 3)
    return;

  //stk_classic::CommAll comm_all(eMesh.get_bulk_data()->parallel());
  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size == 2) return; // FIXME
  if (p_size == 2)
  {
    unsigned elem_num=(12/p_size)*3*3;

    const stk_classic::mesh::Entity* element_1_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num);
    dw() << "P["<<p_rank<<"] element_1_p = " << element_1_p << DWENDL;
    stk_classic::mesh::EntityRank stk_mesh_Edge = 1;

    NeededEntityType entity_rank(stk_mesh_Edge, 1u);
    unsigned iSubDimOrd = 0u;

    const stk_classic::mesh::Entity& element_1 = *eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num);

    //nodeRegistry.getSubDimEntity(subDimEntity, element_1, entity_rank, iSubDimOrd);
    nodeRegistry.registerNeedNewNode(element_1, entity_rank, iSubDimOrd, true);

    // must be a ghost on one processor
    int n_check = 0;
    if (eMesh.isGhostElement(element_1))
    {
      std::cout  << "P["<<p_rank<<"] element no= " << elem_num
                 << " TEST::nodeRegistry::createAddNodes_serial: is ghost on rank= " << p_rank << std::endl;
      STKUNIT_EXPECT_EQ(nodeRegistry.local_size(), 1u);
      ++n_check;
    }
    else
    {
      std::cout  << "P["<<p_rank<<"] element no= " << elem_num
                 << " TEST::nodeRegistry::createAddNodes_serial: is not ghost on rank= " << p_rank << std::endl;
      STKUNIT_EXPECT_EQ(nodeRegistry.local_size(), 1u);
      ++n_check;
    }
    STKUNIT_EXPECT_EQ(n_check, 1);
    //Util::pause(true, "tmp");
    //exit(1);
  }
  else if (p_size == 1)
  {
    unsigned elem_num=1;

    const stk_classic::mesh::Entity* element_1_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num);
    std::cout << "P["<<p_rank<<"] element_1_p = " << element_1_p << std::endl;
    stk_classic::mesh::EntityRank stk_mesh_Edge = 1;

    NeededEntityType entity_rank(stk_mesh_Edge, 1u);
    unsigned iSubDimOrd = 0u;

    const stk_classic::mesh::Entity& element_1 = *eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num);

    nodeRegistry.registerNeedNewNode(element_1, entity_rank, iSubDimOrd, true);

    // must be a ghost on one processor
    if (eMesh.isGhostElement(element_1))
    {
      std::cout << "TEST::nodeRegistry::createAddNodes_serial: is ghost on rank= " << p_rank << std::endl;
      STKUNIT_EXPECT_EQ(nodeRegistry.local_size(), 0u);
    }
    else
    {
      std::cout << "TEST::nodeRegistry::createAddNodes_serial: is not ghost on rank= " << p_rank << std::endl;
      STKUNIT_EXPECT_EQ(nodeRegistry.local_size(), 1u);
    }
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(nodeRegistry, test_parallel_0)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  /* output snippets from this test - note that edge 4 of elem 54 is same as edge 0 of elem 63 (used in next test)

     num_elements_in_bucket = 9 element ids =
     46 47 48 49 50 51 52 53 54
     num_elements_in_bucket = 54 element ids =
     55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74
     75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94
     95 96 97 98 99 100 101 102 103 104 105 106 107 108

     P[1] element_local = Elem: 63 nodes: 107 108 112 111 123 124 128 127
     P[1] element_ghost = Elem: 54 nodes: 91 92 96 95 107 108 112 111
     P[1] element_local edge 0 = 107 108
     P[1] element_local edge 1 = 108 112
     P[1] element_local edge 2 = 111 112
     P[1] element_local edge 3 = 107 111
     P[1] element_local edge 4 = 123 124
     P[1] element_local edge 5 = 124 128
     P[1] element_local edge 6 = 127 128
     P[1] element_local edge 7 = 123 127
     P[1] element_local edge 8 = 107 123
     P[1] element_local edge 9 = 108 124
     P[1] element_local edge 10 = 112 128
     P[1] element_local edge 11 = 111 127

     P[0] element_local = Elem: 54 nodes: 91 92 96 95 107 108 112 111
     P[0] element_ghost = Elem: 63 nodes: 107 108 112 111 123 124 128 127
     P[0] element_local edge 0 = 91 92
     P[0] element_local edge 1 = 92 96
     P[0] element_local edge 2 = 95 96
     P[0] element_local edge 3 = 91 95
     P[0] element_local edge 4 = 107 108
     P[0] element_local edge 5 = 108 112
     P[0] element_local edge 6 = 111 112
     P[0] element_local edge 7 = 107 111
     P[0] element_local edge 8 = 91 107
     P[0] element_local edge 9 = 92 108
     P[0] element_local edge 10 = 96 112
     P[0] element_local edge 11 = 95 111

  */

  // start_demo_nodeRegistry_test_parallel_0
  percept::PerceptMesh eMesh(3u);
  eMesh.new_mesh(percept::GMeshSpec("3x3x12|bbox:0,0,0,1,1,1"));  // create a 3x3x12 hex mesh in the unit cube
  eMesh.commit();
  eMesh.print_info();

  unsigned p_size = eMesh.get_parallel_size();
  unsigned p_rank = eMesh.get_rank();
  Util::setRank(eMesh.get_rank());

  if (p_size != 2) // FIXME
    return;

  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size == 2)
  {
    // pick an element on the processor boundary
    unsigned elem_num_local = (12/p_size)*3*3;
    unsigned elem_num_ghost = elem_num_local+(3*3);

    const stk_classic::mesh::Entity* element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
    const stk_classic::mesh::Entity* element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);

    if (p_rank == 1)
    {
      element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
      element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
    }

    dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
    dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

    const stk_classic::mesh::Entity& element_local = *element_local_p;
    const stk_classic::mesh::Entity& element_ghost = *element_ghost_p;

    std::cout << "P["<<p_rank<<"] element_local = " << element_local << std::endl;
    std::cout << "P["<<p_rank<<"] element_ghost = " << element_ghost << std::endl;

    stk_classic::mesh::EntityRank stk_mesh_Edge = 1;

    stk_classic::mesh::EntityRank needed_entity_rank = stk_mesh_Edge;
    for (unsigned iSubDimOrd = 0; iSubDimOrd < 12; iSubDimOrd++)
    {
      SubDimCell_SDSEntityType subDimEntity;
      nodeRegistry.getSubDimEntity(subDimEntity, element_local, needed_entity_rank, iSubDimOrd);
      std::cout << "P[" << p_rank << "] element_local edge " << iSubDimOrd << " = " << subDimEntity << std::endl;
    }

  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(nodeRegistry, test_parallel_1)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  // start_demo_nodeRegistry_test_parallel_1

  percept::PerceptMesh eMesh(3u);
  eMesh.new_mesh(percept::GMeshSpec("3x3x12|bbox:0,0,0,1,1,1"));  // create a 3x3x12 hex mesh in the unit cube
  eMesh.commit();
  eMesh.print_info();
  eMesh.save_as("./cube_hex9-orig.e");

  unsigned p_size = eMesh.get_parallel_size();
  unsigned p_rank = eMesh.get_rank();
  Util::setRank(eMesh.get_rank());

  if (p_size != 2) // FIXME
    return;

  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size == 2)
  {
    // pick an element on the processor boundary
    unsigned elem_num_local = (12/p_size)*3*3;
    unsigned elem_num_ghost = elem_num_local+(3*3);

    stk_classic::mesh::Entity* element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
    stk_classic::mesh::Entity* element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
    if (p_rank == 1)
    {
      element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
      element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
    }

    dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
    dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

    stk_classic::mesh::Entity& element_local = *element_local_p;
    stk_classic::mesh::Entity& element_ghost = *element_ghost_p;

    std::cout << "P["<<p_rank<<"] element_local = " << element_local << std::endl;
    std::cout << "P["<<p_rank<<"] element_ghost = " << element_ghost << std::endl;

    // choose edges to be used for new node locations (i.e., this would model a serendipity-like element with only edge Lagrange nodes)

    stk_classic::mesh::EntityRank stk_mesh_Edge = 1;
    NeededEntityType needed_entity_rank( stk_mesh_Edge, 1u);
    std::vector<NeededEntityType> needed_entity_ranks(1, needed_entity_rank);

    /*
     * 1st of three steps to create and associate new nodes - register need for new nodes, then check if node is remote, then get
     *   from remote proc if necessary; finally, the local node database is ready to be queried
     *
     * The pattern is to begin the step, loop over all elements (including ghosts) and invoke the local operation
     * The method doForAllSubEntities is a utility for performing the operation on all the sub entities.
     * If more granularity is desired, the member functions can be invoked directly for a particular sub-entity.
     */
    nodeRegistry.beginRegistration();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_ghost, needed_entity_ranks);
    nodeRegistry.endRegistration();

    std::cout << "P["<<p_rank<<"] nodeRegistry size  = " << nodeRegistry.total_size() << std::endl;
    std::cout << "P["<<p_rank<<"] nodeRegistry lsize = " << nodeRegistry.local_size() << std::endl;

    dw() << "P["<<p_rank<<"] nodeRegistry size       = " << nodeRegistry.total_size() << DWENDL;
    dw() << "P["<<p_rank<<"] nodeRegistry lsize      = " << nodeRegistry.local_size() << DWENDL;

    // could do local create of elements here
    nodeRegistry.beginLocalMeshMods();
    nodeRegistry.endLocalMeshMods();

    // check if the newly requested nodes are local or remote
    nodeRegistry.beginCheckForRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_ghost, needed_entity_ranks);
    nodeRegistry.endCheckForRemote();

    // get the new nodes from other procs if they are nonlocal
    nodeRegistry.beginGetFromRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_ghost, needed_entity_ranks);
    nodeRegistry.endGetFromRemote();

    // now we can get the new node's id and entity
    unsigned iSubDimOrd = 4u;
    if (p_rank)
    {
      iSubDimOrd = 0u;
    }
    NodeIdsOnSubDimEntityType& nodeIds_onSE = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_rank.first, iSubDimOrd));
    stk_classic::mesh::Entity*  node   = eMesh.get_bulk_data()->get_entity(stk_classic::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE[0]->identifier());

    // should be the same node on each proc
    std::cout << "P[" << p_rank << "] nodeId = " << nodeIds_onSE << " node= " << node << std::endl;

    // end_demo

    // change element to be a strange element with 9 nodes
    eMesh.get_bulk_data()->modification_begin();
    eMesh.get_bulk_data()->declare_relation(element_local, *node, 8);
    eMesh.get_bulk_data()->modification_end();

    eMesh.save_as("./cube_hex9.e");
    //exit(1);
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(nodeRegistry, test_parallel_1_0)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  // start_demo_nodeRegistry_test_parallel_1

  percept::PerceptMesh eMesh(3u);

  unsigned p_size = eMesh.get_parallel_size();
  unsigned p_rank = eMesh.get_rank();
  Util::setRank(eMesh.get_rank());

  eMesh.new_mesh(percept::GMeshSpec(std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,1")));

  // prepare for adding some quadratic elements
  mesh::Part& block_hex_20 = eMesh.get_fem_meta_data()->declare_part("block_hex_20", eMesh.element_rank());
  /// set cell topology for the part block_hex_20
  mesh::fem::set_cell_topology< shards::Hexahedron<20>  >( block_hex_20 );
  stk_classic::io::put_io_part_attribute(block_hex_20);

  eMesh.commit();
  eMesh.print_info();
  eMesh.save_as("./cube1x1x2_hex-20-orig.e");

  mesh::Part* block_hex_8 = const_cast<mesh::Part *>(eMesh.getPart("block_1"));

  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size <= 2)
  {
    // pick an element on the processor boundary
    unsigned elem_num_local = 1;
    unsigned elem_num_ghost = 2;
    if (p_size == 1)
      elem_num_ghost = 1;

    stk_classic::mesh::Entity* element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
    stk_classic::mesh::Entity* element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
    if (p_rank == 1)
    {
      element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
      element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
    }

    dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
    dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

    stk_classic::mesh::Entity& element_local = *element_local_p;
    stk_classic::mesh::Entity& element_ghost = *element_ghost_p;

    std::cout << "P["<<p_rank<<"] element_local = " << element_local << std::endl;
    std::cout << "P["<<p_rank<<"] element_ghost = " << element_ghost << std::endl;

    // choose edges to be used for new node locations (i.e., this would model a serendipity-like element with only edge Lagrange nodes)
    stk_classic::mesh::EntityRank stk_mesh_Edge = 1;
    NeededEntityType needed_entity_rank( stk_mesh_Edge, 1u);
    std::vector<NeededEntityType> needed_entity_ranks(1, needed_entity_rank);

    /*
     * 1st of three steps to create and associate new nodes - register need for new nodes, then check if node is remote, then get
     *   from remote proc if necessary; finally, the local node database is ready to be queried
     *
     * The pattern is to begin the step, loop over all elements (including ghosts) and invoke the local operation
     * The method doForAllSubEntities is a utility for performing the operation on all the sub entities.
     * If more granularity is desired, the member functions can be invoked directly for a particular sub-entity.
     */
    nodeRegistry.beginRegistration();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_ghost, needed_entity_ranks);
    nodeRegistry.endRegistration();

    std::cout << "P["<<p_rank<<"] nodeRegistry size  = " << nodeRegistry.total_size() << std::endl;
    std::cout << "P["<<p_rank<<"] nodeRegistry lsize = " << nodeRegistry.local_size() << std::endl;

    dw() << "P["<<p_rank<<"] nodeRegistry size       = " << nodeRegistry.total_size() << DWENDL;
    dw() << "P["<<p_rank<<"] nodeRegistry lsize      = " << nodeRegistry.local_size() << DWENDL;

    // could do local create of elements here
    nodeRegistry.beginLocalMeshMods();
    nodeRegistry.endLocalMeshMods();

    // check if the newly requested nodes are local or remote
    nodeRegistry.beginCheckForRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_ghost, needed_entity_ranks);
    nodeRegistry.endCheckForRemote();

    // get the new nodes from other procs if they are nonlocal
    nodeRegistry.beginGetFromRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_ghost, needed_entity_ranks);
    nodeRegistry.endGetFromRemote();

    // now we can get the new node's id and entity
    unsigned iSubDimOrd = 4u;
    if (p_rank)
    {
      iSubDimOrd = 0u;
    }
    NodeIdsOnSubDimEntityType& nodeIds_onSE_0 = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_rank.first, iSubDimOrd));
    stk_classic::mesh::Entity*  node_0   = eMesh.get_bulk_data()->get_entity(stk_classic::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE_0[0]->identifier());

    // should be the same node on each proc
    std::cout << "P[" << p_rank << "] nodeId_0 = " << nodeIds_onSE_0 << " node_0= " << node_0 << std::endl;

    // end_demo

#if STK_ADAPT_HAVE_YAML_CPP
    if (p_size == 1)
      {
        if (1) {
          YAML::Emitter out;
          out << YAML::Anchor("NodeRegistry::map");
          out << YAML::BeginMap;
          out << YAML::Key << YAML::BeginSeq << 1 << 2 << YAML::EndSeq << YAML::Value << YAML::BeginSeq << -1 << -2 << YAML::EndSeq;
          out << YAML::Key << 1;
          out << YAML::Value << 2;
          out << YAML::Key << 3;
          out << YAML::Value << 4;
          out << YAML::EndMap;
          //std::cout << "out=\n" << out.c_str() << "\n=out" << std::endl;
          std::string expected_result = "&NodeRegistry::map\n?\n  - 1\n  - 2\n:\n  - -1\n  - -2\n1: 2\n3: 4";
          //std::cout << "out2=\n" << expected_result << std::endl;
          STKUNIT_EXPECT_TRUE(expected_result == std::string(out.c_str()));
        }

        YAML::Emitter yaml;
        std::cout << "\nnodeRegistry.serialize_write(yaml)" << std::endl;
        SerializeNodeRegistry::serialize_write(nodeRegistry, yaml, 0);
        //std::cout << yaml.c_str() << std::endl;
        if (!yaml.good())
          {
            std::cout << "Emitter error: " << yaml.good() << " " <<yaml.GetLastError() << "\n";
            STKUNIT_EXPECT_TRUE(false);
          }
        std::ofstream file1("out.yaml");
        file1 << yaml.c_str();
        file1.close();
        std::ifstream file2("out.yaml");
        YAML::Parser parser(file2);
        YAML::Node doc;

        try {
          while(parser.GetNextDocument(doc)) {
            std::cout << "\n read doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (doc.Type() == YAML::NodeType::Map)
              {
                for(YAML::Iterator it=doc.begin();it!=doc.end();++it) {
                  int key, value;
                  std::cout << "read it.first().Type() = " << it.first().Type() << " it.first().Tag()= " << it.first().Tag() << std::endl;
                  std::cout << "read it.second().Type() = " << it.second().Type() << " it.second().Tag()= " << it.second().Tag() << std::endl;
                  const YAML::Node& keySeq = it.first();
                  for(YAML::Iterator itk=keySeq.begin();itk!=keySeq.end();++itk) {
                    *itk >> key;
                    std::cout << "read key= " << key << std::endl;
                  }
              
                  const YAML::Node& valSeq = it.second();
                  for(YAML::Iterator itv=valSeq.begin();itv!=valSeq.end();++itv) {
                    *itv >> value;
                    std::cout << "read value= " << value << std::endl;
                  }
              
                }
              }
          }
        }
        catch(YAML::ParserException& e) {
          std::cout << e.what() << "\n";
          STKUNIT_EXPECT_TRUE(false);
        }

        file2.close();
        std::ifstream file3("out.yaml");
        NodeRegistry nrNew(eMesh);
        SerializeNodeRegistry::serialize_read(nrNew, file3);
        YAML::Emitter yaml3;
        std::cout << "\nnrNew.serialize_write(yaml3)" << std::endl;
        SerializeNodeRegistry::serialize_write(nrNew, yaml3, 0);
        std::cout << yaml3.c_str() << std::endl;
        
        //exit(1);
      }
#endif
    // start_demo_nodeRegistry_test_parallel_1_quadratic_elem

    // change element to be a serendipity quadratic element
    eMesh.get_bulk_data()->modification_begin();

    //getCellTopologyData< shards::Node  >()
    const CellTopologyData *const cell_topo_data =stk_classic::percept::PerceptMesh::get_cell_topology(block_hex_20);
    CellTopology cell_topo(cell_topo_data);

    for (unsigned isd = 0; isd < 12; isd++)
    {
      nodeRegistry.makeCentroidCoords(element_local, needed_entity_rank.first, isd);
      NodeIdsOnSubDimEntityType& nodeIds_onSE_0_loc = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_rank.first, isd));

      stk_classic::mesh::Entity*  node   = eMesh.get_bulk_data()->get_entity(stk_classic::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE_0_loc[0]->identifier());

      unsigned edge_ord = 8u + isd;
      //unsigned n_edge_ord = cell_topo_data->edge[isd].topology->node_count;
      //std::cout << "n_edge_ord = " << n_edge_ord << std::endl;
      edge_ord = cell_topo_data->edge[isd].node[2];
      eMesh.get_bulk_data()->declare_relation(element_local, *node, edge_ord);
    }

    std::vector<stk_classic::mesh::Part*> add_parts(1, &block_hex_20);
    std::vector<stk_classic::mesh::Part*> remove_parts(1, block_hex_8);
    eMesh.get_bulk_data()->change_entity_parts( element_local, add_parts, remove_parts );

    eMesh.get_bulk_data()->modification_end();
    eMesh.print_info("After quadratic");

    eMesh.save_as("./cube1x1x2_hex-20.e");
    //exit(1);
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

#ifdef STK_BUILT_IN_SIERRA

STKUNIT_UNIT_TEST(nodeRegistry, test_serial_hex8_tet4_24_1)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  // start_demo_nodeRegistry_test_serial_hex8_tet4_24_1

  percept::PerceptMesh eMesh(3u);

  unsigned p_size = eMesh.get_parallel_size();
  unsigned p_rank = eMesh.get_rank();
  Util::setRank(eMesh.get_rank());

  std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,1");
  eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));

  // prepare for breaking into tet elements
  mesh::Part& block_tet_4 = eMesh.get_fem_meta_data()->declare_part("block_tet_4", eMesh.element_rank());
  /// set cell topology for the part block_tet_4
  mesh::fem::set_cell_topology< shards::Tetrahedron<4>  >( block_tet_4 );
  stk_classic::io::put_io_part_attribute(block_tet_4);

  eMesh.commit();
  eMesh.print_info();
  eMesh.save_as(std::string("./")+std::string("cube1x1x")+toString(p_size)+std::string("-orig.e"));

  //mesh::Part* block_hex_8 = const_cast<mesh::Part *>(eMesh.getPart("block_1"));

  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size <= 2)
  {
    // pick an element on the processor boundary
    unsigned elem_num_local = 1;
    unsigned elem_num_ghost = 2;
    if (p_size == 1)
      elem_num_ghost = 1;

    stk_classic::mesh::Entity* element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
    stk_classic::mesh::Entity* element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
    if (p_rank == 1)
    {
      element_local_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_ghost);
      element_ghost_p = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_num_local);
    }

    dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
    dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

    stk_classic::mesh::Entity& element_local = *element_local_p;
    stk_classic::mesh::Entity& element_ghost = *element_ghost_p;

    std::cout << "P["<<p_rank<<"] element_local = " << element_local << std::endl;
    std::cout << "P["<<p_rank<<"] element_ghost = " << element_ghost << std::endl;

    // choose faces to be used for new node locations

    std::vector<NeededEntityType> needed_entity_ranks(2);
    needed_entity_ranks[0] = NeededEntityType(eMesh.face_rank(), 1u);
    needed_entity_ranks[1] = NeededEntityType(eMesh.element_rank(), 1u);

    /*
     * 1st of three steps to create and associate new nodes - register need for new nodes, then check if node is remote, then get
     *   from remote proc if necessary; finally, the local node database is ready to be queried
     *
     * The pattern is to begin the step, loop over all elements (including ghosts) and invoke the local operation
     * The method doForAllSubEntities is a utility for performing the operation on all the sub entities.
     * If more granularity is desired, the member functions can be invoked directly for a particular sub-entity.
     */
    nodeRegistry.beginRegistration();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_ghost, needed_entity_ranks);
    nodeRegistry.endRegistration();

    std::cout << "P["<<p_rank<<"] nodeRegistry size  = " << nodeRegistry.total_size() << std::endl;
    std::cout << "P["<<p_rank<<"] nodeRegistry lsize = " << nodeRegistry.local_size() << std::endl;

    dw() << "P["<<p_rank<<"] nodeRegistry size       = " << nodeRegistry.total_size() << DWENDL;
    dw() << "P["<<p_rank<<"] nodeRegistry lsize      = " << nodeRegistry.local_size() << DWENDL;

    // could do local create of elements here
    nodeRegistry.beginLocalMeshMods();
    nodeRegistry.endLocalMeshMods();

    // check if the newly requested nodes are local or remote
    nodeRegistry.beginCheckForRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_ghost, needed_entity_ranks);
    nodeRegistry.endCheckForRemote();

    // get the new nodes from other procs if they are nonlocal
    nodeRegistry.beginGetFromRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_local, needed_entity_ranks);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_ghost, needed_entity_ranks);
    nodeRegistry.endGetFromRemote();

    // now we can get the new node's id and entity
    unsigned iSubDimOrd = 5u;  // 6'th face of element 1 is shared with 5th of element 2 (top, bot resp.)
    if (p_rank)
    {
      iSubDimOrd = 4u;
    }
    NodeIdsOnSubDimEntityType& nodeIds_onSE_0 = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_ranks[0].first, iSubDimOrd));
    stk_classic::mesh::Entity*  node_0   = eMesh.get_bulk_data()->get_entity(stk_classic::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE_0[0]->identifier());

    // should be the same node on each proc
    std::cout << "P[" << p_rank << "] nodeId_0 = " << nodeIds_onSE_0 << " node_0= " << node_0 << std::endl;
    if (p_size == 2)
    {
      STKUNIT_EXPECT_EQ(nodeIds_onSE_0[0]->identifier(), 20u);
    }

    // end_demo

    // start_demo_nodeRegistry_test_serial_hex8_tet4_24_1_a

    // check for center node
    if (p_rank == 0)
    {
      NodeIdsOnSubDimEntityType& nodeIds_onSE_1 = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_ranks[1].first, 0u));
      stk_classic::mesh::Entity*  node_1   = eMesh.get_bulk_data()->get_entity(stk_classic::mesh::fem::FEMMetaData::NODE_RANK, nodeIds_onSE_1[0]->identifier());
      std::cout << "P[" << p_rank << "] nodeId_1 = " << nodeIds_onSE_1 << " node_1= " << node_1 << std::endl;
      if (p_size == 2)
      {
        STKUNIT_EXPECT_EQ(nodeIds_onSE_1[0]->identifier(), 19u);
      }
    }
    //exit(1);
    // end_demo
  }
}
#endif

} // namespace unit_tests
} // namespace adapt
} // namespace stk_classic

