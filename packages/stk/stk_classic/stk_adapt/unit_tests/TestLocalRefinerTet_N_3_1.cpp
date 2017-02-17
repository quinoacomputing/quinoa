
#include <unit_tests/TestLocalRefinerTet_N_3_1.hpp>

namespace stk_classic {
  namespace adapt {

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_3_1 in UnitTestLocalRefiner.cpp)

    TestLocalRefinerTet_N_3_1::TestLocalRefinerTet_N_3_1(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk_classic::mesh::FieldBase *proc_rank_field, unsigned edge_mark_bitcode) : 
      Refiner(eMesh, bp, proc_rank_field), m_edge_mark_bitcode(edge_mark_bitcode)
    {
    }


    void TestLocalRefinerTet_N_3_1::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks)
    {
      //static int n_seq = 400;

      const CellTopologyData * const cell_topo_data = stk_classic::percept::PerceptMesh::get_cell_topology(element);
                
      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

      //VectorFieldType* coordField = m_eMesh.get_coordinates_field();

      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;
          stk_classic::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_rank == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
            }
          else if (needed_entity_rank == m_eMesh.element_rank())
            {
              numSubDimNeededEntities = 1;
            }

          // see how many edges are already marked
          int num_marked=0;
          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
                  if (!is_empty) ++num_marked;
                }
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
              //SubDimCell_SDSEntityType subDimEntity;
              //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
              //bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
              //if(1||!is_empty)

              if (needed_entity_rank == m_eMesh.edge_rank())
                {
#if 0
                  stk_classic::mesh::Entity & node0 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk_classic::mesh::Entity & node1 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = stk_classic::mesh::field_data( *coordField , node0 );
                  double * const coord1 = stk_classic::mesh::field_data( *coordField , node1 );
                  
                  // vertical line position
                  const double vx = 0.21;

                  // horizontal line position
                  const double vy = 1.21;

                  // choose to refine or not 
                  if (
                      ( std::fabs(coord0[0]-coord1[0]) > 1.e-3 &&
                        ( (coord0[0] < vx && vx < coord1[0]) || (coord1[0] < vx && vx < coord0[0]) )
                        )
                      ||
                      ( std::fabs(coord0[1]-coord1[1]) > 1.e-3 &&
                        ( (coord0[1] < vy && vy < coord1[1]) || (coord1[1] < vy && vy < coord0[1]) )
                        )
                      )
                    {
                      (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true);
                    }

#endif
                  // mark only the first element


                  if ( ((1 << iSubDimOrd) & m_edge_mark_bitcode ) && 1 == element.identifier()  )
                    {
                      if (1)
                        {
                          std::cout << "tmp TestLocalRefinerTet_N_3_1 element.identifier() = " << element.identifier() 
                                    << " edge_mark_bitcode = " << m_edge_mark_bitcode << "  iSubDimOrd= " << iSubDimOrd << std::endl;
                        }
                      (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true);
                    }

                }

            } // iSubDimOrd
        } // ineed_ent
    }



  }
}

