#ifndef stk_adapt_TestLocalRefinerTet_N_4_hpp
#define stk_adapt_TestLocalRefinerTet_N_4_hpp

#include <stk_adapt/Refiner.hpp>

namespace stk_classic {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that marks some edges randomly to test RefinerPattern_Tri3_Tri3_N_4
     */
    class TestLocalRefinerTet_N_4 : public Refiner
    {
    public:
      TestLocalRefinerTet_N_4(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0, unsigned edge_mark_bitcode=1);

      ElementUnrefineCollection  buildTestUnrefineList();

    protected:


      virtual void 
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
            vector<NeededEntityType>& needed_entity_ranks);


      unsigned m_edge_mark_bitcode;
    };

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_4 in UnitTestLocalRefiner.cpp)

    TestLocalRefinerTet_N_4::TestLocalRefinerTet_N_4(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk_classic::mesh::FieldBase *proc_rank_field, unsigned edge_mark_bitcode) : 
      Refiner(eMesh, bp, proc_rank_field), m_edge_mark_bitcode(edge_mark_bitcode)
    {
    }


    void TestLocalRefinerTet_N_4::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks)
    {
      //static int n_seq = 400;

      const CellTopologyData * const cell_topo_data = stk_classic::percept::PerceptMesh::get_cell_topology(element);
                
      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

      VectorFieldType* coordField = m_eMesh.get_coordinates_field();

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
                  stk_classic::mesh::Entity & node0 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk_classic::mesh::Entity & node1 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = stk_classic::mesh::field_data( *coordField , node0 );
                  double * const coord1 = stk_classic::mesh::field_data( *coordField , node1 );
                  
                  // vertical plane position
                  const double vx = 0.21;

                  // horizontal plane position
                  const double vy = 1.21;

                  // in plane position
                  const double vz = 1.21;

                  // choose to refine or not 
                  bool do_refine = false;
                  if (
                      ( std::fabs(coord0[0]-coord1[0]) > 1.e-3 &&
                        ( (coord0[0] < vx && vx < coord1[0]) || (coord1[0] < vx && vx < coord0[0]) )
                        )
#if 1
                      ||
                      ( std::fabs(coord0[1]-coord1[1]) > 1.e-3 &&
                        ( (coord0[1] < vy && vy < coord1[1]) || (coord1[1] < vy && vy < coord0[1]) )
                        )
                      ||
                      ( std::fabs(coord0[2]-coord1[2]) > 1.e-3 &&
                        ( (coord0[2] < vz && vz < coord1[2]) || (coord1[2] < vz && vz < coord0[2]) )
                        )
#endif
                      )
                    {
                      do_refine = true;
                    }

                  (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, do_refine);
                }

            } // iSubDimOrd
        } // ineed_ent
    }

    ElementUnrefineCollection TestLocalRefinerTet_N_4::buildTestUnrefineList()
    {
      ElementUnrefineCollection elements_to_unref;

      const vector<stk_classic::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

      for ( vector<stk_classic::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          //if (removePartSelector(**k)) 
          {
            stk_classic::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk_classic::mesh::Entity& element = bucket[ientity];

                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error is isParentElement)
                const bool check_for_family_tree = false;  
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
              
                if (isParent)
                  continue;
                
                const mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

                if (elem_nodes.size() && m_eMesh.isChildWithoutNieces(element, false) )
                  {
                    bool found = true;
                    for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                      {
                        stk_classic::mesh::Entity *node = elem_nodes[inode].entity();
                        double *coord = stk_classic::mesh::field_data( *m_eMesh.get_coordinates_field(), *node );
                        if (coord[0] > 2.1 || coord[1] > 2.1 || coord[2] > 2.1)
                          {
                            found = false;
                            break;
                          }
                      }
                    if (found)
                      {
                        elements_to_unref.insert(&element);
                        //std::cout << "tmp element id= " << element.identifier() << " ";
                        //m_eMesh.print_entity(std::cout, element);
                      }
                  }
              }
          }
        }

      return elements_to_unref;
    }


  }
}
#endif
