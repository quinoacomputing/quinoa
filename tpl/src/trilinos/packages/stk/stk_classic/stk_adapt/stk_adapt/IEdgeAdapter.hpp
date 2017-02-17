#ifndef stk_adapt_IEdgeAdapter_hpp
#define stk_adapt_IEdgeAdapter_hpp

#include <stk_adapt/IAdapter.hpp>

namespace stk_classic {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * An IEdgeAdapter is an abstract base class for derived classes that are required to overload the mark method,
     *    which provides info such as the element the edge belongs to, which edge ordinal it is, the nodes of the edge
     *    and the edge coordinates.  
     */
    class IEdgeAdapter : public IAdapter
    {
    public:
      IEdgeAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0);

      /// can be overriden
      virtual ElementUnrefineCollection  buildUnrefineList();

    protected:

      /// Client supplies these methods - given an element, which edge, and the nodes on the edge, return instruction on what to do to the edge,
      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE

      virtual int mark(const stk_classic::mesh::Entity& element, unsigned which_edge, stk_classic::mesh::Entity & node0, stk_classic::mesh::Entity & node1,
                           double *coord0, double *coord1, std::vector<int>* existing_edge_marks) = 0;

      /// This convenience method calls mark and if all edges are marked for unrefine, it returns -1 to unrefine the element.
      /// This method can be overriden to allow for an "element-based" determination that doesn't need to visit edges.
      virtual int markUnrefine(const stk_classic::mesh::Entity& element);

      virtual void 
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
                                              vector<NeededEntityType>& needed_entity_ranks);


    };

    IEdgeAdapter::IEdgeAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk_classic::mesh::FieldBase *proc_rank_field) : 
      IAdapter(eMesh, bp, proc_rank_field)
    {
    }

    void IEdgeAdapter::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element, 
                                            vector<NeededEntityType>& needed_entity_ranks)
    {
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
              throw std::runtime_error("IEdgeAdapter::apply can't use IEdgeAdapter for RefinerPatterns that require face nodes");
            }
          else if (needed_entity_rank == m_eMesh.element_rank())
            {
              numSubDimNeededEntities = 1;
              throw std::runtime_error("IEdgeAdapter::apply can't use IEdgeAdapter for RefinerPatterns that require volume nodes");
            }

          // see how many edges are already marked
          int num_marked=0;
          std::vector<int> edge_marks(numSubDimNeededEntities,0);
          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
                  if (!is_empty) 
                    {
                      edge_marks[iSubDimOrd] = 1;
                      ++num_marked;
                    }
                }
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              if (needed_entity_rank == m_eMesh.edge_rank())
                {
                  stk_classic::mesh::Entity & node0 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk_classic::mesh::Entity & node1 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 = stk_classic::mesh::field_data( *coordField , node0 );
                  double * const coord1 = stk_classic::mesh::field_data( *coordField , node1 );
                  

                  int markInfo = mark(element, iSubDimOrd, node0, node1, coord0, coord1, &edge_marks);

                  bool needNodes = (DO_REFINE & markInfo);
                    {
                      (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, needNodes);
                    }
                }

            } // iSubDimOrd
        } // ineed_ent
    }

    int IEdgeAdapter::markUnrefine(const stk_classic::mesh::Entity& element)
    {
      const CellTopologyData * const cell_topo_data = stk_classic::percept::PerceptMesh::get_cell_topology(element);
                
      CellTopology cell_topo(cell_topo_data);
      const mesh::PairIterRelation elem_nodes = element.relations(stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

      VectorFieldType* coordField = m_eMesh.get_coordinates_field();

      unsigned numSubDimNeededEntities = 0;
      numSubDimNeededEntities = cell_topo_data->edge_count;

      bool unrefAllEdges = true;
      for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
        {
          stk_classic::mesh::Entity & node0 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
          stk_classic::mesh::Entity & node1 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
          double * const coord0 = stk_classic::mesh::field_data( *coordField , node0 );
          double * const coord1 = stk_classic::mesh::field_data( *coordField , node1 );
                  
          int markInfo = mark(element, iSubDimOrd, node0, node1, coord0, coord1, 0);
          bool do_unref = markInfo & DO_UNREFINE;
          if (!do_unref)
            {
              unrefAllEdges = false;
              break;
            }
        }
      if (unrefAllEdges)
        return -1;
      else
        return 0;
    }

    ElementUnrefineCollection IEdgeAdapter::buildUnrefineList()
    {
      ElementUnrefineCollection elements_to_unref;

      const vector<stk_classic::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

      for ( vector<stk_classic::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
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

                if (elem_nodes.size() && m_eMesh.isChildWithoutNieces(element, false))
                  {
                    int markInfo = markUnrefine(element);
                    if (markInfo & DO_UNREFINE)
                      {
                        elements_to_unref.insert(&element);
                        //std::cout << "tmp unref element id= " << element.identifier() << std::endl;
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
