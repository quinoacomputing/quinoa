#ifndef stk_adapt_Refiner_hpp
#define stk_adapt_Refiner_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <utility>
#include <math.h>
#include <map>
#include <set>
#include <vector>


#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include <stk_percept/stk_mesh.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/ProgressMeter.hpp>
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_adapt/Colorer.hpp>

#include <stk_adapt/NodeRegistry.hpp>

#include <stk_adapt/SubDimCell.hpp>

#include <stk_adapt/RefinementInfoByType.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <stk_percept/mesh/geometry/kernel/MeshGeometry.hpp>
#endif

#define UNIFORM_REF_REMOVE_OLD_STD_SET 1
#define UNIFORM_REF_REMOVE_OLD_STD_VECTOR 0
#define UNIFORM_REF_REMOVE_OLD_BOOST_SET 0

#if UNIFORM_REF_REMOVE_OLD_BOOST_SET
#include <boost/unordered_set.hpp>
#endif


namespace stk_classic {
  namespace adapt {

    typedef std::set<stk_classic::mesh::Entity *> ElementUnrefineCollection;
    //typedef std::map<stk_classic::mesh::Part*, stk_classic::mesh::Part*> SidePartMap;
    typedef std::map<std::string, std::string> SidePartMap;

    using std::vector;
    using std::map;
    using std::set;


#if UNIFORM_REF_REMOVE_OLD_STD_SET
    typedef std::set<stk_classic::mesh::Entity *> elements_to_be_destroyed_type;
#endif
#if UNIFORM_REF_REMOVE_OLD_STD_VECTOR
    typedef std::vector<stk_classic::mesh::Entity *> elements_to_be_destroyed_type;
#endif
#if UNIFORM_REF_REMOVE_OLD_BOOST_SET
    typedef boost::unordered_set<stk_classic::mesh::Entity *> elements_to_be_destroyed_type;
#endif


    /// e.g. UniformRefiner<shards::Hex<8>, shards::Tet<4> >
    //template<typename FromTopology, typename ToTopology>
#if 0
    class ParallelMeshModAlgorithm
    {
    public:
      virtual void planActions()=0;
      virtual void performActions()=0;
    protected:
      void helperFunction1();
      void helperFunction2();
      //...
    };
#endif

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    //template<class UniformRefinerPattern>
    class Refiner
#ifndef SWIG //NLM    
  : public stk_classic::percept::Observable<ProgressMeterData>
#endif
    {
    public:
      Refiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk_classic::mesh::FieldBase *proc_rank_field=0);
      Refiner(percept::PerceptMesh& eMesh, Pattern refine_type, stk_classic::mesh::FieldBase *proc_rank_field=0);
      //Refiner(percept::PerceptMesh& eMesh, std::vector<UniformRefinerPatternBase *>&  bp, stk_classic::mesh::FieldBase *proc_rank_field=0);
      
      virtual ~Refiner();
      
      void
      doBreak();

      void
      setRemoveOldElements(bool do_remove);
      bool
      getRemoveOldElements();

      /* for future
      void
      setIOSaveInactiveElements(bool do_save) { m_doIOSaveInactiveElements = do_save; }
      bool
      getIOSaveInactiveElements() { return m_doIOSaveInactiveElements; }
      */

      void
      setGeometryFile(std::string file_name);

      void
      setSmoothGeometry(bool do_smooth) { m_doSmoothGeometry = do_smooth; }
      bool
      getSmoothGeometry() { return m_doSmoothGeometry; }

      void
      setRemoveGeometryBlocks(bool do_remove) { m_removeGeometryBlocks = do_remove; }
      bool
      getRemoveGeometryBlocks() { return m_removeGeometryBlocks; }

      void
      setIgnoreSideSets(bool ignore_sidesets) ;

      bool
      getIgnoreSideSets();

      std::vector< RefinementInfoByType >&
      getRefinementInfoByType();

      void
      setQueryPassOnly(bool doQueryOnly);

      void
      setDoProgressMeter(bool do_progress);
      bool
      getDoProgressMeter();


      // ================================ unrefine

      typedef std::set<stk_classic::mesh::Entity *> NodeSetType;
      typedef std::set<stk_classic::mesh::Entity *> SetOfEntities;


      void
      unrefineTheseElements(ElementUnrefineCollection& elements_to_unref);

      void
      unrefineAll();

      void
      setAlwaysInitializeNodeRegistry(bool do_init) { m_alwaysInitNodeRegistry = do_init; }

      bool
      getAlwaysInitializeNodeRegistry() { return m_alwaysInitNodeRegistry; }

#if defined( STK_PERCEPT_HAS_MESQUITE ) && defined(STK_PERCEPT_HAS_GEOMETRY)

      enum SMOOTHING_OPTIONS {
        // snaps and discards original coord field, tries to smooth
        SNAP_PLUS_SMOOTH,    
        // keeps original and snapped states; does line search between; tries to keep always-valid mesh
        USE_LINE_SEARCH_WITH_MULTIPLE_STATES
        //,END_OPTIONS
      };

      void smoothGeometry(MeshGeometry& mesh_geometry, SMOOTHING_OPTIONS option);
#endif

      void snapAndSmooth(bool geomSnap, std::string geomFile);

      void deleteParentElements();

      void check_db(std::string msg="") ;

      void check_sidesets(std::string msg="");
      void check_sidesets_1(std::string msg);
      void check_sidesets_2(std::string msg);
      void fix_side_sets_1();
      void fix_side_sets_2();
      void fix_side_sets_3(bool checkParentChild, SidePartMap& side_part_map);

      /// determine side part to elem part relations
      void get_side_part_relations(bool checkParentChild, SidePartMap& side_part_map);

      bool connectSides(stk_classic::mesh::Entity *element, stk_classic::mesh::Entity *side_elem, SidePartMap* side_part_map=0);
      void fixElementSides2();
      void fixSides(stk_classic::mesh::Entity *parent);

      NodeRegistry& getNodeRegistry() { return *m_nodeRegistry; }
      percept::PerceptMesh& getMesh() { return m_eMesh; }
    protected:

      //============= unrefine
      void 
      filterUnrefSet(ElementUnrefineCollection& elements_to_unref);

      void
      getKeptNodes(NodeSetType& kept_nodes, ElementUnrefineCollection& elements_to_unref);

      void
      getDeletedNodes(NodeSetType& deleted_nodes, const NodeSetType& kept_nodes, ElementUnrefineCollection& elements_to_unref);


      void 
      removeDeletedNodes(NodeSetType& deleted_nodes);

      void getChildrenToBeRemoved(ElementUnrefineCollection& elements_to_unref,
                                  SetOfEntities& children_to_be_removed, SetOfEntities& children_to_be_removed_with_ghosts, 
                                  SetOfEntities& copied_children_to_be_removed,
                                  SetOfEntities& family_trees_to_be_removed, 
                                  SetOfEntities& parent_elements);

      void removeFamilyTrees(SetOfEntities& family_trees_to_be_removed);


      void getSideElemsToBeRemoved(SetOfEntities& children_to_be_removed, SetOfEntities& side_elem_set_to_be_removed, SetOfEntities& family_trees_to_be_removed, SetOfEntities& parent_side_elements);

      void removeChildElements(SetOfEntities& children_to_be_removed);

      void removeSideElements(SetOfEntities& side_elem_set_to_be_removed, SetOfEntities& elements_to_be_deleted);

      void remesh(SetOfEntities& parent_elements);
      //============= unrefine end

      void check_db_ownership_consistency(std::string msg="");
      void check_db_hanging_nodes();
      void check_db_entities_exist(std::string msg="");

      /**  Overrides start =======>
       */

      /** Overrides
       *   m_nodeRegistry data member.  The same loop should be executed every time this method is called (i.e. the same elements should be visited).  Also, there
       *   is policy associated with the @param only_count and @param doAllElements inputs.  If doAllElements==true, then ghost elements should not be skipped.
       *   If only_count==true, then *do not* call the supplied @param function, rather just count the number of elements to be visited and return that number.
       *   @return Always return the number of elements visited (modulo the doAllElements parameter).
       * 
       *   Note: client code can just use the supplied implementation in this base class - it is not necessary to override
       *
       *   @see UniformRefiner implementation of this method for an example.
       */

      virtual unsigned
      doForAllElements(unsigned irank, std::string function_info,
                       stk_classic::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function,
                       vector< ColorerSetType >& elementColors, unsigned elementType,
                       vector<NeededEntityType>& needed_entity_ranks,  
                       bool only_count=false, bool doAllElements=true) ;

      /** Create a list of nodes from the new nodes that can be easily deciphered by the UniformRefinerPattern.
       *
       *  This is a helper function that gets all the nodes ready for use in creating new elements.  For each rank of those needed by
       *  the algorithm (supplied information from the break pattern), it fills the 3D array new_sub_entity_nodes with the new nodes.
       *
       *  Returns the 3D array new_sub_entity_nodes[entity_rank_of_entity_or_sub_dim_entity][ordinal_of_sub_dim_entity][ordinal_of_node_on_sub_dim_entity]
       *
       */
      virtual bool
      createNewNeededNodeIds(const CellTopologyData * const cell_topo_data,
                             const stk_classic::mesh::Entity& element, vector<NeededEntityType>& needed_entity_ranks, NewSubEntityNodesType& nodes) ;

      /** Method that actually creates new elements by first calling createNewNeededNodeIds then calls the break pattern's createNewElements method.
       *
       *  A sample implementation is shown in @see UniformRefiner
       */
      virtual void
      createElementsAndNodesAndConnectLocal(unsigned irank,  stk_classic::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern,
                                            vector< ColorerSetType >& elementColors,   vector<NeededEntityType>& needed_entity_ranks,
                                            vector<stk_classic::mesh::Entity *>& new_elements_pool) ;

      /** This is a helper method that loops over all sub-dimensional entities whose rank matches on of those in @param needed_entity_ranks
       *    and registers that sub-dimensional entity as needing a new node, or whatever other function NodeRegistry requires (getFromRemote(), etc)
       *  Override it to only apply the @param function to the desired sub-entities (e.g. for non-uniform/local refinement)
       *
       *  It is a copy of NodeRegistry's doForAllSubEntities method.  Provided here so it can be overridden.
       * 
       *  Note: this is the minimal function that needs to be overridden to get different marking/refining behavior
       */

      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk_classic::mesh::Entity& element,
                                              vector<NeededEntityType>& needed_entity_ranks);

      /// =========>  Overrides  end



      void
      removeFamilyTrees();

      void
      removeOldElements(unsigned irank, stk_classic::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern );

      void
      removeElements( elements_to_be_destroyed_type& elements_to_be_destroyed, unsigned irank=0);

      void
      removeEmptyElements();

      // empty nodes (nodes not referred to by any elements) are possibly created during refine, this method removes them
      void
      removeDanglingNodes();

      void
      addOldElementsToPart(stk_classic::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern, unsigned *elementType = 0u);

      void 
      removeFromOldPart(stk_classic::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern);

      void
      renameNewParts(stk_classic::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern);

      void
      fixSurfaceAndEdgeSetNames(stk_classic::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern);

      void
      fixElementSides();

      void
      fixElementSides1();

      void
      fixElementSides(stk_classic::mesh::EntityRank side_rank);

      void
      fixElementSides1(stk_classic::mesh::EntityRank side_rank);

      void
      checkFixElementSides(stk_classic::mesh::EntityRank side_rank, stk_classic::mesh::EntityRank elem_rank);

      void
      buildElementSideDB(SubDimCellToDataMap& cell_2_data_map);

      void
      trace_print(std::string msg);

      void
      checkBreakPatternValidityAndBuildRanks(std::vector<stk_classic::mesh::EntityRank>& ranks);

      void
      set_active_part();


    protected:
      percept::PerceptMesh& m_eMesh;

      //UniformRefinerPatternBase & m_breakPattern;
      std::vector<UniformRefinerPatternBase *> m_breakPattern;

      NodeRegistry* m_nodeRegistry;
      stk_classic::mesh::FieldBase *m_proc_rank_field;
      bool m_doRemove;

      // for future: 
      // bool m_doIOSaveInactiveElements; // default false

      std::vector<stk_classic::mesh::EntityRank> m_ranks;
      bool m_ignoreSideSets;
      std::string m_geomFile;
      bool m_geomSnap;

      std::vector< RefinementInfoByType > m_refinementInfoByType;
      bool m_doQueryOnly;

      int m_progress_meter_frequency;
      bool m_doProgress;
      
      bool m_alwaysInitNodeRegistry;
      bool m_doSmoothGeometry;
      bool m_allocated;

      bool m_removeGeometryBlocks;

    };



  }
}
#endif
