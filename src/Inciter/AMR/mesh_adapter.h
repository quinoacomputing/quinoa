#ifndef QUINOA_MESH_ADAPTER_H
#define QUINOA_MESH_ADAPTER_H

#include "DerivedData.h"

#include "types.h"
#include "util.h"
#include "id_generator.h"

#include "marked_refinements_store.h"
#include "tet_store.h"

#include "node_connectivity.h"
#include "node_store.h"

#include "refinement.h"
//#include "derefinement.h"

#include "Refinement_State.h"

namespace AMR {
    class mesh_adapter_t {

        public:

            //! Default constructor for migration
            mesh_adapter_t() {}

            //! Constructor taking a mesh graph
            mesh_adapter_t( const std::vector< std::size_t >& inpoel ) :
                node_connectivity( tk::npoin(inpoel) )
            {
                consume_tets( inpoel );
                tet_store.generate_edges();
            }

            //AMR::refinement_t init(const std::vector<size_t>& tetinpoel, size_t num_nodes);
            void init_node_store(coord_type* m_x, coord_type* m_y, coord_type* m_z, size_t graph_size);

            // TODO: Set these in a better way
            const real_t derefinement_cut_off = 0.2;
            const real_t refinement_cut_off = 0.9;

            AMR::tet_store_t tet_store;
            AMR::node_connectivity_t node_connectivity;

            // for coord type stuff
            AMR::node_store_t node_store;

            AMR::refinement_t refiner;

            void consume_tets(const std::vector<std::size_t>& tetinpoel );

            void evaluate_error_estimate();
            void uniform_refinement();
            void error_refinement( const std::vector< edge_t >& edge,
                    const std::vector< real_t >& criteria );

            int detect_compatibility(
                    int num_locked_edges,
                    int num_intermediate_edges,
                    AMR::Refinement_Case refinement_case,
                    int normal=0
            );

            void lock_intermediates();

            void mark_refinement();
            void perform_refinement();

            void refinement_class_one(int num_to_refine, size_t tet_id);
            void refinement_class_two(edge_list_t edge_list, size_t tet_id);
            void refinement_class_three(size_t tet_id);

            void lock_tet_edges(size_t tet_id);
            void deactivate_tet_edges(size_t tet_id);
            bool check_valid_refinement_case(size_t child_id);

            void check_derefinement();
            //std::vector< std::size_t >& get_active_inpoel();

            void print_tets();

            void reset_intermediate_edges();
            void update_tet_edges_lock_type(size_t tet_id, AMR::Edge_Lock_Case check, AMR::Edge_Lock_Case new_case);
            void remove_edge_locks(int intermediate = 0);
            void remove_normals();

    };
}

#endif //QUINOA_MESH_ADAPTER_H
