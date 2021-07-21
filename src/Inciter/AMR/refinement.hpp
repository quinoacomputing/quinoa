#ifndef AMR_refinement_h
#define AMR_refinement_h

#include <algorithm>

#include "Macro.hpp"
#include "tet_store.hpp"
#include "node_connectivity.hpp"

// TODO: make this have a base class to support multiple generator schemes
// using the policy design pattern

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

namespace AMR {

    class refinement_t {
        private:

            size_t DEFAULT_REFINEMENT_LEVEL = 0; //TODO: Is this in the right place?
            size_t MIN_REFINEMENT_LEVEL = DEFAULT_REFINEMENT_LEVEL;
            // list of "intermediate" edges to be deleted
            std::vector< edge_t > delete_list;

        public:

            size_t MAX_REFINEMENT_LEVEL = 3;

            // TODO: Document this
            child_id_list_t generate_child_ids( tet_store_t& tet_store, size_t parent_id, size_t count = MAX_CHILDREN)
            {
                //return morton_id_generator_t::get_children_ids(parent_id);
                return tet_store.generate_child_ids(parent_id, count);
            }

            /**
             * @brief function to detect when an invalid refinement is
             * invoked
             *
             * @param tet_store Tet store to use
             * @param tet_id Id the of the tet which will be refined
             *
             * @return A bool stating if the tet can be validly refined
             */
            bool check_allowed_refinement( tet_store_t& tet_store, size_t tet_id)
            {
                Refinement_State& master_element = tet_store.data(tet_id);

                // These asserts mean we never actually try refine a 1:2 or 1:4
                assert( master_element.refinement_case !=
                        Refinement_Case::one_to_two);
                assert( master_element.refinement_case !=
                        Refinement_Case::one_to_four);

                // cppcheck-suppress assertWithSideEffect
                assert( tet_store.is_active(tet_id) );

                // Check this won't take us past the max refinement level
                if (master_element.refinement_level >= MAX_REFINEMENT_LEVEL)
                {
                    return false;
                }

                // If we got here, we didn't detect anything which tells us not
                // to refine
                return true;
            }

            /**
             * @brief Method which takes a tet id, and deduces the other
             * parameters needed to perform a 1:2
             *
             * @param tet_store Tet store to use
             * @param node_connectivity Mesh node connectivity (graph)
             * @param tet_id The id to refine 1:2
             */
            void refine_one_to_two( tet_store_t& tet_store, node_connectivity_t& node_connectivity, size_t tet_id)
            {
                edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);
                node_pair_t nodes = find_single_refinement_nodes(tet_store,edge_list);
                refine_one_to_two( tet_store, node_connectivity, tet_id, nodes[0], nodes[1]);
            }

/*
            //! @brief Method which takes a tet id, and transforms arguments
            //!   into the form needed for the main 1:2 refinement method
            //! @param tet_id The id to refine 1:2
            void refine_one_to_two(
                    size_t tet_id,
                    std::string edge_key
            )
            {
                std::vector<std::string> nodes = util::split(edge_key,KEY_DELIM);
                size_t edge_node_A_id =  std::stoul (nodes[0],nullptr,0);
                size_t edge_node_B_id =  std::stoul (nodes[1],nullptr,0);
                refine_one_to_two( tet_id, edge_node_A_id, edge_node_B_id);
            }
*/
            /**
             * @brief Refine a given tet id into 2 children.
             * NOTE: Does not do any validity checking (currently?)
             *
             * @param tet_store Tet store to use
             * @param node_connectivity Mesh node connectivity (graph)
             * @param tet_id Id of tet to refine
             * @param edge_node_A_id The first node of id of the edge which
             * will be split
             * @param edge_node_B_id The second node of id of the
             * edge which will be split
             */
            void refine_one_to_two(
                    tet_store_t& tet_store,
                    node_connectivity_t& node_connectivity,
                    size_t tet_id,
                    size_t edge_node_A_id,
                    size_t edge_node_B_id
            )
            {

                trace_out << "refine_one_to_two" << std::endl;
                if (!check_allowed_refinement(tet_store,tet_id)) return;

                tet_t original_tet = tet_store.get(tet_id);

                //coordinate_t original_tet_c = node_connectivity->id_to_coordinate(id);

                size_t new_node_id = node_connectivity.add( edge_node_A_id, edge_node_B_id );

                /// Split existing tet into two new tets

                // The two new tets will be the same, but for each an edge will
                // be cut, losing an edge  replaced by E

                tet_t new_tet1;
                tet_t new_tet2;

                // Create a new tet that is based on the original
                copy_tet(&new_tet1, &original_tet);

                // Replace all node ids in tet that were pointing to A with new_node_id
                replace_node(&new_tet1, edge_node_A_id, new_node_id);

                // Create a new tet that is based on the original
                copy_tet(&new_tet2, &original_tet);

                // Replace all node ids in tet that were pointing to B with new_node_id
                replace_node(&new_tet2, edge_node_B_id, new_node_id);

                // Now, update the edge list

                // Generate edges for split
                tet_store.edge_store.split(edge_node_A_id, edge_node_B_id, new_node_id,
                        Edge_Lock_Case::intermediate);

                child_id_list_t child_list = generate_child_ids(tet_store,tet_id, 2);

                size_t first_child_id = child_list[0];
                size_t second_child_id = child_list[1];

                // Add the two new tets to the system
                size_t new_tet_id = first_child_id;
                tet_store.add(
                        first_child_id,
                        new_tet1,
                        Refinement_Case::one_to_two,
                        tet_id
                );

                //size_t new_tet_id2 = second_child_id;
                tet_store.add(
                        second_child_id,
                        new_tet2,
                        Refinement_Case::one_to_two,
                        tet_id
                );

                //trace_out << "1:2 DOING REFINE OF " << tet_id << ". Adding " << child_list[0] << " and " << child_list[1] << std::endl;

                // This call is only needed to add a single edge, from the new
                // node to the node on the normal to that face, but avoids
                // directly calculating which nodes that is
                tet_store.generate_edges(new_tet_id);

                // Currently we lock one per tet, around the split node. We
                // also need to lock the two "arms" which come out from it
                //lock_edges_from_node(new_tet_id, new_node_id, Edge_Lock_Case::intermediate);
                //lock_edges_from_node(new_tet_id2, new_node_id, Edge_Lock_Case::intermediate);

                // Deactivate parent tet?
                tet_store.deactivate(tet_id);
                //lock_edges_from_node(new_node_id, Edge_Lock_Case::intermediate);
                trace_out << "Adding " << new_node_id << " to intermediate list " << std::endl;
                tet_store.intermediate_list.insert(new_node_id);
            }

            /**
             * @brief Method which takes a tet id, and deduces the other
             * parameters needed to perform a 1:4
             *
             * @param tet_store Tet store to use
             * @param node_connectivity Mesh node connectivity (graph)
             * @param tet_id The id to refine 1:4
            */
            void refine_one_to_four( tet_store_t& tet_store,
                    node_connectivity_t& node_connectivity, size_t tet_id)
            {
                trace_out << "do refine 1:4 " << std::endl;
                //bool face_refine = false;
                size_t face_refine_id = 0; // FIXME: Does this need a better default
                face_list_t face_list = tet_store.generate_face_lists(tet_id);

                // Iterate over each face
                for (size_t face = 0; face < NUM_TET_FACES; face++)
                {
                    int num_face_refine_edges = 0;

                    face_ids_t face_ids = face_list[face];
                    trace_out << "face ids " <<
                        face_ids[0] << ", " <<
                        face_ids[1] << ", " <<
                        face_ids[2] << ", " <<
                        std::endl;

                    edge_list_t face_edge_list = AMR::edge_store_t::generate_keys_from_face_ids(face_ids);
                    // For this face list, see which ones need refining
                    trace_out << "Looping to " << NUM_FACE_NODES << std::endl;
                    for (size_t k = 0; k < NUM_FACE_NODES; k++)
                    {
                        trace_out << "nodes " << k << std::endl;

                        edge_t edge = face_edge_list[k];
                        if (tet_store.edge_store.get(edge).needs_refining == 1)
                        {
                            num_face_refine_edges++;
                            trace_out << "Ref " << edge << " Num face => " << num_face_refine_edges << std::endl;
                        }

                        // Check for locked edges
                            // This case only cares about faces with no locks
                        if (tet_store.edge_store.lock_case(edge) != Edge_Lock_Case::unlocked)
                        {
                            // Abort this face
                            trace_out << "Face has lock it's not this one " << face << std::endl;
                            num_face_refine_edges = 0;
                            break;
                        }
                        trace_out << "Num face => " << num_face_refine_edges << std::endl;
                    }
                    if (num_face_refine_edges >= 2)
                    {
                        assert(num_face_refine_edges < 4);
                        //face_refine = true;
                        trace_out << "Accepting face " << face << std::endl;
                        face_refine_id = face;
                        break;
                    }
                }

                tet_t tet = tet_store.get(tet_id);
                size_t opposite_offset = AMR::node_connectivity_t::face_list_opposite(face_list, face_refine_id);
                size_t opposite_id = tet[opposite_offset];

                trace_out << "1:4 tet mark id " << tet_id << std::endl;
                trace_out << "opposite offset " << opposite_offset << std::endl;
                trace_out << "opposite id " << opposite_id << std::endl;
                trace_out << "face refine id " << face_refine_id << std::endl;
                trace_out << "face list 0 " << face_list[face_refine_id][0] << std::endl;
                trace_out << "face list 1 " << face_list[face_refine_id][1] << std::endl;
                trace_out << "face list 2 " << face_list[face_refine_id][2] << std::endl;

                refine_one_to_four(tet_store, node_connectivity, tet_id, face_list[face_refine_id], opposite_id);
            }

            /**
             * @brief Method which takes a tet id, and deduces the other
             * parameters needed to perform a 1:4, as a part of an 8:4 deref
             *
             * @param tet_store Tet store to use
             * @param node_connectivity Mesh node connectivity (graph)
             * @param tet_id The id to refine 1:4
            */
            void deref_refine_one_to_four( tet_store_t& tet_store,
                    node_connectivity_t& node_connectivity, size_t tet_id)
            {
                trace_out << "do refine 1:4 " << std::endl;
                //bool face_refine = false;
                size_t face_refine_id = 0; // FIXME: Does this need a better default
                face_list_t face_list = tet_store.generate_face_lists(tet_id);

                // Iterate over each face
                for (size_t face = 0; face < NUM_TET_FACES; face++)
                {
                    int num_face_refine_edges = 0;

                    face_ids_t face_ids = face_list[face];
                    trace_out << "face ids " <<
                        face_ids[0] << ", " <<
                        face_ids[1] << ", " <<
                        face_ids[2] << ", " <<
                        std::endl;

                    edge_list_t face_edge_list = AMR::edge_store_t::generate_keys_from_face_ids(face_ids);
                    // For this face list, see which ones need refining
                    trace_out << "Looping to " << NUM_FACE_NODES << std::endl;
                    for (size_t k = 0; k < NUM_FACE_NODES; k++)
                    {
                        trace_out << "nodes " << k << std::endl;

                        edge_t edge = face_edge_list[k];
                        if (tet_store.edge_store.get(edge).needs_refining == 2)
                        {
                            num_face_refine_edges++;
                            trace_out << "Ref " << edge << " Num face => " << num_face_refine_edges << std::endl;
                        }

                        // Check for locked edges
                            // This case only cares about faces with no locks
                        if (tet_store.edge_store.lock_case(edge) != Edge_Lock_Case::unlocked)
                        {
                            // Abort this face
                            trace_out << "Face has lock it's not this one " << face << std::endl;
                            num_face_refine_edges = 0;
                            break;
                        }
                        trace_out << "Num face => " << num_face_refine_edges << std::endl;
                    }
                    if (num_face_refine_edges >= 2)
                    {
                        assert(num_face_refine_edges < 4);
                        //face_refine = true;
                        trace_out << "Accepting face " << face << std::endl;
                        face_refine_id = face;
                        break;
                    }
                }

                tet_t tet = tet_store.get(tet_id);
                size_t opposite_offset = AMR::node_connectivity_t::face_list_opposite(face_list, face_refine_id);
                size_t opposite_id = tet[opposite_offset];

                trace_out << "1:4 tet mark id " << tet_id << std::endl;
                trace_out << "opposite offset " << opposite_offset << std::endl;
                trace_out << "opposite id " << opposite_id << std::endl;
                trace_out << "face refine id " << face_refine_id << std::endl;
                trace_out << "face list 0 " << face_list[face_refine_id][0] << std::endl;
                trace_out << "face list 1 " << face_list[face_refine_id][1] << std::endl;
                trace_out << "face list 2 " << face_list[face_refine_id][2] << std::endl;

                refine_one_to_four(tet_store, node_connectivity, tet_id, face_list[face_refine_id], opposite_id);
            }

            /**
             * @brief Refine a given tet id into 4 children.
             * NOTE: Does not do any validity checking (currently?)
             *
             * @param tet_store Tet store to use
             * @param node_connectivity Mesh node connectivity (graph)
             * @param tet_id The id of the tet to refine
             * @param face_ids The ids which make the face to be split
             * @param opposite_id The remaining id which is "opposite" the
             * split face
             */
            void refine_one_to_four(
                    tet_store_t& tet_store,
                    node_connectivity_t& node_connectivity,
                    size_t tet_id,
                    std::array<size_t, NUM_FACE_NODES> face_ids,
                    size_t opposite_id
            )
            {

                trace_out << "refine_one_to_four" << std::endl;
                if (!check_allowed_refinement(tet_store,tet_id)) return;

                trace_out << "Refining tet_id " << tet_id <<
                    " 1:4 opposite edge " << opposite_id << std::endl;

                tet_t t = tet_store.get(tet_id);
                trace_out  << "Tet has nodes " <<
                    t[0] << ", " <<
                    t[1] << ", " <<
                    t[2] << ", " <<
                    t[3] << ", " <<
                    std::endl;

                trace_out << "face_ids " <<
                    face_ids[0] << ", " <<
                    face_ids[1] << ", " <<
                    face_ids[2] << ", " <<
                    std::endl;

                size_t A = face_ids[0];
                size_t B = face_ids[1];
                size_t C = face_ids[2];
                size_t D = opposite_id;

                trace_out <<
                    " A " << A <<
                    " B " << B <<
                    " C " << C <<
                    " D " << D <<
                    std::endl;

                // Make new nodes
                //coordinate_t AB_mid = node_connectivity->find_mid_point(A, B);
                size_t AB = node_connectivity.add(A,B);

                //coordinate_t AC_mid = node_connectivity->find_mid_point(A, C);
                size_t AC = node_connectivity.add(A,C);

                //coordinate_t BC_mid = node_connectivity->find_mid_point(B, C);
                size_t BC = node_connectivity.add(B,C);

                // Use nodes to update edges
                // All added edges will be locked due to containing intermediate points
                // Split Outer face  edges
                tet_store.edge_store.split(A, C, AC, Edge_Lock_Case::intermediate);
                tet_store.edge_store.split(A, B, AB, Edge_Lock_Case::intermediate);
                tet_store.edge_store.split(B, C, BC, Edge_Lock_Case::intermediate);

                // Connect D to intermediate points
                tet_store.edge_store.generate(D, AC, Edge_Lock_Case::intermediate);
                tet_store.edge_store.generate(D, BC, Edge_Lock_Case::intermediate);
                tet_store.edge_store.generate(D, AB, Edge_Lock_Case::intermediate);
                // Connect inner edges
                tet_store.edge_store.generate(AC, BC, Edge_Lock_Case::intermediate);
                tet_store.edge_store.generate(AC, AB, Edge_Lock_Case::intermediate);
                tet_store.edge_store.generate(AB, BC, Edge_Lock_Case::intermediate);

                // Make new Tets
                //  This is just the node opposite the face plus each pair
                //  of the news nodes, and the old corner
                //  FIXME: How to find that near corner programatically?

                // Hard coded solution
                // A AC AB D
                // AC AB BC D
                // AC BC C D
                // AB B BC D

                size_t num_children = 4;
                child_id_list_t child = generate_child_ids(tet_store,tet_id, num_children);

                // Outsides
                tet_store.add(child[0], A,  AB, AC, D, Refinement_Case::one_to_four, tet_id);
                tet_store.add(child[2], AC, BC, C,  D, Refinement_Case::one_to_four, tet_id);
                tet_store.add(child[3], AB, B,  BC, D, Refinement_Case::one_to_four, tet_id);

                // Center
                size_t center_id = child[1]; // 1 to preserve Jacobian order
                tet_store.add(center_id, AC, AB, BC, D, Refinement_Case::one_to_four, tet_id);


                // TODO: replace this with a more concise way to lock the correct edges

                tet_store.add_center(center_id);
                /*
                lock_edges_from_node(child[0], AB, Edge_Lock_Case::intermediate);
                lock_edges_from_node(child[0], AC, Edge_Lock_Case::intermediate);
                lock_edges_from_node(child[2], AC, Edge_Lock_Case::intermediate);
                lock_edges_from_node(child[2], BC, Edge_Lock_Case::intermediate);
                lock_edges_from_node(child[3], AB, Edge_Lock_Case::intermediate);
                lock_edges_from_node(child[3], BC, Edge_Lock_Case::intermediate);
                lock_edges_from_node(center_id, AC, Edge_Lock_Case::intermediate);
                lock_edges_from_node(center_id, AB, Edge_Lock_Case::intermediate);
                lock_edges_from_node(center_id, BC, Edge_Lock_Case::intermediate);
                */


                tet_store.deactivate(tet_id);

                //trace_out << "1:4 DOING REFINE OF " << tet_id << ". Adding "
                    // << child[0] << ", "
                    // << child[1] << ", "
                    // << child[2] << ", "
                    // << child[3]
                    // << std::endl;

                /*
                lock_edges_from_node(AB, Edge_Lock_Case::intermediate);
                lock_edges_from_node(AC, Edge_Lock_Case::intermediate);
                lock_edges_from_node(BC, Edge_Lock_Case::intermediate);
                */

                trace_out << "Adding " << AB << " to intermediate list " << std::endl;
                tet_store.intermediate_list.insert(AB);
                trace_out << "Adding " << AC << " to intermediate list " << std::endl;
                tet_store.intermediate_list.insert(AC);
                trace_out << "Adding " << BC << " to intermediate list " << std::endl;
                tet_store.intermediate_list.insert(BC);

            }

            /**
             * @brief Refine a given tet id into 8 children.
             * NOTE: Does not do any validity checking (currently?)
             *
             * @param tet_store Tet store to use
             * @param node_connectivity Mesh node connectivity (graph)
             * @param tet_id Id of tet to refine
             */
            void refine_one_to_eight( tet_store_t& tet_store,
                    node_connectivity_t& node_connectivity, size_t tet_id)
            {

                trace_out << "refine_one_to_eight" << std::endl;
                if (!check_allowed_refinement(tet_store,tet_id)) return;

                // Split every edge into two
                // Makes 4 tets out of the old corners and 3 near mid-points
                // Make 4 out of the midpoints

                // For Tet {ABCD} need to know all (non-repeating) node pairs
                // {AB} {AC} {AD} {BC} {BD} {CD}
                // This can either be hard coded, or generated with a 2d loop
                // The loop would just be i=0..4, j=i..4
                //


                tet_t tet = tet_store.get(tet_id);

                size_t A = tet[0];
                size_t B = tet[1];
                size_t C = tet[2];
                size_t D = tet[3];

                trace_out << "A " << A << " B " << B << " C " << C << " D " << D
                    << std::endl;

                // Generate pairs of nodes (i.e edges)
                // Hard coding for now, can swap out for loop
                //coordinate_t AB_mid = node_connectivity->find_mid_point(A,B);
                size_t AB = node_connectivity.add(A,B);

                //coordinate_t AC_mid = node_connectivity->find_mid_point(A,C);
                size_t AC = node_connectivity.add(A,C);

                //coordinate_t AD_mid = node_connectivity->find_mid_point(A,D);
                size_t AD = node_connectivity.add(A,D);

                //coordinate_t BC_mid = node_connectivity->find_mid_point(B,C);
                size_t BC = node_connectivity.add(B,C);

                //coordinate_t BD_mid = node_connectivity->find_mid_point(B,D);
                size_t BD = node_connectivity.add(B,D);

                //coordinate_t CD_mid = node_connectivity->find_mid_point(C,D);
                size_t CD = node_connectivity.add(C,D);

                // Update edges

                tet_store.edge_store.split(A, C, AC, Edge_Lock_Case::unlocked);
                tet_store.edge_store.split(A, B, AB, Edge_Lock_Case::unlocked);
                tet_store.edge_store.split(A, D, AD, Edge_Lock_Case::unlocked);
                tet_store.edge_store.split(B, C, BC, Edge_Lock_Case::unlocked);
                tet_store.edge_store.split(B, D, BD, Edge_Lock_Case::unlocked);
                tet_store.edge_store.split(C, D, CD, Edge_Lock_Case::unlocked);


                // Outside edges for face ABC
                tet_store.edge_store.generate(AC, BC, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(AC, AB, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(AB, BC, Edge_Lock_Case::unlocked);

                // Outside edges for face ACD
           	tet_store.edge_store.generate(AC, AD, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(AD, CD, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(AC, CD, Edge_Lock_Case::unlocked);

                // Outside edges for face BCD
                tet_store.edge_store.generate(BD, CD, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(BD, BC, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(CD, BC, Edge_Lock_Case::unlocked);

                // Outside edges for face ABD
                 tet_store.edge_store.generate(AD, BD, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(AB, AD, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(AB, BD, Edge_Lock_Case::unlocked);

                // Interior Edges
                   tet_store.edge_store.generate(AC, BD, Edge_Lock_Case::unlocked);
                tet_store.edge_store.generate(CD, AD, Edge_Lock_Case::unlocked);

                // Add the new tets
                //
                // External
                // A AB AC AD - A
                // B BA BC BD - B
                // C CA CB CD - C
                // D DA DB DC - D
                // -
                // Internal (for a face BDC, it's the intermediate and mid opposite)
                // BC CD DB AC - BDC
                // AB BD AD AC - ABD
                // AB AC BC BD - ABC
                // AC AD CD BD - ACD
                //

                // TODO: This is actually generating IDs not trying to get them
                child_id_list_t child = generate_child_ids(tet_store,tet_id);

                // This order should give a positive Jacobian
                tet_store.add(child[0], A, AB, AC, AD, Refinement_Case::one_to_eight, tet_id);
                tet_store.add(child[1], B, BC, AB, BD, Refinement_Case::one_to_eight, tet_id);
                tet_store.add(child[2], C, AC, BC, CD, Refinement_Case::one_to_eight, tet_id);
                tet_store.add(child[3], D, AD, CD, BD, Refinement_Case::one_to_eight, tet_id);

                tet_store.add(child[4], BC, CD, AC, BD, Refinement_Case::one_to_eight, tet_id);
                tet_store.add(child[5], AB, BD, AC, AD, Refinement_Case::one_to_eight, tet_id);
                tet_store.add(child[6], AB, BC, AC, BD, Refinement_Case::one_to_eight, tet_id);
                tet_store.add(child[7], AC, BD, CD, AD, Refinement_Case::one_to_eight, tet_id);

                tet_store.deactivate(tet_id);

                //trace_out << "1:8 DOING REFINE OF " << tet_id << ". "
                    // << child[0] << ", "
                    // << child[1] << ", "
                    // << child[2] << ", "
                    // << child[3] << ", "
                    // << child[4] << ", "
                    // << child[5] << ", "
                    // << child[6] << ", "
                    // << child[7]
                    // << std::endl;

            }

            // This is just a simple assignment, but I wanted to abstract it
            // for if we change the underlying type to something which a simple
            // assignment is no longer safe
            /**
             * @brief Function to duplicate (deep copy) a tet. Useful for when
             * you want to make a tet that's very similar to an existing one
             *
             * @param out The tet to store the copy
             * @param original The tet to copy the data from
             */
            void copy_tet(tet_t* out, tet_t* original)
            {
                // NOTE: This will do a deep copy, so is safer than it may look
                *out = *original;
            }

            // This is just a std::replace, but may need to be more complicated
            // in the future?
            /**
             * @brief function to take an existing list of tet ids and
             * replace one. This can be useful for when you want to build very
             * similar tets which share nodes
             *
             * @param tet Tet to perform operation on
             * @param remove Element to be replaced
             * @param add Element to replace with
             */
            void replace_node(tet_t* tet, size_t remove, size_t add)
            {
                std::replace(tet->begin(), tet->end(), remove, add);
            }

            /**
             * @brief Function to find out slot in the x,y,z data arrays a tet lives
             *
             * // NOTE: this is _currently_ trivial, but will be nice if we for
             * example swap data stores to a map
             *
             * @param tet_store Tet store to use
             * @param tet tet of the tet to look for
             * @param element offset into that tet to look at
             *
             * @return tet into data arrays the tet lives
             */
            // TODO: Move this (or rename?)
            size_t tet_id_to_node_id( tet_store_t& tet_store, size_t tet, size_t element) {
                return tet_store.get(tet)[element];
            }

            /**
             * @brief Function to find the nodes which make up the
             * single (or first?º edge which needs to be refined in an given
             * edge_list
             *
             * @param tet_store Tet store to use
             * @param edge_list The edge list to search for a refinement edge
             *
             * @return The node pair which represent the edge which needs
             * refining
             */
            node_pair_t find_single_refinement_nodes( tet_store_t& tet_store, edge_list_t edge_list)
            {
                node_pair_t returned_nodes;
                bool found_break = false;
                for (size_t k = 0; k < NUM_TET_EDGES; k++)
                {
                    edge_t edge = edge_list[k];

                    if (tet_store.edge_store.get(edge).needs_refining == 1)
                    {
                        returned_nodes[0] = edge.first();
                        returned_nodes[1] = edge.second();

                        trace_out << "1:2 needs to be split on " <<
                            returned_nodes[0] << " and " <<
                            returned_nodes[1] << std::endl;

                        found_break = true;
                        break;
                    }
                }

                assert(found_break);

                return returned_nodes;
            }

            void lock_intermediates(
                    tet_store_t& tet_store,
                    std::unordered_set<size_t> intermediate_list,
                    Edge_Lock_Case lock_case
                )
            {
                // Loop over all edges
                // If the edge is in the intermediate_list, deal with it
                for (const auto& p : tet_store.edge_store.edges)
                {
                    auto e = p.first;
                    size_t k1 = e.first();
                    size_t k2 = e.second();
                    // Can we make this double search cheaper?
                    if (
                            (intermediate_list.count(k1)) ||
                            (intermediate_list.count(k2))
                       )
                    {
                        trace_out << "Locking intermediate " << e << " from " << k1 << " and " << k2 << std::endl;
                        tet_store.edge_store.get(e).lock_case = lock_case;
                        tet_store.edge_store.get(e).needs_refining = 0;
                    }
                }

            }

            // TODO: remove this, it's horrible and not efficient.
            // WARNING: THIS GOES OVER ALL TETS!!!!
            void lock_edges_from_node(
                    tet_store_t& tet_store,
                    size_t node_id,
                    Edge_Lock_Case lock_case
            )
            {
                // Iterate over edges of ALL tet
                for (const auto& kv : tet_store.tets)
                {
                    size_t tet_id = kv.first;
                    edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);
                    for (size_t k = 0; k < NUM_TET_EDGES; k++)
                    {
                        // If it contains that node id, mark it using lock_case
                        edge_t edge = edge_list[k];

                        size_t edge_node_A_id = edge.first();
                        size_t edge_node_B_id = edge.second();

                        if ((edge_node_A_id == node_id) || (edge_node_B_id == node_id)) {
                            trace_out << " found node in " << edge_node_A_id << " - " << edge_node_B_id << " set to " << lock_case << std::endl;
                            tet_store.edge_store.get(edge).lock_case = lock_case;
                            tet_store.edge_store.get(edge).needs_refining = 0;
                        }
                    }
                }
            }
//            void lock_edges_from_node(
//                    tet_store_t& tet_store,
//                    size_t tet_id,
//                    size_t node_id,
//                    Edge_Lock_Case lock_case
//            )
//            {
//                // Iterate over edges of of tet
//                edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);
//                for (size_t k = 0; k < NUM_TET_EDGES; k++)
//                {
//                    // If it contains that node id, mark it using lock_case
//                    edge_t edge = edge_list[k];
//
//                    size_t edge_node_A_id = edge.first();
//                    size_t edge_node_B_id = edge.second();
//
//                    if ((edge_node_A_id == node_id) || (edge_node_B_id == node_id)) {
//                        tet_store.edge_store.get(edge).lock_case = lock_case;
//                    }
//                }
//            }


            ///// DEREFINEMENT STARTS HERE /////
            /**
             * @brief Function to iterate over children and remove them
             *
             * @param tet_store Tet store to use
             * @param parent_id Id of the parent for whom you will delete the
             * children
             */
            void derefine_children(tet_store_t& tet_store, size_t parent_id)
            {
                // For a given tet_id, find and delete its children
                Refinement_State& parent = tet_store.data(parent_id);
                for (auto c : parent.children)
                {
                    tet_store.erase(c);
                    //tet_store.deactivate(c);

                    /*
                    auto children = tet_store.data(c).children;
                    // Debug printing
                    std::cout << "tet " << c << "has ";
                    for (auto child : children)
                    {
                        std::cout << " _child " << child;
                    }
                    std::cout << std::endl;
                    */
                }
                parent.children.clear();
            }

            /**
             * @brief Common code for derefinement. Deactives the children and
             * actives the parent
             *
             * @param tet_store Tet store to use
             * @param parent_id The id of the parent
             */
            void generic_derefine(tet_store_t& tet_store, size_t parent_id)
            {
                derefine_children(tet_store,parent_id);
                tet_store.activate(parent_id);
            }

            /**
             * @brief Perform 2->1 derefinement on tet
             *
             * @param tet_store Tet store to use
             * @param parent_id The id of the parent
             */
            void derefine_two_to_one(tet_store_t& tet_store, node_connectivity_t&, size_t parent_id)
            {
                //if (!check_allowed_derefinement(tet_store,parent_id)) return;
                // build a delete-list of edges/intermediates first, mesh_adapter
                // deletes edges from this list later
                determine_deletelist_of_intermediates(tet_store, parent_id);
                generic_derefine(tet_store,parent_id);
            }

            /**
             * @brief Perform 4->1 derefinement on tet
             *
             * @param tet_store Tet store to use
             * @param parent_id The id of the parent
             */
            void derefine_four_to_one(tet_store_t& tet_store, node_connectivity_t&, size_t parent_id)
            {
                //if (!check_allowed_derefinement(tet_store,parent_id)) return;
                // build a delete-list of edges/intermediates first, mesh_adapter
                // deletes edges from this list later
                determine_deletelist_of_intermediates(tet_store, parent_id);
                generic_derefine(tet_store,parent_id);
            }

            /**
             * @brief Perform 8->1 derefinement on tet
             *
             * @param tet_store Tet store to use
             * @param parent_id The id of the parent
             */
            void derefine_eight_to_one(tet_store_t& tet_store, node_connectivity_t&, size_t parent_id)
            {
                //if (!check_allowed_derefinement(tet_store,parent_id)) return;

                generic_derefine(tet_store,parent_id);

                // TODO: Do we delete the nodes? Do we even have nodes?

                // Delete the center edges
                    // If edge isn't in the parent, delete it? Is there a better way?
                edge_list_t parent_edges = tet_store.generate_edge_keys(parent_id);

                Refinement_State& parent = tet_store.data(parent_id);
                for (auto c : parent.children)
                {
                    edge_list_t child_edges = tet_store.generate_edge_keys(c);
                    // build a delete-list of non-matching edges first, then
                    // mesh_adapter deletes edges from this list later
                    determine_deletelist_of_non_matching_edges(child_edges, parent_edges);
                }
            }

            // TODO: Document This.
            void derefine_four_to_two(tet_store_t& tet_store, node_connectivity_t& node_connectivity, size_t parent_id)
            {
                //if (!check_allowed_derefinement(tet_store,parent_id)) return;
                auto edge = find_edge_not_derefined(tet_store,
                  node_connectivity, parent_id);
                derefine_four_to_one(tet_store, node_connectivity, parent_id);
                refine_one_to_two( tet_store, node_connectivity, parent_id,
                  edge.first(), edge.second() );
            }

            // TODO: Document This.
            void derefine_eight_to_two(tet_store_t& tet_store, node_connectivity_t& node_connectivity, size_t parent_id)
            {
                //if (!check_allowed_derefinement(tet_store,parent_id)) return;

                auto edge = find_edge_not_derefined(tet_store,
                  node_connectivity, parent_id);
                derefine_eight_to_one(tet_store, node_connectivity, parent_id);
                refine_one_to_two( tet_store, node_connectivity, parent_id,
                  edge.first(), edge.second() );
            }

            // TODO: Document This.
            void derefine_eight_to_four(tet_store_t& tet_store, node_connectivity_t& node_connectivity, size_t parent_id)
            {
                //if (!check_allowed_derefinement(tet_store,parent_id)) return;
                // TODO: think about if the logic for these derefs are right
                derefine_eight_to_one(tet_store, node_connectivity, parent_id);
                deref_refine_one_to_four( tet_store, node_connectivity, parent_id);
            }

            /**
             * @brief Loop over children and determine delete-list of all intermediate edges
             *
             * @param tet_store Tet store to use
             * @param parent_id Id of parent
             */
            void determine_deletelist_of_intermediates(tet_store_t& tet_store, size_t parent_id)
            {
                Refinement_State& parent = tet_store.data(parent_id);
                auto parent_edges = tet_store.generate_edge_keys(parent_id);
                std::set< edge_t > parent_edge_set;
                for (auto pe:parent_edges) parent_edge_set.insert(pe);
                for (auto c : parent.children)
                {
                    edge_list_t edge_list = tet_store.generate_edge_keys(c);
                    for (size_t k = 0; k < NUM_TET_EDGES; k++)
                    {
                        edge_t edge = edge_list[k];
                        // accept this code may try delete an edge which has already gone
                        if (tet_store.edge_store.exists(edge)) {
                            if (parent_edge_set.count(edge) == 0)
                            {
                                trace_out << "child " << c << " adding to delete list: "
                                  << edge.first() << " - " << edge.second() << std::endl;
                                delete_list.push_back(edge);
                            }
                        }
                    }
                }
            }

            /**
             * @brief Deletes the intermediate edge in the delete list for derefinement
             *
             * @param tet_store Tet store to use
             */
            void delete_intermediates_of_children(tet_store_t& tet_store)
            {
              for (const auto& edge : delete_list) {
                tet_store.edge_store.erase(edge);
              }

              delete_list.clear();
            }

            /**
             * @brief If edge in candidate is not present in basis, add edge
             * (candidate) to delete list
             *
             * @param candidate The edge list which is to be searched and deleted
             * @param basis The edge list to check against
             */
            void determine_deletelist_of_non_matching_edges(edge_list_t candidate,
              edge_list_t basis)
            {
                trace_out << "Looking for edges to delete" << std::endl;

                // TODO: Sanity check this now we changed to edge_t

                // Loop over the edges in each child. Look over the basis and
                // if we can't find it, delete it
                for (size_t k = 0; k < NUM_TET_EDGES; k++)
                {
                    edge_t search_key = candidate[k];

                    // Search the basis for it
                    bool found_it = false;

                    for (size_t l = 0; l < NUM_TET_EDGES; l++)
                    {
                        edge_t key = basis[l];
                        if (search_key == key)
                        {
                            found_it = true;
                        }
                    }

                    // If we didn't find it, delete it
                    if (!found_it)
                    {
                        // Delete it
                        //tet_store.edge_store.erase(search_key);
                        delete_list.push_back(search_key);
                    }
                }
            }


            /**
             * @brief function to detect when an invalid derefinement is
             * invoked
             *
             * @param tet_store Tet store to use
             * @param tet_id Id the of the tet which will be de-refined
             *
             * @return A bool stating if the tet can be validly de-refined
             */
            bool check_allowed_derefinement( tet_store_t& tet_store, size_t tet_id)
            {
                Refinement_State& master_element = tet_store.data(tet_id);

                // Check this won't take us past the max refinement level
                if (master_element.refinement_level <= MIN_REFINEMENT_LEVEL)
                {
                    return false;
                }

                // If we got here, we didn't detect anything which tells us not
                // to
                return true;
            }


            // HERE BE DRAGONS! THIS IS DANGEROUS IF YOU USE IT WRONG
            // For every child of parent_id, set his children to our won
            // TODO: set a flag for the curious user to know we trashed the children
            void overwrite_children(
                    tet_store_t& tet_store,
                    const child_id_list_t& to_be_replaced,
                    const child_id_list_t& replace_with
            )
            {
                for (auto c : to_be_replaced)
                {
                    tet_store.data(c).children = replace_with;
                }
            }

            /**
             * @brief function to detect which edge should not get derefined
             *
             * @param tet_store Tet store to use
             * @param node_connectivity Node connectivity to use
             * @param tet_id Id the of the tet which will be de-refined
             *
             * @return Array of size two containing nodes of required edge
             */
            edge_t find_edge_not_derefined(
              tet_store_t& tet_store,
              node_connectivity_t& node_connectivity,
              size_t tet_id)
            {
              // 2 nonparent nodes set to derefine for a 4:2
              // 5 nonparent nodes set to derefine for an 8:2
              // will have 2 or 5 nonparent nodes set to deref; Figure out which
              // edge is the one that is not set to deref

              auto derefine_node_set = find_derefine_node_set(tet_store, tet_id);

              //// Do number of points
              //std::unordered_set<size_t> derefine_node_set;

              // Find the set of nodes which are not in the parent
              std::unordered_set<size_t> non_parent_nodes =
                child_exclusive_nodes(tet_store, tet_id);

              // from the above non_parent_nodes set and derefine_node_set,
              // figureout which node should be removed
              std::size_t ed_A(0), ed_B(0);
              for (auto npn:non_parent_nodes) {
                if (derefine_node_set.count(npn) == 0) {
                  // we've found the node that should not be removed, now we
                  // need to find the edge it belongs to
                  auto nonderef_edge = node_connectivity.get(npn);
                  ed_A = nonderef_edge[0];
                  ed_B = nonderef_edge[1];
                  //std::cout << "do-not-deref-APAN " << "A " << nd_edge[0]
                  //        << " B " << nd_edge[1] << std::endl;
                }
              }

              assert(ed_A!=ed_B);
              edge_t nd_edge(ed_A, ed_B);
              return nd_edge;
            }

            /**
             * @brief function to detect what intermediate/non-parent nodes are
             *   marked for derefinement
             *
             * @param tet_store Tet store to use
             * @param tet_id Id of the tet which will be de-refined
             *
             * @return Set of nodes of marked for derefinement
             */
            std::unordered_set< size_t > find_derefine_node_set(
              tet_store_t& tet_store,
              size_t tet_id)
            {
              // Set of nodes which are not in the parent
              std::unordered_set<size_t> non_parent_nodes =
                child_exclusive_nodes(tet_store, tet_id);
              std::unordered_set<size_t> derefine_node_set, unmarked_deref_node_set,
                final_deref_node_set;

              child_id_list_t children = tet_store.data(tet_id).children;

              // Look at children
              trace_out << tet_id << " Looping over " << children.size() << "children" << std::endl;
              for (size_t i = 0; i < children.size(); i++)
              {
                  trace_out << "child: " << children[i] << std::endl;
                  // TODO: Is this in element or tet ids?
                  edge_list_t edge_list = tet_store.generate_edge_keys(children[i]);
                  for (size_t k = 0; k < NUM_TET_EDGES; k++)
                  {
                      edge_t edge = edge_list[k];
                      // TODO: where do we makr the edges that need to be derefed? parent of child?
                      // Check each node, see if its an intermediate
                      size_t A = edge.first();
                      size_t B = edge.second();
                      trace_out << "checking edge for deref " << A << " - " << B << std::endl;

                      //if (tet_store.is_intermediate(A))
                      if (non_parent_nodes.count(A) )
                      {
                        if (tet_store.edge_store.get(edge).needs_derefining) {
                          trace_out << "Adding " << A << std::endl;
                          derefine_node_set.insert(A);
                        }
                        else {
                          unmarked_deref_node_set.insert(A);
                          //trace_out << "NOT added " << A << std::endl;
                        }
                      }

                      //if (tet_store.is_intermediate(B))
                      if (non_parent_nodes.count(B))
                      {
                        if (tet_store.edge_store.get(edge).needs_derefining) {
                          trace_out << "Adding " << B << std::endl;
                          derefine_node_set.insert(B);
                        }
                        else {
                          unmarked_deref_node_set.insert(B);
                          //trace_out << "NOT added " << B << std::endl;
                        }
                      }
                  }
              }

              //trace_out << "marked for deref: " << derefine_node_set.size() << std::endl;
              //trace_out << "NOT marked for deref: " << unmarked_deref_node_set.size() << std::endl;

              // remove nodes that are unmarked for derefinement
              for (auto drnode : derefine_node_set) {
                if (unmarked_deref_node_set.count(drnode) == 0)
                  final_deref_node_set.insert(drnode);
              }
              derefine_node_set = final_deref_node_set;
              return derefine_node_set;
            }


            std::unordered_set<size_t> child_exclusive_nodes(tet_store_t& tet_store,
              size_t tet_id)
            {
              std::unordered_set<size_t> non_parent_nodes;

              // array
              auto parent_tet = tet_store.get(tet_id);

              // convert to set
              std::unordered_set<size_t> parent_set(begin(parent_tet), end(parent_tet));

              child_id_list_t children = tet_store.data(tet_id).children;
              for (size_t i = 0; i < children.size(); i++)
              {
                      auto child_tet = tet_store.get( children[i] );

                      // Look at nodes, if not present add to set
                      for (std::size_t j = 0; j < NUM_TET_NODES; j++)
                      {
                              auto node = child_tet[j];
                              if (parent_set.count(node) == 0)
                              {
                                      non_parent_nodes.insert(node);
                              }
                      }
              }

              trace_out <<" Found " << non_parent_nodes.size() << " non parent nodes " << std::endl;
              return non_parent_nodes;

            }
    };
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // guard
