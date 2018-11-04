#ifndef AMR_refinement_h
#define AMR_refinement_h

#include <algorithm>

#include "Macro.h"
#include "tet_store.h"
#include "node_connectivity.h"

// TODO: make this have a base class to support multiple generator schemes
// using the policy design pattern

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

namespace AMR {

    class refinement_t {
        private:

            const size_t DEFAULT_REFINEMENT_LEVEL = 0; //TODO: Is this in the right place?
            const size_t MIN_REFINEMENT_LEVEL = DEFAULT_REFINEMENT_LEVEL;

        public:

            const size_t MAX_REFINEMENT_LEVEL = 4;

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
             * @param tet_id The id to refine 1:2
             */
            void refine_one_to_two( tet_store_t& tet_store, node_connectivity_t& node_connectivity, size_t tet_id)
            {
                edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);
                node_pair_t nodes = find_single_refinement_nodes(tet_store,edge_list);
                refine_one_to_two( tet_store, node_connectivity, tet_id, nodes[0], nodes[1]);
            }

            /**
             * @brief Method which takes a tet id, and transforms arguments
             * into the form needed for the main 1:2 refinement method
             *
             * @param tet_id The id to refine 1:2
             */
/*
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
                tet_store.intermediate_list.insert(new_node_id);
            }

            /**
             * @brief Method which takes a tet id, and deduces the other
             * parameters needed to perform a 1:4
             *
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
                    for (size_t k = 0; k < NUM_FACE_NODES; k++)
                    {
                        edge_t edge = face_edge_list[k];
                        if (tet_store.edge_store.get(edge).needs_refining == true)
                        {
                            num_face_refine_edges++;
                        }

                        // Check for locked edges
                            // This case only cares about faces with no locks
                        if (tet_store.edge_store.lock_case(edge) != Edge_Lock_Case::unlocked)
                        {
                            // Abort this face
                            num_face_refine_edges = 0;
                            continue;
                        }
                    }
                    if (num_face_refine_edges >= 2)
                    {
                        assert(num_face_refine_edges < 4);
                        //face_refine = true;
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

                tet_store.intermediate_list.insert(AB);
                tet_store.intermediate_list.insert(AC);
                tet_store.intermediate_list.insert(BC);

            }

            /**
             * @brief Refine a given tet id into 8 children.
             * NOTE: Does not do any validity checking (currently?)
             *
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
                // Internal (for a face BDC, it's the intermediat and mid opposite)
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

                    if (tet_store.edge_store.get(edge).needs_refining == true)
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
                        tet_store.edge_store.get(e).lock_case = lock_case;
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
                        }
                    }
                }
            }
            void lock_edges_from_node(
                    tet_store_t& tet_store,
                    size_t tet_id,
                    size_t node_id,
                    Edge_Lock_Case lock_case
            )
            {
                // Iterate over edges of of tet
                edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);
                for (size_t k = 0; k < NUM_TET_EDGES; k++)
                {
                    // If it contains that node id, mark it using lock_case
                    edge_t edge = edge_list[k];

                    size_t edge_node_A_id = edge.first();
                    size_t edge_node_B_id = edge.second();

                    if ((edge_node_A_id == node_id) || (edge_node_B_id == node_id)) {
                        tet_store.edge_store.get(edge).lock_case = lock_case;
                    }
                }
            }


            ///// DEREFINEMENT STARTS HERE /////
            /**
             * @brief Function to iterate over children and remove them
             *
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
                    parent.num_children--; // Could directly set to 0
                }
                parent.children.clear();
            }

            /**
             * @brief Common code for derefinement. Deactives the children and
             * actives the parent
             *
             * @param parent_id The id of the parent
             */
            void generic_derefine(tet_store_t& tet_store, size_t parent_id)
            {
                derefine_children(tet_store,parent_id);
                tet_store.activate(parent_id);
            }

            // TODO: Document This.
            void derefine_two_to_one(tet_store_t& tet_store, size_t parent_id)
            {
                delete_intermediates_of_children( tet_store, parent_id);
                generic_derefine(tet_store,parent_id);
            }

            // TODO: Document This.
            void derefine_four_to_one(tet_store_t& tet_store, size_t parent_id)
            {
                delete_intermediates_of_children(tet_store, parent_id);
                generic_derefine(tet_store,parent_id);
            }

            // TODO: Document This.
            void derefine_eight_to_one(tet_store_t& tet_store, size_t parent_id)
            {
                generic_derefine(tet_store,parent_id);

                // Delete the center edges
                    // If edge isn't in the parent, delete it? Is there a better way?
                edge_list_t parent_edges = tet_store.generate_edge_keys(parent_id);

                Refinement_State& parent = tet_store.data(parent_id);
                for (auto c : parent.children)
                {
                    edge_list_t child_edges = tet_store.generate_edge_keys(c);
                    delete_non_matching_edges( tet_store, child_edges, parent_edges);
                }
            }

            // TODO: Document This.
            /*
            void derefine_four_to_two(size_t parent_id)
            {
                assert(0);
            }

            // TODO: Document This.
            void derefine_eight_to_two(size_t parent_id)
            {
                assert(0);
            }

            // TODO: Document This.
            void derefine_eight_to_four(size_t parent_id)
            {
                assert(0);
            }
            */

            /**
             * @brief Loop over children and delete all intermediate edges
             *
             * @param parent_id Id of parent
             */
            void delete_intermediates_of_children( tet_store_t& tet_store, size_t parent_id)
            {
                Refinement_State& parent = tet_store.data(parent_id);
                for (auto c : parent.children)
                {
                    delete_intermediates(tet_store,c);
                }
            }

            // TODO: Document this
            void delete_intermediates( tet_store_t& tet_store, size_t tet_id)
            {
                edge_list_t edge_list = tet_store.generate_edge_keys(tet_id);
                for (size_t k = 0; k < NUM_TET_EDGES; k++)
                {
                    edge_t edge = edge_list[k];
                    // accept this code may try delete an edge which has already gone
                    if (tet_store.edge_store.exists(edge)) {
                        if (tet_store.edge_store.get(edge).lock_case == Edge_Lock_Case::intermediate)
                        {
                            tet_store.edge_store.erase(edge);
                        }
                    }
                }
            }

            /**
             * @brief If edge in candidate is not present in basis, delete the
             * edge (candidate) from the main edge store
             *
             * @param candidate The edge list which is to be searched and deleted
             * @param basis The edge list to check against
             */
            void delete_non_matching_edges( tet_store_t& tet_store, edge_list_t candidate, edge_list_t basis)
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
                        tet_store.edge_store.erase(search_key);
                    }
                }
            }


            /**
             * @brief function to detect when an invalid derefinement is
             * invoked
             *
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



    };
}

#if defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // guard
