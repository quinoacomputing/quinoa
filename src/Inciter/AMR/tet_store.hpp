#ifndef AMR_tet_store_h
#define AMR_tet_store_h

#include <unordered_set>
#include <vector>

#include "AMR_types.hpp"
#include "active_element_store.hpp"
#include "master_element_store.hpp"
#include "marked_refinements_store.hpp"
#include "edge_store.hpp"
#include "util.hpp"
#include "id_generator.hpp"

namespace AMR {

    class tet_store_t {
        public:
            // FIXME: Remove this (center_tets) data structure!
                // This is a horrendous code abuse, and I'm sorry. I'm fairly
                // certain we'll be re-writing how this detection is done and just
                // wanted a quick-fix so I could move on :(
            std::set<size_t> center_tets; // Store for 1:4 centers

            std::set<size_t> delete_list; // For marking deletions in deref

            AMR::active_element_store_t active_elements;
            AMR::master_element_store_t master_elements;

            std::vector< std::size_t > active_tetinpoel;
            std::set< std::size_t > active_nodes;

            AMR::id_generator_t id_generator;

            std::unordered_set<size_t> intermediate_list;
            // Public so it can be trivially grabbed for looping over.
            std::vector< std::size_t > active_id_mapping;

            tet_list_t tets;
            AMR::edge_store_t edge_store;

            // TODO: Make this (and others) private at some point
            AMR::marked_refinements_store_t<AMR::Refinement_Case> marked_refinements;
            AMR::marked_refinements_store_t<AMR::Derefinement_Case> marked_derefinements;

            /**
             * @brief function to return the number of tets stored
             *
             * @return Num of tets
             */
            size_t size() {
                return tets.size();
            }

            /**
             * @brief Helper to check if an given tet is active
             *
             * @param id id of the tool to check
             *
             * @return active status of tet
             */
            bool is_active(size_t id)
            {
                return active_elements.exists(id);
            }

            /**
             * @brief Return refinement case for a given id
             *
             * @param id id to get case for
             *
             * @return Refinement case for id
             */
            Refinement_Case get_refinement_case(size_t id)
            {
                return data(id).refinement_case;
            }

            /**
             * @brief Set value of normal
             *
             * @param id Id to set
             * @param val true/false to set
             */
            void set_normal(size_t id, bool val)
            {
                data(id).normal = val;
            }

            /**
             * @brief Get normal value for given id
             *
             * @param id Id of tet to check
             *
             * @return true/false of normal value
             */
            bool is_normal(size_t id)
            {
                // should the underlying type be bool?
                return data(id).normal;
            }

            /**
             * @brief set a tet as normal
             *
             * @param id id to set
             */
            void mark_normal(size_t id)
            {
                set_normal(id, true);
            }

            /**
             * @brief get data for a tet from master element
             *
             * @param id id of tet to get
             *
             * @return state of the tet
             */
            Refinement_State& data(size_t id)
            {
                return master_elements.get(id);
            }

            const Refinement_State& data(size_t id) const
            {
                return master_elements.get(id);
            }

            /**
             * @brief Method to insert tet into the tet store, so the
             * underlying data structure doesn't have to be interfaced with
             * directly
             *
             * @param id Id of the added tet
             * @param t The tet element
             */
            void insert(size_t id, tet_t t)
            {
                // cppcheck-suppress assertWithSideEffect
                assert( !exists(id) );
                tets.insert( std::pair<size_t, tet_t>(id, t));
            }

            /**
             * @brief Getter for tet element
             *
             * @param id Id of tet to get
             *
             * @return Copy of the tet
             */
            tet_t get( size_t id )
            {
                // cppcheck-suppress assertWithSideEffect
                assert( exists(id) );
                return tets.at(id);
            }

            /**
             * @brief Function to check if a tet exists. Useful for debugging
             * access to invalid tets, or trying to re-create a tet which
             * already exists
             *
             * @param id Id of the tet to check
             *
             * @return Bool stating if the tet already exists
             */
            bool exists(size_t id)
            {
                auto f = tets.find(id);
                if (f != tets.end())
                {
                    //trace_out << "tet " << id << " exists." << std::endl;
                    return true;
                }
                return false;
            }

            /**
             * @brief Function to store a tet from a list of nodes
             *
             * @param id The ID of the tetrahedron to insert
             * @param nodes The node ids which make up the tet
             */
            void store_tet(size_t id, tet_t nodes) {

                insert(id, nodes);

                // Sanity check the storage ids
                // (this is probably better in a function/2d loop)
                assert( nodes[0] != nodes[1] );
                assert( nodes[0] != nodes[2] );
                assert( nodes[0] != nodes[3] );
                assert( nodes[1] != nodes[2] );
                assert( nodes[1] != nodes[3] );
                assert( nodes[2] != nodes[3] );
            }

            /**
             * @brief Convenience function to store a tet without first building
             * a list
             *
             * @param id The ID of the tetrahedron to store
             * @param first First Node
             * @param second Second Node
             * @param third Third Node
             * @param forth Forth Node
             */
            void store_tet(
                    size_t id,
                    size_t first,
                    size_t second,
                    size_t third,
                    size_t forth
                    )
            {
                store_tet( id, { {first, second, third, forth} } );
            }

            void add(
                    size_t id,
                    const tet_t& nodes,
                    Refinement_Case refinement_case,
                    size_t parent_id
            )
            {

                //std::cout << "id " << id << " parent " << parent_id << std::endl;

                add_to_master(id, nodes, refinement_case, parent_id, true);

                master_elements.get(id).refinement_level =
                    master_elements.get(parent_id).refinement_level+1;

                // Deal with updating parent
                master_elements.add_child(parent_id, id);

                trace_out << "Added child " << id << std::endl;
            }

            /**
             * @brief Convenience function to add a tet to the master_elements
             * and active_elements store
             *
             * @param id The ID of the tetrahedron to add
             * @param nodes A list of the nodes which form th etet
             * @param refinement_case The refinement case which caused this tet
             * to be generated
             * @param parent_id The ID of the parent tetrahedron
             * @param has_parent True if element has a parent
            */
            void add_to_master(size_t id, const tet_t& nodes,
              Refinement_Case refinement_case, size_t parent_id=0,
              bool has_parent=false)
            {
                store_tet(id, nodes);

                size_t refinement_level = 0;

                // Add to master list
                master_elements.add(id, refinement_case, refinement_level, parent_id, has_parent);

                // The new master element should start as active
                active_elements.add(id);
            }

            /**
             * @brief Interface to add a tet from the original mesh (no parent)
            */
            void add(const tet_t& nodes, Refinement_Case refinement_case)
            {
                size_t id = id_generator.get_next_tet_id();
                add_to_master(id, nodes, refinement_case);
            }

            void add(
                    size_t id,
                    size_t first,
                    size_t second,
                    size_t third,
                    size_t forth,
                    Refinement_Case refinement_case,
                    size_t parent_id
                    )
            {
                return add(
                        id,
                        { {first, second, third, forth} },
                        refinement_case,
                        parent_id
                        );

            }

            // NOTE: this does *not* deal with edges
            /**
             * @brief Function to delete a tet from the tet store (useful in
             * derefinement)
             *
             * @param id id of the tet to delete
             */
            void erase(size_t id)
            {
                deactivate(id);
                master_elements.erase(id);
                tets.erase(id);
                // TODO: Should this update the number of children here rather than at the call site?
            }

            /**
             * @brief Function to remove a tet from the active tet list
             *
             * @param id The id of the tet to deactivate
             */
            void deactivate(size_t id) {
                active_elements.erase(id);
                // TODO: For safety, should we also mark it's edges as not
                // needing to be refined?
            }

            void activate(size_t id) {
                if (!is_active(id) )
                {
                    active_elements.add(id);
                }
            }

            /**
             * @brief Function to add a tet to a list which maintains what is a
             * center tet. (The need to maintain a list could be replaced by a
             * geometric check on the tet itself)
             *
             * @param id Id of the tet to add
             */
            void add_center(size_t id)
            {
                // cppcheck-suppress assertWithSideEffect
                assert( !is_center(id) );
                center_tets.insert(id);
            }

            /**
             * @brief function to check if a tet is a center tet in a 1:4
             *
             * @param id Id of the tet to check
             *
             * @return Bool stating if it's a center tet or not
             */
            bool is_center(size_t id)
            {
                if (center_tets.find(id) != center_tets.end())
                {
                    return true;
                }
                return false;
            }

            /**
             * @brief Function to get a list of refinement levels, useful for
             * vis
             *
             * @return Vector containing refinement levels of tets
             */
            std::vector< real_t > get_refinement_level_list() const
            {
                std::vector<real_t> refinement_level_list;

                for (const auto& kv : tets)
                {
                    size_t element_id = kv.first;
                    if (active_elements.exists( element_id  )) {
                        real_t val = static_cast< tk::real >(
                          master_elements.get(element_id).refinement_level );
                        refinement_level_list.push_back(val);
                    }
                }

                trace_out << "Made refinement level list of len " << refinement_level_list.size() << std::endl;
                return refinement_level_list;
            }

            /**
             * @brief Function to return a list of cell types, useful when
             * invoking the vis to do coloring by cell type
             *
             * @return Vector listening the types of cells
             */
            std::vector< real_t > get_cell_type_list() const
            {
                std::vector<real_t> cell_type_list;

                for (const auto& kv : tets)
                {
                    size_t element_id = kv.first;

                    if (active_elements.exists( element_id  )) {

                        real_t val = 0.0;

                        // Be a good citizen, make this enum human readable
                        switch (master_elements.get(element_id).refinement_case)
                        {
                            case Refinement_Case::one_to_two:
                                val = 2.0;
                                break;
                            case Refinement_Case::one_to_four:
                                val = 4.0;
                                break;
                            case Refinement_Case::one_to_eight:
                                val = 8.0;
                                break;
                            case Refinement_Case::initial_grid:
                                val = 1.0;
                                break;
                            // TODO: this will never actually happen, as a 2:8 currently views
                            // itself as a 1:8 (as it did a 2:1, and a 1:8)
                            case Refinement_Case::two_to_eight:
                                val = 2.8;
                                break;
                            case Refinement_Case::four_to_eight:
                                val = 4.8;
                                break;
                            case Refinement_Case::none:
                                val = 0.0;
                                break;
                        }
                        cell_type_list.push_back(val);
                    }
                }

                trace_out << "Made cell type list of len " << cell_type_list.size() << std::endl;
                return cell_type_list;
            }

            /**
             * @brief The function gives a way to go back from active_inpoel to
             * real AMR id
             *
             * @return A vector which hold the AMR ids of the active inpoel
             */
            std::vector< std::size_t >& get_active_id_mapping()
            {
                active_id_mapping.clear();

                for (const auto& kv : tets)
                {
                    size_t element_id = kv.first;

                    if (is_active( element_id )) {
                        active_id_mapping.push_back( element_id );
                    }
                }

                return active_id_mapping;
            }

            /**
             * @brief Function to extract only the active elements from tetinpeol
             *
             * @return List of only active elements
             */
            // TODO: need equiv for m_x/m_y/m_z? Otherwise there will be
            // useless nodes in m_x etc? (Although that's functionally fine)
            std::vector< std::size_t >& get_active_inpoel()
            {
                active_tetinpoel.clear();
                active_nodes.clear();

                for (const auto& kv : tets)
                {

                    size_t element_id = kv.first;
                    auto t = kv.second;

                    if (is_active( element_id )) {
                        active_tetinpoel.push_back( t[0] );
                        active_tetinpoel.push_back( t[1] );
                        active_tetinpoel.push_back( t[2] );
                        active_tetinpoel.push_back( t[3] );

                        active_nodes.insert( t[0] );
                        active_nodes.insert( t[1] );
                        active_nodes.insert( t[2] );
                        active_nodes.insert( t[3] );
                    }
                }

                return active_tetinpoel;
            }

            // TODO: These mark methods can probably be a single one to which a
                // Refinement_Case is passed, depending on how extra edges are marked
            /**
             * @brief Function to mark a given tet as needing a 1:2 refinement
             *
             * @param tet_id The tet to mark
             */
            void mark_one_to_two(size_t tet_id)
            {
                // TODO: If none of these methods need extra markings, then
                // change them into a single method
                trace_out << "Mark " << tet_id << " as 1:2" << std::endl;
                marked_refinements.add(tet_id, Refinement_Case::one_to_two);
            }

            /**
             * @brief Mark a given tet as needing a 1:4 refinement
             *
             * @param tet_id The tet to mark
             */
            void mark_one_to_four(size_t tet_id)
            {
                trace_out << "Mark " << tet_id << " as 1:4" << std::endl;
                marked_refinements.add(tet_id, Refinement_Case::one_to_four);
            }

            /**
             * @brief Mark a given tet as needing a 2:8 refinement
             *
             * @param tet_id id of tet to mark
             */
            void mark_two_to_eight(size_t tet_id)
            {
                trace_out << "Mark " << tet_id << " as 2:8" << std::endl;
                marked_refinements.add(tet_id, Refinement_Case::two_to_eight);
            }

            /**
             * @brief Mark a given tet as needing a 4:8 refinement
             *
             * @param tet_id id of tet to mark
             */
            void mark_four_to_eight(size_t tet_id)
            {
                trace_out << "Mark " << tet_id << " as 4:8" << std::endl;
                marked_refinements.add(tet_id, Refinement_Case::four_to_eight);
            }

            /**
             * @brief Function to mark a given tet as needing a 1:8 refinement
             *
             * @param tet_id The tet to mark
             */
            void mark_one_to_eight(size_t tet_id)
            {
                trace_out << "Mark " << tet_id << " as 1:8" << std::endl;
                marked_refinements.add(tet_id, Refinement_Case::one_to_eight);
            }

            bool has_refinement_decision(size_t id)
            {
                return marked_refinements.exists(id);
            }

            /**
             * @brief Helper debug function to print tet information
             */
            void print_tets() {
                for (const auto& kv : tets)
                {
                    tet_t tet = kv.second;

                    trace_out << "Tet " << kv.first << " has edges :" <<
                        tet[0] << ", " <<
                        tet[1] << ", " <<
                        tet[2] << ", " <<
                        tet[3] << ", " <<
                    std::endl;
                }
            }

            void print_edges()
            {
                edge_store.print();
            }

            void print_node_types()
            {

                int initial_grid = 0;
                int one_to_two = 0;
                int one_to_four = 0;
                int one_to_eight = 0;
                int other = 0;

                for (const auto& kv : tets)
                {
                    size_t tet_id = kv.first;
                    if (is_active(tet_id))
                    {
                        switch(get_refinement_case(tet_id))
                        {
                            case Refinement_Case::one_to_two:
                                one_to_two++;
                                break;
                            case Refinement_Case::one_to_four:
                                one_to_four++;
                                break;
                            case Refinement_Case::one_to_eight:
                                one_to_eight++;
                                break;
                            case Refinement_Case::initial_grid:
                                initial_grid++;
                                break;
                            case Refinement_Case::two_to_eight:
                                // Don't care (yet)
                                other++;
                                break;
                            case Refinement_Case::four_to_eight:
                                // Don't care (yet)
                                other++;
                                break;
                            case Refinement_Case::none:
                                // Don't care (yet)
                                other++;
                                break;
                        }

                    }
                }

                //std::cout << "Active Totals:" << std::endl;
                //std::cout << "  --> Initial = " << initial_grid << std::endl;
                //std::cout << "  --> 1:2 = " << one_to_two << std::endl;
                //std::cout << "  --> 1:4 = " << one_to_four << std::endl;
                //std::cout << "  --> 1:8 = " << one_to_eight << std::endl;
            }

            edge_list_t generate_edge_keys(size_t tet_id)
            {
                tet_t tet = get(tet_id);
                return edge_store.generate_keys(tet);
            }

            void generate_edges(size_t i) {
                //  For tet ABCD, edges are:
                //  AB, AC, AD, BC, BD, CD
                //
                edge_list_t edge_list = generate_edge_keys(i);

                for (size_t j = 0; j < NUM_TET_EDGES; j++)
                {
                    edge_t edge = edge_list[j];

                    size_t A = edge.first();
                    size_t B = edge.second();

                    Edge_Refinement er = Edge_Refinement(A, B, false,
                            false, Edge_Lock_Case::unlocked);

                    edge_store.add(edge, er);
                }
            }

            /**
             * @brief function to take a tet_id, finds it's nodes, and
             * expresses them as faces
             *
             * Take tet ABCD, generate faces {ABC, ABD, ACD, BCD}
             *
             * @param tet_id The tet to generate faces for
             *
             * @return A list of faces making this tet
             */
            face_list_t generate_face_lists(size_t tet_id)
            {
                // Hard code this for now...
                tet_t tet = get(tet_id);

                trace_out  << "Tet has nodes " <<
                    tet[0] << ", " <<
                    tet[1] << ", " <<
                    tet[2] << ", " <<
                    tet[3] << ", " <<
                    std::endl;

                face_list_t face_list;

                // ABC
                face_list[0][0] = tet[0];
                face_list[0][1] = tet[1];
                face_list[0][2] = tet[2];

                // ABD
                face_list[1][0] = tet[0];
                face_list[1][1] = tet[3];
                face_list[1][2] = tet[1];

                // ACD
                face_list[2][0] = tet[0];
                face_list[2][1] = tet[2];
                face_list[2][2] = tet[3];

                // BCD
                face_list[3][0] = tet[1];
                face_list[3][1] = tet[3];
                face_list[3][2] = tet[2];

                return face_list;
            }

            /**
             * @brief Function which marks all edges in a given tet as needing
             * to be refined
             *
             * @param tet_id ID of the tet to mark
             */
            void mark_edges_for_refinement(size_t tet_id)
            {
                edge_list_t edge_list = generate_edge_keys(tet_id);
                for (size_t k = 0; k < NUM_TET_EDGES; k++)
                {
                    edge_t edge = edge_list[k];
                    edge_store.mark_for_refinement(edge);
                    trace_out << "Marking edge " << edge << " for refine " << std::endl;
                }
            }

            /**
             * @brief Delete existing edges. Iterate over the tets, and add them to
             * the edge store.
             */
            // FIXME: Better name for this?
            void generate_edges() {

                // Go over tets, and generate all known edges
                edge_store.edges.clear();

                // Jump over tets
                for (const auto& kv : tets)
                {
                    size_t tet_id = kv.first;
                    generate_edges(tet_id);
                }
            }

            void unset_marked_children(size_t parent_id)
            {
                Refinement_State& parent = data(parent_id);
                for (auto c : parent.children)
                {
                    marked_refinements.erase(c);
                }
            }

            child_id_list_t generate_child_ids(size_t parent_id, size_t count = MAX_CHILDREN)
            {
                return id_generator.generate_child_ids(parent_id, count);
            }
            size_t get_child_id(size_t parent_id, size_t offset) const
            {
                return master_elements.get_child_id(parent_id, offset);
            }

            size_t get_parent_id(size_t id) const
            {
                return master_elements.get_parent(id);
            }

            // Deref
            void process_delete_list()
            {
                trace_out << "process_delete_list " << delete_list.size() << std::endl;
                size_t original_size = size();
                for(auto f : delete_list) {
                    erase(f);
                }

                size_t end_size = size();
                trace_out << "Deleted " << original_size-end_size << std::endl;

                delete_list.clear();
            }
            /**
             * @brief Function to mark a given tet as a specific Derefinement_Case
             *
             * @param tet_id The tet to mark
             * @param decision The Derefinement_case to set
             */
            void mark_derefinement_decision(size_t tet_id, AMR::Derefinement_Case decision)
            {
                trace_out << "MARKING_DEREF_DECISION" << std::endl;
                marked_derefinements.add(tet_id, decision);
            }
            bool has_derefinement_decision(size_t id)
            {
                return marked_derefinements.exists(id);
            }

    };
}

#endif // guard
