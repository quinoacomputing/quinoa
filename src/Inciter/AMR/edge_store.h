#ifndef AMR_edge_store_h
#define AMR_edge_store_h

#include "Base/Exception.h"

namespace AMR {

    class edge_store_t {
        public:
            edges_t edges;

            size_t size()
            {
                return edges.size();
            }

            /**
             * @brief Function to create new edge between two nodes with an
             * intermediate. Given nodes A, B, and AB makes edge A->AB and AB->B
             *
             * @param A First end node
             * @param B Second end node
             * @param AB Intermediate node
             * @param lockcase Lock case for the new edges
             */
            void split(size_t A, size_t B, size_t AB, Edge_Lock_Case lockcase)
            {
                generate(A, AB, lockcase);
                generate(B, AB, lockcase);

                //children.insert( std::pair<edge_t, size_t>(edge_t(A,B), AB));
                // Generate pertinent keys
                //edge_t keyAB = nodes_to_key(A, B);

                // NOTE: This isn't explicitly needed in the paper, and may be
                    // implicitly dealt with somewhere?
                //mark_edge_for_refinement(keyAB);
            }

            /**
             * @brief Given nodes A and B, generate an edge between them
             *
             * @param A First node
             * @param B Second node
             * @param lockcase Lock case for new edge
             */
            void generate(size_t A, size_t B, Edge_Lock_Case lockcase)
            {
                Assert(A != B, "Trying to add edge between duplicate IDs");
                // Generate key
                edge_t keyAB = nodes_to_key(A, B);
                //Create refined edge
                Edge_Refinement edgeAB = Edge_Refinement(A, B, 0.00, false,
                        false, false, lockcase);
                // Add edge to store
                add(keyAB, edgeAB);
            }

            bool exists(edge_t key)
            {
                if (edges.find(key) != edges.end())
                {
                    return true;
                }
                return false;
            }

            /**
             * @brief Function to retrieve an edge from the edge store
             *
             * @param key Key of the edge to get
             *
             * @return A reference to the fetched edge
             */
            Edge_Refinement& get(edge_t key)
            {
                Assert( exists(key), "Key does not exist" );
                return edges[key];
            }

            Edge_Lock_Case lock_case(edge_t key)
            {
                return get(key).lock_case;
            }

            void erase(edge_t key)
            {
                edges.erase(key);
            }

            /**
             * @brief Function to add edge to ede store
             *
             * @param key The key for the given edge
             * @param e The edge data
             *
             * Note: This tolerate the addition of duplicate edges
             */
            void add(edge_t key, Edge_Refinement e)
            {
                // Add edge if it doesn't exist (default behavior of insert)
                edges.insert( std::pair<edge_t, Edge_Refinement>(key, e));

                // TODO: It may be worth adding a check here to ensure if we're
                // trying to add a new edge that exists it should contain the
                // same data
            }

            static edge_t nodes_to_key(size_t A, size_t B)
            {
                return edge_t(A,B);
            }

            /**
             * @brief Function to build a  string key from two node ids
             * NOTE: Regardless of order of arguments, the same key will be generated
             */
            //static std::string nodes_to_key(size_t A, size_t B)
            //{
                //return std::to_string(std::min(A,B)) + KEY_DELIM + std::to_string(std::max(A,B));
            //}

            /**
             * @brief function to take the nodes representing a face
             * and to build the possible edges based on that
             *
             * For a given face {ABC}, generate the edge pairs {AB, AC, BC}
             *
             * @param face_ids The ids of the face to generate this for
             *
             * @return A (partially filled) list of all edges present on the
             * face
             */
            // FIXME: Is it OK that it leaves some of the array blank?
            static edge_list_t generate_keys_from_face_ids(face_ids_t face_ids)
            {
                edge_list_t key_list;
                size_t A = face_ids[0];
                size_t B = face_ids[1];
                size_t C = face_ids[2];

                // TODO: Investigate if AB BC CA is a much better jacobian ordering
                edge_t key = nodes_to_key(A,B);
                key_list[0] = key; // TODO: Is it OK to use copy assignment here?

                key = nodes_to_key(A,C);
                key_list[1] = key;

                key = nodes_to_key(B,C);
                key_list[2] = key;

                return key_list;
            }

            /**
             * @brief function to take a list of edge and mark them all
             * as needing to be refined
             *
             * @param ids List of ids to mark for refinement
             */
            void mark_edges_for_refinement(std::vector<node_pair_t> ids) {
                for (const auto& id : ids)
                {
                    edge_t key = nodes_to_key(id[0], id[1]);

                    mark_for_refinement(key);
                }
            }


            /**
             * @brief function to mark a single edge as needing
             * refinement (provides a nice abstraction from messing with the
             * struct directly).
             *
             * @param key The edge key to mark as refinement
             */
            void mark_for_refinement(edge_t key)
            {
                    Assert( exists(key), "Key does not exist" );
                    get(key).needs_refining = true;
            }

            /**
             * @brief Function to unmark and edge as needing refinement
             *
             * @param key The key representing the edge to unmark
             */
            void unmark_for_refinement(edge_t key)
            {
                    Assert( exists(key), "Key does not exist");
                    get(key).needs_refining = false;
            }

            // TODO: Document this (and implement!)
            void mark_edges_for_derefinement(std::vector<node_pair_t> ids) {
                for (const auto& id : ids)
                {
                    edge_t key = nodes_to_key(id[0], id[1]);

                    mark_edge_for_derefinement(key);
                }
            }
            void mark_edge_for_derefinement(edge_t key) {
                    get(key).needs_derefining = true;
            }


            /**
             * @brief Function to generate a list of edge keys from a tet
             *
             * @param tet The tet to generate edge pairs for
             *
             * @return A list (array) of edge keys which can be separated out to
             * name the two composing node ids
             */
            // TODO: Should this return a pointer/reference?
            edge_list_t generate_keys(tet_t tet)
            {
                // NOTE: Generate these with a (2d) loop and not hard code them?
                edge_list_t key_list;

                size_t A = tet[0];
                size_t B = tet[1];
                size_t C = tet[2];
                size_t D = tet[3];

                edge_t key;

                key = nodes_to_key(A,B);
                key_list[0] = key;

                key = nodes_to_key(A,C);
                key_list[1] = key;

                key = nodes_to_key(A,D);
                key_list[2] = key;

                key = nodes_to_key(B,C);
                key_list[3] = key;

                key = nodes_to_key(B,D);
                key_list[4] = key;

                key = nodes_to_key(C,D);
                key_list[5] = key;

                return key_list;
            }

            /**
             * @brief Helper debug function to print edge information
             */
            void print() {
                for (const auto& kv : edges)
                {
                    std::cout << "edge " << kv.first << " between " <<
                        kv.second.A << " and " << kv.second.B << " val " <<
                        kv.second.refinement_criteria <<
                    std::endl;
                }
            }

            /*
            // This function is not needed, as we're no longer required to
            // update ids inplace
            void replace(size_t old_id, size_t new_id)
            {
                for (const auto& kv : edges)
                {
                    // Find all slots which have the id in the key
                    if (false) //  TODO: Implement this
                    {
                        // Cache value
                        auto value = kv.second;
                        edge_t key = kv.first;

                        // Find the bit of the old key we need to keep
                        size_t remainder_key = 0; // TODO:

                        // Delete them
                        erase(key);

                        // Build new key
                        edge_t new_key = nodes_to_key(new_id, remainder_key);

                        // Re add with new key
                        add(new_key, value);
                    }

                }
            }
            */

            /*
            // TODO: Document this
            // If this returns 0, it couldn't find it. 0 is a really unlikely
            // added node id..
            size_t find_intermediate_nodes(size_t A, size_t B)
            {
                size_t id = 0;
                edge_t key = nodes_to_key(A, B);

                if (children.find(key) != children.end())
                {
                    id = children[key];
                }

                return id;
            }
            */


    };
}

#endif // guard
