#ifndef AMR_node_connectivity_h
#define AMR_node_connectivity_h

#include <unordered_map>
#include <UnsMesh.hpp>

namespace AMR {

    /**
     * @brief This class stores the connectivity of the node. Simply what this
     * means is that it just a mapping of node ids to a unique index.  The
     * value of the map is the index, and the key is the two nodes the new node
     * joins
     */
    class node_connectivity_t {

        private:
            //using Hash = tk::UnsMesh::Hash;
            //using Eq = tk::UnsMesh::Eq;

            using node_list_key_t = node_pair_t;
            using node_list_value_t = size_t;
            using node_list_t = std::unordered_map<node_list_key_t, node_list_value_t,  tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2>>;
            using inv_node_list_t = std::unordered_map<node_list_value_t, node_list_key_t>;

            node_list_t nodes;
            inv_node_list_t inv_nodes;

        public:

            size_t empty_node_count = 0;


            node_connectivity_t() { } // default cons

            //! Non-const-ref accessor to state
            //! \return  All node pairs
            node_list_t& data() { return nodes; }
            inv_node_list_t& inv_data() { return inv_nodes; }

            /**
             * @brief Method to add initial nodes to the store
             *
             * @param initial_size Size of the list to fill to
             */
            explicit node_connectivity_t(size_t initial_size)
            {
                for (size_t i = 0; i < initial_size; i++)
                {
                    // These can initially be 0 as initial nodes don't join any
                    // two others.. this could be updated to track
                    // intermediates, but this currently tracks "added" nodes
                    // nicely
                    add(0,0);
                }
            }

            /**
             * @brief Return size of node container -- the number of nodes
             *
             * @return Number of nodes
             */
            size_t size()
            {
                return nodes.size();
            }

            /**
             * @brief Getter into node storage *VALUE*
             *
             * @param id VALUE of the node to get
             *
             * @return The node_pair at the given id
             */
            node_pair_t get(size_t id)
            {
                //trace_out << "PROBLEM FINDING ID " << id << std::endl;

                // Ban getting of a node whos parents are {0,0}
                //assert(id > empty_node_count-1); //[0..empty_node_counts)

                auto iter = inv_nodes.find(id);

                /* old linear search code
                auto it = nodes.begin();
                for (; it != nodes.end(); ++it) {
                    if (it->second == id) break;
                }
                */

                //assert(iter != inv_nodes.end());
                //return iter->second;
                return (iter != inv_nodes.end() ? iter->second : node_pair_t{{id,id}});
            }

            /**
             * @brief function to calculate which node is opposite a
             * tet face
             *
             * @param face_list A list of faces on the tet
             * @param opposite_index The index for the face you want to know
             * the opposite node for
             *
             * @return An index (0-3) to tell you if A, B, C, or D
             * (respectively) is opposite the given face_list_t
             *
             * This function is tightly coupled (too coupled) to generate_face_lists
             *
             * generate_face_lists generates the faces {ABC, ABD, ACD, BCD} in
             * a fixed order. Opposite_index says the face from a face list we
             * care about. This function returns a number in the range {0,3} to
             * tell you  which node is missing from that face.
             *
             * I.e If opposite_index is 1, Node C is missing => 2.
             */
            static size_t face_list_opposite(face_list_t face_list, size_t opposite_index)
            {
                // FIXME: make this actually inspect the face_list and be much
                    // more robust...
                size_t result = face_list[0][0];
                switch(opposite_index)
                {
                    case 0:  // ABC
                        result = 3;
                        break;
                    case 1:  // ABD
                        result = 2;
                        break;
                    case 2:  // ACD
                        result = 1;
                        break;
                    case 3:  // BCD
                        result = 0;
                        break;
                    default: // something went horribly wrong..
                        assert(0);
                        break;
                }
                return result;
            }

            /**
             * @brief Add connectivity, unless it already exists
             *
             * @param A First node
             * @param B Second node
             *
             * @return Id/unique identifier of the node
             */
            node_list_value_t add(size_t A, size_t B)
            {
                if ((A == 0) && (B == 0))
                {
                    trace_out << "empty nodes = " << empty_node_count << std::endl;
                    node_list_value_t value = nodes.size() + empty_node_count;
                    empty_node_count++;
                    return value;
                }
                else {
                    assert(A != B);
                }

                node_list_key_t key = {{std::min(A,B), std::max(A,B)}};
                auto iter = nodes.find(key);

                trace_out << "A " << A << " B " << B << std::endl;

                // return the corresponding value if we find the key in the map
                if(iter != nodes.end()) {
                    trace_out << "Reuse " << iter->second << std::endl;
                    return iter->second;
                }
                else {
                    // if not in map
                    node_list_value_t value = nodes.size() + empty_node_count;
                    nodes[key] = value;
                    inv_nodes[value] = key;
                    trace_out << "Made new node " << value << std::endl;
                    return value;
                }
            }

            /**
             * @brief Print connectivity a id: a-b
             */
            void print()
            {
//                 std::cout << "Connectivity" << std::endl;
//                 for (size_t i = 0; i < size(); i ++)
//                 {
//                     std::cout << i << ": A " << get(i)[0] << " B " << get(i)[1] << std::endl;
//                 }
            }

    };
}

#endif // guard
