#ifndef AMR_node_store_h
#define AMR_node_store_h

#include <cmath>

#include "types.h"
#include "tet_store.h"

// TODO: make this have a base class to support multiple generator schemes
// using the policy design pattern
namespace AMR {

    class node_store_t
    {

        private:
            coord_type* m_x;
            coord_type* m_y;
            coord_type* m_z;

            // We really don't want people to pass this by value..
            // (Because we store refs in here, which are consts..)
            // NonCopyable & operator=(const NonCopyable&) = delete;
            node_store_t(const node_store_t& c) = delete;
            node_store_t& operator=(const node_store_t&) = delete;

        public:
            node_store_t() { } // default cons

            size_t* m_graphsize;

            void set_x(coord_type* x_in) { m_x = x_in; }
            void set_y(coord_type* y_in) { m_y = y_in; }
            void set_z(coord_type* z_in) { m_z = z_in; }

            /**
             * @brief Function to add z coordinate data
             *
             * @param x Data to add
             */
            void add_x(real_t x) { m_x->push_back(x); }

            /**
             * @brief Function to add z coordinate data
             *
             * @param y Data to add
             */
            void add_y(real_t y) { m_y->push_back(y); }

            /**
             * @brief Function to add z coordinate data
             *
             * @param z data to add
             */
            void add_z(real_t z) { m_z->push_back(z); }
            real_t x(size_t id)
            {
                return (*m_x)[id];
            }
            real_t y(size_t id)
            {
                return (*m_y)[id];
            }
            real_t z(size_t id)
            {
                return (*m_z)[id];
            }

            real_t size()
            {
                return m_x->size();
            }

            // TODO: Document this
            void print()
            {
								for (size_t i = 0; i < size(); i++)
								{
										trace_out << "Node " << i << " has coords :" <<
												x(i) << ", " <<
												y(i) << ", " <<
												z(i) << ", " <<
												std::endl;
								}
						}


            /**
             * @brief Function to add a new node
             *
             * @param x x val of node
             * @param y y val of node
             * @param z z val of node
             *
             * @return id of node added
             */
            size_t add(real_t x, real_t y, real_t z) {

                // Need to: Add to {x,y,z} Add any connectivity?

                // Check if the node already exists
                int already_exists = check_node_exists(x,y,z);

                if (already_exists == -1) {
                    size_t return_node_id = add_coordinates(x,y,z);
                    (*m_graphsize)++; // TODO: how to deal with this?
                    trace_out << "--> Made " << return_node_id << std::endl;
                    return return_node_id;
                }
                else {
                    trace_out << "--> Reusing " << already_exists << std::endl;
                    return static_cast<size_t>(already_exists);
                }

            }

            /**
             * @brief Function to add a new node
             *
             * @param coord_tuple The coordinate data to add for the node
             *
             * @return Id of node added
             */
            size_t add(coordinate_t coord_tuple)
            {
                return add( coord_tuple[0], coord_tuple[1], coord_tuple[2]);
            }

            /**
             * @brief Function to add a new point/coordinates
             *
             * @param x x val
             * @param y y val
             * @param z z val
             *
             * @return id of coordinate added
             */
            size_t add_coordinates(real_t x, real_t y, real_t z) {
                add_x(x);
                add_y(y);
                add_z(z);
                return size()-1; // -1 because of the 0 index
            }

            /**
             * @brief Function to add a new node
             *
             * @param x x val of node
             * @param y y val of node
             * @param z z val of node
             *
             * @return id of node added
             */
            size_t add_node(real_t x, real_t y, real_t z) {

                // Need to: Add to {x,y,z} Add any connectivity?

                // Check if the node already exists
                int already_exists = check_node_exists(x,y,z);

                if (already_exists == -1) {
                    size_t return_node_id = add_coordinates(x,y,z);
                    (*m_graphsize)++; // TODO: how to deal with this?
                    trace_out << "--> Made " << return_node_id << std::endl;
                    return return_node_id;
                }
                else {
                    trace_out << "--> Reusing " << already_exists << std::endl;
                    return static_cast<size_t>(already_exists);
                }

            }

            // TODO: Remove all calls to this as it's fairly expensive...
            // TODO: Find a more cost effective way to implement this
                // Most likely change data structure for a faster search
                // This is also going to be a potential problem in async parallel
            /**
             * @brief Helper function to check if a node already exists at
             * coords {x,y,z} to avoid it being duplicated
             *
             * @param x_in X coord to check
             * @param y_in Y coord to check
             * @param z_in Z coord to check
             *
             * @return The id of the node if it exists, -1 if it doesn't.
             */
            int check_node_exists(real_t x_in, real_t y_in, real_t z_in)
            {
                const real_t eps = 1e-7;
                for (size_t i = 0; i < size(); i++)
                {
                    if (
                        std::abs( x(i) - x_in) < eps &&
                        std::abs( y(i) - y_in) < eps &&
                        std::abs( z(i) - z_in) < eps
                    )
                    {
                        trace_out << "!!!! x " << x_in << " y " << y_in <<
                            " z " << z_in << " exits " << std::endl;
                            return static_cast<int>(i);
                    }
                }

                return -1;
            }

            /**
             * @brief function to find the mid point between two points (nodes)
             * based on ids
             *
             * @param id1 id of the first node/point
             * @param id2 id of the second node/point
             *
             * @return The mid point
             */
            coordinate_t find_mid_point(size_t id1, size_t id2)
            {
                coordinate_t edge_node_A = id_to_coordinate(id1);
                coordinate_t edge_node_B = id_to_coordinate(id2);
                return AMR::util::find_mid_point(edge_node_A, edge_node_B);
            }

            /**
             * @brief Function to gather the {x,y,z} coordinates of a node from
             * index id
             *
             * @param id Id of node to gather (direct index into m_x[])
             *
             * @return List (array) of coordinate data
             */
            coordinate_t id_to_coordinate(size_t id)
            {
                assert( id < size());

                // Note: extra braces are to appease Clangs warning generator.
                //   (It's probably ok to remove them....)
                coordinate_t c = { {x(id), y(id), z(id) } };
                return c;
            }

    }; // end class
}

#endif  // guard
