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
            coord_type m_x;
            coord_type m_y;
            coord_type m_z;

            // We really don't want people to pass this by value..
            // (Because we store refs in here, which are consts..)
            // NonCopyable & operator=(const NonCopyable&) = delete;
            //node_store_t(const node_store_t& c) = delete;
            //node_store_t& operator=(const node_store_t&) = delete;

        public:
            node_store_t() { } // default cons

            size_t m_graphsize;

            void set_x(coord_type x_in) { m_x = x_in; }
            void set_y(coord_type y_in) { m_y = y_in; }
            void set_z(coord_type z_in) { m_z = z_in; }

            /**
             * @brief Function to add x coordinate data
             *
             * @param xc Data to add
             */
            void add_x(real_t xc) { m_x.push_back(xc); }

            /**
             * @brief Function to add y coordinate data
             *
             * @param yc Data to add
             */
            void add_y(real_t yc) { m_y.push_back(yc); }

            /**
             * @brief Function to add z coordinate data
             *
             * @param zc data to add
             */
            void add_z(real_t zc) { m_z.push_back(zc); }

            coord_type get_x_array() {
                return m_x;
            }
            coord_type get_y_array() {
                return m_y;
            }
            coord_type get_z_array() {
                return m_z;
            }

            real_t x(size_t id)
            {
                return m_x[id];
            }
            real_t y(size_t id)
            {
                return m_y[id];
            }
            real_t z(size_t id)
            {
                return m_z[id];
            }

            size_t size()
            {
                return m_x.size();
            }

            // TODO: Document this
            void print()
            {
                for (size_t i = 0; i < size(); i++)
                {
                    std::cout << "Node " << i << " has coords :" <<
                        x(i) << ", " <<
                        y(i) << ", " <<
                        z(i) << ", " <<
                        std::endl;
                }
            }


            /**
             * @brief Function to add a new node
             *
             * @param xc x val of node
             * @param yc y val of node
             * @param zc z val of node
             *
             * @return id of node added
             */
            size_t add(real_t xc, real_t yc, real_t zc) {

                // Need to: Add to {x,y,z} Add any connectivity?

                // Check if the node already exists
                int already_exists = check_node_exists(xc,yc,zc);

                if (already_exists == -1) {
                    size_t return_node_id = add_coordinates(xc,yc,zc);
                    m_graphsize++; // TODO: how best to deal with this?
                    trace_out << "Made new node " << return_node_id << std::endl;
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
             * @param xc x val
             * @param yc y val
             * @param zc z val
             *
             * @return id of coordinate added
             */
            size_t add_coordinates(real_t xc, real_t yc, real_t zc) {
                add_x(xc);
                add_y(yc);
                add_z(zc);
                return size()-1; // -1 because of the 0 index
            }

            /**
             * @brief Function to add a new node
             *
             * @param xc x val of node
             * @param yc y val of node
             * @param zc z val of node
             *
             * @return id of node added
             */
// TODO: Why is there 2 add functions
            size_t add_node(real_t xc, real_t yc, real_t zc) {

                // Need to: Add to {x,y,z} Add any connectivity?

                // Check if the node already exists
                int already_exists = check_node_exists(xc,yc,zc);

                if (already_exists == -1) {
                    size_t return_node_id = add_coordinates(xc,yc,zc);
                    m_graphsize++; // TODO: how to deal with this?
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
