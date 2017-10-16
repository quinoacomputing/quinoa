#include "util.h"

namespace AMR {
    namespace util {

        /**
         * @brief Split a string based on a given delimiter
         *
         * @param s The string to split
         * @param delim The delimiter to use
         * @param elems The vector which will be filled
         */
        void split(const std::string &s, char delim,
                std::vector<std::string> &elems)
        {
            std::stringstream ss;
            ss.str(s);
            std::string item;
            while (std::getline(ss, item, delim)) {
                elems.push_back(item);
            }
        }

        /**
         * @brief Build a vector of split strings based on a delmiter
         *
         * @param s The string to split
         * @param delim The delimiter to use
         *
         * @return A vector split by delim. Does not handle empty tokens
         */
        std::vector<std::string> split(const std::string &s, char delim) {
            std::vector<std::string> elems;
            split(s, delim, elems);
            return elems;
        }


        /**
         * @brief function to find the mid point betwen two points (nodes)
         *
         * @param edge_node_A First node/point
         * @param edge_node_B Second node/point
         *
         * @return  The mid point between the two nodes
         */
        coordinate_t find_mid_point(
                coordinate_t edge_node_A,
                coordinate_t edge_node_B
                )
        {
            coordinate_t mid_point;

            for(size_t i = 0; i < DIMENSION; i++)
            {
                mid_point[i] = (edge_node_A[i]+edge_node_B[i])/2.0;
            }

            return mid_point;
        }

        /**
         * @brief Function to find the mid point between two points (nodes)
         *
         * @param x1 x coord of first point
         * @param y1 y coord of first point
         * @param z1 z coord of first point
         * @param x2 x coord of second point
         * @param y2 y coord of second point
         * @param z2 z coord of second point
         *
         * @return The mid point between the two nodes
         */
        coordinate_t find_mid_point(real_t x1, real_t y1,
                real_t z1, real_t x2, real_t y2, real_t z2)
        {

            coordinate_t mid_point;

            mid_point[0] = (x1+x2)/2.0;
            mid_point[1] = (y1+y2)/2.0;
            mid_point[2] = (z1+z2)/2.0;

            return mid_point;
        }

        /*
           real_t jacobian(size_t a, size_t b, size_t c, size_t d)
           {

           std::cout << "A:  x " << m_x[a] << " y " << m_y[a] << " z " << m_z[a] << std::endl;
           std::cout << "b:  x " << m_x[b] << " y " << m_y[b] << " z " << m_z[b] << std::endl;
           std::cout << "c:  x " << m_x[c] << " y " << m_y[c] << " z " << m_z[c] << std::endl;
           std::cout << "d:  x " << m_x[d] << " y " << m_y[d] << " z " << m_z[d] << std::endl;

           const std::array< tk::real, 3 > ba = {{ m_x[b]-m_x[a], m_y[b]-m_y[a], m_z[b]-m_z[a] }};
           const std::array< tk::real, 3 > ca = {{ m_x[c]-m_x[a], m_y[c]-m_y[a], m_z[c]-m_z[a] }};
           const std::array< tk::real, 3 > da = {{ m_x[d]-m_x[a], m_y[d]-m_y[a], m_z[d]-m_z[a] }};

           const auto J = tk::triple( ba, ca, da );

           return J;

           }*/

    }
}

