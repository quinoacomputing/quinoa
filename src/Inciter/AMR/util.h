#ifndef AMR_util_h
#define AMR_util_h

#include <iostream>
#include <vector>
#include <sstream>

#include "AMR_types.h"

namespace AMR {
    namespace util {

        // Prototypes
        void split(const std::string &s, char delim, std::vector<std::string> &elems);
        std::vector<std::string> split(const std::string &s, char delim);
        coordinate_t find_mid_point(coordinate_t edge_node_A, coordinate_t edge_node_B);
        coordinate_t find_mid_point(real_t x1, real_t y1, real_t z1, real_t x2, real_t y2, real_t z2);

    }
}
#endif // AMR_util
