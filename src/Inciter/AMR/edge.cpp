#include "edge.hpp"

namespace AMR {

// Implementation of the friend function, extending osstream and not an edge_t
// member
std::ostream& operator<<(std::ostream& os, const AMR::edge_t& e)
{
    os << e.get_data()[0] << "-" << e.get_data()[1];
    return os;
}

} // AMR::
