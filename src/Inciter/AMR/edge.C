#include "edge.h"

// Implementation of the friend function, extending osstream and not an edge_t
// member
std::ostream& operator<<(std::ostream& os, const edge_t& e)
{
    os << e.data.first << "-" << e.data.second;
    return os;
}
