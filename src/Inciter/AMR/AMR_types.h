#ifndef AMR_types_h
#define AMR_types_h

#include <array>
#include <vector>
#include <map>

#include "../Base/Types.h"
#include "edge.h"

// TODO: Do we need to merge this with Base/Types.h?

namespace AMR {

const int DIMENSION = 3;
const size_t NUM_TET_NODES = 4;
const size_t NUM_FACE_NODES = 3;
const size_t NUM_TET_EDGES = 6;
const size_t NUM_TET_FACES = 4;
const size_t ID_SHIFT = 3;
const size_t MAX_CHILDREN = 8;
const char KEY_DELIM = '-';

using real_t = tk::real;
using coordinate_t = std::array<real_t, DIMENSION>;
using tet_t = std::array<size_t, NUM_TET_NODES>;

//using node_pair_t  = std::array<std::size_t, 2>;
using node_pair_t  = std::pair<std::size_t, std::size_t>;

using face_ids_t = std::array<std::size_t, NUM_FACE_NODES>;
using face_list_t  = std::array< face_ids_t, NUM_TET_FACES>;

//using child_id_list_t = std::array<size_t, MAX_CHILDREN>;
using child_id_list_t = std::vector<size_t>;

using tet_list_t = std::map<size_t, tet_t>;

using inpoel_t = std::vector< std::size_t >;     //!< Tetrahedron connectivity
using node_list_t = std::vector<real_t>;

enum Edge_Lock_Case {unlocked = 0, locked, intermediate, temporary};

 // TODO: Make these class enums? (breaks printing)
struct Edge_Refinement {
    size_t A;
    size_t B;
    real_t refinement_criteria;
    bool needs_refining; // TODO: This could possibly be deduced implicitly
    bool needs_derefining; // TODO: Marge this with needs_refining
    bool is_dead;
    Edge_Lock_Case lock_case; // TODO: Refactor this to match _ style?

    // Explicit Empty Constructor
    Edge_Refinement() :
        A(0),
        B(0),
        refinement_criteria(0.0),
        needs_refining(false),
        needs_derefining(false),
        is_dead(false),
        lock_case(Edge_Lock_Case::unlocked)
    {
        // Empty
    }

    // This abstraction is hardly any better than using an explicit initialisation
    // list but it makes it easier if we decide to add/remove a parameter
    Edge_Refinement(
            size_t A_in,
            size_t B_in,
            real_t refinement_criteria_in,
            bool needs_refining_in,
            bool needs_derefining_in,
            bool is_dead_in,
            Edge_Lock_Case lock_case_in
            ) :
        A(A_in),
        B(B_in),
        refinement_criteria(refinement_criteria_in),
        needs_refining(needs_refining_in),
        needs_derefining(needs_derefining_in),
        is_dead(is_dead_in),
        lock_case(lock_case_in)
    {
        // Empty, all implicit.
        // Could add logic here to reconcile needs_refining and needs_derefining
    }
};
// Complex types
struct Edge_Refinement; // forward declare
using edges_t = std::map<edge_t, Edge_Refinement>;
using edge_list_t  = std::array<edge_t, NUM_TET_EDGES>;
using edge_list_ids_t  = std::array<std::size_t, NUM_TET_EDGES>;

using coord_type = std::vector< tk::real >;

}  // AMR::

#endif
