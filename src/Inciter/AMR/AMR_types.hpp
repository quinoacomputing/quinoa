#ifndef AMR_types_h
#define AMR_types_h

#include <array>
#include <vector>
#include <map>

#include "../Base/Types.hpp"
#include "edge.hpp"
#include "UnsMesh.hpp"

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

using node_pair_t  = std::array<std::size_t, 2>;

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
    int needs_refining; // value of 1= refinement; 2= deref-ref (as in for 8:4)
    bool needs_derefining; // TODO: Marge this with needs_refining
    Edge_Lock_Case lock_case; // TODO: Refactor this to match _ style?

    // Explicit Empty Constructor
    Edge_Refinement() :
        A(0),
        B(0),
        needs_refining(0),
        needs_derefining(false),
        lock_case(Edge_Lock_Case::unlocked)
    {
        // Empty
    }

   // bool operator==( const Edge_Refinement& r ) const {
   //   return A == r.A &&
   //          B == r.B &&
   //          //std::abs(refinement_criteria-r.refinement_criteria) < 1.0e-12 &&
   //          needs_refining == r.needs_refining &&
   //          needs_derefining == r.needs_derefining &&
   //          is_dead == r.is_dead &&
   //          lock_case == r.lock_case;
   // }

    // This abstraction is hardly any better than using an explicit initialisation
    // list but it makes it easier if we decide to add/remove a parameter
    Edge_Refinement(
            size_t A_in,
            size_t B_in,
            int needs_refining_in,
            bool needs_derefining_in,
            Edge_Lock_Case lock_case_in
            ) :
        A(A_in),
        B(B_in),
        needs_refining(needs_refining_in),
        needs_derefining(needs_derefining_in),
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

//! \brief Needs refinement and edge lock case associated to an edge given by global
//!    parent IDs
using EdgeData =
   std::unordered_map< tk::UnsMesh::Edge,
                       std::tuple< int, int, Edge_Lock_Case >,
                       tk::UnsMesh::Hash<2>,
                       tk::UnsMesh::Eq<2> >;

//! Enum used to tag an edge for refinement or derefinement
enum class edge_tag : uint8_t { REFINE, DEREFINE };

}  // AMR::

#endif
