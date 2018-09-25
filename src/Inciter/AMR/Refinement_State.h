//
// Created by Robert Francis Bird on 10/27/16.
//

#ifndef TETAMR_AMR_DATA_H
#define TETAMR_AMR_DATA_H

#include "NoWarning/charm++.h"

#include <vector>
#include <cstddef>

#include "types.h"

namespace AMR {
    // Types
    using real_t = tk::real;

    const size_t DEFUALT_CHILD_NUMBER = 0; // Used to mark "doesn't have children"
    const size_t DEFAULT_NUM_CHILDREN = 0; // Default to no children

    // TODO: Populate this enum with the available refinement cases
    enum Refinement_Case { initial_grid = 0, one_to_two, one_to_four, one_to_eight,
        two_to_eight, four_to_eight, none };

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

        /** @name Charm++ pack/unpack serializer member functions */
        ///@{
        //! \brief Pack/Unpack serialize member function
        //! \param[in,out] p Charm++'s PUP::er serializer object reference
        void pup( PUP::er &p ) {
          p | A;
          p | B;
          p | refinement_criteria;
          p | needs_refining;
          p | needs_derefining;
          p | is_dead;
          p | lock_case;
        }
        //! \brief Pack/Unpack serialize operator|
        //! \param[in,out] p Charm++'s PUP::er serializer object reference
        //! \param[in,out] e Edge_Refinement object reference
        friend void operator|( PUP::er& p, Edge_Refinement& e ) { e.pup(p); }
        //@}

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

    //stop it being copied around?

    class Refinement_State {

        public:

            /// Common to active and master elements
            size_t active_element_number; // TODO: Some of these can be removed?
            Refinement_Case refinement_case;
            size_t num_children; // TODO: this could be replace with children.size()?
            child_id_list_t children;
            size_t refinement_level;
            size_t child_number;
            size_t parent_id;
            int normal; // TODO: make this a bool?

            // Only needed for active
            //size_t master_element_number; // TODO: Some of these can be removed?

            Refinement_State() {} // TODO: Remove?


            /**
             * @brief Constructor which allows for all data fields to be explicitly
             * specified
             *
             * @param active_element_number_in The active element id
             * @param refinement_case_in The refinement case
             * @param num_children_in The number of children
             * @param children_in The children ids
             * @param refinement_level_in The level of refinement
             * @param child_number_in ??  // TODO: What is this?
             * @param parent_id_in Id of parent element
            */
            Refinement_State(
                    size_t active_element_number_in,
                    Refinement_Case refinement_case_in,
                    size_t num_children_in,
                    child_id_list_t children_in,
                    size_t refinement_level_in,
                    size_t child_number_in,
                    size_t parent_id_in
            ) :
                    active_element_number(active_element_number_in),
                    refinement_case(refinement_case_in),
                    num_children(num_children_in),
                    children(children_in),
                    refinement_level(refinement_level_in),
                    child_number(child_number_in),
                    parent_id(parent_id_in),
                    normal(0)
            {
                // Empty
            }

            /**
             * @brief Constructor which assumes sensible Defaults for new nodes
             *
             * @param active_element_number_in The active element id
             * @param refinement_case_in The refinement case
             * @param refinement_level_in The level of refinement
             * @param parent_id_in Id of parent element
             */
            Refinement_State(
                    size_t active_element_number_in,
                    Refinement_Case refinement_case_in,
                    size_t refinement_level_in,
                    size_t parent_id_in
            ) :
                    active_element_number(active_element_number_in),
                    refinement_case(refinement_case_in),
                    num_children(DEFAULT_NUM_CHILDREN), // No children by default
                    refinement_level(refinement_level_in),
                    child_number(DEFUALT_CHILD_NUMBER), // Give it default child id
                    parent_id(parent_id_in),
                    normal(0)
            {
                // Set default size of children to be sensible
                children.reserve(MAX_CHILDREN);
            }
    };

}

#endif //TETAMR_AMR_DATA_H

