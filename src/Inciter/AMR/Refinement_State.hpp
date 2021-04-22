//
// Created by Robert Francis Bird on 10/27/16.
//

#ifndef TETAMR_AMR_DATA_H
#define TETAMR_AMR_DATA_H

#include <cstddef>

#include "Types.hpp"
#include "AMR_types.hpp"

namespace AMR {
    // Types
    using real_t = tk::real;

    const size_t DEFUALT_CHILD_NUMBER = 0; // Used to mark "doesn't have children"

    // TODO: Populate this enum with the available refinement cases
    enum Refinement_Case { initial_grid = 0, one_to_two, one_to_four, one_to_eight,
        two_to_eight, four_to_eight, none };

    enum Derefinement_Case { two_to_one = 0, four_to_one, four_to_two,
        eight_to_one, eight_to_two, eight_to_four, skip };

    class Refinement_State {

        public:
            /// Common to active and master elements
            size_t active_element_number; // TODO: Some of these can be removed?
            Refinement_Case refinement_case;
            child_id_list_t children;
            size_t refinement_level;
            size_t child_number;
            size_t parent_id;
            bool normal;
            bool has_parent;

            Refinement_State() {}

            /**
             * @brief Constructor which allows for all data fields to be explicitly
             * specified
             *
             * @param active_element_number_in The active element id
             * @param refinement_case_in The refinement case
             * @param children_in The children ids
             * @param refinement_level_in The level of refinement
             * @param child_number_in ??  // TODO: What is this?
             * @param parent_id_in Id of parent element
             * @param has_parent_in True if element has a parent, default is true
            */
            Refinement_State(
                    size_t active_element_number_in,
                    Refinement_Case refinement_case_in,
                    const child_id_list_t& children_in,
                    size_t refinement_level_in,
                    size_t child_number_in,
                    size_t parent_id_in,
                    bool has_parent_in=true
            ) :
                    active_element_number(active_element_number_in),
                    refinement_case(refinement_case_in),
                    children(children_in),
                    refinement_level(refinement_level_in),
                    child_number(child_number_in),
                    parent_id(parent_id_in),
                    normal(0),
                    has_parent(has_parent_in)
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
             * @param has_parent_in True if element has a parent, default is true
             */
            Refinement_State(
                    size_t active_element_number_in,
                    Refinement_Case refinement_case_in,
                    size_t refinement_level_in,
                    size_t parent_id_in,
                    bool has_parent_in=true
            ) :
                    active_element_number(active_element_number_in),
                    refinement_case(refinement_case_in),
                    refinement_level(refinement_level_in),
                    child_number(DEFUALT_CHILD_NUMBER), // Give it default child id
                    parent_id(parent_id_in),
                    normal(0),
                    has_parent(has_parent_in)
            {
                // Set default size of children to be sensible
                children.reserve(MAX_CHILDREN);
            }
    };

}

#endif //TETAMR_AMR_DATA_H

