#ifndef AMR_master_element_store_h
#define AMR_master_element_store_h

#include <map>
#include <algorithm>
#include <cassert>

#include "Refinement_State.hpp"
#include "AMR/Loggers.hpp"                   // for trace_out

namespace AMR {

    class master_element_store_t {
        private:
            std::map<size_t, Refinement_State> master_elements;
        public:
            //! Non-const-ref access to state
            std::map<size_t, Refinement_State>& data() {
              return master_elements;
            }

            /**
             * @brief Add an element to the master element list.
             *
             * @param element_number The element number to add (currently
             * somewhat redundant?)
             * @param refinement_case The refinement_case which gave rise to
             * this element
             * @param refinement_level The level of refinement
             * @param parent_id The id of the parent element
             * @param has_parent True if element has a parent, default is true
             *
             * @return The id of the added element
            */
            // TODO: These add methods could probably be much tidier. The
            // current none reliance on default input values is nice, as it
            // means that we don't have to do any checks on if a valid value
            // was passed for the base case..
            size_t add(
                 size_t element_number,
                 Refinement_Case refinement_case,
                 size_t refinement_level,
                 size_t parent_id,
                 bool has_parent=true
            )
            {
                Refinement_State d = Refinement_State(
                        element_number,
                        refinement_case,
                        refinement_level,
                        parent_id,
                        has_parent
                );

                master_elements.insert( std::pair<size_t, Refinement_State>(element_number, d));

                return element_number;
            }

            /**
             * @brief Accessor method to retrieve master element by element id
             *
             * @param id The element_id of the master_element to fetch
             *
             * @return The master_element which represents the corresponding tet
             */
            Refinement_State& get(size_t id)
            {
                // cppcheck-suppress assertWithSideEffect
                assert( exists(id) );
                return master_elements.at(id);
            }

            const Refinement_State& get(size_t id) const
            {
                assert( exists(id) );
                return master_elements.at(id);
            }

            /**
             * @brief Function to check if master element entry exists. Useful
             * for debugging access to invalid elements, or trying to re-create
             * an element which already exists
             *
             * @param id Id to check
             *
             * @return Bool stating if the element already exists
             */
            bool exists(size_t id) const
            {
                auto f = master_elements.find(id);
                if (f != master_elements.end())
                {
                    //trace_out << "master_elements " << id << " exists." << std::endl;
                    return true;
                }
                return false;
            }

            /**
             * @brief Getting function to access parent of a given ID
             *
             * @param id ID to find parent of
             *
             * @return ID of parent
             */
            size_t get_parent(size_t id) const
            {
                return get(id).parent_id;
            }

            /**
             * @brief Getting function to find child id for a given parent
             *
             * @param parent_id ID of parent
             * @param offset Offset into child list (i.e child number)
             *
             * @return ID of the child
             */
            size_t get_child_id(size_t parent_id, size_t offset) const
            {
                // cppcheck-suppress assertWithSideEffect
                assert(offset < get(parent_id).children.size());
                return get(parent_id).children[offset];
            }

            /**
             * @brief Wrapper function to calculate number of master elements
             *
             * @return total number of master elements
             */
            size_t size() {
                return master_elements.size();
            }

            /**
             * @brief Add child to a given node, to track the relationship and
             * allow graph traversal
             *
             * @param parent_id ID of the parent
             * @param child_id ID of the child
             */
            void add_child(size_t parent_id, size_t child_id)
            {
                get(parent_id).children.push_back(child_id);

                // cppcheck-suppress assertWithSideEffect
                assert( get(parent_id).children.size() <= 8);
            }

            /**
             * @brief  Remove id from master_element (relevant for
             * de-refinement)
             *
             * @param id ID of element to remove
             */
            void erase(size_t id) {
                master_elements.erase(id);
            }


    };
}

#endif // guard
