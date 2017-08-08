#ifndef AMR_master_element_store_h
#define AMR_master_element_store_h

#include <map>
#include <algorithm>

#include "Base/Exception.h"
#include "Refinement_State.h"

namespace AMR {

    class master_element_store_t {
        private:
            std::map<size_t, Refinement_State> master_elements;
        public:
            /**
             * @brief Add an element to the master element list.
             *
             * @param element_number The element number to add (currently
             * somewhat redundant?)
             * @param refinement_case The refinement_case which gave rise to
             * this element
             * @param refinement_level The level of refinement
             * @param parent_id The id of the parent element
             *
             * @return The id of the added element
            */
            // TODO: These add methods could probably be much tidier. The
            // current none reliance on default input values is nice, as it
            // means that we don't have to do any checks on if a valid value
            // was passed for the base case..
            size_t add
            (
                 size_t element_number,
                 Refinement_Case refinement_case,
                 size_t refinement_level,
                 size_t parent_id
            )
            {
                Refinement_State data = Refinement_State(
                        element_number,
                        refinement_case,
                        refinement_level,
                        parent_id
                );

                master_elements.insert( std::pair<size_t, Refinement_State>(element_number, data));

                return element_number;
            }

            /**
             * @brief Add a master element, with a default parent_id (0) and
             * default refinement level (0)
             *
             * @param element_number The element number to add
             * @param refinement_case The refinement_case which gave rise to
             * this element
             *
             * @return The id of the added element
            */
            size_t add(
                 size_t element_number,
                 Refinement_Case refinement_case
            )
            {
                add(element_number, refinement_case, 0, 0);
                return element_number;
            }

            /**
             * @brief Add a master element, with a specified parent_id which is
             * used to udpate the refinement_level
             *
             * @param element_number The element number to add
             * @param refinement_case The refinement_case which gave rise to
             * this element
             * @param parent_id The id of the parent element
             *
             * @return The id of the added element
             */
            size_t add(
                 size_t element_number,
                 Refinement_Case refinement_case,
                 size_t parent_id
            )
            {
                size_t refinement_level =
                    get(parent_id).refinement_level + 1;

                add(element_number, refinement_case,
                        refinement_level, parent_id);

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
                Assert( exists(id), "ID does not exist" );
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
            bool exists(size_t id)
            {
                auto f = master_elements.find(id);
                if (f != master_elements.end())
                {
                    return true;
                }
                return false;
            }

            // TODO: document this
            void erase(size_t id) {
                master_elements.erase(id);
            }

            size_t get_parent(size_t id)
            {
                return get(id).parent_id;
            }

            /**
             * @brief Wrapper function to calculate number of master elements
             *
             * @return total number of master elements
             */
            size_t size() {
                return master_elements.size();
            }

            // TODO: Document
            void add_child(size_t parent_id, size_t child_id)
            {
                get(parent_id).children.push_back(child_id);
                get(parent_id).num_children++;

                Assert(
                    get(parent_id).num_children <= 8,
                    "Addition would violate max child rules"
                );
            }

            size_t get_child_id(size_t parent_id, size_t offset)
            {
                Assert(offset < get(parent_id).children.size(), "Child ID out of range");
                return get(parent_id).children[offset];
            }

    };
}

#endif // guard
