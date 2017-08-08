#ifndef AMR_marked_refinements_store_h
#define AMR_marked_refinements_store_h

#include <map>
#include "Refinement_State.h"

namespace AMR {

    // TODO: Should we more careful about deleting things from here?
        // Or is it OK to just set them to none?
    class marked_refinements_store_t {
        private:
            // TODO: Make this unordered
            std::map<size_t, Refinement_Case> marked_refinements;

            // TODO: This probably isn't the right place for this
            // We will use this variable to check if anything has changed
            // during a round of refinement marking
            bool refinement_state_changed = false;

        public:
            /**
             * @brief function to see if a given id has already had a
             * refinement decision made
             *
             * @param id The tet_id to check
             *
             * @return A bool stating if a refinement decision/marking exists
             */
            bool exists(size_t id)
            {
                auto f = marked_refinements.find(id);
                if (f != marked_refinements.end())
                {
                    return true;
                }
                return false;
            }

            /**
             * @brief Accessor function to get a marked_refinement from the store
             *
             * @param id The id of the tet to fetch the marked_refinement for
             *
             * @return The marked refinement case for the given tet
             */
            Refinement_Case& get(size_t id)
            {
                return marked_refinements.at(id);
            }

            // TODO: This should probably actually delete the element, not mark it as none?
            // TODO: Document
            void erase(size_t id)
            {
                marked_refinements[id] = Refinement_Case::none;
            }

            /**
             * @brief Function to add a marked refinement to the store, which
             * is context aware based on existing values
             *
             * @param id The id of the tet to mark
             * @param r The refinement decision for the tet
             */
            void add(size_t id, Refinement_Case r)
            {
                // TODO: This is doing a nice double find, once in exist, and
                    // once when looking/setting the value. We can avoid this

                // Check if that active element already exists
                if (exists(id))
                {
                    if (marked_refinements[id] != r)
                    {
                        marked_refinements[id] = r;

                        // TODO: Find a better way to handle/update this global
                        refinement_state_changed = true;
                    }
                }
                else {
                    marked_refinements.insert( std::pair<size_t, Refinement_Case>(id, r));
                    refinement_state_changed = true;
                }
            }

            // TODO: document this
            bool get_state_changed()
            {
                return refinement_state_changed;
            }

            // TODO: document this
            void set_state_changed(bool t)
            {
                refinement_state_changed = t;
            }

    };
}

#endif // guard
