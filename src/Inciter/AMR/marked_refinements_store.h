#ifndef AMR_marked_refinements_store_h
#define AMR_marked_refinements_store_h

#include <map>
#include "Refinement_State.h"

namespace AMR {

    class marked_refinements_store_t {
        private:
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
                    trace_out << "Marked element " << id << " has value " <<
                        f->second << std::endl;

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
                // Check if that active element already exists
                if (exists(id))
                {
                    if (marked_refinements[id] != r)
                    {
                        trace_out << "Updating marked value to " << r <<
                            " was " << marked_refinements[id] << std::endl;

                        marked_refinements[id] = r;

                        // TODO :Find a better way to handle/update this global
                        refinement_state_changed = true;
                    }
                    else {
                        trace_out << "Not setting marked refinement val as same val"<< std::endl;
                    }
                }
                else {
                    trace_out << "Adding new marked value " << id << std::endl;
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

            void replace(size_t old_id, size_t new_id)
            {
                // Swap id out in map
                auto i = marked_refinements.find(old_id);
                auto value = i->second;
                marked_refinements.erase(i);
                marked_refinements[new_id] = value;
            }
    };
}

#endif // guard
