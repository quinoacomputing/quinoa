#ifndef AMR_marked_refinements_store_h
#define AMR_marked_refinements_store_h

#include <unordered_map>
#include "Refinement_State.hpp"

namespace AMR {

    /**
     * @brief This class stores the decisions made during the iterative AMR
     * algorithm.
     *
     * The template parameter determines the internal enum type.
     */
    template<class case_t> class marked_refinements_store_t {
        private:
            // TODO: make this have a more generic name
            std::unordered_map<size_t, case_t> marked_refinements;

            // TODO: This may not be right place for this
            // We will use this variable to check if anything has changed
            // during a round of refinement marking
            bool state_changed = false;

        public:
            //! Const-ref access to number of tets
            //! \return Map of marked refinements
            std::size_t size() const {
              return marked_refinements.size();
            }

            //! Non-const-ref access to state
            //! \return Map of marked refinements
            std::unordered_map<size_t, case_t>& data() {
              return marked_refinements;
            }

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
            case_t& get(size_t id)
            {
                // TODO: is there any performance hit for at?
                return marked_refinements.at(id);
            }

            void erase(size_t id)
            {
                //marked_refinements[id] = case_t::none;

                // Changing to actually erase
                marked_refinements.erase(id);
            }

            /**
             * @brief Function to add a marked refinement to the store, which
             * is context aware based on existing values
             *
             * @param id The id of the tet to mark
             * @param r The refinement decision for the tet
             */
            void add(size_t id, case_t r)
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
                        state_changed = true;
                    }
                    else {
                        trace_out << "Not setting marked refinement val as same val"<< std::endl;
                    }
                }
                else {
                    trace_out << "Adding new marked value " << id << " = " << r << std::endl;
                    marked_refinements.insert( std::pair<size_t, case_t>(id, r));
                    state_changed = true;
                }
            }

            /**
             * @brief Accessor for variable which tracks state change
             *
             * @return Bool stating is the state has changed
             */
            bool& get_state_changed()
            {
                return state_changed;
            }
    };
}

#endif // guard
