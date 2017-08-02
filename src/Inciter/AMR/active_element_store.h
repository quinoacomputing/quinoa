#ifndef AMR_active_element_store_h
#define AMR_active_element_store_h

#include <set>

namespace AMR {

    class active_element_store_t {
        private:
            std::set<size_t> active_elements;
        public:
            /**
             * @brief Function to add active elements
             *
             * @param id The master_id this should map to
             */
            void add(size_t id)
            {
                // Check if that active element already exists
                assert( !exists(id) );
                active_elements.insert(id);
            }

            void erase(size_t id)
            {
                active_elements.erase(id);
            }

            /**
             * @brief function to check if an active element exists in the
             * active_elements store
             *
             * @param id The id of the element to check
             *
             * @return A bool for if the element was found or no
             */
            bool exists(size_t id)
            {
                if (active_elements.find(id) != active_elements.end())
                {
                    return true;
                }
                return false;
            }

            void replace(size_t old_id, size_t new_id)
            {
                erase(old_id);
                add(new_id);
            }
    };
}

#endif // guard
