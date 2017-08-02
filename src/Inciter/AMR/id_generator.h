#ifndef AMR_id_generator_h
#define AMR_id_generator_h

#include <limits>
#include <cassert>

// TODO: make this have a base class to support multiple generator schemes
// using the policy design pattern
namespace AMR {

    class id_generator_t {
        protected:
            size_t START_TET_ID = 0;
            size_t start_id;

            // Used to track which tet_id to give the next parent
            // The range between this and get_child_id(this) determinate how
            // many tets can be on the first level
            size_t next_tet_id;

        public:
            // Constructor
            id_generator_t()
            {
                start_id = START_TET_ID;
                next_tet_id = start_id;
            }

            /**
             * @brief Helper function to generate tet ids
             *
             * @return Id of the next tet
             */
            size_t get_next_tet_id()
            {
                return next_tet_id++;
            }

            /**
             * @brief Helper function to get the all the child ids for a given
             * parent
             *
             * WARNING: If you don't use all the children you ask for, you may have a bad time...
             *
             * @param parent_id The id of the parent
             *
             * @return The list of children ids
             */
            child_id_list_t generate_child_ids(size_t parent_id, size_t count = MAX_CHILDREN)
            {
                child_id_list_t c;
                c.resize(count);
                c[0] = parent_id; // TODO: Remove this. Horrific hack to suppress warning
                for (size_t i = 0; i < count; i++)
                {
                    c[i] = get_next_tet_id();
                }
                return c;
            }
    };

    class morton_id_generator_t : public id_generator_t
    {
        public:
            // This can be calculated as something like floor( num_bits - log_2(START_TET_ID) ) / log_2(MAX_CHILDREN)
            // So for a simulation which supports an initial grid of ~1,048,576 elements in 3D:
            // (64 - 20) / 3  = 14

            // This basically says the number of tets which can be in an initial grid
            // A sensible value is 2^20 (1,048,576) for big simulations, and anything
            // smaller for toy problems
            size_t START_TET_ID = 1024;

            // Constructor to reset START_TET_ID on the new value
            morton_id_generator_t() : id_generator_t() {
                // TODO: Is there a nice way to remove this code duplication
                // between this and the base class?
                id_generator_t::start_id = START_TET_ID;
                next_tet_id = start_id;
            }

            /**
             * @brief Helper function to convert from parent id to (base) child
             * id
             *
             * @param parent_id Id of parent
             *
             * @return Id of child
             */
            static size_t get_child_id(size_t parent_id)
            {
                size_t child_id = (parent_id << ID_SHIFT);
                return child_id;
            }

            /**
             * @brief Helper function to get the all the child ids for a given
             * parent
             *
             * @param parent_id The id of the parent
             *
             * @return The list of children ids
             */
            static child_id_list_t generate_child_ids(size_t parent_id, size_t count = MAX_CHILDREN)
            {
                child_id_list_t c;
                c.resize(count);
                for (size_t i = 0; i < count; i++)
                {
                    c[i] = get_child_id(parent_id, i);
                }
                return c;
            }

            /**
             * @brief Helper function to get the child id of a specific child
             * of a parent
             *
             * @param parent_id id of the parent
             * @param offset offset into the parent list (i.e number of child)
             *
             * @return the id of the given child
             */
            static size_t get_child_id(size_t parent_id, size_t offset)
            {
                // Try detect overflow
                assert( parent_id <= get_parent_id(std::numeric_limits<size_t>::max()));
                return get_child_id(parent_id) + offset;
            }

            /**
             * @brief Helper function to get the parent id of a child
             *
             * @param child_id Child for whom we should find the parent
             *
             * @return The parent id
             */
            static size_t get_parent_id(size_t child_id)
            {
                return (child_id >> ID_SHIFT);
            }

    };
}

#endif // guard
