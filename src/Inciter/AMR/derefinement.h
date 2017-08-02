#ifndef AMR_derefinement_h
#define AMR_derefinement_h

#include "tet_store.h"

// TODO: make this have a base class to support multiple generator schemes
// using the policy design pattern
namespace AMR {

    class derefinement_t {
        private:

            // TODO: Tidy up edge store references here 
            tet_store_t* tet_store;

        public:

            const size_t DEFAULT_REFINEMENT_LEVEL = 0; //TODO: Is this in the right place?
            const size_t MIN_REFINEMENT_LEVEL = DEFAULT_REFINEMENT_LEVEL;

            derefinement_t(tet_store_t* ts) : tet_store(ts)
            {
                assert(0);
            }

    };
}

#endif  // guard
