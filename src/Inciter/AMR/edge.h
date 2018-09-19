#ifndef AMR_edge_t_h
#define AMR_edge_t_h

#include "NoWarning/charm++.h"

#include <iostream>
#include <ostream>
#include <array>

class edge_t {
    using edge_ = std::array< std::size_t, 2 >;
    private:
        edge_ data;
        friend std::ostream& operator<<(std::ostream&, const edge_t&);

    public:
        /** @name Charm++ pack/unpack serializer member functions */
        ///@{
        //! \brief Pack/Unpack serialize member function
        //! \param[in,out] p Charm++'s PUP::er serializer object reference
        void pup( PUP::er &p ) {
          p | data;
        }
        //! \brief Pack/Unpack serialize operator|
        //! \param[in,out] p Charm++'s PUP::er serializer object reference
        //! \param[in,out] e edge_t object reference
        friend void operator|( PUP::er& p, edge_t& e ) { e.pup(p); }
        //@}

        edge_ get_data() const {
            return data;
        }

        // Constructors
        edge_t()
        {
        }
        edge_t(size_t A, size_t B)
        {
            data = {{ std::min(A,B), std::max(A,B) }};
        }

        // Operators
            // Piggy back underlying edge_ type where possible
        bool operator==(const edge_t& rhs) const
        {
            return (this->data==rhs.get_data());
        }
        bool operator>(const edge_t& rhs) const
        {
          return (data > rhs.get_data());
        }
        bool operator<(const edge_t& rhs) const
        {
          return (data < rhs.get_data());
        }

        size_t first() const
        {
            return data[0];
        }
        size_t second() const
        {
            return data[1];
        }

        void replace(size_t new_id, size_t old_id)
        {
            if (data[0] == old_id)
            {
                data[0] = new_id;
            }

            if (data[1] == old_id)
            {
                data[1] = new_id;
            }
        }

};

#endif
