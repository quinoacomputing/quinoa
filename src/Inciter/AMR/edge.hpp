#ifndef AMR_edge_t_h
#define AMR_edge_t_h

#include <iostream>
#include <stddef.h>
#include <array>
#include <algorithm>

namespace AMR {

class edge_t {
    using edge_ = std::array< std::size_t, 2 >;
    private:
        edge_ data;
        friend std::ostream& operator<<(std::ostream&, const edge_t&);

    public:
        edge_& get_data() {
            return data;
        }

        const edge_& get_data() const {
            return data;
        }

        // Constructors
        edge_t()
        {
        }

        edge_t(size_t A, size_t B) : data( {{std::min(A,B), std::max(A,B)}} )
        {
        }

        explicit edge_t( edge_ e ) : data( std::move(e) ) {}

        // Operators
            // Piggy back underlying edge_ type where possible
        bool operator==(const edge_t& rhs) const
        {
          // ensure entries of rhs and this are in ascending order
          auto this_copy = this->data;
          this_copy[0] = std::min(this->data[0], this->data[1]);
          this_copy[1] = std::max(this->data[0], this->data[1]);
          std::array< std::size_t, 2 > rhs_copy;
          rhs_copy[0] = std::min(rhs.get_data()[0], rhs.get_data()[1]);
          rhs_copy[1] = std::max(rhs.get_data()[0], rhs.get_data()[1]);

          if (this_copy[0] == rhs_copy[0] && this_copy[1] == rhs_copy[1])
            return true;
          else
            return false;
        }
        //bool operator>(const edge_t& rhs) const
        //{
        //  return (data > rhs.get_data());
        //}
        bool operator<(const edge_t& rhs) const
        {
          // ensure entries of rhs and this are in ascending order
          auto this_copy = this->data;
          this_copy[0] = std::min(this->data[0], this->data[1]);
          this_copy[1] = std::max(this->data[0], this->data[1]);
          std::array< std::size_t, 2 > rhs_copy;
          rhs_copy[0] = std::min(rhs.get_data()[0], rhs.get_data()[1]);
          rhs_copy[1] = std::max(rhs.get_data()[0], rhs.get_data()[1]);

          if (this_copy[0] < rhs_copy[0])
            return true;
          else if (this_copy[0] == rhs_copy[0] && this_copy[1] < rhs_copy[1])
            return true;
          else
            return false;
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

}  // AMR::

#endif
