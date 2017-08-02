#ifndef AMR_edge_t_h
#define AMR_edge_t_h

#include <iostream>
#include <ostream>

class edge_t {
    using edge_ = std::pair<size_t, size_t>;
    private:
        edge_ data;
        friend std::ostream& operator<<(std::ostream&, const edge_t&);

    public:

        edge_ get_data() const {
            return data;
        }

        // Constructors
        edge_t()
        {
        }
        edge_t(size_t A, size_t B)
        {
            data = std::make_pair(std::min(A,B), std::max(A,B));
        }

        // Copy constructor
        edge_t(const edge_t &e)
        {
            data = e.data;
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

        size_t first()
        {
            return data.first;
        }
        size_t second()
        {
            return data.second;
        }

        void replace(size_t new_id, size_t old_id)
        {
            if (data.first == old_id)
            {
                data.first = new_id;
            }

            if (data.second == old_id)
            {
                data.second = new_id;
            }
        }

};

#endif
