// *****************************************************************************
/*!
  \file      src/Mesh/Around.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Helper class for iterating through linked lists of derived data
  \details   Helper class for iterating through every item in a linked list data
    structure derived from unstructured mesh connectivity.
  \see src/Mesh/DerivedData.h
*/
// *****************************************************************************
#ifndef Around_h
#define Around_h

#include <vector>
#include <map>
#include <utility>
#include <cstddef>
#include "Types.h"
#include "Fields.h"
#include "UnsMesh.h"

namespace tk {

//! \brief Helper class simplifying client code for iterating on entries
//!   surrounding entries via linked lists derived from unstructured mesh
//!   connectivity
class Around {

  private:
    //! Linked list type: T1: item list, T2: index list
    using List =
      std::pair< std::vector< std::size_t >, std::vector< std::size_t > >;

    //! Const iterator to the item list in List::T1
    using const_iterator = List::first_type::const_iterator;

    //! Difference type for iterator/pointer arithmetics
    using diff_type = List::first_type::difference_type;

  public:
    //! Constructor
    //! \param[in] list Linked list (vectors) storing derived data
    //! \param[in] idx Index of entry whose surrounding items we will iterate on
    explicit Around( const List& list, std::size_t idx ) :
      m_list( list ), m_idx( idx ) {}

    //! Iterator to the beginning of the entries of surrounding entries
    //! \return Iterator to the beginning of the entries of surrounding entries
    const_iterator begin() {
      return m_list.first.begin() +
             static_cast< diff_type >( m_list.second[m_idx] + 1 );
    }

    //! Iterator to the entry after the last of the surrounding entries
    //! \return Iterator to the entry after the last of the surrounding entries
    const_iterator end() {
      return m_list.first.begin() +
             static_cast< diff_type >( m_list.second[m_idx+1] + 1 );
    }

  private:
    const List& m_list; //!< Linked list to iterate in
    std::size_t m_idx;  //!< Index whose surrounding entries to iterate on
};

} // tk::

#endif // Around_h
