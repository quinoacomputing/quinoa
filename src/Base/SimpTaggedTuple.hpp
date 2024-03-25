// *****************************************************************************
/*!
  \file      src/Base/SimpTaggedTuple.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.,
             2022-2024 J. Bakosi
             All rights reserved. See the LICENSE file for details.
  \brief     Tagged tuple allowing tag-based access
  \details   Tagged tuple allowing tag-based access. This is very much like
    [std::tuple](http://en.cppreference.com/w/cpp/utility/tuple), but instead of
    having to index the elements by integers, it allows access by a tag, which
    can be an empty struct with a unique name. Credit goes to
    ecatmur_at_stackoverflow.com, for more details, see
    http://stackoverflow.com/questions/13065166/c11-tagged-tuple.
*/
// *****************************************************************************
#pragma once

#include <type_traits>
#include <tuple>

#include <brigand/adapted/tuple.hpp>

#include "NoWarning/any.hpp"
#include "NoWarning/partition.hpp"
#include "NoWarning/index_of.hpp"

#include "PUPUtil.hpp"
#include "Exception.hpp"

namespace tag {
//! Printable tag for SimpTaggedTuple that returns its name
#define DEFTAG(n) struct n { static const char* name() { return #n; } }
} // tag::

namespace tk {

//! \brief Tagged tuple, allowing tag-based access
//! \details "Tag" here is any type, but mostly an empty struct with a good name
//!   for the data member
//! \tparam List Type list as brigand::list
//! \see https://stackoverflow.com/a/42988897
//! \see https://gist.github.com/underdoeg/4c5c616c1ad4cbb718f787eefcab902d
template< class List >
class SimpTaggedTuple {

  private:
    //! Generate index for every 2nd type of a type list
    template< typename T >
    using is_odd = brigand::size_t< (brigand::index_of<List,T>::value%2) != 0 >;

    //! Partition a type list into two lists with the even and the odd types
    using Pair = brigand::partition< List, brigand::bind<is_odd,brigand::_1> >;

    //! List of member types
    using Data = typename Pair::first_type;

    //! Tuple of member types
    using Tuple = brigand::as_tuple< Data >;

    //! False-overload for detecting if T is a tagged tuple
    template< typename T, typename = std::void_t<> >
    struct is_tagged_tuple_t : std::false_type {};

    //! True-overload for detecting if T is a tagged tuple
    template< typename T >
    struct is_tagged_tuple_t< T, std::void_t< typename T::i_am_tagged_tuple > >
     : std::true_type {};

    //! Member data as a tuple
    Tuple m_members;

  public:
    //! List of key-value pairs
    using PairList = List;

    //! List of keys
    using Keys = typename Pair::second_type;

    //! Typedef defining self for identifying self
    using i_am_tagged_tuple = void;

    //! Acces type in tuple behind tag
    template< typename Tag >
    using TupleElement =
      std::tuple_element_t< brigand::index_of<Keys,Tag>::value, Tuple >;

    //! Query if the type behind Tag is a SimpTaggedTuple
    //! Usage: if constexpr( is_tagged_tuple<Tag>::value ) { ... }
    template< typename Tag >
    using is_tagged_tuple =
      is_tagged_tuple_t< std::decay_t< TupleElement<Tag> > >;

    //! Default constructor
    explicit SimpTaggedTuple() = default;
    //! Initializer constructor
    explicit SimpTaggedTuple( Tuple&& tuple ) : m_members( std::move(tuple) ) {}

    //! Const-ref access to member tuple
    const Tuple& tuple() const { return m_members; }

    //! Const-reference data member accessor of field of tagged tuple at depth
    template< typename Tag, typename... Tags >
    const auto& get() const noexcept {
      constexpr std::size_t idx = brigand::index_of< Keys, Tag >::value;
      if constexpr( is_tagged_tuple<Tag>::value and sizeof...(Tags) != 0 )
        return std::get< idx >( m_members ).template get< Tags... >();
      else
        return std::get< idx >( m_members );
    }

    //! Reference data member accessor of field of tagged tuple at depth
    template< typename Tag, typename... Tags >
    auto& get() noexcept {
      constexpr std::size_t idx = brigand::index_of< Keys, Tag >::value;
      if constexpr( is_tagged_tuple<Tag>::value and sizeof...(Tags) != 0 )
        return std::get< idx >( m_members ).template get< Tags... >();
      else
        return std::get< idx >( m_members );
    }

    //! Operator == between two SimpTaggedTuple objects
    //! \tparam L Type list as brigand::list for other SimpTaggedTuple
    //! \return True if the lhs and rhs equal
    template< typename L >
    bool operator== ( const SimpTaggedTuple< L >& t ) const {
      static_assert( std::is_same_v< L, List >, "Invoking operator== on "
        "SimpTaggedTuple objects with different typelists" );
      static_assert( !brigand::any< List,
        std::is_floating_point<brigand::_1> >::value, "Invoking operator== on "
        "SimpTaggedTuple objects containing a floating point type is unreliable" );
      return m_members == t.tuple();
    }

    //! Operator < between two SimpTaggedTuple objects
    //! \tparam L Type list as brigand::list for other SimpTaggedTuple
    //! \return True if lhs < rhs
    template< typename L >
    bool operator< ( const SimpTaggedTuple< L >& t ) const {
      static_assert( std::is_same_v< L, List >, "Invoking operator< on "
        "SimpTaggedTuple objects with different typelists" );
      return m_members < t.tuple();
    }

    //! Return number of tuple entries
    static constexpr std::size_t size() { return std::tuple_size_v< Tuple >; }

    //! Pack/Unpack
    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    // cppcheck-suppress constParameter
    void pup( PUP::er& p ) { p | m_members; }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] t SimpTaggedTuple object reference
    friend void operator|( PUP::er& p, SimpTaggedTuple<List>& t ) { t.pup(p); }
    //@}
};

} // tk::
