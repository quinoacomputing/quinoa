// *****************************************************************************
/*!
  \file      src/Base/TaggedTuple.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Tagged tuple allowing tag-based access
  \details   Tagged tuple allowing tag-based access. This is very much like
    [std::tuple](http://en.cppreference.com/w/cpp/utility/tuple), but instead of
    having to index the elements by integers, it allows access by a tag, which
    can be an empty struct with a unique name. Credit goes to
    ecatmur_at_stackoverflow.com, for more details, see
    http://stackoverflow.com/questions/13065166/c11-tagged-tuple. For tags, see
    Control/Tags.h. Tagged tuples are extensively used for transferring data
    from the parser to an internal data structure in a type-save manner, which
    is a tagged tuple containing a hierarchy of various containers. As an
    example on how tagged tuples are used for parsing an input file, see
    Control/Walker/InputDeck/InputDeck.h. Another way to use a tagged tuple is a
    compile-time associated container between tags and an arbitrary type. As an
    example, see rngtest::TestU01Stack::runner.
*/
// *****************************************************************************
#ifndef TaggedTuple_h
#define TaggedTuple_h

#include <type_traits>
#include <tuple>
#include <sstream>
#include <string>

#include <brigand/adapted/tuple.hpp>
#include "NoWarning/any.h"
#include "NoWarning/partition.h"
#include "NoWarning/index_of.h"

#include "PUPUtil.h"
#include "Exception.h"

namespace tk {

//! \brief Tagged tuple, allowing tag-based access
//! \details "Tag" here is any type, but mostly an empty struct with a good name
//!   for the data member
//! \tparam List Type list as brigand::list
//! \see https://stackoverflow.com/a/42988897
//! \see https://gist.github.com/underdoeg/4c5c616c1ad4cbb718f787eefcab902d
template< class List >
class TaggedTuple{

  private:
    //! Generate index for every 2nd type of a type list
    template< typename T >
    using is_odd = brigand::size_t< (brigand::index_of<List,T>::value%2) != 0 >;

    //! Partition a type list into lists of even and odd types
    using Pair = brigand::partition< List, brigand::bind<is_odd,brigand::_1> >;

    //! List of member types
    using Data = typename Pair::first_type;

    //! List of keys
    using Keys = typename Pair::second_type;

    //! Tuple of member types
    using Tuple = brigand::as_tuple< Data >;

    //! Member data as a tuple
    Tuple m_members;

  public:
    //! Acces type in tuple behing tag
    template< typename Tag >
    using TupleElement =
      std::tuple_element_t< brigand::index_of<Keys,Tag>::value, Tuple >;

    //! Default constructor
    explicit TaggedTuple() = default;
    //! Initializer constructor
    explicit TaggedTuple( Tuple&& tuple ) : m_members( std::move(tuple) ) {}

    //! Const-ref access to member tuple
    const Tuple& tuple() const { return m_members; }
 
    //! Const-ref data member accessor
    template< typename Tag >
    const auto& get() const noexcept {
      return std::get< brigand::index_of<Keys,Tag>::value >( m_members );
    }
    //! Non-const-ref data member accessor
    template< typename Tag >
    auto& get() noexcept {
      return std::get< brigand::index_of<Keys,Tag>::value >( m_members );
    }

    //! Const-ref data member access behind nested tuple
    template< typename Tag1, typename Tag2 >
    auto& get() noexcept {
      return get< Tag1 >().template get< Tag2 >();
    }
    //! Non-const-ref data member access behind nested tuple
    template< typename Tag1, typename Tag2 >
    const auto& get() const noexcept {
      return get< Tag1 >().template get< Tag2 >();
    }

    //! Const-ref data member access behind nested tuple
    template< typename Tag1, typename Tag2, typename Tag3 >
    auto& get() noexcept {
      return get< Tag1 >().template get< Tag2 >().template get< Tag3 >();
    }
    //! Non-const-ref data member access behind nested tuple
    template< typename Tag1, typename Tag2, typename Tag3 >
    const auto& get() const noexcept {
      return get< Tag1 >().template get< Tag2 >().template get< Tag3 >();
    }

    //! Push back value to vector
    //! \param[in] value Value to push back
    template< typename Tag >
    void push_back( const typename TupleElement< Tag >::
                          value_type& value = {} )
    {
      get< Tag >().push_back( value );
    }
    //! Push back value to vector
    //! \param[in] value Value to push back
    template< typename Tag, typename Tag2 >
    void push_back( const typename TupleElement< Tag >::
                          template TupleElement< Tag2 >::
                          value_type& value = {} )
    {
      get< Tag, Tag2 >().push_back( value );
    }
    //! Push back value to vector
    //! \param[in] value Value to push back
    template< typename Tag, typename Tag2, typename Tag3 >
    void push_back( const typename TupleElement< Tag >::
                          template TupleElement< Tag2 >::
                          template TupleElement< Tag3 >::
                          value_type& value = {} )
    {
      get< Tag, Tag2, Tag3 >().push_back( value );
    }

    //! Push back value to back of nested vector
    //! \param[in] value Value to push back
    template< typename Tag >
    void push_back_back( const typename TupleElement< Tag >::
                               value_type::value_type& value = {} )
    {
      get< Tag >().back().push_back( value );
    }
    //! Push back value to back of nested vector
    //! \param[in] value Value to push back
    template< typename Tag, typename Tag2 >
    void push_back_back( const typename TupleElement< Tag >::
                               template TupleElement< Tag2 >::
                               value_type::value_type& value = {} )
    {
      get< Tag, Tag2 >().back().push_back( value );
    }
    //! Push back value to back of nested vector
    //! \param[in] value Value to push back
    template< typename Tag, typename Tag2, typename Tag3 >
    void push_back_back( const typename TupleElement< Tag >::
                               template TupleElement< Tag2 >::
                               template TupleElement< Tag3 >::
                               value_type::value_type& value = {} )
    {
      get< Tag, Tag2, Tag3 >().back().push_back( value );
    }

    //! Convert and store value converting from string
    //! \param[in] value Value to convert and store
    template< typename Tag >
    void store( const std::string& value ) noexcept {
      get< Tag >() = convert< TupleElement<Tag> >( value );
    }
    //! Convert and store value converting from string
    //! \param[in] value Value to convert and store
    template< typename Tag, typename Tag2 >
    void store( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::template TupleElement<Tag2>;
      get< Tag, Tag2 >() = convert< T >( value );
    }
    //! Convert and store value converting from string
    //! \param[in] value Value to convert and store
    template< typename Tag, typename Tag2, typename Tag3 >
    void store( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::
                template TupleElement<Tag2>::
                template TupleElement<Tag3>;
      get< Tag, Tag2, Tag3 >() = convert< T >( value );
    }

    //! Convert and push back value, converting from string, to vector
    //! \param[in] value Value to convert and store
    template< typename Tag >
    void store_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::value_type;
      get< Tag >().push_back( convert< T >( value ) );
    }
    //! Convert and push back value, converting from string, to vector
    //! \param[in] value Value to convert and store
    template< typename Tag, typename Tag2 >
    void store_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::
                template TupleElement<Tag2>::value_type;
      get< Tag, Tag2 >().push_back( convert< T >( value ) );
    }
    //! Convert and push back value, converting from string, to vector
    //! \param[in] value Value to convert and store
    template< typename Tag, typename Tag2, typename Tag3 >
    void store_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::
                template TupleElement<Tag2>::
                template TupleElement<Tag3>::value_type;
      get< Tag, Tag2, Tag3 >().push_back( convert< T >( value ) );
    }

    //! \brief Convert and push back value, converting from string, to back of
    //!   a nested vector
    //! \param[in] value Value to convert and store
    template< typename Tag >
    void store_back_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::value_type::value_type;
      get< Tag >().back().push_back( convert< T >( value ) );
    }
    //! \brief Convert and push back value, converting from string, to back of
    //!   a nested vector
    //! \param[in] value Value to convert and store
    template< typename Tag, typename Tag2 >
    void store_back_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::
                template TupleElement<Tag2>::value_type::value_type;
      get< Tag, Tag2 >().back().push_back( convert< T >( value ) );
    }
    //! \brief Convert and push back value, converting from string, to back of
    //!   a nested vector
    //! \param[in] value Value to convert and store
    template< typename Tag, typename Tag2, typename Tag3 >
    void store_back_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::
                template TupleElement<Tag2>::
                template TupleElement<Tag3>::value_type::value_type;
      get< Tag, Tag2, Tag3 >().back().push_back( convert< T >( value ) );
    }

    //! \brief Convert and push back value, converting from string, to back of
    //!   a doubly nested vector
    //! \param[in] value Value to convert and store
    template< typename Tag >
    void store_back_back_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::value_type::value_type::value_type;
      get< Tag >().back().back().push_back( convert< T >( value ) );
    }
    //! \brief Convert and push back value, converting from string, to back of
    //!   a doubly nested vector
    //! \param[in] value Value to convert and store
    template< typename Tag, typename Tag2 >
    void store_back_back_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::
                template TupleElement<Tag2>::value_type::value_type::value_type;
      get< Tag, Tag2 >().back().back().push_back( convert< T >( value ) );
    }
    //! \brief Convert and push back value, converting from string, to back of
    //!   a doubly nested vector
    //! \param[in] value Value to convert and store
    template< typename Tag, typename Tag2, typename Tag3 >
    void store_back_back_back( const std::string& value ) noexcept {
      using T = typename TupleElement<Tag>::
                template TupleElement<Tag2>::
                template TupleElement<Tag3>::value_type::value_type::value_type;
      get< Tag, Tag2, Tag3 >().back().back().push_back( convert< T >( value ) );
    }

    //! Insert key-value pair to map of nested TaggedTuple
    template< typename FieldTag, typename FieldType,
              typename Tag, typename... Tags, typename Key >
    void insert_field( const Key& key, const FieldType& value ) {
      get< Tag, Tags... >()[ key ].template get< FieldTag >() = value;
    }

    //! \brief Insert key-value pair, converting value from string, to map of
    //!   nested TaggedTuple
    template< typename FieldTag, typename FieldType,
              typename Tag, typename... Tags, typename Key >
    void insert_field( const Key& key, const std::string& value ) {
      get< Tag, Tags... >()[ key ].template get< FieldTag >() =
        convert< FieldType >( value );
    }

    //! Insert key-value pair, converting value from string, to map
    template< typename Tag, typename... Tags, typename Key, typename Value >
    void insert( const Key& key, const Value& value ) {
      get< Tag, Tags... >()[ key ] = value;
    }

    //! Operator == between two TaggedTuple objects
    //! \tparam L Type list as brigand::list for other TaggedTuple
    //! \return True if the lhs and rhs equal
    template< typename L >
    bool operator== ( const TaggedTuple< L >& t ) const {
      static_assert( std::is_same_v< L, List >, "Invoking operator== on "
        "TaggedTuple objects with different typelists" );
      static_assert( !brigand::any< List,
        std::is_floating_point<brigand::_1> >::value, "Invoking operator== on "
        "TaggedTuple objects containing a floating point type is unreliable" );
      return m_members == t.tuple();
    }

    //! Operator < between two TaggedTuple objects
    //! \tparam L Type list as brigand::list for other TaggedTuple
    //! \return True if lhs < rhs
    template< typename L >
    bool operator< ( const TaggedTuple< L >& t ) const {
      static_assert( std::is_same_v< L, List >, "Invoking operator< on "
        "TaggedTuple objects with different typelists" );
      return m_members < t.tuple();
    }

    //! Return number of tuple entries
    static constexpr std::size_t size() { return std::tuple_size_v< Tuple >; }

    //! Pack/Unpack
    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { p | m_members; }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] t TaggedTuple object reference
    friend void operator|( PUP::er& p, TaggedTuple<List>& t ) { t.pup(p); }
    //@}

    //! \brief Convert string to a type given by the template argument using
    //!   std::stringstream
    //! \param[in] str String to convert
    //! \return A value of type given by the template argument
    template< typename type >
    type convert( const std::string& str ) {
      std::stringstream ss( str );
      type num;
      ss >> std::boolalpha >> num;
      if (ss.fail())
        Throw( "Failed to convert '" + str +
               "' to typeid " + typeid(num).name() );
      return num;
    }
};

} // tk::

#endif // TaggedTuple_h
