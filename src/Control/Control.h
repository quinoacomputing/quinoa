// *****************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Control base contains generic accessors to tagged tuple elements
  \details   Control is a slightly more specialized level of a tagged tuple,
    implementing still very generic accessors to tuple elements at various
    depths at the tuple hierarchy. At this time, max 3 levels are implemented,
    but it would be nice to replace the triple-overloads with a single generic
    one that works at all depths.
*/
// *****************************************************************************
#ifndef Control_h
#define Control_h

#include <string>
#include <sstream>

#include "TaggedTuple.h"
#include "Exception.h"

namespace tk {

//! Control is a slightly more specialized level of a tagged tuple, implementing
//! still very generic accessors to tuple elements at various depths at the
//! tuple hierarchy. At this time, max 3 levels are implemented, but it would be
//! nice to replace the triple-overloads with a single generic one that works at
//! all depths. For an example specialization, i.e., client-side code, see
//! walker::ctr::InputDeck.
//! \author J. Bakosi
template<typename... Ts>
class Control : public tuple::tagged_tuple<Ts...> {

  private:
    //! Short-hand to innherited tagged tuple
    using Tuple = tuple::tagged_tuple<Ts...>;

  public:
    /** @name Const-ref accessors at three different depths */
    ///@{
    //! \brief Const-ref accessor to single element at 1st level
    //! \return A constant reference behind a tag given by the template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    constexpr const typename Tuple::template nT<tag>&
    get() const noexcept {
      return Tuple::template get<tag>();
    }
    //! \brief Const-ref accessor to single element at 2nd level
    //! \return A constant reference behind a tag, subtag given by the template
    //!   arguments
    // TODO Combine the three overloads into a single variadic one
    //! \author J. Bakosi
    template< typename tag, typename subtag >
    constexpr const typename Tuple::template nT<tag>::template nT<subtag>&
    get() const noexcept {
      return Tuple::template get<tag>().template get<subtag>();
    }
    //! \brief Const-ref accessor to single element at 3rd level
    //! \return A constant reference behind a tag, subtag, subsubtag given by
    //! the template arguments
    // TODO Combine the three overloads into a single variadic one
    //! \author J. Bakosi
    template< typename tag, typename subtag, typename subsubtag >
    constexpr const typename Tuple::template nT<tag>
                                  ::template nT<subtag>
                                  ::template nT<subsubtag>&
    get() const noexcept {
      return Tuple::template get<tag>().
                    template get<subtag>().
                    template get<subsubtag>();
    }
    ///@}

    /** @name Rvalue accessors at three different depths */
    ///@{
    //! \brief Rvalue accessor to single element at 1st level
    //! \return An rvalue reference behind a tag given by the template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    typename Tuple::template nT<tag>& get() noexcept {
      return Tuple::template get<tag>();
    }
    //! \brief Rvalue accessor to single element at 2nd level
    //! \return An rvalue reference behind a tag, subtag given by the template
    //!   arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    typename Tuple::template nT<tag>::template nT<subtag>& get() noexcept {
      return Tuple::template get<tag>().template get<subtag>();
    }
    //! \brief Rvalue accessor to single element at 3rd level
    //! \return An rvalue reference behind a tag, subtag, subsubtag given by the
    //!   template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    typename Tuple::template nT<tag>
                  ::template nT<subtag>
                  ::template nT<subsubtag>& get() noexcept {
      return Tuple::template get<tag>().
                    template get<subtag>().
                    template get<subsubtag>();
    }
    ///@}

    /** @name Set value at three different depths */
    ///@{
    //! \brief Set value at slot at tag at 1st level
    //! \param[in] value Value to set behind tag given by the template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    void set(const typename Tuple::template nT<tag>& value) noexcept {
      Tuple::template get<tag>() = value;
    }
    //! \brief Set value at slot at tag at 2nd level
    //! \param[in] value Value to set behind tag and subtag given by the
    //!   template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    void set(const typename Tuple::template nT<tag>
                                 ::template nT<subtag>& value) noexcept {
      Tuple::template get<tag>().template get<subtag>() = value;
    }
    //! \brief Set value at slot at tag at 3rd level
    //! \param[in] value Value to set behind tag, subtag, and subsubtag given by
    //!   the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    void set(const typename Tuple::template nT<tag>
                                 ::template nT<subtag>
                                 ::template nT<subsubtag>& value) noexcept {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>() = value;
    }
    ///@}

    /** @name Convert and set value at tag at three different depths */
    ///@{
    //! \brief Convert and set value at tag at 1st level
    //! \param[in] value Value to convert and set behind tag given by the
    //!   template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    void store(const std::string& value) noexcept {
      Tuple::template get<tag>() =
        convert<typename Tuple::template nT<tag>>( value );
    }
    //! \brief Convert and set value at tag at 2nd level
    //! \param[in] value Value to convert and set behind tag and subtag given by
    //!   the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    void store(const std::string& value) noexcept {
      Tuple::template get<tag>().
             template get<subtag>() =
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>>( value );
    }
    //! \brief Convert and set value at tag at 3rd level
    //! \param[in] value Value to convert and set behind tag, subtag, and
    //    subsubtag given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    void store(const std::string& value) noexcept {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>() =
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::template nT<subsubtag>>( value );
    }
    ///@}

    /** @name Push back value to vector at tag at three different depths */
    ///@{
    //! \brief Push back value to vector at tag at 1st level without conversion
    //! \param[in] value Value to push back behind tag given by the template
    //!   argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    void push_back( const typename Tuple::template nT<tag>::value_type& value =
                          typename Tuple::template nT<tag>::value_type() ) {
      Tuple::template get<tag>().push_back( value );
    }
    //! \brief Push back value to vector at tag at 2nd level without conversion
    //! \param[in] value Value to push back behind tag and subtag given by the
    //!   template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    void push_back( const typename Tuple::template nT<tag>
                                        ::template nT<subtag>
                                        ::value_type& value =
                          typename Tuple::template nT<tag>
                                        ::template nT<subtag>
                                        ::value_type() ) {
      Tuple::template get<tag>().
             template get<subtag>().push_back( value );
    }
    //! \brief Push back value to vector at tag at 3rd level without conversion
    //! \param[in] value Value to push back behind tag, subtag, and subsubtag
    //!   given by the template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    void push_back( const typename Tuple::template nT<tag>
                                        ::template nT<subtag>
                                        ::template nT<subsubtag>
                                        ::value_type& value =
                          typename Tuple::template nT<tag>
                                        ::template nT<subtag>
                                        ::template nT<subsubtag>
                                        ::value_type() ) {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>().push_back( value );
    }
    ///@}

    /** @name Push back value to vector of back of vector at tag at three different depths */
    ///@{
    //! \brief Push back value to vector of back of vector at tag at 1st level
    //!   without conversion.
    //! \details This is similar to store_back_back but performes no conversion.
    //! \param[in] value Value to push back behind tag given by the template
    //!   argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    void push_back_back( const typename Tuple::template nT<tag>
                                             ::value_type::value_type& value =
                          typename Tuple::template nT<tag>
                                         ::value_type::value_type() ) {
      Tuple::template get<tag>().back().push_back( value );
    }
    //! \brief Push back value to vector of back of vector at tag at 2nd level
    //!   without conversion
    //! \details This is similar to store_back_back but performes no conversion.
    //!   no conversion.
    //! \param[in] value Value to push back behind tag and subtag given by the
    //!   template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    void push_back_back( const typename Tuple::template nT<tag>
                                             ::template nT<subtag>
                                             ::value_type::value_type& value =
                          typename Tuple::template nT<tag>
                                        ::template nT<subtag>
                                        ::value_type::value_type() ) {
      Tuple::template get<tag>().
             template get<subtag>().back().push_back( value );
    }
    //! \brief Push back value to vector of back of vector at tag at 3rd level
    //!   without conversion
    //! \details This is similar to store_back_back but performes no conversion.
    //!   no conversion.
    //! \param[in] value Value to push back behind tag, subtag, and subsubtag
    //!   given by the template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    void push_back_back( const typename Tuple::template nT<tag>
                                             ::template nT<subtag>
                                             ::template nT<subsubtag>
                                             ::value_type::value_type& value =
                          typename Tuple::template nT<tag>
                                        ::template nT<subtag>
                                        ::template nT<subsubtag>
                                        ::value_type::value_type() ) {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>().back().push_back( value );
    }
    ///@}

    /** @name Convert and push back value to vector at tag at three different
      * depths */
    ///@{
    //! \brief Convert and push back value to vector at tag at 1st level
    //! \param[in] value Value to convert and push back behind tag given by the
    //!   template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    void store_back(const std::string& value) {
      Tuple::template get<tag>().push_back(
        convert<typename Tuple::template nT<tag>
                              ::value_type>( value ));
    }
    //! \brief Convert and push back value to slot at tag at 2nd level
    //! \param[in] value Value to convert and push back behind tag and subtag
    //!   given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    void store_back(const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().push_back(
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::value_type>( value ));
    }
    //! \brief Convert and push back value to slot at tag at 3rd level
    //! \param[in] value Value to convert and push back behind tag, subtag, and
    //!   subsubtag given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    void store_back(const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>().push_back(
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::template nT<subsubtag>
                              ::value_type>( value ));
    }
    ///@}

    /** @name Convert and push back value to vector of back of vector at tag at three different depths */
    ///@{
    //! \brief Convert and push back value to vector of back of vector at tag at
    //!   1st level
    //! \param[in] value Value to convert and push back behind tag given by the
    //!   template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    void store_back_back(const std::string& value) {
      Tuple::template get<tag>().back().push_back(
        convert<typename Tuple::template nT<tag>
                              ::value_type::value_type>( value ) );
    }
    //! \brief Convert and push back value to vector of back of vector at tag at
    //!   2nd level
    //! \param[in] value Value to convert and push back behind tag and subtag
    //!   given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    void store_back_back(const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().back().push_back(
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::value_type::value_type>( value ) );
    }
    //! \brief Convert and push back value to vector of back of vector at tag at
    //!   3rd level
    //! \param[in] value Value to convert and push back behind tag, subtag, and
    //!   subsubtag given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    void store_back_back(const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>().back().push_back(
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::template nT<subsubtag>
                              ::value_type::value_type>( value ) );
    }
    ///@}

    /** @name Convert and push back value to vector of back of vector of back of vector at tag at three different depths */
    ///@{
    //! \brief Convert and push back value to vector of back of vector of back
    //!    of vector at tag at 1st level
    //! \param[in] value Value to convert and push back behind tag given by the
    //!   template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    void store_back_back_back(const std::string& value) {
      Tuple::template get<tag>().back().back().push_back(
        convert<typename Tuple::template nT<tag>
                              ::value_type::value_type::value_type>( value ) );
    }
    //! \brief Convert and push back value to vector of back of vector of back
    //!    of vector at tag at 2nd level
    //! \param[in] value Value to convert and push back behind tag and subtag
    //!   given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    void store_back_back_back(const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().back().back().push_back(
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::value_type::value_type::value_type>( value ) );
    }
    //! \brief Convert and push back value to vector of back of vector of back
    //!    of vector at tag at 3rd level
    //! \param[in] value Value to convert and push back behind tag, subtag, and
    //!   subsubtag given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    void store_back_back_back(const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>().back().back().push_back(
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::template nT<subsubtag>
                              ::value_type::value_type::value_type>( value ) );
    }
    ///@}

    /** @name Insert key-value pair to map at tag at three different depths */
    ///@{
    //! \brief Insert key-value pair to map at tag at 1st level using
    //!   std::map::operator[]
    //! \param[in] key Key to insert to std::map behind tag given by the
    //!   template argument
    //! \param[in] value Value to insert to std::map behind tag given by the
    //!   template argument; optional argument, if not given, use default
    //!    constructor
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag >
    void insert( const typename Tuple::template nT<tag>
                                     ::key_type& key,
                 const typename Tuple::template nT<tag>
                                     ::mapped_type& value =
                       typename Tuple::template nT<tag>
                                     ::mapped_type() )
    {
      Tuple::template get<tag>()[ key ] = value;
    }
    //! \brief Insert key-value pair to map at tag at 2nd level using
    //!   std::map::operator[]
    //! \param[in] key Key to insert to std::map behind tag and subtag given by
    //!   the template arguments
    //! \param[in] value Value to insert to std::map behind tag and subtag given
    //!   by the template arguments; optional argument, if not given, use default
    //!    constructor
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag >
    void insert( const typename Tuple::template nT<tag>
                                     ::template nT<subtag>
                                     ::key_type& key,
                 const typename Tuple::template nT<tag>
                                     ::template nT<subtag>
                                     ::mapped_type& value =
                       typename Tuple::template nT<tag>
                                     ::template nT<subtag>
                                     ::mapped_type() )
    {
      Tuple::template get<tag>().
             template get<subtag>()[ key ] = value;
    }
    //! \brief Insert key-value pair to map at tag at 3rd level using
    //!   std::map::operator[]
    //! \param[in] key Key to insert to std::map behind tag, subtag, and
    //!   subsubtag given by the template arguments
    //! \param[in] value Value to insert to std::map behind tag, subtag, and
    //!   subsubstag given by the template arguments; optional argument, if not
    //!    given, use default constructor
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename tag, typename subtag, typename subsubtag >
    void insert( const typename Tuple::template nT<tag>
                                     ::template nT<subtag>
                                     ::template nT<subsubtag>
                                     ::key_type& key,
                 const typename Tuple::template nT<tag>
                                     ::template nT<subtag>
                                     ::template nT<subsubtag>
                                     ::mapped_type& value =
                       typename Tuple::template nT<tag>
                                     ::template nT<subtag>
                                     ::template nT<subsubtag>
                                     ::mapped_type() )
    {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>()[ key ] = value;
    }
    ///@}

    /** @name Insert key-value pair with converting value to map at tag at three different depths */
    ///@{
    //! \brief Insert key-value pair with converting value to map at tag at 1st
    //!   level using std::map::operator[]
    //! \details This member function is used to set a value behind a field
    //!   given by the field template argument of a tagged tuple that exist as a
    //!   value of a std::map behind a key at the 1st level of Control object
    //!   given by the tag template argument. The assumed hierarchy is: Control
    //!   (this object) -> tag -> std::map< key_type, tagged_tuple > -> field =
    //!   value. This is similar to insert_opt, but performs conversion.
    //! \param[in] key Key to insert to std::map behind tag given by the
    //!   template argument
    //! \param[in] value Value to insert to std::map behind tag given by the
    //!   template argument
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename key_type, typename field, typename tag >
    void insert_field(const key_type& key, const std::string& value) {
      Tuple::template get<tag>()[ key ].template get<field>() =
        convert<typename Tuple::template nT<tag>
                              ::mapped_type
                              ::template nT<field>>( value );
    }
    //! \brief Insert key-value pair with converting value to map at tag at 2nd
    //!   level using std::map::operator[]
    //! \details This member function is used to set a value behind a field
    //!   given by the field template argument of a tagged tuple that exist as a
    //!   value of a std::map behind a key at the 2nd level of Control object
    //!   given by the tag and subtag template arguments. The assumed hierarchy
    //!   is: Control (this object) -> tag -> subtag -> std::map< key_type,
    //!   tagged_tuple > -> field = value. This is similar to insert_opt, but
    //!   performs conversion.
    //! \param[in] key Key to insert to std::map behind tag and subtag given by
    //!   the template arguments
    //! \param[in] value Value to insert to std::map behind tag and subtag given
    //!   by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename key_type, typename field, typename tag, typename subtag >
    void insert_field(const key_type& key, const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>()[ key ].template get<field>() =
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::mapped_type
                              ::template nT<field>>( value );
    }
    //! \brief Insert key-value pair with converting value to map at tag at 3rd
    //!   level using std::map::operator[]
    //! \details This member function is used to set a value behind a field
    //!   given by the field template argument of a tagged tuple that exist as a
    //!   value of a std::map behind a key at the 3rd level of Control object
    //!   given by the tag, subtag, and subsubtag template arguments. The
    //!   assumed hierarchy is: Control (this object) -> tag -> subtag ->
    //!   subsubtag -> std::map< key_type, tagged_tuple > -> field = value. This
    //!   is similar to insert_opt, but performs conversion.
    //! \param[in] key Key to insert to std::map behind tag, subtag, and
    //!   subsubtag given by the template arguments
    //! \param[in] value Value to insert to std::map behind tag, subtag, and
    //!   subsubtag given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename key_type, typename field, typename tag, typename subtag,
              typename subsubtag >
    void insert_field(const key_type key, const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>()[ key ].template get<field>() =
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::template nT<subsubtag>
                              ::mapped_type
                              ::template nT<field>>( value );
    }
    ///@}

    /** @name Insert key-value pair without conversion of value to map at tag at three different depths */
    ///@{
    //! \brief Insert value to field of tagged tuple behind a key of a map at
    //!   tag at 1st level.
    //! \details This member function is used to set a value behind a field
    //!   given by the field template argument of a tagged tuple that exist as a
    //!   value of a std::map behind a key at the 1st level of Control object
    //!   given by the tag template argument. The assumed hierarchy is: Control
    //!   (this object) -> tag -> subtag -> std::map< key_type, tagged_tuple >
    //!   -> field = value. This is similar to insert_field, but performs no
    //!   conversion.
    //! \param[in] key Key used to access the std::map value using
    //!   std::map::operator[], behind which a type that defines the get()
    //!   member function (e.g., a tagged_tuple) is assumed to exist
    //! \param[in] value Value to insert
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename key_type, typename field, typename field_type,
              typename tag >
    void insert_opt( const key_type& key, const field_type& value ) {
      Tuple::template get<tag>()[ key ].template get<field>() = value;
    }

    //! \brief Insert value to field of tagged tuple behind a key of a map at
    //!   tag at 2nd level.
    //! \details This member function is used to set a value behind a field
    //!   given by the field template argument of a tagged tuple that exist as a
    //!   value of a std::map behind a key at the 2nd level of Control object
    //!   given by the tag and subtag template arguments. The assumed hierarchy
    //!   is: Control (this object) -> tag -> subtag -> std::map< key_type,
    //!   tagged_tuple > -> field = value. This is similar to insert_field, but
    //!   performs no conversion.
    //! \param[in] key Key used to access the std::map value using
    //!   std::map::operator[], behind which a type that defines the get()
    //!   member function (e.g., a tagged_tuple) is assumed to exist
    //! \param[in] value Value to insert
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename key_type, typename field, typename field_type,
              typename tag, typename subtag >
    void insert_opt( const key_type& key, const field_type& value ) {
      Tuple::template get<tag>().
             template get<subtag>()[ key ].template get<field>() = value;
    }

    //! \brief Insert key-value pair with converting value to map at tag at 3rd
    //!   level using std::map::operator[]
    //! \details This member function is used to set a value behind a field
    //!   given by the field template argument of a tagged tuple that exist as a
    //!   value of a std::map behind a key at the 3rd level of Control object
    //!   given by the tag, subtag, and subsubtag template arguments. The
    //!   assumed hierarchy is: Control (this object) -> tag -> subtag ->
    //!   subsubtag -> std::map< key_type, tagged_tuple > -> field = value. This
    //!   is similar to insert_field, but performs no conversion.
    //! \param[in] key Key to insert to std::map behind tag, subtag, and
    //!   subsubtag given by the template arguments
    //! \param[in] value Value to insert to std::map behind tag, subtag, and
    //!   subsubtag given by the template arguments
    //! \author J. Bakosi
    // TODO Combine the three overloads into a single variadic one
    template< typename key_type, typename field, typename field_type,
              typename tag, typename subtag, typename subsubtag >
    void insert_opt( const key_type& key, const field_type& value ) {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>()[ key ].template get<field>() = value;
    }
    ///@}

    //! \brief Convert string to a type given by the template argument using
    //!   std::stringstream
    //! \param[in] str String to convert
    //! \return A value of type given by the template argument
    //! \author J. Bakosi
    template< typename type >
    type convert( const std::string& str ) {
      std::stringstream ss( str );
      type num;
      ss >> num;
      if (ss.fail())
        Throw( "Failed to convert '" + str +
               "' to typeid " + typeid(num).name() );
      return num;
    }

    //! \brief Convert value of type given by the template argument to
    //!   std::string using std::stringstream
    //! \param[in] val Value of type given by the template argument
    //! \return std::string of value converted
    //! \author J. Bakosi
    template< typename type >
    std::string convert( const type& val ) {
      std::stringstream ss;
      ss << val;
      if (ss.fail())
        Throw( "Failed to convert value to string" );
      return ss.str();
    }

    /** @name Pack/Unpack: Serialize Control object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \author J. Bakosi
    void pup( PUP::er& p ) { Tuple::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c Control object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, Control<Ts...>& c ) { c.pup(p); }
    ///@}
};

} // tk::

#endif // Control_h
