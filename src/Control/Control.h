//******************************************************************************
/*!
  \file      src/Control/Control.h
  \author    J. Bakosi
  \date      Sat 16 Aug 2014 09:28:15 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Control base
  \details   Control base
*/
//******************************************************************************
#ifndef Control_h
#define Control_h

#include <string>
#include <sstream>

#include <TaggedTuple.h>

namespace tk {

//! Control : tagged_tuple
template<typename... Ts>
class Control : public tuple::tagged_tuple<Ts...> {

  private:
    //! Short-hand to innherited tagged tuple
    using Tuple = tuple::tagged_tuple<Ts...>;

  public:
    //! Const-ref accessor
    //! TODO: Replace the overloads below with a variadic one
    //! Const-ref accessor to single element at 1st level
    template< typename tag >
    constexpr const typename Tuple::template nT<tag>&
    get() const noexcept {
      return Tuple::template get<tag>();
    }
    //! Const-ref accessor to single element at 2nd level
    template< typename tag, typename subtag >
    constexpr const typename Tuple::template nT<tag>::template nT<subtag>&
    get() const noexcept {
      return Tuple::template get<tag>().template get<subtag>();
    }
    //! Const-ref accessor to single element at 3rd level
    template< typename tag, typename subtag, typename subsubtag >
    constexpr const typename Tuple::template nT<tag>
                                  ::template nT<subtag>
                                  ::template nT<subsubtag>&
    get() const noexcept {
      return Tuple::template get<tag>().
                    template get<subtag>().
                    template get<subsubtag>();
    }

    //! Rvalue accessor
    //! TODO: Replace the overloads below with a variadic one
    //! Rvalue accessor to single element at 1st level
    template< typename tag >
    typename Tuple::template nT<tag>& get() noexcept {
      return Tuple::template get<tag>();
    }
    //! Rvalue accessor to single element at 2nd level
    template< typename tag, typename subtag >
    typename Tuple::template nT<tag>::template nT<subtag>& get() noexcept {
      return Tuple::template get<tag>().template get<subtag>();
    }
    //! Rvalue accessor to single element at 3rd level
    template< typename tag, typename subtag, typename subsubtag >
    typename Tuple::template nT<tag>
                  ::template nT<subtag>
                  ::template nT<subsubtag>& get() noexcept {
      return Tuple::template get<tag>().
                    template get<subtag>().
                    template get<subsubtag>();
    }

    //! Set value
    //! TODO: Replace the overloads below with a variadic one
    //! Set value at slot at tag at 1st level
    template< typename tag >
    void set(const typename Tuple::template nT<tag>& value) noexcept {
      Tuple::template get<tag>() = value;
    }
    //! Set value at slot at tag at 2nd level
    template< typename tag, typename subtag >
    void set(const typename Tuple::template nT<tag>
                                 ::template nT<subtag>& value) noexcept {
      Tuple::template get<tag>().template get<subtag>() = value;
    }
    //! Set value at slot at tag at 3rd level
    template< typename tag, typename subtag, typename subsubtag >
    void set(const typename Tuple::template nT<tag>
                                 ::template nT<subtag>
                                 ::template nT<subsubtag>& value) noexcept {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>() = value;
    }

    //! Convert and move value to slot
    //! TODO: Replace the overloads below with a variadic one
    //! Convert and move value to slot at tag at 1st level
    template< typename tag >
    void store(const std::string& value) noexcept {
      Tuple::template get<tag>() =
        convert<typename Tuple::template nT<tag>>( value );
    }
    //! Convert and move value to slot at tag at 2nd level
    template< typename tag, typename subtag >
    void store(const std::string& value) noexcept {
      Tuple::template get<tag>().
             template get<subtag>() =
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>>( value );
    }
    //! Convert and move value to slot at tag at 3rd level
    template< typename tag, typename subtag, typename subsubtag >
    void store(const std::string& value) noexcept {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>() =
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::template nT<subsubtag>>( value );
    }

    //! Push back value to vector at slot
    //! TODO: Replace the overloads below with a variadic one
    //! Push back value to vector at tag at 1st level without conversion
    template< typename tag >
    void push_back(const typename Tuple::template nT<tag>::value_type& value) {
      Tuple::template get<tag>().push_back( value );
    }
    //! Push back value to vector at tag at 2nd level without conversion
    template< typename tag, typename subtag >
    void push_back(const typename Tuple::template nT<tag>
                                       ::template nT<subtag>
                                       ::value_type& value) {
      Tuple::template get<tag>().
             template get<subtag>().push_back( value );
    }
    //! Push back value to vector at tag at 3rd level without conversion
    template< typename tag, typename subtag, typename subsubtag >
    void push_back(const typename Tuple::template nT<tag>
                                       ::template nT<subtag>
                                       ::template nT<subsubtag>
                                       ::value_type& value) {
      Tuple::template get<tag>().
             template get<subtag>().
             template get<subsubtag>().push_back( value );
    }

    //! Convert and push back value to vector at slot
    //! TODO: Replace the overloads below with a variadic one
    //! Convert and push back value to vector at tag at 1st level
    template< typename tag >
    void store_back(const std::string& value) {
      Tuple::template get<tag>().push_back(
        convert<typename Tuple::template nT<tag>
                              ::value_type>( value ));
    }
    //! Convert and move value to slot at tag at 2nd level
    template< typename tag, typename subtag >
    void store_back(const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().push_back(
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::value_type>( value ));
    }
    //! Convert and move value to slot at tag at 3rd level
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

    //! Convert and push back value to vector of back of vector at slot
    //! TODO: Replace the overloads below with a variadic one
    //! Convert and push back value to vector of back of vector at tag at 1st
    //! level
    template< typename tag >
    void store_back_back(const std::string& value) {
      Tuple::template get<tag>().back().push_back(
        convert<typename Tuple::template nT<tag>
                              ::value_type>( value ) );
    }
    //! Convert and push back value to vector of back of vector at tag at 2nd
    //! level
    template< typename tag, typename subtag >
    void store_back_back(const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>().back().push_back(
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::value_type>( value ) );
    }
    //! Convert and push back value to vector of back of vector at tag at 3rd
    //! level
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

    //! Convert and insert value to field of map at slot
    //! TODO: Replace the overloads below with a variadic one
    //! Convert and insert value to field of map at tag at 1st level
    template< typename key_type, typename field, typename tag >
    void insert_field(const key_type& key, const std::string& value) {
      Tuple::template get<tag>()[ key ].template get<field>() =
        convert<typename Tuple::template nT<tag>
                              ::mapped_type
                              ::template nT<field>>( value );
    }
    //! Convert and insert value to field of map at tag at 2nd level
    template< typename key_type, typename field, typename tag, typename subtag >
    void insert_field(const key_type& key, const std::string& value) {
      Tuple::template get<tag>().
             template get<subtag>()[ key ].template get<field>() =
        convert<typename Tuple::template nT<tag>
                              ::template nT<subtag>
                              ::mapped_type
                              ::template nT<field>>( value );
    }
    //! Convert and insert value to field of map at tag at 3rd level
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

    //! Convert and insert option value to field of map at tag at 2nd level
    template< typename key_type, typename field, typename field_type,
              typename tag, typename subtag >
    void insert_opt(const key_type& key, const field_type& value) {
      Tuple::template get<tag>().
             template get<subtag>()[ key ].template get<field>() = value;
    }

    // convert string to 'type' via std::stringstream
    template< typename type >
    type convert(const std::string& str) {
      std::stringstream ss(str);
      type num;
      ss >> num;
      return num;
    }

    // convert 'type' to string via std::stringstream
    template< typename type >
    std::string convert(const type& val) {
      std::stringstream ss;
      ss << val;
      return ss.str();
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) { Tuple::pup(p); }
    friend void operator|( PUP::er& p, Control<Ts...>& c ) { c.pup(p); }
};

} // tk::

#endif // Control_h
