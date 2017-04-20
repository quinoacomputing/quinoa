// *****************************************************************************
/*!
  \file      src/Base/TaggedTuple.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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

#include <tuple>
#include "PUPUtil.h"

namespace tk {
//! Tagged tuple allowing tag-based access to tuple members
namespace tuple {

template<typename... Ts> struct typelist {
  template<typename T> using prepend = typelist<T, Ts...>;
};

template<typename T, typename... Ts> struct index;
template<typename T, typename... Ts> struct index<T, T, Ts...> :
  std::integral_constant<int, 0> {};
template<typename T, typename U, typename... Ts> struct index<T, U, Ts...> :
  std::integral_constant<int, index<T, Ts...>::value + 1> {};

template<int n, typename... Ts> struct nth_impl;
template<typename T, typename... Ts> struct nth_impl<0, T, Ts...> {
  using type = T; };
template<int n, typename T, typename... Ts> struct nth_impl<n, T, Ts...> {
  using type = typename nth_impl<n - 1, Ts...>::type; };
template<int n, typename... Ts> using nth = typename nth_impl<n, Ts...>::type;

template<int n, int m, typename... Ts> struct extract_impl;
template<int n, int m, typename T, typename... Ts>
struct extract_impl<n, m, T, Ts...> : extract_impl<n, m - 1, Ts...> {};
template<int n, typename T, typename... Ts>
struct extract_impl<n, 0, T, Ts...> { using types = typename
  extract_impl<n, n - 1, Ts...>::types::template prepend<T>; };
template<int n, int m> struct extract_impl<n, m> {
  using types = typelist<>; };
template<int n, int m, typename... Ts>
  using extract = typename extract_impl<n, m, Ts...>::types;

template<typename S, typename T> struct tt_impl;
template<typename... Ss, typename... Ts>
struct tt_impl<typelist<Ss...>, typelist<Ts...>> : public std::tuple<Ts...> {
  //! Accessor to type of nth element
  template<typename S>
  using nT = nth<index<S, Ss...>::value, Ts...>;
  //! Constructor
  template<typename... Args> tt_impl(Args &&...args) :
    std::tuple<Ts...>(std::forward<Args>(args)...) {}
  //! Const-ref accessor
  template<typename S> constexpr const nT<S>& get() const {
    return std::get<index<S, Ss...>::value>(*this); }
  //! Rvalue accessor
  template<typename S> nT<S>& get() {
    return std::get<index<S, Ss...>::value>(*this); }
  //! Set value by copying source (for lvalues)
  template<typename S> void set(const nT<S>& value) {
    std::get<index<S, Ss...>::value>(*this) = value; }
  //! Set value by moving source (for rvalues)
  template<typename S> void set(nT<S>&& value) {
    std::get<index<S, Ss...>::value>(*this) = std::forward<nT<S>>(value); }
  //! Pack/Unpack
  void pup( PUP::er& p ) { PUP::pup( p, *this ); }
  friend void
  operator|( PUP::er& p, tt_impl<typelist<Ss...>, typelist<Ts...>>& t )
  { t.pup(p); }
};

//! Tagged tuple. Client-side interface. Tagged tuple allowing tag-based access.
//! This is very much like
//! [std::tuple](http://en.cppreference.com/w/cpp/utility/tuple), but instead of
//! having to index the elements by integers, it allows access by a tag, which
//! can be an empty struct with a unique name. Credit goes to
//! ecatmur_at_stackoverflow.com, for more details, see
//! http://stackoverflow.com/questions/13065166/c11-tagged-tuple. For tags, see
//! Control/Tags.h. Tagged tuples are extensively used for transferring data
//! from the parser to an internal data structure in a type-save manner, which
//! is a tagged tuple containing a hierarchy of various containers. As an
//! example on how tagged tuples are used for parsing an input file, see
//! Control/Walker/InputDeck/InputDeck.h. Another way to use a tagged tuple is a
//! compile-time associated container between tags and an arbitrary type. As an
//! example, see rngtest::TestU01Stack::runner.
template<typename... Ts> struct tagged_tuple :
  tt_impl<extract<2, 0, Ts...>, extract<2, 1, Ts...>> {
  //! Constructor
  template<typename... Args> tagged_tuple(Args &&...args) :
    tt_impl<extract<2, 0, Ts...>, extract<2, 1, Ts...>>(
      std::forward<Args>(args)...) {}
  //! Pack/Unpack
  void pup( PUP::er& p ) {
    tt_impl<extract<2, 0, Ts...>, extract<2, 1, Ts...>>::pup(p); }
  friend void operator|( PUP::er& p, tagged_tuple<Ts...>& t ) { t.pup(p); }
};

//! tagged_tuple_size = std::tuple_size/2
template<typename _Tp>
  struct tagged_tuple_size;

template<typename _Tp>
  struct tagged_tuple_size<const _Tp> :
    public std::integral_constant<
           typename std::remove_cv<decltype(std::tuple_size<_Tp>::value)>::type,
           std::tuple_size<_Tp>::value> { };

template<typename _Tp>
  struct tagged_tuple_size<volatile _Tp> :
    public std::integral_constant<
           typename std::remove_cv<decltype(std::tuple_size<_Tp>::value)>::type,
           std::tuple_size<_Tp>::value> { };

template<typename _Tp>
  struct tagged_tuple_size<const volatile _Tp> :
    public std::integral_constant<
           typename std::remove_cv<decltype(std::tuple_size<_Tp>::value)>::type,
           std::tuple_size<_Tp>::value> { };

template<typename... _Elements>
  struct tagged_tuple_size<tagged_tuple<_Elements...>> :
    public std::integral_constant<std::size_t, sizeof...(_Elements)/2> { };

} // tuple::
} // tk::

#endif // TaggedTuple_h
