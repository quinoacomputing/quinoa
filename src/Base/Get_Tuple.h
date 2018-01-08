// *****************************************************************************
/*!
  \file      src/Base/Get_Tuple.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Define std::get<T>(std::tuple) until C++14
  \details   Define std::get<T>(std::tuple) until C++14. When we switch to
    C++14, this can go away. The code below is from
    https://stackoverflow.com/a/16707966.
*/
// *****************************************************************************
#ifndef Get_Tuple_h
#define Get_Tuple_h

#include <type_traits>
#include <tuple>

namespace tk {

#if not defined (_LIBCPP_STD_VER) || (_LIBCPP_STD_VER <= 11)

namespace detail {

template <class T, std::size_t N, class... Args>
struct get_number_of_element_from_tuple_by_type_impl {
  static constexpr auto value = N;
};

template <class T, std::size_t N, class... Args>
struct get_number_of_element_from_tuple_by_type_impl<T, N, T, Args...> {
  static constexpr auto value = N;
};

template <class T, std::size_t N, class U, class... Args>
struct get_number_of_element_from_tuple_by_type_impl<T, N, U, Args...> {
  static constexpr auto value =
    get_number_of_element_from_tuple_by_type_impl<T, N + 1, Args...>::value;
};

} // detail::

template <class T, class... Args>
T get(const std::tuple<Args...>& t) {
  return std::get<detail::
    get_number_of_element_from_tuple_by_type_impl<T, 0, Args...>::value>(t);
}

template <class T, class... Args>
T& get(std::tuple<Args...>& t) {
  return std::get<detail::
    get_number_of_element_from_tuple_by_type_impl<T, 0, Args...>::value>(t);
}

#else

using std::get;

#endif

} // tk::

#endif // Get_Tuple_h
