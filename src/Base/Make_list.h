// *****************************************************************************
/*!
  \file      src/Base/Make_list.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Convert a variadic template argument pack to boost::mpl::list
  \details   Convert a variadic template argument pack to boost::mpl::list. For
    more information on the Boost MetaProgramming Library (MPL), see
    http://www.boost.org/doc/libs/release/libs/mpl, in particular, be sure to
    check out mpl::list at
    http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/list.html.
    Credit goes to Dmytro Shandyba: http://blog.shandyba.com/2009/12/17/converting-variadic-template-arguments-pack-to-boost-mpl-sequence.
*/
// *****************************************************************************
#ifndef Make_list_h
#define Make_list_h

#include <boost/mpl/list.hpp>
#include <boost/mpl/push_front.hpp>

namespace tk {

// General definition of the helper class
template< typename ...Args > struct make_list;

// This specialization does the actual job: it splits the whole pack into two
// parts: one single type T and the rest of types Args... As soon as it is done
// T is added to an mpl::list. "bottom--up" recursion is used to fetch all
// types.
template< class T, typename ...Args >
struct make_list< T, Args... > {
  typedef typename
  boost::mpl::push_front<typename make_list<Args...>::type, T>::type type;
};

// This is a specialization for the case when only one type is passed and also
// finishes recursive descent
template< class T >
struct make_list< T > {
  typedef boost::mpl::list< T > type;
};

// This one handles the case when no types where passed at all
template<>
struct make_list<> {
  typedef boost::mpl::list<> type;
};

} // tk::

#endif // Make_list_h
