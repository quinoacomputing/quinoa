//******************************************************************************
/*!
  \file      src/Base/make_list.h
  \author    J. Bakosi
  \date      Thu 13 Feb 2014 07:37:49 PM CET
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Convert a variadic template argument pack to boost::mpl::list
  \details   Convert a variadic template argument pack to boost::mpl::list
*/
//******************************************************************************
#ifndef make_list_h
#define make_list_h

#include <boost/mpl/list.hpp>
#include <boost/mpl/push_front.hpp>

namespace tk {

// Credit goes to Dmytro Shandyba: http://blog.shandyba.com/2009/12/17/
// converting-variadic-template-arguments-pack-to-boost-mpl-sequence/

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

#endif // make_list_h
