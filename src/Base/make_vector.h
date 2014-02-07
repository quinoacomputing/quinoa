//******************************************************************************
/*!
  \file      src/Base/make_vector.h
  \author    J. Bakosi
  \date      Fri 07 Feb 2014 07:44:27 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Convert a variadic template argument pack to boost::mpl::vector
  \details   Convert a variadic template argument pack to boost::mpl::vector
*/
//******************************************************************************
#ifndef make_vector_h
#define make_vector_h

#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_front.hpp>

namespace tk {

// Credit goes to Dmytro Shandyba: http://blog.shandyba.com/2009/12/17/
// converting-variadic-template-arguments-pack-to-boost-mpl-sequence/

//General definition of the helper class
template< typename ...Args > struct make_vector;

// This specialization does the actual job: it splits the whole pack into 2
// parts: one single type T and the rest of types Args... As soon as it is done
// T is added to an mpl::vector. "bottom--up" recursion is used to fetch all
// types
template< class T, typename ...Args >
struct make_vector< T, Args... > {
  typedef typename
  boost::mpl::push_front<typename make_vector<Args...>::type, T>::type type;
};

// This is a specialization for the case when only one type is passed and also
// finishes recursive descent
template< class T >
struct make_vector< T > {
  typedef boost::mpl::vector< T > type;
};

// This one handles the case when no types where passed at all
template<>
struct make_vector<> {
  typedef boost::mpl::vector<> type;
};

} // tk::

#endif // make_vector_h
