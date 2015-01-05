//******************************************************************************
/*!
  \file      src/Base/if.h
  \author    J. Bakosi
  \date      Sat 03 Jan 2015 10:01:28 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Compile-time type selection
  \details   Compile-time type selection
*/
//******************************************************************************
#ifndef if_h
#define if_h

namespace tk {

//! \brief Type selection: if_< Condition, Then, Else >::type
//! \details Selectively defines type 'Then' or 'Else' based on the value of
//!   'Condition'
//! \author David Rodriguez at stackoverflow.com
//! \see http://stackoverflow.com/a/11814074
template < bool Condition, typename Then, typename Else = void >
struct if_ {
   using type = Then;
};
template < typename Then, typename Else >
struct if_< false, Then, Else > {
   using type = Else;
};

} // tk::

#endif // if_h
