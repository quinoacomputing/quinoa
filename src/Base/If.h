//******************************************************************************
/*!
  \file      src/Base/If.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:01:46 PM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Compile-time type selection
  \details   Compile-time type selection
*/
//******************************************************************************
#ifndef If_h
#define If_h

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

#endif // If_h
