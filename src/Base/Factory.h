//******************************************************************************
/*!
  \file      src/Base/Factory.h
  \author    J. Bakosi
  \date      Tue Jan 14 09:03:29 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Factory utils
  \details   Factory utils
*/
//******************************************************************************
#ifndef Factory_h
#define Factory_h

#include <list>

namespace tk {

//! Register into factory
//! \param[in] C       Type of the (derived) class constructor
//! \param[in] F       Type of factory to add to
//! \param[in] O       Type of option to add
//! \param[in] E       Type of enum to add
//! \param[in] Args... Types of variable number of arguments to constructor
//! \param[in] f       Factory instance to add to
//! \param[in] reg     List of enums to add enum to
//! \param[in] o       Type of option to add
//! \param[in] e       Enum key to factory's std::map
//! \param[in] args    Variable number of arguments to constructor
template< class C, class F, class O, typename E, typename... Args >
void regist( F& f, std::list<E>& reg, const O& o, E e, const Args&... args ) {
  reg.push_back( o.template add<C>( f, e, std::move(args)... ) );
}

} // tk::

#endif // Factory_h
