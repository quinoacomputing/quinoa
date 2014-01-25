//******************************************************************************
/*!
  \file      src/Base/Factory.h
  \author    J. Bakosi
  \date      Sat 25 Jan 2014 03:49:12 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Factory utils
  \details   Factory utils
*/
//******************************************************************************
#ifndef Factory_h
#define Factory_h

#include <list>
#include <functional>

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
void regist( F& f, std::list<E>& reg, const O& o, E e, Args&&... args ) {
  reg.push_back( o.template add<C>( f, e, std::forward<Args>(args)... ) );
}

} // tk::

#endif // Factory_h
