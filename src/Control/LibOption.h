//******************************************************************************
/*!
  \file      src/Control/LibOption.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:13:08 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Library option
  \details   Library option
*/
//******************************************************************************
#ifndef LibOption_h
#define LibOption_h

#include <Option.h>

namespace tk {
namespace ctr {

//! Generic library option interface templated on option 'Type'
template< class Type >
class LibOption : public Option<Type> {

  public:
    //! Destructor
    ~LibOption() noexcept override = default;

    //! Lookup library option parameter
    const typename Type::ParamType& param(typename Type::EnumType value) const {
      return Option<Type>::m_option.param(value);
    }

    //! Lookup library option library
    typename Type::LibType lib(typename Type::EnumType value) const {
      return Option<Type>::m_option.lib(value);
    }

  private:
    //! Permit compiler to generate copy constructor
    //! Permit compiler to generate move constructor
    //! Don't permit copy assigment
    LibOption& operator=(const LibOption&) = delete;
    //! Don't permit move assigment
    LibOption& operator=(LibOption&&) = delete;
};

} // ctr::
} // tk::

#endif // LibOption_h
