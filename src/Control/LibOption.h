//******************************************************************************
/*!
  \file      src/Control/LibOption.h
  \author    J. Bakosi
  \date      Sat 24 Aug 2013 10:15:37 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Library option
  \details   Library option
*/
//******************************************************************************
#ifndef LibOption_h
#define LibOption_h

#include <string>

#include <Option.h>

namespace Quinoa {

namespace control {

//! Generic library option interface templated on option 'Type'
template< class Type >
class LibOption : public Option<Type> {

  public:
    //! Destructor
    virtual ~LibOption() noexcept = default;

    //! Lookup library option parameter
    const typename Type::ParamType& param(typename Type::EnumType value) const {
      return Option<Type>::m_option.param(value);
    }

    //! Lookup library option library
    typename Type::LibType lib(typename Type::EnumType value) const {
      return Option<Type>::m_option.lib(value);
    }

  private:
    //! Permit compiler to generate private copy constructor
    //! Permit compiler to generate private move constructor
    //! Don't permit copy assigment
    LibOption& operator=(const LibOption&) = delete;
    //! Don't permit move assigment
    LibOption& operator=(LibOption&&) = delete;
};

} // namespace control

} // namespace Quinoa

#endif // LibOption_h
