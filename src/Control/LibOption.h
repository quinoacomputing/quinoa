//******************************************************************************
/*!
  \file      src/Control/LibOption.h
  \author    J. Bakosi
  \date      Sat 26 Oct 2013 08:21:20 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Library option
  \details   Library option
*/
//******************************************************************************
#ifndef LibOption_h
#define LibOption_h

#include <Option.h>

namespace tk {

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

} // tk::

#endif // LibOption_h
