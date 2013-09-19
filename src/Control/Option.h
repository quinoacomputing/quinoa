//******************************************************************************
/*!
  \file      src/Control/Option.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:24:44 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Option base
  \details   Option base
*/
//******************************************************************************
#ifndef Option_h
#define Option_h

#include <string>

namespace quinoa {
namespace ctr {

//! Generic option interface templated on option 'Type'
template< class Type >
class Option {

  public:
    //! Destructor
    virtual ~Option() noexcept = default;

    //! Lookup option value
    typename Type::EnumType value(const std::string& keyword) const {
      return m_option.value(keyword);
    }

    //! Lookup option name
    const std::string& name(typename Type::EnumType value) const {
      return m_option.name(value);
    }

  protected:
    Type m_option;                         //!< Option

  private:
    //! Permit compiler to generate (public) copy constructor
    //! Permit compiler to generate (public) move constructor
    //! Don't permit copy assigment
    Option& operator=(const Option&) = delete;
    //! Don't permit move assigment
    Option& operator=(Option&&) = delete;
};

} // ctr::
} // quinoa::

#endif // Option_h
