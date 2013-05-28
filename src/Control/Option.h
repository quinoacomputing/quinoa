//******************************************************************************
/*!
  \file      src/Control/Option.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 07:11:45 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Option base
  \details   Option base
*/
//******************************************************************************
#ifndef Option_h
#define Option_h

namespace Quinoa {

namespace control {

//! Generic option templated on option 'Type'
template< class Type >
class Option {

  public:
    //! Lookup option value
    typename Type::EnumType option(const std::string& keyword) const {
      return optionType.value(keyword);
    }

    //! Lookup option name
    const std::string& name(typename Type::EnumType value) const {
      return optionType.name(value);
    }

  private:
    //! Permit compiler to generate private copy constructor
    //! Permit compiler to generate private move constructor
    //! Don't permit copy assigment
    Option& operator=(const Option&) = delete;
    //! Don't permit move assigment
    Option& operator=(Option&&) = delete;

    Type optionType;                         //!< Option type class
};

} // namespace control

} // namespace Quinoa

#endif // Option_h
