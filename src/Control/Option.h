//******************************************************************************
/*!
  \file      src/Control/Option.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 11:58:49 AM MDT
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
    //! Constructor
    explicit Option() = default;
    //! Destructor
    ~Option() noexcept = default;

    //! Lookup option value
    typename Type::Enum option(const std::string& keyword) const {
      return optionType.value(keyword);
    }

    //! Lookup option name
    const std::string& name(typename Type::Enum value) const {
      return optionType.name(value);
    }

  private:
    //! Don't permit copy constructor
    Option(const Option&) = delete;
    //! Don't permit copy assigment
    Option& operator=(const Option&) = delete;
    //! Don't permit move constructor
    Option(Option&&) = delete;
    //! Don't permit move assigment
    Option& operator=(Option&&) = delete;

    Type optionType;                         //!< Option type class
};

} // namespace control

} // namespace Quinoa

#endif // Option_h
