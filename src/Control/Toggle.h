//******************************************************************************
/*!
  \file      src/Control/Toggle.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 08:04:07 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics options and associations
  \details   Physics options and associations
*/
//******************************************************************************
#ifndef Toggle_h
#define Toggle_h

#include <map>
#include <sstream>

#include <Exception.h>

namespace Quinoa {

namespace select {

template< typename Enum >
class Toggle {

  public:
    using EnumType = Enum;     //! Used to access template typename from outside

    //! Constructor
    Toggle(const std::map<Enum, std::string>& n,
           const std::map<std::string, Enum>& v) :
      names(n), values(v) {}

    //! Destructor
    virtual ~Toggle() noexcept = default;

    //! Lookup Enum value based on keyword
    Enum value(const std::string keyword) const {
      auto it = values.find(keyword);
      Assert(it != values.end(), FATAL,
            "Cannot find value for keyword \"" + keyword + "\"");
      return it->second;
    }

    //! Lookup option name based on Enum
    const std::string& name(Enum value) const {
      auto it = names.find(value);
      Assert(it != names.end(), FATAL,
             std::string("Cannot find name for value \"") + value + "\"");
      return it->second;
    }

  private:
    //! Don't permit copy constructor
    Toggle(const Toggle&) = delete;
    //! Don't permit copy assigment
    Toggle& operator=(const Toggle&) = delete;
    //! Don't permit move constructor
    Toggle(Toggle&&) = delete;
    //! Don't permit move assigment
    Toggle& operator=(Toggle&&) = delete;

    const std::map<Enum, std::string>& names;
    const std::map<std::string, Enum>& values;
};

// Operators defined outside of class (still in namespace select) to equate
// operator scope with that of enums

//! Operator + for adding Enum to a std::string
template< typename Enum >
std::string operator+ (const std::string& lhs, Enum e) {
  std::stringstream ss;
  ss << lhs << e ;
  std::string rhs = ss.str();
  return rhs;
}

//! Operator << for writing Enum to output streams
template< typename Enum >
std::ostream& operator<< (std::ostream& os, const Enum& e) {
  os << static_cast<typename std::underlying_type<Enum>::type>(e);
  return os;
}

} // namespace select

} // namespace Quinoa

#endif // Toggle_h
