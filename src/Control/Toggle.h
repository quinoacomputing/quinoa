//******************************************************************************
/*!
  \file      src/Control/Toggle.h
  \author    J. Bakosi
  \date      Fri Aug  2 13:02:04 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Options and associations
  \details   Options and associations
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
           const std::map<std::string, Enum>& v,
           const std::map<Enum, int>* const p = nullptr) :
      names(n), values(v), params(p) {}

    //! Destructor
    virtual ~Toggle() noexcept = default;

    //! Lookup Enum value based on keyword
    Enum value(const std::string keyword) const {
      auto it = values.find(keyword);
      Assert(it != values.end(), ExceptType::FATAL,
            "Cannot find value for keyword \"" + keyword + "\"");
      return it->second;
    }

    //! Lookup option name based on Enum
    const std::string& name(Enum value) const {
      auto it = names.find(value);
      Assert(it != names.end(), ExceptType::FATAL,
             std::string("Cannot find name for value \"") + value + "\"");
      return it->second;
    }

    //! Lookup parameter based on Enum
    int param(Enum value) const {
      Assert(params != nullptr, ExceptType::FATAL,
             "Attempt to search undefined parameter array");
      auto it = params->find(value);
      Assert(it != params->end(), ExceptType::FATAL,
             std::string("Cannot find parameter for value \"") + value + "\"");
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
    const std::map<Enum, int>* const params;
};

// Operators defined outside of class (still in namespace select) to equate
// operator scope with that of enums

//! Operator + for adding Enum to a std::string
template< typename Enum >
std::string operator+ (const std::string& lhs, Enum e) {
  std::stringstream ss;
  // Explicit operator call instead of 'ss << lhs', to avoid gcc's 'ambiguous
  // overload'
  operator<<(ss,lhs) << e;
  std::string rhs = ss.str();
  return rhs;
}

//! Operator << for writing Enum to output streams
template< typename Enum >
std::ostream& operator<< (std::ostream& os, const Enum& e) {
  os << static_cast<unsigned int>(e);
  return os;
}

} // namespace select

} // namespace Quinoa

#endif // Toggle_h
