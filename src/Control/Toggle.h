//******************************************************************************
/*!
  \file      src/Control/Toggle.h
  \author    J. Bakosi
  \date      Sun 15 Dec 2013 03:57:03 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Options and associations
  \details   Options and associations
*/
//******************************************************************************
#ifndef Toggle_h
#define Toggle_h

#include <map>
#include <sstream>

#include <boost/functional/factory.hpp>

#include <Exception.h>

namespace tk {

// Operators defined outside of class (still in namespace tk) to equate operator
// scope with that of enums

//! Operator << for writing Enum to output streams
template< typename Enum, typename Ch, typename Tr >
std::basic_ostream<Ch,Tr>& operator<< (std::basic_ostream<Ch,Tr>& os, Enum e) {
  os << static_cast<unsigned int>(e);
  return os;
}

//! Operator + for adding Enum to a std::string
template< typename Enum, typename Ch, typename Tr >
std::basic_string<Ch,Tr> operator+ (const std::basic_string<Ch,Tr>& lhs,
                                    Enum e) {
  std::stringstream ss;
  ss << lhs << e;
  return ss.str();
}

template< typename Enum >
class Toggle {

  public:
    using EnumType = Enum;     //! Used to access template typename from outside

    //! Lookup group name
    const std::string& group() const {
      return groupname;
    }

    //! Lookup Enum value based on keyword, Enum must exist
    Enum value(const std::string keyword) const {
      auto it = values.find(keyword);
      Assert(it != values.end(), ExceptType::FATAL,
            "Cannot find value for keyword \"" + keyword + "\"");
      return it->second;
    }

    //! Check if keyword exists
    bool exist(const std::string keyword) const {
      auto it = values.find(keyword);
      if (it != values.end()) return true; else return false;
    }

    //! Lookup option name based on Enum
    const std::string& name(Enum value) const {
      auto it = names.find(value);
      Assert(it != names.end(), ExceptType::FATAL,
             std::string("Cannot find name for value \"") + value + "\"");
      return it->second;
    }

    //! Register a Toggle into a factory
    //! \param[in] C       Type of the (derived) class constructor
    //! \param[in] F       Type of factory to add to
    //! \param[in] Args... Types of variable number of arguments to constructor
    //! \param[in] factory Factory instance
    //! \param[in] e       Enum key to factory's std::map
    //! \param[in] args    Variable number of arguments to constructor
    //! \return    The enum of option as returned
    template< typename C, typename F, typename... Args >
    Enum add( F& factory, Enum e, const Args&... args ) const {
      factory[e] = std::bind( boost::factory<C*>(), std::move(args)... );
      return e;
    }

  protected:
    //! Constructor protected: designed to be used only as a base class
    Toggle(const std::string& g,
           const std::map<Enum, std::string>& n,
           const std::map<std::string, Enum>& v) :
      groupname(g), names(n), values(v) {}

    //! Destructor protected: designed to be used only as a base class
    virtual ~Toggle() noexcept = default;

  private:
    //! Don't permit copy constructor
    Toggle(const Toggle&) = delete;
    //! Don't permit copy assigment
    Toggle& operator=(const Toggle&) = delete;
    //! Don't permit move constructor
    Toggle(Toggle&&) = delete;
    //! Don't permit move assigment
    Toggle& operator=(Toggle&&) = delete;

    const std::string groupname;
    const std::map<Enum, std::string>& names;
    const std::map<std::string, Enum>& values;
};

} // tk::

#endif // Toggle_h
