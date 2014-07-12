//******************************************************************************
/*!
  \file      src/Control/Toggle.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:12:51 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Options and associations
  \details   Options and associations
*/
//******************************************************************************
#ifndef Toggle_h
#define Toggle_h

#include <map>
#include <sstream>

#include <Exception.h>

namespace tk {

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
      Assert( it != values.end(),
             "Cannot find value for keyword \"" + keyword + "\"" );
      return it->second;
    }

    //! Check if keyword exists
    bool exist(const std::string keyword) const {
      auto it = values.find(keyword);
      if (it != values.end()) return true; else return false;
    }

    //! Lookup option name based on Enum
    const std::string& name(Enum val) const {
      auto it = names.find(val);
      Assert( it != names.end(),
              std::string("Cannot find name for value \"") + val + "\"" );
      return it->second;
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
