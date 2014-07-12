//******************************************************************************
/*!
  \file      src/Control/Option.h
  \author    J. Bakosi
  \date      Wed Mar 19 15:55:13 2014
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Option base
  \details   Option base
*/
//******************************************************************************
#ifndef Option_h
#define Option_h

#include <string>

namespace tk {

//! Generic option interface templated on option 'Type'
template< class Type >
class Option {

  public:
    //! Destructor
    virtual ~Option() noexcept = default;

    //! Lookup group name
    const std::string& group() const {
      return m_option.group();
    }

    //! Lookup option value
    typename Type::EnumType value(const std::string& keyword) const {
      return m_option.value(keyword);
    }

    //! Check if keyword exist
    bool exist(const std::string& keyword) const {
      return m_option.exist(keyword);
    }

    //! Lookup option name
    const std::string& name(typename Type::EnumType val) const {
      return m_option.name(val);
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

} // tk::

#endif // Option_h
