// *****************************************************************************
/*!
  \file      src/Control/Toggle.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Toggle is the base for an Option, doing generic searches
  \details   Toggle is the base for an Option, doing generic searches.
*/
// *****************************************************************************
#ifndef Toggle_h
#define Toggle_h

#include <map>

#include "Exception.h"
#include "StrConvUtil.h"

namespace tk {

//! \brief Toggle is the base for an Option, doing generic searches
//! \details Toggle is templated on an enum type (a strongly typed enum), whose
//!   values are used as keys in maps of associated option values.
//! \see Control/Option for client code examples
//! \author J. Bakosi
template< typename Enum >
class Toggle {

  public:
    using EnumType = Enum;     //! Used to access template typename from outside

    //! Lookup group name
    const std::string& group() const { return groupname; }

    //! Lookup Enum value based on keyword, Enum must exist
    //! \param[in] keyword Keyword to search for
    //! \return Strongly-typed enum value associated to keyword if found
    Enum value( const std::string& keyword ) const {
      auto it = values.find( keyword );
      Assert( it != end(values),
              std::string("Cannot find value for keyword \"") + keyword + "\"" );
      return it->second;
    }

    //! Check if keyword exists
    //! \param[in] keyword Keyword to search for
    //! \return True if keyword if found
    bool exist( const std::string& keyword ) const {
      return values.find( keyword ) != end( values ) ? true : false;
    }

    //! Lookup option name based on Enum
    //! \param[in] val Enum value to search for
    //! \return Keyword if enum value is found
    const std::string& name( Enum val ) const {
      auto it = names.find( val );
      Assert( it != end(names),
              std::string("Cannot find name for value \"") << val << "\"" );
      return it->second;
    }

  protected:
    //! Constructor protected: designed to be used only as a base class
    //! \param[in] g Group name of the option set
    //! \param[in] n Map associating enum values to option names
    //! \param[in] v Map associating keyword strings to enum values
    //! \details Note that all arguments are rvalues, thus the objects gutted
    //!    and this is the only constructor provided.
    explicit Toggle( std::string&& g,
                     std::map< Enum, std::string >&& n,
                     std::map< std::string, Enum >&& v ) :
      groupname( std::move(g) ), names( std::move(n) ), values( std::move(v) )
    { Assert( names.size() == values.size(), "map sizes differ in Toggle" ); }

  private:
    std::string groupname;                     //!< Group name of the option set
    std::map< Enum, std::string > names;       //!< Map of enums -> names
    std::map< std::string, Enum > values;      //!< Map of keywords -> enums
};

} // tk::

#endif // Toggle_h
