//******************************************************************************
/*!
  \file      src/Control/Toggle.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 08:13:20 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Options and associations
  \details   Options and associations
*/
//******************************************************************************
#ifndef Toggle_h
#define Toggle_h

#include <map>

#include <Exception.h>

namespace tk {

template< typename Enum >
class Toggle {

  public:
    using EnumType = Enum;     //! Used to access template typename from outside

    //! Lookup group name
    const std::string& group() const { return groupname; }

    //! Lookup Enum value based on keyword, Enum must exist
    Enum value( const std::string& keyword ) const {
      auto it = values.find( keyword );
      Assert( it != end(values),
              std::string("Cannot find value for keyword \"") + keyword + "\"" );
      return it->second;
    }

    //! Check if keyword exists
    bool exist( const std::string& keyword ) const {
      return values.find( keyword ) != end( values ) ? true : false;
    }

    //! Lookup option name based on Enum
    const std::string& name( Enum val ) const {
      auto it = names.find( val );
      Assert( it != end(names),
              std::string("Cannot find name for value \"") << val << "\"" );
      return it->second;
    }

  protected:
    //! Constructor protected: designed to be used only as a base class
    explicit Toggle( std::string&& g,
                     std::map< Enum, std::string >&& n,
                     std::map< std::string, Enum >&& v ) :
      groupname(std::move(g)), names(std::move(n)), values(std::move(v))
    { Assert( names.size() == values.size(), "map sizes differ in Toggle" ); }

  private:
    std::string groupname;
    std::map< Enum, std::string > names;
    std::map< std::string, Enum > values;
};

} // tk::

#endif // Toggle_h
