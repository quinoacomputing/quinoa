//******************************************************************************
/*!
  \file      src/Base/Print.h
  \author    J. Bakosi
  \date      Fri 07 Feb 2014 10:07:22 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Print
  \details   Print
*/
//******************************************************************************
#ifndef Print_h
#define Print_h

#include <iostream>
#include <iomanip>
#include <list>

#include <boost/format.hpp>

#include <tkTypes.h>
#include <Option.h>
#include <Options/RNG.h>

namespace tk {

//! Print base
class Print {

  public:
    //! Constructor
    explicit Print( std::ostream& stream = std::cout ) : m_stream( stream ) {}

    //! Destructor
    virtual ~Print() = default;

    //! Print header
    void header(const std::string& title) const;

    //! Print part header: title
    void part(const std::string& title) const;

    //! Print section header: title
    void section(const std::string& title) const {
      m_stream << m_section_title_fmt % m_section_indent
                                      % m_section_bullet
                                      % title;
      m_stream << m_section_underline_fmt
                  % m_section_indent
                  % std::string(m_section_indent_size + 2 + title.size(),'-');
    }
    //! Print section header: title : value
    template<typename T>
    void section(const std::string& name, const T& value) const {
      m_stream << m_section_title_value_fmt % m_section_indent
                                            % m_section_bullet
                                            % name
                                            % value;
      m_stream << m_section_underline_fmt
                  % m_section_indent
                  % std::string(m_section_indent_size + 3 + name.size() +
                                value.size(), '-');
    }

    //! Print subsection header: title
    void subsection(const std::string& title) const {
      m_stream << m_subsection_title_fmt % m_subsection_indent
                                         % m_subsection_bullet
                                         % title;
    }

    //! Print item: name
    void item(const std::string& name) const {
      m_stream << m_item_name_fmt % m_item_indent % name;
    }

    //! Print item: name : value
    template< typename T >
    void item(const std::string& name, const T& value) const {
      m_stream << m_item_name_value_fmt % m_item_indent % name % value;
    }

    //! Print list: name: entries...
    void list(const std::string& name,
              const std::list< std::string >& entries) const {
      if (!entries.empty()) section( name );
      for (auto& e : entries) {
        m_stream << m_list_item_fmt % m_item_indent % e;
      }
    }

    //! Print list: name: option names...
    template< typename Enum, typename OptionType >
    void list(const std::string& name,
              const OptionType& opt,
              const std::list< Enum >& entries) const {
      if (!entries.empty()) section( name );
      for (auto& e : entries) {
        m_stream << m_list_item_fmt % m_item_indent % opt.name(e);
      }
    }

    //! Print note
    void note( const std::string& msg ) const {
      m_stream << m_note_fmt % m_section_indent % msg;
    }

    //! Print end of part
    void endpart() const { m_stream << "\n\n"; }

    //! Print end of subsection
    void endsubsection() const { m_stream << "\n"; }

    //! Print raw
    template< typename T >
    void raw( const T& raw ) const { m_stream << raw; }

    //! Raw stream access
    std::ostream& stream() const noexcept { return m_stream; }

    //! Print all fields of MKL RNG parameters
    void MKLParams( const std::vector< ctr::RNGType >& vec,
                    const ctr::MKLRNGParameters& map ) const;

    //! Print all fields of RNGSSE parameters
    void RNGSSEParams( const std::vector< ctr::RNGType >& vec,
                       const ctr::RNGSSEParameters& map ) const;

  protected:
    //! bullets
    const char m_section_bullet = '*';
    const char m_subsection_bullet = '<';
    //! indents
    const std::string m_section_indent = " ";
    const std::string m_subsection_indent =
      std::operator+(m_section_indent,"  ");
    const std::string m_item_indent = std::operator+(m_subsection_indent,"  ");
    //! indent sizes
    const std::string::size_type m_section_indent_size =
      m_section_indent.size();

    //! Format strings
    using format = boost::format;
    mutable format m_header_fmt = format("%|=80|\n");
    mutable format m_part_fmt = format("\n%|=80|\n");
    mutable format m_section_title_fmt = format("\n%s%c %s:\n");
    mutable format m_section_title_value_fmt = format("\n%s%c %s: %s\n");
    mutable format m_subsection_title_fmt = format("%s%c %s >\n");
    mutable format m_list_item_fmt = format("%s%-30s\n");
    mutable format m_note_fmt = format("\n%s%-30s\n");
    mutable format m_item_name_fmt = format("%s%-30s : ");
    mutable format m_item_name_value_fmt = format("%s%-30s : %s\n");
    mutable format m_item_widename_value_fmt = format("%s%-65s : %s\n");
    mutable format m_part_underline_fmt = format("      %|=68|\n");
    mutable format m_section_underline_fmt = format("%s%s\n");

    // Steam object
    std::ostream& m_stream;

  private:
    //! Don't permit copy constructor
    Print(const Print&) = delete;
    //! Don't permit copy assigment
    Print& operator=(const Print&) = delete;
    //! Don't permit move constructor
    Print(Print&&) = delete;
    //! Don't permit move assigment
    Print& operator=(Print&&) = delete;

    //! Echo information on MKL random number generator
    void echoMKLParams( const ctr::MKLRNGParam& p ) const;

    //! Echo information on RNGSSE random number generator
    void echoRNGSSEParams( const ctr::RNGSSEParam& p,
                           const ctr::RNG& rng,
                           const ctr::RNGType& r ) const;
};

} // tk::

#endif // Print_h
