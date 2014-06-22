//******************************************************************************
/*!
  \file      src/Base/Print.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 04:29:36 PM MDT
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

#include <Option.h>

namespace tk {

//! Print base
class Print {

  public:
    //! Constructor
    explicit Print( std::ostream& str = std::cout ) : m_stream( str ) {}

    //! Destructor
    virtual ~Print() = default;

    //! Save pointer to stream
    std::streambuf* save() const { return m_stream.rdbuf(); }

    //! Reset stream to streambuf
    std::streambuf* reset( std::streambuf* buf ) const
    { return m_stream.rdbuf( buf ); }

    //! Reset stream to ostream
    std::streambuf* reset( std::ostream& os ) const
    { return m_stream.rdbuf( os.rdbuf() ); }

    //! Define operator << for any type
    template <typename T>
    friend const Print& operator<<( const Print& os, const T& t )
    { os.m_stream << t; return os; }

    //! Define operator << for function pointer taking ostream returning ostream
    friend const Print& operator<<( const Print& os,
                                    std::ostream& (*pf)(std::ostream&) )
    { os.m_stream << pf; return os; }

    //! Print header
    void header( const std::string& title ) const;

    //! Print part header: title
    void part( const std::string& title ) const;

    //! Print section header: title
    void section( const std::string& title ) const {
      m_stream << m_section_title_fmt % m_section_indent
                                      % m_section_bullet
                                      % title;
      m_stream << m_section_underline_fmt
                  % m_section_indent
                  % std::string(m_section_indent.size() + 2 + title.size(),'-');
    }
    //! Print section header: title : value
    template<typename T>
    void section( const std::string& name, const T& value ) const {
      m_stream << m_section_title_value_fmt % m_section_indent
                                            % m_section_bullet
                                            % name
                                            % value;
      m_stream << m_section_underline_fmt
                  % m_section_indent
                  % std::string(m_section_indent.size() + 3 + name.size() +
                                value.size(), '-');
    }

    //! Print subsection header: title
    void subsection( const std::string& title ) const {
      m_stream << m_subsection_title_fmt % m_subsection_indent
                                         % m_subsection_bullet
                                         % title;
    }

    //! Print item: name
    void item( const std::string& name ) const
    { m_stream << m_item_name_fmt % m_item_indent % name; }

    //! Print item: name : value
    template< typename T >
    void item( const std::string& name, const T& value ) const
    { m_stream << m_item_name_value_fmt % m_item_indent % name % value; }

    //! Print list: name: entries...
    void list( const std::string& name,
               const std::list< std::string >& entries ) const {
      if (!entries.empty()) {
        section( name );
        for (auto& e : entries) m_stream << m_list_item_fmt % m_item_indent % e;
      }
    }

    //! Print list: name: option names...
    template< class Option, class Factory >
    void list( const std::string& title, const Factory& factory ) const {
      if (!factory.empty()) {
        section( title );
        Option option;
        for (const auto& f : factory)
          m_stream << m_list_item_fmt % m_item_indent % option.name( f.first );
      }
    }

    //! Print note
    void note( const std::string& msg ) const
    { m_stream << m_note_fmt % m_section_indent % msg; }

    //! Print end of part
    void endpart() const { m_stream << '\n'; }

    //! Print end of subsection
    void endsubsection() const { m_stream << "\n"; }

    //! Print raw
    template< typename T >
    void raw( const T& r ) const { m_stream << r; }

    //! Raw stream access
    std::ostream& stream() const noexcept { return m_stream; }

    //! Implicit conversion to ostream: allows *this in stl algorithms via an
    //! std::osteam_iterator, e.g., std::ostream_iterator< int >( print, ", " )
    operator std::ostream&() const noexcept { return m_stream; }

  protected:
    //! bullets
    const char m_section_bullet = '*';
    const char m_subsection_bullet = '<';
    //! indents
    const std::string m_section_indent = " ";
    const std::string m_subsection_indent =
      std::operator+(m_section_indent,"  ");
    const std::string m_item_indent = std::operator+(m_subsection_indent,"  ");

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

    // Stream object
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
};

} // tk::

#endif // Print_h
