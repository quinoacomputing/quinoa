//******************************************************************************
/*!
  \file      src/Base/Print.h
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 07:13:56 PM MST
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

namespace tk {

//! Print base
class Print {

  public:
    //! Constructor
    explicit Print() = default;

    //! Destructor
    virtual ~Print() = default;

    //! Print header
    void header(const std::string& title) const {
      std::cout << m_header_fmt % boost::io::group(std::setfill('='), "");
      std::cout << std::endl;
      std::cout << m_header_fmt % title;
      std::cout << std::endl;
      std::cout << m_header_fmt % boost::io::group(std::setfill('='), "");
    }

    //! Print part header: title
    void part(const std::string& title) const {
      using std::operator+;
      std::string::size_type half_length = title.size()/2 + 1;
      std::string s(half_length, '-');
      std::string underline(s + " o " + s);
      std::string upper(title);
      std::transform(title.begin(), title.end(), upper.begin(), ::toupper);
      upper = "< " + upper + " >";
      std::cout << m_part_fmt % upper;
      std::cout << m_part_underline_fmt % underline;
    }

    //! Print section header: title
    void section(const std::string& title) const {
      std::cout << m_section_title_fmt % m_section_indent
                                       % m_section_bullet
                                       % title;
      std::cout << m_section_underline_fmt
                   % m_section_indent
                   % std::string(m_section_indent_size + 2 + title.size(),'-');
    }
    //! Print section header: title : value
    template<typename T>
    void section(const std::string& name, const T& value) const {
      std::cout << m_section_title_value_fmt % m_section_indent
                                             % m_section_bullet
                                             % name
                                             % value;
      std::cout << m_section_underline_fmt
                   % m_section_indent
                   % std::string(m_section_indent_size + 3 + name.size() +
                                 value.size(), '-');
    }

    //! Print subsection header: title
    void subsection(const std::string& title) const {
      std::cout << m_subsection_title_fmt % m_subsection_indent
                                          % m_subsection_bullet
                                          % title;
    }

    //! Print item: name
    void item(const std::string& name) const {
      std::cout << m_item_name_fmt % m_item_indent % name;
    }

    //! Print item: name : value
    template< typename T >
    void item(const std::string& name, const T& value) const {
      std::cout << m_item_name_value_fmt % m_item_indent % name % value;
    }

    //! Print list: name: entries...
    void list(const std::string& name,
              const std::list< std::string >& entries) const {
      section( name );
      for (const auto& e : entries) {
        std::cout << m_list_item_fmt % m_item_indent % e;
      }
    }

    //! Print list: name: option names...
    template< typename Enum, typename OptionType >
    void list(const std::string& name,
              const OptionType& opt,
              const std::list< Enum >& entries) const {
      section( name );
      for (const auto& e : entries) {
        std::cout << m_list_item_fmt % m_item_indent % opt.name(e);
      }
    }

    //! Print end of part
    void endpart() const { std::cout << "\n\n"; }

    //! Print end of subsection
    void endsubsection() const { std::cout << "\n"; }

    //! Print raw
    template<typename T>
    void raw(const T& raw) const { std::cout << raw; }

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
    mutable format m_item_name_fmt = format("%s%-30s : ");
    mutable format m_item_name_value_fmt = format("%s%-30s : %s\n");
    mutable format m_part_underline_fmt = format("      %|=68|\n");
    mutable format m_section_underline_fmt = format("%s%s\n");

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
