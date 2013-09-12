//******************************************************************************
/*!
  \file      src/Base/Printer.h
  \author    J. Bakosi
  \date      Wed 11 Sep 2013 10:30:31 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Printer
  \details   Printer
*/
//******************************************************************************
#ifndef Printer_h
#define Printer_h

#include <string>
#include <iomanip>

#include <boost/format.hpp>

namespace quinoa {

//! Printer base
class Printer {

  public:
    //! Constructor
    explicit Printer() = default;

    //! Destructor
    virtual ~Printer() noexcept = default;

    //! Print title
    void title(const std::string& title) const {
      boost::format fmt("%|=80|\n");
      std::cout << fmt % boost::io::group(std::setfill('='), "");
      std::cout << fmt % title;
      std::cout << fmt % boost::io::group(std::setfill('='), "");
    }

    //! Print section header: title
    void section(const std::string& title) const {
      std::cout << boost::format("\n%s%c %s:\n")
                   % m_section_indent % m_section_bullet % title;
      std::cout << boost::format("%s%s\n")
                   % m_section_indent
                   % std::string(m_section_indent_size + 2 + title.size(),'-');
    }
    //! Print section header: title : value
    template<typename T>
    void section(const std::string& name, const T& value) const {
      std::cout << boost::format("\n%s%c %s: %s\n")
                   % m_section_indent % m_section_bullet % name % value;
      std::cout << boost::format("%s%s\n")
                   % m_section_indent
                   % std::string(m_section_indent_size + 3 + name.size() +
                                 value.size(), '-');
    }

    //! Print subsection header: title
    void subsection(const std::string& title) const {
      std::cout << boost::format("\n%s%c %s:\n")
                   % m_subsection_indent
                   % m_subsection_bullet
                   % title;
    }

    //! Print item: name
//     void item(const std::string& name) const {
//        std::cout << boost::format("%30s\n") % name;
//     }
    //! Print item: name : value
    template<typename T>
    void item(const std::string& name, const T& value) const {
      std::cout << boost::format("%s%-30s : %s\n") %
                   m_item_indent % name % value;
    }

  protected:
    const char m_section_bullet = '*';
    const char m_subsection_bullet = '-';

    const std::string m_section_indent = " ";
    const std::string m_subsection_indent = m_section_indent + "  ";
    const std::string m_item_indent = m_subsection_indent + "  ";

    const std::string::size_type m_section_indent_size = m_section_indent.size();

  private:
    //! Don't permit copy constructor
    Printer(const Printer&) = delete;
    //! Don't permit copy assigment
    Printer& operator=(const Printer&) = delete;
    //! Don't permit move constructor
    Printer(Printer&&) = delete;
    //! Don't permit move assigment
    Printer& operator=(Printer&&) = delete;
};

} // namespace quinoa

#endif // Printer_h
