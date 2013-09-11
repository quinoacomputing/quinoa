//******************************************************************************
/*!
  \file      src/Base/Printer.h
  \author    J. Bakosi
  \date      Wed Sep 11 15:48:38 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Printer
  \details   Printer
*/
//******************************************************************************
#ifndef Printer_h
#define Printer_h

#include <string>

#include <boost/format.hpp>

namespace quinoa {

//! Printer base
class Printer {

  public:
    //! Constructor
    explicit Printer() noexcept = default;

    //! Destructor
    virtual ~Printer() noexcept = default;

    //! Print title
    void title(const std::string& title) const;

    //! Print section header
    void section(const std::string& title) const;

    //! Print item: name
    void item(const std::string& name) const {
       using boost::format;
       using boost::io::group;

       std::cout << format("%26s\n") % name;
    }
    //! Print item: name : value
    template<typename T>
    void item(const std::string& name, const T& value) const {
       using boost::format;
       using boost::io::group;

       std::cout << format("%26s : %s\n") % name % value;
    }

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
