//******************************************************************************
/*!
  \file      src/Base/Printer.h
  \author    J. Bakosi
  \date      Wed Sep 11 12:50:27 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Printer
  \details   Printer
*/
//******************************************************************************
#ifndef Printer_h
#define Printer_h

#include <string>

namespace quinoa {

//! Printer base
class Printer {

  public:
    //! Constructor
    explicit Printer() noexcept = default;

    //! Destructor
    virtual ~Printer() noexcept = default;

    //! Printer header
    void header(const std::string& title) const;

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
