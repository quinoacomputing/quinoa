//******************************************************************************
/*!
  \file      src/Base/QuinoaPrinter.h
  \author    J. Bakosi
  \date      Wed Sep 11 12:47:08 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's printer
  \details   Quinoa's printer
*/
//******************************************************************************
#ifndef QuinoaPrinter_h
#define QuinoaPrinter_h

#include <Printer.h>

namespace quinoa {

//! QuinoaPrinter : Printer
class QuinoaPrinter : public Printer {

  public:
    //! Constructor
    explicit QuinoaPrinter() noexcept = default;

    //! Destructor
    ~QuinoaPrinter() noexcept override = default;

  private:
    //! Don't permit copy constructor
    QuinoaPrinter(const QuinoaPrinter&) = delete;
    //! Don't permit copy assigment
    QuinoaPrinter& operator=(const QuinoaPrinter&) = delete;
    //! Don't permit move constructor
    QuinoaPrinter(QuinoaPrinter&&) = delete;
    //! Don't permit move assigment
    QuinoaPrinter& operator=(QuinoaPrinter&&) = delete;
};

} // namespace quinoa

#endif // QuinoaPrinter_h
