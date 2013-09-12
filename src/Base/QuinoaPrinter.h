//******************************************************************************
/*!
  \file      src/Base/QuinoaPrinter.h
  \author    J. Bakosi
  \date      Thu 12 Sep 2013 06:57:30 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's printer
  \details   Quinoa's printer
*/
//******************************************************************************
#ifndef QuinoaPrinter_h
#define QuinoaPrinter_h

#include <Printer.h>
#include <QuinoaControl.h>

namespace quinoa {

//! QuinoaPrinter : Printer
class QuinoaPrinter : public Printer {

  public:
    //! Constructor
    explicit QuinoaPrinter() = default;

    //! Destructor
    ~QuinoaPrinter() noexcept override {}

    //! Echo vector of vector of element names
    //! Fields of vector<vector< struct{field, name, plot} >> must exist
    //! See src/Control/QuinoaControlTypes.h for the definitions of operator <<
    //! for outputing Term and vector<Term>, and operator <<= for outputing
    //! requested (i.e., plotted) Term
    template<typename... tags>
    void vecvecNames(const QuinoaControl& ctr,
                     const std::string& msg,
                     const bool req = false) const {
      std::cout << m_item_name_fmt % m_item_indent % msg;
      if (req) for (auto& v : ctr.get<tags...>()) std::cout <<= v;
      else for (auto& v : ctr.get<tags...>()) std::cout << v;
      std::cout << "\n";
    }

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
