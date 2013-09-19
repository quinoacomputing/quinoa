//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:28:45 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa control file parser
  \details   Quinoa control file parser
*/
//******************************************************************************
#ifndef QuinoaParser_h
#define QuinoaParser_h

#include <vector>

#include <Parser.h>
#include <QuinoaControl.h>
#include <QuinoaPrint.h>

namespace quinoa {

//! QuinoaParser : Parser
class QuinoaParser : public Parser {

  public:
    //! Constructor
    //! \param[in]    filename  Control file name to read from
    //! \param[in]    print     Pretty printer
    //! \param[inout] control   Control object to put parsed data in
    explicit QuinoaParser(const std::string& filename,
                          const QuinoaPrint& print,
                          QuinoaControl& control)
      : Parser(filename),
        m_print(print),
        m_control(control) {
      m_control.set<ctr::io,ctr::control>(filename);
    }

    //! Destructor
    ~QuinoaParser() noexcept override = default;

    //! Parse quinoa control file
    void parse() override;

    //! Echo problem setup
    void echo() const override;

  private:
    //! Don't permit copy constructor
    QuinoaParser(const QuinoaParser&) = delete;
    //! Don't permit copy assigment
    QuinoaParser& operator=(const QuinoaParser&) = delete;
    //! Don't permit move constructor
    QuinoaParser(QuinoaParser&&) = delete;
    //! Don't permit move assigment
    QuinoaParser& operator=(QuinoaParser&&) = delete;

    //! Make requested statistics unique
    void unique(std::vector<ctr::Product>& statistics);

    const QuinoaPrint& m_print;       //!< Pretty printer

    QuinoaControl& m_control;         //!< Control object
};

} // namespace quinoa

#endif // QuinoaParser_h
