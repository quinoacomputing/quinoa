//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.h
  \author    J. Bakosi
  \date      Sun 08 Sep 2013 05:27:40 PM MDT
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

namespace quinoa {

//! QuinoaParser : Parser
class QuinoaParser : public Parser {

  public:
    //! Constructor
    //! \param[in]  filename  Control file name to read from
    //! \param[in]  control   Control object to put parsed data in
    explicit QuinoaParser(const std::string& filename, QuinoaControl& control)
      : Parser(filename),
        m_control(control) {}

    //! Destructor
    ~QuinoaParser() noexcept override = default;

    //! Parse quinoa control file
    void parse() override;

    //! Echo parsed information from quinoa control
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
    void unique(std::vector<control::Product>& statistics);

    QuinoaControl& m_control;     //!< Control object
};

} // namespace quinoa

#endif // QuinoaParser_h
