//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.h
  \author    J. Bakosi
  \date      Mon 02 Sep 2013 11:35:29 AM MDT
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

    //! Echo parsed data specific to geometry
    void echoGeometry() const;

    //! Echo parsed data specific to physics
    void echoPhysics() const;

    //! Echo parsed data specific to mass model
    void echoMass() const;

    //! Echo parsed data specific to hydrodynamics model
    void echoHydro() const;

    //! Echo parsed data specific to mix model
    void echoMix() const;

    //! Echo parsed data specific to turbulence frequency model
    void echoFrequency() const;

    QuinoaControl& m_control;     //!< Control object
};

} // namespace quinoa

#endif // QuinoaParser_h
