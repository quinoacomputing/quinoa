//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.h
  \author    J. Bakosi
  \date      Wed 28 Aug 2013 07:50:41 PM MDT
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

namespace Quinoa {

//! QuinoaParser : Parser
class QuinoaParser : public Parser {

  public:
    //! Constructor
    //! \param[in]  filename  Control file name to read from
    //! \param[in]  control   Control object to put parsed data in
    explicit QuinoaParser(const std::string& filename,
                          QuinoaControl* const control)
      : Parser(filename),
        m_control(control) {}

    //! Destructor
    virtual ~QuinoaParser() noexcept = default;

    //! Parse quinoa control file
    virtual void parse();

    //! Echo parsed information from quinoa control
    virtual void echo() const;

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

    QuinoaControl* const m_control;     //!< Control object
};

} // namespace Quinoa

#endif // QuinoaParser_h
