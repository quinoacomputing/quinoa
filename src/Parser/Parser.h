//******************************************************************************
/*!
  \file      src/Parser/Parser.h
  \author    J. Bakosi
  \date      Wed 19 Jun 2013 08:49:58 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************
#ifndef Parser_h
#define Parser_h

#include <fstream>

#include <Option.h>
#include <PhysicsOptions.h>
#include <PositionOptions.h>

using namespace std;

namespace Quinoa {

class Control;

//! Parser base
class Parser {

  public:
    //! Constructor
    explicit Parser(const string& filename, Control* const control);

    //! Destructor
    ~Parser() = default;

    //! Parse
    void parse();

    //! Echo information on stuff parsed
    void echo() const;

  private:
    //! Don't permit copy constructor
    Parser(const Parser&) = delete;
    //! Don't permit copy assigment
    Parser& operator=(const Parser&) = delete;
    //! Don't permit move constructor
    Parser(Parser&&) = delete;
    //! Don't permit move assigment
    Parser& operator=(Parser&&) = delete;

    //! Make requested statistics unique
    void unique(vector<control::Product>& statistics);

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

    const string m_filename;                    //!< Name of file to parse
    Control* const m_control;                   //!< Control category
    ifstream m_q;                               //!< Control file input stream
};

} // namespace Quinoa

#endif // Parser_h
