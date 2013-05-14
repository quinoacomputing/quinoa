//******************************************************************************
/*!
  \file      src/Parser/Parser.h
  \author    J. Bakosi
  \date      Mon 13 May 2013 08:39:34 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************
#ifndef Parser_h
#define Parser_h

#include <fstream>

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

    //! Echo parsd data common to all physics
    void echoPhysicsCommon() const;

    //! Echo parsed data specific to mass model
    void echoMass() const;

    //! Echo parsed data specific to hydrodynamics model
    void echoHydro() const;

    //! Echo parsed data specific to mix model
    void echoMix() const;

    const string m_filename;            //!< Name of file to parse
    Control* const m_control;           //!< Main control category
    ifstream m_q;                       //!< Control file input stream
};

} // namespace Quinoa

#endif // Parser_h
