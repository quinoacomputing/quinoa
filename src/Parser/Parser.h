//******************************************************************************
/*!
  \file      src/Parser/Parser.h
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 10:57:35 AM MDT
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
    Parser(const string& filename, Control* const control);

    //! Destructor
    ~Parser() = default;

    //! Parse
    void parse();

    //! Echo information on stuff parsed
    void echo();

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

    //! Echo material mix block
    void echoMix();

    //! Echo hydrodynamics block
    void echoHydro();

    const string m_filename;            //!< Name of file to parse
    Control* const m_control;           //!< Main control category
    ifstream m_q;                       //!< Control file input stream
};

} // namespace Quinoa

#endif // Parser_h
