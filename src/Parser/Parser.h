//******************************************************************************
/*!
  \file      src/Parser/Parser.h
  \author    J. Bakosi
  \date      Tue 29 Jan 2013 09:04:59 PM MST
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

    const string m_filename;            //!< Name of file to parse
    Control* const m_control;           //!< Main control category
    ifstream m_q;                       //!< Control file input stream
};

} // namespace Quinoa

#endif // Parser_h
