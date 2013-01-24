//******************************************************************************
/*!
  \file      src/Parser/Parser.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 09:48:41 PM MST
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

//! Parser base
class Parser {

  public:
    //! Constructor
    Parser(const string& filename);

    //! Destructor
    virtual ~Parser();

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
    ifstream m_q;                       //!< Q (control) file input stream

    // Include token IDs
    #include <Parser.def.h>
};

} // namespace Quinoa

#endif // Parser_h
