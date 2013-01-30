//******************************************************************************
/*!
  \file      src/Parser/Parser.h
  \author    J. Bakosi
  \date      Tue 29 Jan 2013 09:55:31 PM MST
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

//     //! Store problem title
//     void storeTitle(const string& title);
//     //! Convert physics
//     void convertPhysics(const string& physics);
//     //! Convert hydrodynamics model
//     void convertHydro(const string& hydro);
//     //! Convert material mix model
//     void convertMix(const string& mix);
//     //! Convert number of time steps to take
//     void convertNstep(const string& nstep);
//     //! Convert value at which to stop simulation
//     void convertTerm(const string& term);
//     //! Convert size of time step
//     void convertDt(const string& dt);
//     //! Convert number of mixing scalars
//     void convertNscalar(const string& nscalar);
//     //! Convert total number of particles
//     void convertNpar(const string& npar);
//     //! Convert echo interval
//     void convertEcho(const string& echo);

    const string m_filename;            //!< Name of file to parse
    Control* const m_control;           //!< Main control category
    ifstream m_q;                       //!< Control file input stream
};

} // namespace Quinoa

#endif // Parser_h
