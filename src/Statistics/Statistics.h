//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 10:41:47 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics
  \details   Statistics
*/
//******************************************************************************
#ifndef Statistics_h
#define Statistics_h

#include <QuinoaTypes.h>
#include <Memory.h>
#include <Distribution.h>
#include <StatException.h>
#include <ControlTypes.h>

using namespace std;

namespace Quinoa {

class Memory;
class Paradigm;
class Control;
class Mix;

//! Statistics estimator
class Statistics {

  public:
    //! Constructor
    Statistics(Memory* const memory,
               Paradigm* const paradigm,
               Control* const control,
               Mix* const mix);

    //! Destructor
    virtual ~Statistics();

    //! Accumulate statistics
    void accumulate();

    //! Number of ordinary moments accessor
    int nord() const { return m_nord; }

    //! Number of central moments accessor
    int ncen() const { return m_ncen; }

    //! Ordinary moments accessor
    const real* ordinary() const { return m_ordinary.ptr; }

    //! Central moments accessor
    const real* central() const { return m_central.ptr; }

    //! Find out whether ordinary moment is to be plotted
    bool plotOrdinary(const int m) const;

    //! Return the name of ordinary moment
    const string& nameOrdinary(const int m) const;

    //! Return the name of central moment
    const string& nameCentral(const int m) const;

  private:
    //! Don't permit copy constructor
    Statistics(const Statistics&) = delete;
    //! Don't permit copy assigment
    Statistics& operator=(const Statistics&) = delete;
    //! Don't permit move constructor
    Statistics(Statistics&&) = delete;
    //! Don't permit move assigment
    Statistics& operator=(Statistics&&) = delete;

    //! Estimate ordinary moments
    void estimateOrdinary();

    //! Estimate central moments
    void estimateCentral();

    //! Find out whether product only contains ordinary moment terms
    bool ordinary(const vector<control::Term>& product);

    //! Return mean for fluctuation
    int mean(const string name) const;

    //! Convert string to upper case
    string toUpper(const string s) const;

    //! Return true if string is all lower case
    bool isLower(const string s) const;

    Memory* const m_memory;                     //!< Memory object
    const int m_nthread;                        //!< Number of threads
    const int m_npar;                           //!< Number of particles
    Mix* const m_mix;                           //!< Mix model object
    const int m_nscalar;                        //!< Number of mixing scalars
    const vector<control::Product> m_statistics;//!< Requested tatistics

    //! Instantaneous variable pointers for computing ordinary moments
    vector<vector<const real*>> m_instOrd;
    Data<real> m_ordinary;                      //!< Ordinary moments
    vector<bool> m_plotOrdinary;                //!< Whether to plot ord moments
    vector<string> m_nameOrdinary;              //!< Names of ordinary moments
    int m_nord;                                 //!< Number of ordinary moments

    //! Instantaneous variable pointers for computing central moments
    vector<vector<const real*>> m_instCen;
    Data<real> m_central;                       //!< Central moments
    //! Ordinary moments about which to compute central moments
    vector<vector<const real*>> m_center;
    vector<string> m_nameCentral;               //!< Names of central moments
    int m_ncen;                                 //!< Number of central moments
};

} // namespace Quinoa

#endif // Statistics_h
