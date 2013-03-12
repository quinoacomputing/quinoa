//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Mon 11 Mar 2013 06:40:19 PM MDT
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

    //! Ordinary moments accessor
    const real* ordinary() const { return m_ordinary.ptr; }

  private:
    //! Don't permit copy constructor
    Statistics(const Statistics&) = delete;
    //! Don't permit copy assigment
    Statistics& operator=(const Statistics&) = delete;
    //! Don't permit move constructor
    Statistics(Statistics&&) = delete;
    //! Don't permit move assigment
    Statistics& operator=(Statistics&&) = delete;

    Memory* const m_memory;                     //!< Memory object
    Paradigm* const m_paradigm;                 //!< Parallel programming object
    Control* const m_control;                   //!< Control object
    const int m_nthread;                        //!< Number of threads
    const int m_npar;                           //!< Number of particles
    Mix* const m_mix;                           //!< Mix model object
    const vector<control::Product> m_statistics;//!< Requested tatistics

    vector<const real*> m_instantaneous;  //!< Instantaneous variable pointers
    Data<real> m_ordinary;                //!< Ordinary moments
    int m_nord;                           //!< Number of ordinary moments
};

} // namespace Quinoa

#endif // Statistics_h
