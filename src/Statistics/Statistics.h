//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 01:46:07 PM MDT
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

  protected:
    Memory* const m_memory;       //!< Memory object
    Paradigm* const m_paradigm;   //!< Parallel programming object
    Control* const m_control;     //!< Control object
    const int m_nthread;          //!< Number of threads

  private:
    //! Don't permit copy constructor
    Statistics(const Statistics&) = delete;
    //! Don't permit copy assigment
    Statistics& operator=(const Statistics&) = delete;
    //! Don't permit move constructor
    Statistics(Statistics&&) = delete;
    //! Don't permit move assigment
    Statistics& operator=(Statistics&&) = delete;

    Mix* const m_mix;                     //!< Mix model object
    const vector<control::Product> m_statistics;   //!< Requested tatistics

    vector<const real*> m_instantaneous;  //!< Instantaneous variable pointers
    Data<real> m_ordinary;                //!< Ordinary moments
    vector<const real*> m_ordinaryThread; //!< Ordinary moments for threads
    int m_nord;                           //!< Number of ordinary moments
};

} // namespace Quinoa

#endif // Statistics_h
