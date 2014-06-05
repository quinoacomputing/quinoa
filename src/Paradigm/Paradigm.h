//******************************************************************************
/*!
  \file      src/Paradigm/Paradigm.h
  \author    J. Bakosi
  \date      Sun 25 May 2014 06:12:00 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parallel programming paradigms
  \details   Parallel programming paradigms
*/
//******************************************************************************
#ifndef Paradigm_h
#define Paradigm_h

#include <OpenMP.h>
#include <Print.h>

namespace tk {

//! Parallel programming paradigms
class Paradigm {

  public:
    //! Constructor
    explicit Paradigm() = default;
    explicit Paradigm( const Print& print ) { info( print ); }

    //! Destructor
    virtual ~Paradigm() = default;

    //! Output info on compute environment
    void info( const tk::Print& print );

    //! Query if OpenMP is available
    bool ompAvailable() const noexcept { return OpenMP().available(); }

    //! Query if OpenMP is used
    bool ompUsed() const noexcept { return OpenMP().used(); }

    //! Accessor to number of OpenMP threads
    int ompNthreads() const noexcept { return OpenMP().nthreads(); }

  private:
    //! Don't permit copy constructor
    Paradigm(const Paradigm&) = delete;
    //! Don't permit copy assigment
    Paradigm& operator=(const Paradigm&) = delete;
    //! Don't permit move constructor
    Paradigm(Paradigm&&) = delete;
    //! Don't permit move assigment
    Paradigm& operator=(Paradigm&&) = delete;
};

} // tk::

#endif // Paradigm_h
