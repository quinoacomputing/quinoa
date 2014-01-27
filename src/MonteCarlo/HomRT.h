//******************************************************************************
/*!
  \file      src/MonteCarlo/HomRT.h
  \author    J. Bakosi
  \date      Mon 27 Jan 2014 03:43:50 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Rayleigh-Taylor
  \details   Homogeneous Rayleigh-Taylor
*/
//******************************************************************************
#ifndef HomRT_h
#define HomRT_h

#include <Physics.h>

namespace quinoa {

//! HomRT : Physics
class HomRT : public Physics {

  public:
    //! Constructor
    explicit HomRT( const Base& base ) : Physics( base ) {
      ErrChk( mass(), tk::ExceptType::FATAL, "No mass model specified" );
      ErrChk( hydro(), tk::ExceptType::FATAL, "No hydro model specified" );
    }

    //! Destructor
    ~HomRT() override = default;

    //! Run
    void run() override;

  private:
    //! Don't permit copy constructor
    HomRT(const HomRT&) = delete;
    //! Don't permit copy assigment
    HomRT& operator=(const HomRT&) = delete;
    //! Don't permit move constructor
    HomRT(HomRT&&) = delete;
    //! Don't permit move assigment
    HomRT& operator=(HomRT&&) = delete;

    //! One-liner report
    void reportHeader() const;
    void report(const uint64_t it,
                const uint64_t nstep,
                const tk::real t,
                const tk::real dt,
                const bool wroteJpdf,
                const bool wroteGlob,
                const bool wrotePlot);

    //! Advance
    void advance(tk::real dt);

    //! Output joint scalar PDF
    void outJpdf(const tk::real t);
};

} // quinoa::

#endif // HomRT_h
