//******************************************************************************
/*!
  \file      src/MonteCarlo/HomMix.h
  \author    J. Bakosi
  \date      Mon 27 Jan 2014 03:40:42 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <Physics.h>
#include <Base.h>

namespace quinoa {

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    explicit HomMix( const Base& base ) : Physics( base ) {
      ErrChk( mix(), tk::ExceptType::FATAL, "No material mix model specified" );
    }

    //! Destructor
    ~HomMix() override = default;

    //! Run
    void run() override;

  private:
    //! Don't permit copy constructor
    HomMix(const HomMix&) = delete;
    //! Don't permit copy assigment
    HomMix& operator=(const HomMix&) = delete;
    //! Don't permit move constructor
    HomMix(HomMix&&) = delete;
    //! Don't permit move assigment
    HomMix& operator=(HomMix&&) = delete;

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

#endif // HomMix_h
