//******************************************************************************
/*!
  \file      src/MonteCarlo/HomMix.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 05:06:47 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <Physics.h>

namespace quinoa {

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    explicit HomMix() : Physics() {
      ErrChk( mix(), "No material mix model specified" );
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

    //! Advance
    void advance( tk::real dt );

    //! Output joint scalar PDF
    void outJpdf( tk::real t );
};

} // quinoa::

#endif // HomMix_h
