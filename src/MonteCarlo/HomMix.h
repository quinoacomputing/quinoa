//******************************************************************************
/*!
  \file      src/MonteCarlo/HomMix.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:09:41 PM MDT
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
    explicit HomMix( const Base& b ) : Physics( b ) {
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
