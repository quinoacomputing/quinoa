//******************************************************************************
/*!
  \file      src/MonteCarlo/HomRT.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:10:04 PM MDT
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
    explicit HomRT( const Base& b ) : Physics( b ) {
      ErrChk( mass(), "No mass model specified" );
      ErrChk( hydro(), "No hydro model specified" );
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

    //! Advance
    void advance(tk::real dt);

    //! Output joint scalar PDF
    void outJpdf(const tk::real t);
};

} // quinoa::

#endif // HomRT_h
