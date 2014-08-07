//******************************************************************************
/*!
  \file      src/MonteCarlo/HomRT.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 05:17:56 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
    explicit HomRT() {
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
