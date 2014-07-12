//******************************************************************************
/*!
  \file      src/MonteCarlo/HomHydro.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:09:49 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Homogeneous hydrodynamics
  \details   Homogeneous hydrodynamics
*/
//******************************************************************************
#ifndef HomHydro_h
#define HomHydro_h

#include <Physics.h>

namespace quinoa {

//! HomHydro : Physics
class HomHydro : public Physics {

  public:
    //! Constructor
    explicit HomHydro( const Base& b ) : Physics( b ) {
      ErrChk( hydro(), "No hydro model specified" );
    }

    //! Destructor
    ~HomHydro() override = default;

    //! Run
    void run() override;

  private:
    //! Don't permit copy constructor
    HomHydro(const HomHydro&) = delete;
    //! Don't permit copy assigment
    HomHydro& operator=(const HomHydro&) = delete;
    //! Don't permit move constructor
    HomHydro(HomHydro&&) = delete;
    //! Don't permit move assigment
    HomHydro& operator=(HomHydro&&) = delete;
    
    //! Advance
    void advance(tk::real dt);

    //! Output joint PDF
    void outJpdf(const tk::real t);
};

} // quinoa::

#endif // HomHydro_h
