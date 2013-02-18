//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.h
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 10:23:02 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mix model
  \details   Homogeneous material mix model
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <limits>

#include <Physics.h>
#include <Type.h>

namespace Quinoa {

class Memory;
class Paradigm;
class Mix;

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    HomMix(Memory* memory,
           Paradigm* paradigm,
           const string& name,
           const control::MixType mix,
           const int& nscalar,
           const int& npar,
           const real time,
           const int echo = 1,
           const int nstep = numeric_limits<int>::max());

    //! Destructor
    virtual ~HomMix();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomMix(const HomMix&) = delete;
    //! Don't permit copy assigment
    HomMix& operator=(const HomMix&) = delete;
    //! Don't permit move constructor
    HomMix(HomMix&&) = delete;
    //! Don't permit move assigment
    HomMix& operator=(HomMix&&) = delete;

    //! Output joint scalar PDF
    void outJPDF();

    Mix* m_mix;                     //!< Mix model object
};

} // namespace Quinoa

#endif // HomMix_h
