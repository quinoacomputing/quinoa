//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.h
  \author    J. Bakosi
  \date      Wed 20 Feb 2013 09:32:23 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mix model
  \details   Homogeneous material mix model
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <map>

#include <Physics.h>
#include <ControlTypes.h>

namespace Quinoa {

class Memory;
class Paradigm;
class Mix;

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    HomMix(Memory* const memory,
           Paradigm* const paradigm,
           Control* const control);

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
    void outJPDF(const real t);

    Mix* m_mix;                         //!< Mix model object
    const string m_jpdf_filename_base;  //!< Joint PDF filename base
};

} // namespace Quinoa

#endif // HomMix_h
