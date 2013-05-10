//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.h
  \author    J. Bakosi
  \date      Thu May  9 22:23:05 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Standalone-Particle Incompressible Navier-Stokes Flow
  \details   Standalone-Particle Incompressible Navier-Stokes Flow
*/
//******************************************************************************
#ifndef SPINSFlow_h
#define SPINSFlow_h

#include <Physics.h>
#include <Hydro.h>
#include <Timer.h>

namespace Quinoa {

class Memory;
class Paradigm;

//! SPINSFlow : Physics
class SPINSFlow : public Physics {

  public:
    //! Constructor
    explicit SPINSFlow(Memory* const memory,
                       Paradigm* const paradigm,
                       Control* const control,
                       Timer* const timer,
                       const string& filename);

    //! Destructor
    virtual ~SPINSFlow() noexcept = default;

    //! Echo informaion on model
    virtual void echo() const;

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    SPINSFlow(const SPINSFlow&) = delete;
    //! Don't permit copy assigment
    SPINSFlow& operator=(const SPINSFlow&) = delete;
    //! Don't permit move constructor
    SPINSFlow(SPINSFlow&&) = delete;
    //! Don't permit move assigment
    SPINSFlow& operator=(SPINSFlow&&) = delete;

    const string m_filename;        //!< Unstructured mesh file name
};

} // namespace Quinoa

#endif // SPINSFlow_h
