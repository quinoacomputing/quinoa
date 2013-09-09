//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.h
  \author    J. Bakosi
  \date      Mon Sep  9 08:27:41 2013
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

namespace quinoa {

class Memory;
class Paradigm;

//! SPINSFlow : Physics
class SPINSFlow : public Physics {

  public:
    //! Constructor
    explicit SPINSFlow(Memory* const memory,
                       Paradigm* const paradigm,
                       const QuinoaControl& control,
                       Timer* const timer,
                       const std::string& filename);

    //! Destructor
    ~SPINSFlow() noexcept override = default;

    //! Initialize model
    void init() override;

    //! Solve model
    void solve() override;

  private:
    //! Don't permit copy constructor
    SPINSFlow(const SPINSFlow&) = delete;
    //! Don't permit copy assigment
    SPINSFlow& operator=(const SPINSFlow&) = delete;
    //! Don't permit move constructor
    SPINSFlow(SPINSFlow&&) = delete;
    //! Don't permit move assigment
    SPINSFlow& operator=(SPINSFlow&&) = delete;

    //! Echo information on standalone-particle incompressible Navier-Stokes
    //! physics
    void echo();

    const std::string m_filename;        //!< Unstructured mesh file name
};

} // namespace quinoa

#endif // SPINSFlow_h
