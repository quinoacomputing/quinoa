//******************************************************************************
/*!
  \file      src/Main/RNGDriver.h
  \author    J. Bakosi
  \date      Mon 19 May 2014 09:29:32 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver with RNGs
  \details   Driver with RNGs
*/
//******************************************************************************
#ifndef RNGDriver_h
#define RNGDriver_h

#include <list>
#include <map>

#include <Options/RNG.h>
#include <RNG.h>
#include <RNGSSE.h>
#include <Driver.h>
#include <tkTypes.h>

namespace tk {

//! Random number generator factory type
using RNGFactory = std::map< tk::ctr::RNGType, std::function< tk::RNG*() > >;

//! RNGDriver
class RNGDriver : public Driver {

  public:
    //! Constructor
    explicit RNGDriver() = default;

    //! Destructor
    virtual ~RNGDriver() noexcept = default;

    //! Execute
    virtual void execute() const = 0;

    //! Register random number generators into factory
    void initRNGFactory( tk::RNGFactory& factory,
                         int nthreads,
                         #ifdef HAS_MKL
                         const tk::ctr::MKLRNGParameters& mklparam,
                         #endif
                         const tk::ctr::RNGSSEParameters& rngsseparam );

  private:
    //! Don't permit copy constructor
    RNGDriver(const RNGDriver&) = delete;
    //! Don't permit assigment constructor
    RNGDriver& operator=(const RNGDriver&) = delete;
    //! Don't permit move constructor
    RNGDriver(RNGDriver&&) = delete;
    //! Don't permit move assignment
    RNGDriver& operator=(RNGDriver&&) = delete;

   #ifdef HAS_MKL
   //! Register MKL RNGs into factory
   void regMKL( tk::RNGFactory& factory,
                int nthreads,
                const tk::ctr::MKLRNGParameters& mklparam );
   #endif

   //! Register RNGSSE RNGs into factory
   void regRNGSSE( tk::RNGFactory& factory,
                   int nthreads,
                   const tk::ctr::RNGSSEParameters& param );
};

} // namespace tk

#endif // RNGDriver_h
