//******************************************************************************
/*!
  \file      src/Random/MKLCall.h
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 10:50:17 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-call wrappers with error handling
  \details   MKL-call wrappers with error handling
*/
//******************************************************************************
#ifndef MKLCall_h
#define MKLCall_h

#include <mkl_vsl.h>

#include <QuinoaTypes.h>

namespace Quinoa {

//! MKL-call wrappers with error handling
class MKLCall {

  public:
    //! Constructor: Default, compiler-generated
    MKLCall() = default;

    //! Destructor: Default, compiler-generated
    ~MKLCall() = default;

    //! Call MKL's vdRngUniform() and handle error
    void uniform(const int& method,
                 VSLStreamStatePtr& stream,
                 const int& n,
                 real* r,
                 const real& a,
                 const real& b);

    //! Call MKL's vdRngGaussian() and handle error
    void gaussian(const int& method,
                  VSLStreamStatePtr& stream,
                  const int& n,
                  real* r,
                  const real& a,
                  const real& b);

    //! Call MKL's vdRngGamma() and handle error
    void gamma(const int& method,
               VSLStreamStatePtr& stream,
               const int& n,
               real* r,
               const real& alpha,
               const real& a,
               const real& beta);

  protected:
    //! Call MKL's vslNewStream() and handle error
    void newStream(VSLStreamStatePtr* stream,
                   const int& brng,
                   const unsigned int& seed);

    //! Call MKL's vslCopyStream() and handle error
    void copyStream(VSLStreamStatePtr* newstream,
                    const VSLStreamStatePtr& srcstream);

    //! Call MKL's vslSkipaheadStream() and handle error
    void skipAheadStream(VSLStreamStatePtr& stream,
                         const long long int& nskip);

    //! Call MKL's vslLeapfrogStream() and handle error
    void leapfrogStream(VSLStreamStatePtr& stream,
                        const int& k,
                        const int& nstreams);

  private:
    //! Don't permit copy constructor
    MKLCall(const MKLCall&) = delete;
    //! Don't permit copy assigment
    MKLCall& operator=(const MKLCall&) = delete;
    //! Don't permit move constructor
    MKLCall(MKLCall&&) = delete;
    //! Don't permit move assigment
    MKLCall& operator=(MKLCall&&) = delete;
};

} // namespace Quinoa

#endif // MKLCall_h
