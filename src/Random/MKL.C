//******************************************************************************
/*!
  \file      src/Random/MKL.C
  \author    J. Bakosi
  \date      Sat 10 Nov 2012 10:04:53 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-call wrappers with error handling
  \details   MKL-call wrappers with error handling
*/
//******************************************************************************

#include <iostream>

#include <mkl_vsl.h>

#include <MKL.h>
#include <VSLException.h>

using namespace Quinoa;

void
MKL::uniform(const int& method,
             const VSLStreamStatePtr& stream,
             const int& n,
             real* r,
             const real& a,
             const real& b)
//******************************************************************************
//  Call MKL's vdRngUniform() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] r        Vector of n uniform random numbers
//! \param[in]  a        Left bound
//! \param[in]  b        Right bound
//! \author  J. Bakosi
//******************************************************************************
{
  int vslerr = vdRngUniform(method, stream, n, r, a, b);
  if (vslerr != VSL_STATUS_OK) throw VSLException(FATAL, vslerr);
}

void
MKL::gaussian(const int& method,
              const VSLStreamStatePtr& stream,
              const int& n,
              real* r,
              const real& a,
              const real& b)
//******************************************************************************
//  Call MKL's vdRngGaussian() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] r        Vector of n Gaussian random numbers
//! \param[in]  a        mean
//! \param[in]  b        standard deviation
//! \author  J. Bakosi
//******************************************************************************
{
  int vslerr = vdRngGaussian(method, stream, n, r, a, b);
  if (vslerr != VSL_STATUS_OK) throw VSLException(FATAL, vslerr);
}

void
MKL::gamma(const int& method,
           const VSLStreamStatePtr& stream,
           const int& n,
           real* r,
           const real& alpha,
           const real& a,
           const real& beta)
//******************************************************************************
//  Call MKL's vdRngGamma() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] r        Vector of n gamma-distributed random numbers
//! \param[in]  alpha    Shape
//! \param[in]  a        Displacement
//! \param[in]  beta     Scale factor
//! \author  J. Bakosi
//******************************************************************************
{
  int vslerr = vdRngGamma(method, stream, n, r, alpha, a, beta);
  if (vslerr != VSL_STATUS_OK) throw VSLException(FATAL, vslerr);
}

void
MKL::newStream(VSLStreamStatePtr* stream,
               const int& brng,
               const unsigned int& seed)
//******************************************************************************
//  Call MKL's vslNewStream() and handle error
//! \param[out]  stream  VSL stream state descriptor
//! \param[in]   brng    Index of the basic generator to initialize the stream
//! \param[in]   seed    Initial condition of the stream
//! \author  J. Bakosi
//******************************************************************************
{
  int vslerr = vslNewStream(stream, brng, seed);
  if (vslerr != VSL_STATUS_OK) throw VSLException(FATAL, vslerr);
}

void
MKL::copyStream(VSLStreamStatePtr* newstream,
                const VSLStreamStatePtr& srcstream)
//******************************************************************************
//  Call MKL's vslCopyStream() and handle error
//! \param[out]  newstream  Copied random stream descriptor
//! \param[in]   srcstream  Pointer to the stream state structure to be copied
//! \author  J. Bakosi
//******************************************************************************
{
  int vslerr = vslCopyStream(newstream, srcstream);
  if (vslerr != VSL_STATUS_OK) throw VSLException(FATAL, vslerr);
}

void
MKL::skipAheadStream(VSLStreamStatePtr& stream,
                     const long long int& nskip)
//******************************************************************************
//  Call MKL's vslCopyStream() and handle error
//! \param[in]  stream  Pointer to the stream state structure to which the
//!                     block-splitting method is applied
//! \param[in]  nskip   Number of skipped elements
//! \author  J. Bakosi
//******************************************************************************
{
  int vslerr = vslSkipAheadStream(stream, nskip);
  if (vslerr != VSL_STATUS_OK) throw VSLException(FATAL, vslerr);
}

void
MKL::leapfrogStream(VSLStreamStatePtr& stream,
                    const int& k,
                    const int& nstreams)
//******************************************************************************
//  Call MKL's vslCopyStream() and handle error
//! \param[in]  stream   Pointer to the stream state structure to which leapfrog
//                       method is applied
//! \param[in]  k        Index of the computational node, or stream number
//! \param[in]  nstreams Largest number of computational nodes, or stride
//! \author  J. Bakosi
//******************************************************************************
{
  int vslerr = vslLeapfrogStream(stream, k, nstreams);
  if (vslerr != VSL_STATUS_OK) throw VSLException(FATAL, vslerr);
}
