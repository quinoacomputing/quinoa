//******************************************************************************
/*!
  \file      src/Model/Model.C
  \author    J. Bakosi
  \date      Mon 13 May 2013 10:15:31 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************

#include <cstring>

#include <Model.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

using namespace Quinoa;

void
Model::initGaussian(real* const particles,
                    int numvar,
                    MKLRndStream* const rndstr,
                    const VSLStreamStatePtr& str,
                    real mean,
                    real rms)
//******************************************************************************
//  Initialize numvar * npar particles with uncorrelated joint Gaussian
//! \author  J. Bakosi
//******************************************************************************
{
  real r[numvar];

  // Generate uncorrelated joint Gaussian for numvar * npar particles
  for (int p=0; p<m_npar; ++p) {
    rndstr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                     str, numvar, r, mean, rms);
    memcpy(particles + p*numvar, r, numvar*sizeof(real));
  }
}
