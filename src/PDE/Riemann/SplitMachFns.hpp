// *****************************************************************************
/*!
  \file      src/PDE/Riemann/SplitMachFns.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Split-Mach functions for the AUSM+ Riemann solver
  \details   This file implements the split-Mach functions used by the AUSM
             Riemann solver, with the all-speed corrections.
             Ref. Liou, M. S. (2006). A sequel to AUSM, Part II: AUSM+-up for
             all speeds. Journal of computational physics, 214(1), 137-170.
*/
// *****************************************************************************
#ifndef SplitMachFns_h
#define SplitMachFns_h

namespace inciter {

//! Split Mach polynomials for AUSM+ flux
//! \param[in] fa All-speed parameter
//! \param[in] mach Local Mach numner
//! \return Values of the positive and negative split Mach and pressure
//!   polynomials.
//! \details This function returns a vector with positive and negative Mach
//!   and pressure polynomials, as:
//!   ms[0] = M_4(+),
//!   ms[1] = M_4(-),
//!   ms[2] = P_5(+), and
//!   ms[3] = P_5(-).
//!   For more details, ref. Liou, M. S. (2006). A sequel to AUSM, Part II:
//!   AUSM+-up for all speeds. J. Comp. Phys., 214(1), 137-170.
static std::array< tk::real, 4 > splitmach_ausm( tk::real fa,
                                                 tk::real mach )
{
  std::array< tk::real, 4 > ms;

  std::array< tk::real, 3 > msplus, msminus;
  tk::real psplus, psminus;

  msplus[0] = 0.5*(mach + std::fabs(mach));
  msminus[0]= 0.5*(mach - std::fabs(mach));

  msplus[1] = +0.25*(mach + 1.0)*(mach + 1.0);
  msminus[1]= -0.25*(mach - 1.0)*(mach - 1.0);

  auto alph_fa = 0.0; //(3.0/16.0) * (-4.0 + 5.0*fa*fa);
  auto beta = 0.0; //1.0/8.0;

  if (std::fabs(mach) >= 1.0)
  {
      msplus[2] = msplus[0];
      msminus[2]= msminus[0];
      psplus    = msplus[0]/mach;
      psminus   = msminus[0]/mach;
  }
  else
  {
      msplus[2] = msplus[1]* (1.0 - 16.0*beta*msminus[1]);
      msminus[2]= msminus[1]* (1.0 + 16.0*beta*msplus[1]);
      psplus    = msplus[1]*
                  ((+2.0 - mach) - (16.0 * alph_fa)*mach*msminus[1]);
      psminus   = msminus[1]*
                  ((-2.0 - mach) + (16.0 * alph_fa)*mach*msplus[1]);
  }

  ms[0] = msplus[2];
  ms[1] = msminus[2];
  ms[2] = psplus;
  ms[3] = psminus;

  return ms;
}

} // inciter::

#endif // SplitMachFns_h
