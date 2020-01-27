// *****************************************************************************
/*!
  \file      src/DiffEq/Dissipation/DissipationCoeffPolicy.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Dissipation SDE coefficients policies
  \details   This file defines coefficients policy classes for the dissipation
             SDE, defined in DiffEq/Dissipation/Dissipation.h. For general
             requirements on dissipation SDE coefficients policy classes see the
             header file.
*/
// *****************************************************************************

#include "DissipationCoeffPolicy.hpp"

walker::DissipationCoeffConst::DissipationCoeffConst(
  kw::sde_c3::info::expect::type c3_,
  kw::sde_c4::info::expect::type c4_,
  kw::sde_com1::info::expect::type com1_,
  kw::sde_com2::info::expect::type com2_,
  kw::sde_c3::info::expect::type& c3,
  kw::sde_c4::info::expect::type& c4,
  kw::sde_com1::info::expect::type& com1,
  kw::sde_com2::info::expect::type& com2 )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] c3_ Scalar coefficient used to initialize coefficient c3
//! \param[in] c4_ Scalar coefficient used to initialize coefficient c4
//! \param[in] com1_ Scalar coefficient used to initialize coefficient com1
//! \param[in] com2_ Scalar coefficient used to initialize coefficient com2
//! \param[in,out] c3 Scalar coefficient to be initialized
//! \param[in,out] c4 Scalar coefficient to be initialized
//! \param[in,out] com1 Scalar coefficient to be initialized
//! \param[in,out] com2 Scalar coefficient to be initialized
// *****************************************************************************
{
  c3 = c3_;
  c4 = c4_;
  com1 = com1_;
  com2 = com2_;
}

walker::DissipationCoeffStationary::DissipationCoeffStationary(
  kw::sde_c3::info::expect::type c3_,
  kw::sde_c4::info::expect::type c4_,
  kw::sde_com1::info::expect::type com1_,
  kw::sde_com2::info::expect::type com2_,
  kw::sde_c3::info::expect::type& c3,
  kw::sde_c4::info::expect::type& c4,
  kw::sde_com1::info::expect::type& com1,
  kw::sde_com2::info::expect::type& com2 )
// *****************************************************************************
// Constructor: initialize coefficients
//! \param[in] c3_ Scalar coefficient used to initialize coefficient c3
//! \param[in] c4_ Scalar coefficient used to initialize coefficient c4
//! \param[in] com1_ Scalar coefficient used to initialize coefficient com1
//! \param[in] com2_ Scalar coefficient used to initialize coefficient com2
//! \param[in,out] c3 Scalar coefficient to be initialized
//! \param[in,out] c4 Scalar coefficient to be initialized
//! \param[in,out] com1 Scalar coefficient to be initialized
//! \param[in,out] com2 Scalar coefficient to be initialized
// *****************************************************************************
{
  c3 = c3_;
  c4 = c4_;
  com1 = com1_;
  com2 = com2_;
}
