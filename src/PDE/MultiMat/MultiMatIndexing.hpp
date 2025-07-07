// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/MultiMatIndexing.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Multi-material system indexing functions
  \details   This file defines functions that return indices to specific
    equations for the system of multi-material compressible hydrodynamics.
*/
// *****************************************************************************
#ifndef MultiMatIndexing_h
#define MultiMatIndexing_h

#include "Kokkos_Core.hpp"
using execution_space = Kokkos::Serial;
using memory_space = Kokkos::HostSpace;

namespace inciter {

/** @name Functions that compute indices for physics variables for MultiMat */
///@{

// The functions below must follow the signature of MultiMatIdxFn.

//! Get the index of the required material volume fraction
//! \param[in] kmat Index of required material
//! \return Index of the required material volume fraction
KOKKOS_INLINE_FUNCTION std::size_t volfracIdx( std::size_t /*nmat*/, std::size_t kmat )
{ return kmat; }

//! Get the index of the required material continuity equation
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \return Index of the required material continuity equation
KOKKOS_INLINE_FUNCTION std::size_t densityIdx( std::size_t nmat, std::size_t kmat )
{ return (nmat+kmat); }

//! Get the index of the required momentum equation component
//! \param[in] nmat Number of materials
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \return Index of the required momentum equation component
KOKKOS_INLINE_FUNCTION std::size_t momentumIdx( std::size_t nmat, std::size_t idir )
{ return (2*nmat+idir); }

//! Get the index of the required material total energy equation
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \return Index of the required material total energy equation
KOKKOS_INLINE_FUNCTION std::size_t energyIdx( std::size_t nmat, std::size_t kmat )
{ return (2*nmat+3+kmat); }

//! Get the index of the required material deformation gradient equation
//! \param[in] nmat Number of materials
//! \param[in] ksld Index of required solid
//! \param[in] i Row-index of required tensor component
//! \param[in] j Column-index of required tensor component
//! \return Index of the required material deformation gradient equation
KOKKOS_INLINE_FUNCTION std::size_t deformIdx( std::size_t nmat, std::size_t ksld,
  std::size_t i, std::size_t j )
{ return (2*nmat+3+nmat + 9*(ksld-1)+3*i+j); }

//! Get the index of the required velocity component from vector of primitives
//! \param[in] nmat Number of materials
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \return Index of the required velocity component from vector of primitives
KOKKOS_INLINE_FUNCTION std::size_t velocityIdx( std::size_t nmat, std::size_t idir )
{ return nmat+idir; }

//! Get the index of the required material pressure from vector of primitives
//! \param[in] kmat Index of required material
//! \return Index of the required material pressure from vector of primitives
KOKKOS_INLINE_FUNCTION std::size_t pressureIdx( std::size_t /*nmat*/, std::size_t kmat )
{ return kmat; }

//! Get the index of the required material stress component from primitives
//! \param[in] nmat Number of materials
//! \param[in] ksld Index of required solid
//! \param[in] i Index of required stress component
//! \return Index of the required material Cauchy stress component from vector
//!   of primitives
KOKKOS_INLINE_FUNCTION std::size_t stressIdx( std::size_t nmat, std::size_t ksld,
  std::size_t i )
{ return (nmat+3 + 6*(ksld-1)+i); }

//! \brief Get the index of the required DOF of material volume fraction from
//!   the DG solution vector
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required material volume fraction
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
KOKKOS_INLINE_FUNCTION std::size_t volfracDofIdx( std::size_t nmat, std::size_t kmat,
  std::size_t ndof, std::size_t idof )
{ return volfracIdx(nmat, kmat)*ndof+idof; }

//! \brief Get the index of the required DOF of material continuity equation
//!   from the DG solution vector
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required material continuity equation
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
KOKKOS_INLINE_FUNCTION std::size_t densityDofIdx( std::size_t nmat, std::size_t kmat,
  std::size_t ndof, std::size_t idof )
{ return densityIdx(nmat, kmat)*ndof+idof; }

//! \brief Get the index of the required DOF of momentum equation component from
//!   the DG solution vector
//! \param[in] nmat Number of materials
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required momentum equation component
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
KOKKOS_INLINE_FUNCTION std::size_t momentumDofIdx( std::size_t nmat, std::size_t idir,
  std::size_t ndof, std::size_t idof )
{ return momentumIdx(nmat, idir)*ndof+idof; }

//! \brief Get the index of the required DOF of material total energy equation
//!   from the DG solution vector
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required material total energy equation
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
KOKKOS_INLINE_FUNCTION std::size_t energyDofIdx( std::size_t nmat, std::size_t kmat,
  std::size_t ndof, std::size_t idof )
{ return energyIdx(nmat, kmat)*ndof+idof; }

//! \brief Get the index of the required DOF of material deformation gradient
//!   equation from the DG solution vector
//! \param[in] nmat Number of materials
//! \param[in] ksld Index of required solid
//! \param[in] i Row-index of required tensor component
//! \param[in] j Column-index of required tensor component
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required material total energy equation
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
KOKKOS_INLINE_FUNCTION std::size_t deformDofIdx( std::size_t nmat, std::size_t ksld,
  std::size_t i, std::size_t j, std::size_t ndof, std::size_t idof )
{ return deformIdx(nmat, ksld, i, j)*ndof+idof; }

//! \brief Get the index of the required DOF of velocity component from the DG
//!   vector of primitives
//! \param[in] nmat Number of materials
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required velocity component from vector of primitives
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
KOKKOS_INLINE_FUNCTION std::size_t velocityDofIdx( std::size_t nmat, std::size_t idir,
  std::size_t ndof, std::size_t idof )
{ return velocityIdx(nmat, idir)*ndof+idof; }

//! \brief Get the index of the required DOF of material pressure from the DG
//!   vector of primitives
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required material pressure from vector of primitives
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
KOKKOS_INLINE_FUNCTION std::size_t pressureDofIdx( std::size_t nmat, std::size_t kmat,
  std::size_t ndof, std::size_t idof )
{ return pressureIdx(nmat, kmat)*ndof+idof; }

//! \brief Get the index of the required DOF of material stress component from
//!   the DG vector of primitives
//! \param[in] nmat Number of materials
//! \param[in] ksld Index of required solid
//! \param[in] i Index of required stress component
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required material Cauchy stress component from vector
//!   of primitives
//! \details This function is used to get the index of the required DOF in the
//!   primitives vector, which is of type tk::Fields.
KOKKOS_INLINE_FUNCTION std::size_t stressDofIdx( std::size_t nmat, std::size_t ksld,
  std::size_t i, std::size_t ndof, std::size_t idof )
{ return stressIdx(nmat, ksld, i)*ndof+idof; }

KOKKOS_INLINE_FUNCTION bool matExists( tk::real volfrac )
{ return (volfrac > 1e-10) ? true : false; }

KOKKOS_INLINE_FUNCTION tk::real volfracPRelaxLim()
{ return 1.0e-02; }

//! \brief Get the index of the quantity vel[l]*g[i][j] computed inside the
//!   Riemann flux solver.
//! \param[in] kmat Index of required material
//! \param[in] i Row of inverse deformation tensor
//! \param[in] j Column of inverse deformation tensor
//! \param[in] l Velocity component
//! \return Index of the quantity vel[l]*g[i][j] computed inside the
//!   Riemann flux solver.
//! \details This function is used to get the index of the quantity
//!   vel[l]*g[i][j] computed inside the Riemann flux solver.
KOKKOS_INLINE_FUNCTION std::size_t newSolidsAccFn( std::size_t kmat,
  std::size_t i, std::size_t j, std::size_t l)
{ return 3*9*kmat+3*(3*i+j)+l; }

//! \brief Index for Cauchy stress components, since only the 6 independent
//!   components are stored.
const std::array< std::array< std::size_t, 3 >, 3 > stressCmp{{
  {{0, 3, 4}},
  {{3, 1, 5}},
  {{4, 5, 2}} }};

//! Kokkos version of stressCmp since std::array cant be used
const Kokkos::Array<Kokkos::Array<size_t, 3>, 3> stressCmpKokkos = {{
  {{0, 3, 4}},
  {{3, 1, 5}},
  {{4, 5, 2}} }};

//! Get the index of the required material deformation gradient equation
//! in the context of a list where only the g's of solid materials are present.
//! If one needs to access the deformation tensor within the state array one
//! should use deformIdx instead!
//! \param[in] ksld Index of required solid
//! \param[in] i Row-index of required tensor component
//! \param[in] j Column-index of required tensor component
//! \return Index of the required material deformation gradient equation
//! in the context of a list where only the g's of solid materials are present.
KOKKOS_INLINE_FUNCTION std::size_t solidTensorIdx( std::size_t ksld, std::size_t i, std::size_t j )
{ return 9*ksld+(3*i+j); }


//@}

} //inciter::

#endif // MultiMatIndexing_h
