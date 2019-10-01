// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/MultiMatIndexing.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Multi-material system indexing functions
  \details   This file defines functions that return indices to specific
    equations for the system of multi-material compressible hydrodynamics.
*/
// *****************************************************************************
#ifndef MultiMatIndexing_h
#define MultiMatIndexing_h

namespace inciter {

//! Get the index of the required material volume fraction
//! \param[in] kmat Index of required material
//! \return Index of the required material volume fraction
inline std::size_t volfracIdx( std::size_t /*nmat*/, std::size_t kmat )
{ return kmat; }

//! Get the index of the required material continuity equation
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \return Index of the required material continuity equation
inline std::size_t densityIdx( std::size_t nmat, std::size_t kmat )
{ return (nmat+kmat); }

//! Get the index of the required momentum equation component
//! \param[in] nmat Number of materials
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \return Index of the required momentum equation component
inline std::size_t momentumIdx( std::size_t nmat, std::size_t idir )
{ return (2*nmat+idir); }

//! Get the index of the required material total energy equation
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \return Index of the required material total energy equation
inline std::size_t energyIdx( std::size_t nmat, std::size_t kmat )
{ return (2*nmat+3+kmat); }

//! Get the index of the required velocity component from vector of primitives
//! \param[in] nmat Number of materials
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \return Index of the required velocity component from vector of primitives
inline std::size_t velocityIdx( std::size_t nmat, std::size_t idir )
{ return nmat+idir; }

//! Get the index of the required material pressure from vector of primitives
//! \param[in] kmat Index of required material
//! \return Index of the required material pressure from vector of primitives
inline std::size_t pressureIdx( std::size_t /*nmat*/, std::size_t kmat )
{ return kmat; }

//! \brief Get the index of the required DOF of material volume fraction from
//!   the DG solution vector
//! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required material volume fraction
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
inline std::size_t volfracDofIdx( std::size_t nmat, std::size_t kmat,
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
inline std::size_t densityDofIdx( std::size_t nmat, std::size_t kmat,
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
inline std::size_t momentumDofIdx( std::size_t nmat, std::size_t idir,
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
inline std::size_t energyDofIdx( std::size_t nmat, std::size_t kmat,
  std::size_t ndof, std::size_t idof )
{ return energyIdx(nmat, kmat)*ndof+idof; }

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
inline std::size_t velocityDofIdx( std::size_t nmat, std::size_t idir,
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
inline std::size_t pressureDofIdx( std::size_t nmat, std::size_t kmat,
  std::size_t ndof, std::size_t idof )
{ return pressureIdx(nmat, kmat)*ndof+idof; }

} //inciter::

#endif // MultiMatIndexing_h
