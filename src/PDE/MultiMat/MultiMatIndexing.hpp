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

namespace inciter {

/** @name Functions that compute indices for physics variables for MultiMat */
///@{

// The functions below must follow the signature of MultiMatIdxFn.

//! Get the index of the required material volume fraction
//! \param[in] kmat Index of required material
//! \return Index of the required material volume fraction
inline std::size_t volfracIdx( std::size_t /*nmat*/, std::size_t kmat )
{ return (3 + 3*kmat); }

//! Get the index of the required material continuity equation
// //! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \return Index of the required material continuity equation
inline std::size_t densityIdx( std::size_t /*nmat*/, std::size_t kmat )
{ return (3 + 3*kmat+1); }

//! Get the index of the required momentum equation component
// //! \param[in] nmat Number of materials
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \return Index of the required momentum equation component
inline std::size_t momentumIdx( std::size_t /*nmat*/, std::size_t idir )
{ return idir; }

//! Get the index of the required material total energy equation
// //! \param[in] nmat Number of materials
//! \param[in] kmat Index of required material
//! \return Index of the required material total energy equation
inline std::size_t energyIdx( std::size_t /*nmat*/, std::size_t kmat )
{ return (3 + 3*kmat+2); }

//! Get the index of the required material deformation gradient equation
//! \param[in] nmat Number of materials
//! \param[in] ksld Index of required solid
//! \param[in] i Row-index of required tensor component
//! \param[in] j Column-index of required tensor component
//! \return Index of the required material deformation gradient equation
inline std::size_t deformIdx( std::size_t nmat, std::size_t ksld,
  std::size_t i, std::size_t j )
{ return (2*nmat+3+nmat + 9*(ksld-1)+3*i+j); }

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

//! Get the index of the required material stress component from primitives
//! \param[in] nmat Number of materials
//! \param[in] ksld Index of required solid
//! \param[in] i Index of required stress component
//! \return Index of the required material Cauchy stress component from vector
//!   of primitives
inline std::size_t stressIdx( std::size_t nmat, std::size_t ksld,
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
inline std::size_t deformDofIdx( std::size_t nmat, std::size_t ksld,
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
inline std::size_t stressDofIdx( std::size_t nmat, std::size_t ksld,
  std::size_t i, std::size_t ndof, std::size_t idof )
{ return stressIdx(nmat, ksld, i)*ndof+idof; }

inline bool matExists( tk::real volfrac )
{ return (volfrac > 1e-10) ? true : false; }

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
inline std::size_t newSolidsAccFn( std::size_t kmat,
  std::size_t i, std::size_t j, std::size_t l)
{ return 3*9*kmat+3*(3*i+j)+l; }

//@}

} //inciter::

#endif // MultiMatIndexing_h
