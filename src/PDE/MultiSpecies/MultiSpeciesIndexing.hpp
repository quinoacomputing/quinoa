// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/MultiSpeciesIndexing.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Multi-species system indexing functions
  \details   This file defines functions that return indices to specific
    equations for the system of multi-species compressible fluid dynamics.
*/
// *****************************************************************************
#ifndef MultiSpeciesIndexing_h
#define MultiSpeciesIndexing_h

namespace inciter {

namespace multispecies {

/** @name Functions that compute indices for physics variables for MultiSpecies */
///@{

// The functions below must follow the signature of MultiSpeciesIdxFn.

//! Get the index of the required species continuity equation
// //! \param[in] nspec Number of species
//! \param[in] kspec Index of required species
//! \return Index of the required species continuity equation
inline std::size_t densityIdx( std::size_t /*nspec*/, std::size_t kspec )
{ return (kspec); }

//! Get the index of the required momentum equation component
//! \param[in] nspec Number of species
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \return Index of the required momentum equation component
inline std::size_t momentumIdx( std::size_t nspec, std::size_t idir )
{ return (nspec+idir); }

//! Get the index of the required mode of the total energy equation
//! \param[in] nspec Number of species
//! \param[in] kmode Index of required energy mode;
//!   0: translational-rotational-vibrational
//! \return Index of the required species total energy equation
inline std::size_t energyIdx( std::size_t nspec, std::size_t kmode )
{ return (nspec+3+kmode); }

//! Get the index of the required mode of the temperature from primitives vector
// //! \param[in] nspec Number of species
//! \param[in] kmode Index of required energy mode;
//!   0: translational-rotational-vibrational
//! \return Index of the required mode of temperature
inline std::size_t temperatureIdx( std::size_t /*nspec*/, std::size_t kmode )
{ return (kmode); }

//! \brief Get the index of the required DOF of species continuity equation
//!   from the DG solution vector
//! \param[in] nspec Number of species
//! \param[in] kspec Index of required species
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required species continuity equation
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
inline std::size_t densityDofIdx( std::size_t nspec, std::size_t kspec,
  std::size_t ndof, std::size_t idof )
{ return densityIdx(nspec, kspec)*ndof+idof; }

//! \brief Get the index of the required DOF of momentum equation component from
//!   the DG solution vector
//! \param[in] nspec Number of species
//! \param[in] idir Required component direction;
//!   0: X-component,
//!   1: Y-component,
//!   2: Z-component.
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required momentum equation component
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
inline std::size_t momentumDofIdx( std::size_t nspec, std::size_t idir,
  std::size_t ndof, std::size_t idof )
{ return momentumIdx(nspec, idir)*ndof+idof; }

//! \brief Get the index of the required DOF of total energy equation
//!   from the DG solution vector
//! \param[in] nspec Number of species
//! \param[in] kmode Index of required species
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required species total energy equation
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
inline std::size_t energyDofIdx( std::size_t nspec, std::size_t kmode,
  std::size_t ndof, std::size_t idof )
{ return energyIdx(nspec, kmode)*ndof+idof; }

//! \brief Get the index of the required DOF of temperature from the DG vector
//!   of primitives
//! \param[in] nspec Number of species
//! \param[in] kmode Index of required species
//! \param[in] ndof Number of solution DOFs stored in DG solution vector
//! \param[in] idof Index of required solution DOF from DG solution vector
//! \return Index of the required species total energy equation
//! \details This function is used to get the index of the required DOF in the
//!   solution vector, which is of type tk::Fields.
inline std::size_t temperatureDofIdx( std::size_t nspec, std::size_t kmode,
  std::size_t ndof, std::size_t idof )
{ return temperatureIdx(nspec, kmode)*ndof+idof; }

//@}

} //multispecies::

} //inciter::

#endif // MultiSpeciesIndexing_h
