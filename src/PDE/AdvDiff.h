//******************************************************************************
/*!
  \file      src/PDE/AdvDiff.h
  \author    J. Bakosi
  \date      Wed 03 Feb 2016 03:59:45 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Advection-diffusion equation of a transported scalar
  \details   This file implements the time integration of the
    advection-diffusion equation of a single scalar.
*/
//******************************************************************************
#ifndef AdvDiff_h
#define AdvDiff_h

#include "AdvDiffProblem.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief Advection-diffusion equation used polymorphically with tk::PDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Problem - problem configuration, see PDE/AdvDiffProblem.h
template< class Problem >
class AdvDiff {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! \brief Constructor
    //! \author J. Bakosi
    explicit AdvDiff( ncomp_t c ) :
      m_c( c ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::advdiff >(c) )
    {}

    //! \brief Initalize the advection-diffusion equations, prepare for time
    //!   integration
    //! \param[inout] unk Array of unknowns
    //! \author J. Bakosi
    void initialize( tk::MeshNodes& unk ) {
      //! Set initial conditions using initialization policy
      //Problem::template init< tag::advdiff >( g_inputdeck, unk, m_offset );
    }

    //! \brief Advance unknowns according to the Euler equations
    //! \param[inout] unk Array of unknowns
    //! \param[in] dt Time step size
    //! \param[in] t Physical time
    //! \author J. Bakosi
    void advance( tk::MeshNodes& unk, tk::real dt, tk::real t ) {}

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const ncomp_t m_offset;             //!< Offset this PDE operates from
};

} // inciter::

#endif // AdvDiff_h
