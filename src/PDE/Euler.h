//******************************************************************************
/*!
  \file      src/PDE/Euler.h
  \author    J. Bakosi
  \date      Wed 03 Feb 2016 04:17:29 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Euler equations describing compressible flow
  \details   This file implements the time integration of the Euler equations
    governing compressible fluid flow.
*/
//******************************************************************************
#ifndef Euler_h
#define Euler_h

#include "EulerProblem.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief Euler equations used polymorphically with tk::PDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Problem - problem configuration, see PDE/EulerProblem.h
template< class Problem >
class Euler {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! \brief Constructor
    //! \author J. Bakosi
    explicit Euler( ncomp_t c ) :
      m_c( c ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::euler >(c) )
    {}

    //! Initalize the Euler equations, prepare for time integration
    //! \param[inout] unk Array of unknowns
    //! \author J. Bakosi
    void initialize( tk::MeshNodes& unk ) {
      //! Set initial conditions using initialization policy
      //Problem::template init< tag::euler >( g_inputdeck, unk, m_offset );
    }

    //! \brief Advance unknowns according to the Euler equations
    //! \param[inout] unk Array of unknowns
    //! \param[in] dt Time step size
    //! \param[in] t Physical time
    //! \author J. Bakosi
    void advance( tk::MeshNodes& unk, tk::real dt, tk::real t ) {}

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const ncomp_t m_offset;             //!< Offset PDE operates from
};

} // inciter::

#endif // Euler_h
