// *****************************************************************************
/*!
  \file      src/LinearSolver/ConjugateGradients.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare array for distributed conjugate gradients
  \details   Charm++ chare array for asynchronous distributed
    conjugate gradients linear solver.

    There are a potentially large number of ConjugateGradients Charm++ chares.
    Each ConjugateGradient chare gets a chunk of the full load (due to partiting
    the mesh on which the solve is performed.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality.
*/
// *****************************************************************************
#ifndef ConjugateGradients_h
#define ConjugateGradients_h

#include "Types.hpp"
#include "CSR.hpp"

#include "NoWarning/conjugategradients.decl.h"

namespace tk {

//! \brief ConjugateGradients Charm++ chare array used to perform a distributed
//!   linear solve with the conjugate gradients algorithm
class ConjugateGradients : public CBase_ConjugateGradients {

  public:
    // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
    // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
    //ConjugateGradients_SDAG_CODE

    //! Constructor
    explicit ConjugateGradients(
     std::size_t size,
     std::size_t dof,
     const std::pair< std::vector< std::size_t >,
                      std::vector< std::size_t > >& psup,
     const std::vector< tk::real >& b );

    //! Migrate constructor
    //explicit ConjugateGradients( CkMigrateMessage* ) {}

    /** @name Pack/unpack (Charm++ serialization) routines */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_A;
      p | m_x;
      p | m_r;
      p | m_p;
      p | m_q;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c ConjugateGradients object reference
    friend void operator|( PUP::er& p, ConjugateGradients& c ) { c.pup(p); }
    ///@}

  private:
   CSR m_A;                             //!< Sparse matrix
   std::vector< tk::real > m_x;         //!< Solution/unknown
   std::vector< tk::real > m_r;         //!< Auxiliary vector for CG solve
   std::vector< tk::real > m_p;         //!< Auxiliary vector for CG solve
   std::vector< tk::real > m_q;         //!< Auxiliary vector for CG solve
};

} // tk::

#endif // ConjugateGradients_h
