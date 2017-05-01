// *****************************************************************************
/*!
  \file      src/Inciter/AuxSolver.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Auxiliary solvers used in tk::LinSysMerger by inciter.
  \details   This file defines policy classes for tk::LinSysMerger's auxiliary
    solver. They can be used as the AuxSolver template argument as long as the
    correct member functions are defined with their correct return types and
    arguments.
*/
// *****************************************************************************
#ifndef AuxSolver_h
#define AuxSolver_h

#include <functional>

#include "Macro.h"
#include "Fields.h"
#include "ContainerUtil.h"

namespace inciter {

//! AuxSolverLumpMassDiff: Auxiliary linear solver: lumped mass + diffusion
//! \details This class defines member functions that implement
//!   assembling/merging the left and right hand sides of a linear system
//!   resulting from a high order discretization obtained by lumping the
//!   left hand side consistent mass matrix and adding mass diffusion to the
//!   right hand side vector. Used in flux-corrected transport for transport
//!   equations as the low order system. This is only the functionality that is
//!   needed for the parallel merging of the linear system. Used by and
//!   abstracted away from tk::LinSysMerger as its auxiliary solver, injected in
//!   via the template argument AuxSolver.
//! \see src/LinSys/LinSysMerger.h
//! \note All member functions are static as they are called without an object.
class AuxSolverLumpMassDiff {

  public:
    //! Chares contribute to the lhs of the low order linear system
    //! \param[in] hostproxy The Charm++ proxy of the class we are used as a
    //!   template argument of
    //! \param[in] host Pointer to the object instance we interoperate with
    //! \param[in] lower Lower index of the global linear system rows on my PE
    //! \param[in] upper Upper index of the global linear system rows on my PE
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] mass Portion of the lhs vector contributed
    //! \param[in,out] auxlhsimport Import map associating a list of global row
    //!   ids to a worker chare id during the communication of the low order lhs
    //!   vector
    //! \param[in] auxlhs Part of the low order (auxiliary) system lhs vector
    //!   owned by my PE. These are the vector of values (for each scalar
    //!   equation solved) associated to global mesh point row ids. This vector
    //!   collects the nonzero values of the low order (auxiliary) linear system
    //!   lhs "matrix" solution, in the form of a vector, since the matrix of
    //!   this particular low order system is lumped, i.e., stored as a single
    //!   vector of diagonal elements of the lhs matrix.
    //! \note Static member function as it is called without an object.
    //! \author J. Bakosi
    template< class HostProxy, class Host > static
    void chareauxlhs( const HostProxy& hostproxy,
                      Host* const host,
                      std::size_t lower,
                      std::size_t upper,
                      int fromch,
                      const std::vector< std::size_t >& gid,
                      const tk::Fields& mass,
                      std::map< int, std::vector< std::size_t > >& auxlhsimport,
                      std::map< std::size_t, std::vector< tk::real > >& auxlhs )
    {
      using tk::operator+=;
      Assert( gid.size() == mass.nunk(),
              "Size of mass diffusion lhs and row ID vectors must equal" );
      // Store+add vector of nonzero values owned and pack those to be exported
      std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= lower && gid[i] < upper) {  // if own
          auxlhsimport[ fromch ].push_back( gid[i] );
          auxlhs[ gid[i] ] += mass[i];
        } else
          exp[ host->pe(gid[i]) ][ gid[i] ] = mass[i];
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        hostproxy[ tope ].addauxlhs( fromch, p.second );
      }
    }

    //! \brief Receive and add lhs vector to the auxiliary system from fellow
    //!   group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] mass Portion of the lhs vector contributed to the auxiliary
    //!   (low order) linear system, containing global row indices and values
    //! \param[in,out] auxlhsimport Import map associating a list of global row
    //!   ids to a worker chare id during the communication of the low order lhs
    //!   vector
    //! \param[in,out] auxlhs Part of the low order (auxiliary) system lhs
    //!   vector owned by my PE. These are the vector of values (for each scalar
    //!   equation solved) associated to global mesh point row ids. This vector
    //!   collects the nonzero values of the low order (auxiliary) linear system
    //!   lhs "matrix" solution, in the form of a vector, since the matrix of
    //!   this particular low order system is lumped, i.e., stored as a single
    //!   vector of diagonal elements of the lhs matrix.
    //! \note Static member function as it is called without an object.
    //! \author J. Bakosi
    static void
    addauxlhs( int fromch,
               const std::map< std::size_t, std::vector< tk::real > >& mass,
               std::map< int, std::vector< std::size_t > >& auxlhsimport,
               std::map< std::size_t, std::vector< tk::real > >& auxlhs )
    {
      using tk::operator+=;
      for (const auto& r : mass) {
        auxlhsimport[ fromch ].push_back( r.first );
        auxlhs[ r.first ] += r.second;
      }
    }

    //! Chares contribute their mass diffusion rhs to the low order system
    //! \param[in] hostproxy The Charm++ proxy of the class we are used as a
    //!   template argument of
    //! \param[in] host Pointer to the object instance we interoperate with
    //! \param[in] lower Lower index of the global linear system rows on my PE
    //! \param[in] upper Upper index of the global linear system rows on my PE
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] gid Global row indices of the vector contributed
    //! \param[in] diff Portion of the rhs vector contributed
    //! \param[in,out] auxrhsimport Import map associating a list of global row
    //!   ids to a worker chare id during the communication of the mass
    //!   diffusion rhs vector vector.
    //! \param[in,out] auxrhs Part of the mass diffusion right-hand side vector
    //!   owned by my PE. Vector of values (for each scalar equation solved)
    //!   associated to global mesh point row ids. This vector collects the mass
    //!   diffusion terms to be added to the right-hand side for the low order
    //!   solution (for flux-corrected transport).
    //! \note Static member function as it is called without an object.
    //! \author J. Bakosi
    template< class HostProxy, class Host > static
    void chareauxrhs( const HostProxy& hostproxy,
                      Host* const host,
                      std::size_t lower,
                      std::size_t upper,
                      int fromch,
                      const std::vector< std::size_t >& gid,
                      const tk::Fields& diff,
                      std::map< int, std::vector< std::size_t > >& auxrhsimport,
                      std::map< std::size_t, std::vector< tk::real > >& auxrhs )
    {
      using tk::operator+=;
      Assert( gid.size() == diff.nunk(),
              "Size of mass diffusion rhs and row ID vectors must equal" );
      // Store+add vector of nonzero values owned and pack those to be exported
      std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;
      for (std::size_t i=0; i<gid.size(); ++i)
        if (gid[i] >= lower && gid[i] < upper) {  // if own
          auxrhsimport[ fromch ].push_back( gid[i] );
          auxrhs[ gid[i] ] += diff[i];
        } else
          exp[ host->pe(gid[i]) ][ gid[i] ] = diff[i];
      // Export non-owned vector values to fellow branches that own them
      for (const auto& p : exp) {
        auto tope = static_cast< int >( p.first );
        hostproxy[ tope ].addauxrhs( fromch, p.second );
      }
    }

    //! \brief Receive and add rhs vector to the auxiliary system from fellow
    //!   group branches
    //! \param[in] fromch Charm chare array index contribution coming from
    //! \param[in] diff Portion of the right-hand side vector contributed,
    //!   containing global row indices and values
    //! \param[in,out] auxrhsimport Import map associating a list of global row
    //!   ids to a worker chare id during the communication of the mass
    //!   diffusion rhs vector vector.
    //! \param[in,out] auxrhs Part of the mass diffusion right-hand side vector
    //!   owned by my PE. Vector of values (for each scalar equation solved)
    //!   associated to global mesh point row ids. This vector collects the mass
    //!   diffusion terms to be added to the right-hand side for the low order
    //!   solution (for flux-corrected transport).
    //! \note Static member function as it is called without an object.
    //! \author J. Bakosi
    static void
    addauxrhs( int fromch,
               const std::map< std::size_t, std::vector< tk::real > >& diff,
               std::map< int, std::vector< std::size_t > >& auxrhsimport,
               std::map< std::size_t, std::vector< tk::real > >& auxrhs )
    {
      using tk::operator+=;
      for (const auto& r : diff) {
        auxrhsimport[ fromch ].push_back( r.first );
        auxrhs[ r.first ] += r.second;
      }
    }

    //! Solve low order system
    //! \param[in] host Pointer to the object instance we interoperate with
    //! \param[in] ncomp Number of scalar components per unknown
    //! \param[in] rhs Part of right-hand side vector owned by my PE. Vector of
    //!   values (for each scalar equation solved) associated to global mesh
    //!   point row ids. In the flux-corrected transport algorithm, this is the
    //!   rhs vector of the high order scheme, used to add mass diffusion to
    //!   create the low order rhs.
    //! \param[in] lump Part of auxiliary system left-hand side "matrix"
    //!   (vector) owned by my PE. Vector of values (for each scalar equation
    //!   solved) associated to global mesh point row ids. This vector collects
    //!   the nonzero values of the auxiliary system lhs "matrix" solution. In
    //!   flux-corrected transport this is the lhs forthe low order system
    //!   only storing the nonzero diagonal entries of the lhs matrix.
    //! \param[in,out] diff Part of the mass diffusion right-hand side vector
    //!   owned by my PE. Vector of values (for each scalar equation solved)
    //!   associated to global mesh point row ids. This vector collects the mass
    //!   diffusion terms to be added to the right-hand side for the low order
    //!   solution (for flux-corrected transport).
    //! \note Static member function as it is called without an object.
    //! \author J. Bakosi
    template< class Host > static
    void solve( Host* const host,
                std::size_t ncomp,
                const std::map< std::size_t, std::vector< tk::real > >& rhs,
                const std::map< std::size_t, std::vector< tk::real > >& lump,
                std::map< std::size_t, std::vector< tk::real > >& diff )
    {
      IGNORE(host);
      Assert( host->rhscomplete(),
              "Values of distributed right-hand-side vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
              "order system" );
      Assert( host->auxrhscomplete(),
              "Values of distributed mass diffusion rhs vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
              "order system" );
      Assert( host->auxlhscomplete(),
              "Values of distributed lumped mass lhs vector on PE " +
              std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
              "order system" );
      Assert( tk::keyEqual( rhs, diff ), "Row IDs of rhs and mass "
              "diffusion rhs vector unequal on PE " + std::to_string( CkMyPe() )
              + ": cannot solve low order system" );
      Assert( tk::keyEqual( rhs, lump ), "Row IDs of rhs and lumped mass "
              "lhs vector unequal on PE " + std::to_string( CkMyPe() ) + ": "
              "cannot solve low order system" );
      auto ir = rhs.cbegin();
      auto id = diff.begin();
      auto im = lump.cbegin();
      while (ir != rhs.cend()) {
        const auto& r = ir->second;
        const auto& m = im->second;
        auto& d = id->second;
        Assert( r.size()==ncomp && m.size()==ncomp && d.size()==ncomp,
                "Wrong number of components in solving the low order system" );
        for (std::size_t i=0; i<ncomp; ++i) d[i] = (r[i]+d[i])/m[i];
        ++ir; ++id; ++im;
      }
    }

    //! \brief Group low order solution vector by workers and send each the parts
    //!   back to workers that own them
    //! \param[in] worker Worker proxy to use send update to
    //! \param[in] solimport Import map associating a list of global row ids to a
    //!   worker chare id during the communication of the solution/unknown vector.
    //!   This is used here to find out which solution values associated to which
    //!   global nodes IDs to send the solution updates back to the workers.
    //! \param[in] sol Part of solution vector owned by my PE. Vector of values
    //!   (for each scalar equation solved) associated to global mesh point row
    //!   IDs.
    //! \note Static member function as it is called without an object.
    //! \author J. Bakosi
    template< class WorkerProxy > static
    void update( const WorkerProxy& worker,
                 const std::map< int, std::vector< std::size_t > >& solimport,
                 const std::map< std::size_t, std::vector< tk::real > >& sol )
    {
      for (const auto& w : solimport) {
        std::vector< std::size_t > gid;
        std::vector< tk::real > solution;
        for (auto r : w.second) {
          const auto it = sol.find( r );
          if (it != end(sol)) {
            gid.push_back( it->first );
            solution.insert(end(solution), begin(it->second), end(it->second));
          } else
            Throw( "Can't find global row id " + std::to_string(r) +
                   " to export in low order solution vector" );
        }
        worker[ w.first ].updateLowSol( gid, solution );
      }
    }
};

//! AuxSolverNull: No auxiliary solver (unused for now)
class AuxSolverNull {

};

} // inciter::

#endif // AuxSolver_h
