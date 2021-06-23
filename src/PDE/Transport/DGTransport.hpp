// *****************************************************************************
/*!
  \file      src/PDE/Transport/DGTransport.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Scalar transport using disccontinous Galerkin discretization
  \details   This file implements the physics operators governing transported
     scalars using disccontinuous Galerkin discretization.
*/
// *****************************************************************************
#ifndef DGTransport_h
#define DGTransport_h

#include <vector>
#include <array>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <map>

#include "Macro.hpp"
#include "Exception.hpp"
#include "Vector.hpp"
#include "Inciter/Options/BC.hpp"
#include "UnsMesh.hpp"
#include "Integrate/Basis.hpp"
#include "Integrate/Quadrature.hpp"
#include "Integrate/Initialize.hpp"
#include "Integrate/Mass.hpp"
#include "Integrate/Surface.hpp"
#include "Integrate/Boundary.hpp"
#include "Integrate/Volume.hpp"
#include "Riemann/Upwind.hpp"
#include "Reconstruction.hpp"
#include "Limiter.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief Transport equation used polymorphically with tk::DGPDE
//! \details The template argument(s) specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/Transport/Physics.h
//!   - Problem - problem configuration, see PDE/Transport/Problem.h
//! \note The default physics is DGAdvection, set in
//!    inciter::deck::check_transport()
template< class Physics, class Problem >
class Transport {

  private:
    using eq = tag::transport;

  public:
    //! Constructor
    //! \param[in] c Equation system index (among multiple systems configured)
    explicit Transport( ncomp_t c ) :
      m_physics( Physics() ),
      m_problem( Problem() ),
      m_system( c ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< eq >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< eq >(c) )
    {
      // associate boundary condition configurations with state functions, the
      // order in which the state functions listed matters, see ctr::bc::Keys
      brigand::for_each< ctr::bc::Keys >( ConfigBC< eq >( m_system, m_bc,
        { dirichlet
        , invalidBC  // Symmetry BC not implemented
        , inlet
        , outlet
        , invalidBC  // Characteristic BC not implemented
        , extrapolate } ) );
      m_problem.errchk( m_system, m_ncomp );
    }

    //! Find the number of primitive quantities required for this PDE system
    //! \return The number of primitive quantities required to be stored for
    //!   this PDE system
    std::size_t nprim() const
    {
      // transport does not need/store any primitive quantities currently
      return 0;
    }

    //! Find the number of materials set up for this PDE system
    //! \return The number of materials set up for this PDE system
    std::size_t nmat() const
    {
      return m_ncomp;
    }

    //! Assign number of DOFs per equation in the PDE system
    //! \param[in,out] numEqDof Array storing number of Dofs for each PDE
    //!   equation
    void numEquationDofs(std::vector< std::size_t >& numEqDof) const
    {
      // all equation-dofs initialized to ndofs
      for (std::size_t i=0; i<m_ncomp; ++i) {
        numEqDof.push_back(g_inputdeck.get< tag::discr, tag::ndof >());
      }
    }

    //! Determine elements that lie inside the user-defined IC box
    void IcBoxElems( const tk::Fields&,
      std::size_t,
      std::vector< std::unordered_set< std::size_t > >& ) const
    {}

    //! Initalize the transport equations for DG
    //! \param[in] L Element mass matrix
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    //! \param[in] nielem Number of internal elements
    void
    initialize(
      const tk::Fields& L,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const std::vector< std::unordered_set< std::size_t > >& /*inbox*/,
      tk::Fields& unk,
      tk::real t,
      const std::size_t nielem ) const
    {
      tk::initialize( m_system, m_ncomp, m_offset, L, inpoel, coord,
                      Problem::initialize, unk, t, nielem );
    }

    //! Compute the left hand side mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      tk::mass( m_ncomp, m_offset, ndof, geoElem, l );
    }

    //! Update the primitives for this PDE system
    //! \details This function computes and stores the dofs for primitive
    //!   quantities, which are currently unused for transport.
    void updatePrimitives( const tk::Fields&,
                           const tk::Fields&,
                           const tk::Fields&,
                           tk::Fields&,
                           std::size_t ) const {}

    //! Clean up the state of trace materials for this PDE system
    //! \details This function cleans up the state of materials present in trace
    //!   quantities in each cell. This is currently unused for transport.
    void cleanTraceMaterial( const tk::Fields&,
                             tk::Fields&,
                             tk::Fields&,
                             std::size_t ) const {}

    //! Reconstruct second-order solution from first-order
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] esup Elements-surrounding-nodes connectivity
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] U Solution vector at recent time step
    //! \param[in,out] P Primitive vector at recent time step
    void reconstruct( tk::real t,
                      const tk::Fields& geoFace,
                      const tk::Fields& geoElem,
                      const inciter::FaceData& fd,
                      const std::map< std::size_t, std::vector< std::size_t > >&
                        esup,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      tk::Fields& U,
                      tk::Fields& P,
                      tk::Fields& ) const
    {
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

      // do reconstruction only if P0P1
      if (rdof == 4 && g_inputdeck.get< tag::discr, tag::ndof >() == 1) {
        const auto nelem = fd.Esuel().size()/4;
        const auto intsharp = g_inputdeck.get< tag::param, tag::transport,
          tag::intsharp >()[m_system];

        Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
                "vector must equal "+ std::to_string(rdof*m_ncomp) );
        Assert( fd.Inpofa().size()/3 == fd.Esuf().size()/2,
                "Mismatch in inpofa size" );

        // allocate and initialize matrix and vector for reconstruction
        std::vector< std::array< std::array< tk::real, 3 >, 3 > >
          lhs_ls( nelem, {{ {{0.0, 0.0, 0.0}},
                            {{0.0, 0.0, 0.0}},
                            {{0.0, 0.0, 0.0}} }} );
        // specify how many variables need to be reconstructed
        std::array< std::size_t, 2 > varRange {{0, m_ncomp-1}};

        std::vector< std::vector< std::array< tk::real, 3 > > >
          rhs_ls( nelem, std::vector< std::array< tk::real, 3 > >
            ( m_ncomp,
              {{ 0.0, 0.0, 0.0 }} ) );

        // reconstruct x,y,z-derivatives of unknowns
        // 0. get lhs matrix, which is only geometry dependent
        tk::lhsLeastSq_P0P1(fd, geoElem, geoFace, lhs_ls);

        // 1. internal face contributions
        tk::intLeastSq_P0P1( m_offset, rdof, fd, geoElem, U, rhs_ls, varRange );

        // 2. boundary face contributions
        for (const auto& b : m_bc)
          tk::bndLeastSqConservedVar_P0P1( m_system, m_ncomp, m_offset, rdof,
            b.first, fd, geoFace, geoElem, t, b.second, P, U, rhs_ls, varRange );

        // 3. solve 3x3 least-squares system
        tk::solveLeastSq_P0P1( m_offset, rdof, lhs_ls, rhs_ls, U, varRange );

        for (std::size_t e=0; e<nelem; ++e)
        {
          std::vector< std::size_t > matInt(m_ncomp, 0);
          std::vector< tk::real > alAvg(m_ncomp, 0.0);
          for (std::size_t k=0; k<m_ncomp; ++k)
            alAvg[k] = U(e, k*rdof, m_offset);
          auto intInd = interfaceIndicator(m_ncomp, alAvg, matInt);
          if ((intsharp > 0) && intInd)
          {
            // Reconstruct second-order dofs of volume-fractions in Taylor space
            // using nodal-stencils, for a good interface-normal estimate
            tk::recoLeastSqExtStencil( rdof, m_offset, e, esup, inpoel, geoElem,
              U, varRange );
          }
        }

        // 4. transform reconstructed derivatives to Dubiner dofs
        tk::transform_P0P1( m_offset, rdof, nelem, inpoel, coord, U, varRange );
      }
    }

    //! Limit second-order solution
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] esup Elements surrounding points
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] ndofel Vector of local number of degrees of freedome
    //! \param[in] gid Local->global node id map
    //! \param[in] bid Local chare-boundary node ids (value) associated to
    //!   global node ids (key)
    //! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
    //!   variables
    //! \param[in] pNodalExtrm Chare-boundary nodal extrema for primitive
    //!   variables
    //! \param[in,out] U Solution vector at recent time step
    void limit( [[maybe_unused]] tk::real t,
                [[maybe_unused]] const tk::Fields& geoFace,
                const tk::Fields& geoElem,
                const inciter::FaceData& fd,
                const std::map< std::size_t, std::vector< std::size_t > >& esup,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& ndofel,
                const std::vector< std::size_t >& gid,
                const std::unordered_map< std::size_t, std::size_t >& bid,
                const std::vector< std::vector<tk::real> >& uNodalExtrm,
                [[maybe_unused]] const std::vector< std::vector<tk::real> >&
                  pNodalExtrm,
                tk::Fields& U,
                tk::Fields& ) const
    {
      const auto limiter = g_inputdeck.get< tag::discr, tag::limiter >();
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

      if (limiter == ctr::LimiterType::WENOP1)
        WENO_P1( fd.Esuel(), m_offset, U );
      else if (limiter == ctr::LimiterType::SUPERBEEP1)
        Superbee_P1( fd.Esuel(), inpoel, ndofel, m_offset, coord, U );
      else if (limiter == ctr::LimiterType::VERTEXBASEDP1 && rdof == 4)
        VertexBasedTransport_P1( esup, inpoel, ndofel, fd.Esuel().size()/4,
          m_system, m_offset, coord, gid, bid, uNodalExtrm, U );
      else if (limiter == ctr::LimiterType::VERTEXBASEDP1 && rdof == 10)
        VertexBasedTransport_P2( esup, inpoel, ndofel, fd.Esuel().size()/4,
          m_system, m_offset, geoElem, coord, gid, bid, uNodalExtrm, U );
    }

    //! Compute right hand side
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] U Solution vector at recent time step
    //! \param[in] P Primitive vector at recent time step
    //! \param[in] ndofel Vector of local number of degrees of freedom
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const std::vector< std::unordered_set< std::size_t > >&,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              const tk::Fields& P,
              const tk::Fields& VolFracMax,
              const std::vector< std::size_t >& ndofel,
              tk::Fields& R ) const
    {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
      const auto intsharp = g_inputdeck.get< tag::param, tag::transport,
        tag::intsharp >()[m_system];

      Assert( U.nunk() == P.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( P.nprop() == 0, "Number of components in primitive "
              "vector must equal "+ std::to_string(0) );
      Assert( R.nprop() == ndof*m_ncomp, "Number of components in right-hand "
              "side vector must equal "+ std::to_string(ndof*m_ncomp) );
      Assert( fd.Inpofa().size()/3 == fd.Esuf().size()/2,
              "Mismatch in inpofa size" );

      // set rhs to zero
      R.fill(0.0);

      // empty vector for non-conservative terms. This vector is unused for
      // linear transport since, there are no non-conservative terms in the
      // system of PDEs.
      std::vector< std::vector < tk::real > > riemannDeriv;

      std::vector< std::vector< tk::real > > vriem;
      std::vector< std::vector< tk::real > > riemannLoc;

      // compute internal surface flux integrals
      tk::surfInt( m_system, m_ncomp, m_offset, t, ndof, rdof, inpoel, coord,
                   fd, geoFace, geoElem, Upwind::flux,
                   Problem::prescribedVelocity, U, P, VolFracMax, ndofel, R,
                   vriem, riemannLoc, riemannDeriv, intsharp );

      if(ndof > 1)
        // compute volume integrals
        tk::volInt( m_system, m_ncomp, m_offset, t, ndof, rdof,
                    fd.Esuel().size()/4, inpoel, coord, geoElem, flux,
                    Problem::prescribedVelocity, U, P, ndofel, R, intsharp );

      // compute boundary surface flux integrals
      for (const auto& b : m_bc)
        tk::bndSurfInt( m_system, m_ncomp, m_offset, ndof, rdof, b.first, fd,
          geoFace, geoElem, inpoel, coord, t, Upwind::flux,
          Problem::prescribedVelocity, b.second, U, P, VolFracMax, ndofel, R,
          vriem, riemannLoc, riemannDeriv, intsharp );
    }

    //! Compute the nodal extrema for chare-boundary nodes
    //! \param[in] ncomp Number of conservative variables
    //! \param[in] nprim Number of primitive variables
    //! \param[in] ndof_NodalExtrm Degree of freedom for nodal extrema
    //! \param[in] bndel List of elements contributing to chare-boundary nodes
    //! \param[in] inpoel Element-node connectivity for element e
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] gid Local->global node id map
    //! \param[in] bid Local chare-boundary node ids (value) associated to
    //!   global node ids (key)
    //! \param[in] U Vector of conservative variables
    //! \param[in] P Vector of primitive variables
    //! \param[in,out] uNodalExtrm Chare-boundary nodal extrema for conservative
    //!   variables
    //! \param[in,out] pNodalExtrm Chare-boundary nodal extrema for primitive
    //!   variables
    void evalNodalExtrm( const std::size_t ncomp,
                         [[maybe_unused]] const std::size_t nprim,
                         const std::size_t ndof_NodalExtrm,
                         const std::vector< std::size_t >& bndel,
                         const std::vector< std::size_t >& inpoel,
                         const tk::UnsMesh::Coords& coord,
                         const std::vector< std::size_t >& gid,
                         const std::unordered_map< std::size_t, std::size_t >&
                           bid,
                         const tk::Fields& U,
                         [[maybe_unused]] const tk::Fields& P,
                         std::vector< std::vector<tk::real> >& uNodalExtrm,
                         [[maybe_unused]] std::vector< std::vector<tk::real> >&
                           pNodalExtrm ) const
    {
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

      for (auto e : bndel)
      {
        // access node IDs
        const std::vector<std::size_t> N
          { inpoel[e*4+0], inpoel[e*4+1], inpoel[e*4+2], inpoel[e*4+3] };

        // Loop over nodes of element e
        for(std::size_t ip=0; ip<4; ++ip)
        {
          auto i = bid.find( gid[N[ip]] );
          if (i != end(bid))      // If ip is the chare boundary point
          {
            // Find the nodal extrema of conservative variables
            for (std::size_t c=0; c<ncomp; ++c)
            {
              auto max_mark = c * ndof_NodalExtrm;
              auto min_mark = max_mark + ncomp * ndof_NodalExtrm;
              uNodalExtrm[i->second][max_mark] =
                std::max(uNodalExtrm[i->second][max_mark], U(e,c*rdof,0));
              uNodalExtrm[i->second][min_mark] =
                std::min(uNodalExtrm[i->second][min_mark], U(e,c*rdof,0));
            }

            // If DG(P2) is applied, find the nodal extrema of the gradients of
            // conservative/primitive variables in the physical domain
            if(ndof_NodalExtrm > 1)
            {
              // Vector used to store the first order derivatives
              std::vector< tk::real > grad(3, 0.0);

              const auto& cx = coord[0];
              const auto& cy = coord[1];
              const auto& cz = coord[2];

              std::array< std::array< tk::real, 3>, 4 > coordel {{
                {{ cx[ N[0] ], cy[ N[0] ], cz[ N[0] ] }},
                {{ cx[ N[1] ], cy[ N[1] ], cz[ N[1] ] }},
                {{ cx[ N[2] ], cy[ N[2] ], cz[ N[2] ] }},
                {{ cx[ N[3] ], cy[ N[3] ], cz[ N[3] ] }}
              }};

              auto jacInv = tk::inverseJacobian( coordel[0], coordel[1],
                coordel[2], coordel[3] );

              // Compute the derivatives of basis functions
              auto dBdx = tk::eval_dBdx_p1( rdof, jacInv );

              std::vector< tk::real > center{0.25, 0.25, 0.25};
              tk::evaldBdx_p2(center, jacInv, dBdx);

              // Evaluate the first order derivative in physical domain
              for(std::size_t icomp = 0; icomp < ncomp; icomp++)
              {
                auto mark = icomp * rdof;
                for(std::size_t idir = 0; idir < 3; idir++)
                {
                  for(std::size_t idof = 1; idof < rdof; idof++)
                    grad[idir] += U(e, mark+idof, 0) * dBdx[idir][idof];
                }
              }

              // Store the extrema for the gradients
              for (std::size_t c=0; c<ncomp; ++c)
              {
                for (std::size_t idof = 1; idof < ndof_NodalExtrm; idof++)
                {
                  auto max_mark = c * ndof_NodalExtrm + idof;
                  auto min_mark = max_mark + ncomp * ndof_NodalExtrm;
                  uNodalExtrm[i->second][max_mark] =
                    std::max(uNodalExtrm[i->second][max_mark], grad[idof]);
                  uNodalExtrm[i->second][min_mark] =
                    std::min(uNodalExtrm[i->second][min_mark], grad[idof]);
                }
              }
            }
          }
        }
      }
    }

    //! Compute the minimum time step size
//     //! \param[in] U Solution vector at recent time step
//     //! \param[in] coord Mesh node coordinates
//     //! \param[in] inpoel Mesh element connectivity
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& /*coord*/,
                 const std::vector< std::size_t >& /*inpoel*/,
                 const inciter::FaceData& /*fd*/,
                 const tk::Fields& /*geoFace*/,
                 const tk::Fields& /*geoElem*/,
                 const std::vector< std::size_t >& /*ndofel*/,
                 const tk::Fields& /*U*/,
                 const tk::Fields&,
                 const std::size_t /*nielem*/ ) const
    {
      tk::real mindt = std::numeric_limits< tk::real >::max();
      return mindt;
    }

    //! Return analytic field names to be output to file
    //! \return Vector of strings labelling analytic fields output in file
    std::vector< std::string > analyticFieldNames() const {
      std::vector< std::string > n;
      auto depvar = g_inputdeck.get< tag::param, eq, tag::depvar >()[m_system];
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_analytic" );
      return n;
    }

    //! Return surface field output going to file
    std::vector< std::vector< tk::real > >
    surfOutput( const std::map< int, std::vector< std::size_t > >&,
                tk::Fields& ) const
    {
      std::vector< std::vector< tk::real > > s; // punt for now
      return s;
    }

    //! Return time history field names to be output to file
    //! \return Vector of strings labelling time history fields output in file
    std::vector< std::string > histNames() const {
      std::vector< std::string > s; // punt for now
      return s;
    }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] nunk Number of unknowns to extract
    //! \param[in] vol Volumes associated to elements (or nodes)
    //! \param[in] coord Coordinates at which to evaluate the solution
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    //! \details This functions should be written in conjunction with names(),
    //!   which provides the vector of field names
    //! \note U is overwritten
    std::vector< std::vector< tk::real > >
    fieldOutput( tk::real t,
                 tk::real,
                 std::size_t nunk,
                 std::size_t rdof,
                 const std::vector< tk::real >& vol,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 const tk::Fields& U,
                 [[maybe_unused]] const tk::Fields& = tk::Fields() ) const
    {
      Assert( U.nunk() >= nunk, "Size mismatch" );
      std::vector< std::vector< tk::real > > out;

      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c*rdof, m_offset ) );

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // evaluate analytic solution at time t
      auto E = U;
      for (std::size_t i=0; i<nunk; ++i)
      {
        auto s =
          Problem::solution( m_system, m_ncomp, x[i], y[i], z[i], t );
        for (ncomp_t c=0; c<m_ncomp; ++c)
          E( i, c*rdof, m_offset ) = s[c];
      }

      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( E.extract( c*rdof, m_offset ) );

      // will output error for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) {
        auto mark = c*rdof;
        auto u = U.extract( mark, m_offset );
        auto e = E.extract( mark, m_offset );
        for (std::size_t i=0; i<nunk; ++i)
          e[i] = std::pow( e[i] - u[i], 2.0 ) * vol[i];
        out.push_back( e );
      }

      // will output material indicator function
      std::vector< tk::real > matInd(nunk, 0.0);
      for (std::size_t i=0; i<nunk; ++i) {
        for (std::size_t c=0; c<m_ncomp; ++c) {
          matInd[i] += U(i,c*rdof,m_offset) * static_cast< tk::real >(c+1);
        }
      }
      out.push_back(matInd);

      return out;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const {
      std::vector< std::string > n;
      const auto& depvar =
      g_inputdeck.get< tag::param, eq, tag::depvar >().at(m_system);
      // construct the name of the numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) );
      return n;
    }

    //! Return analytic solution (if defined by Problem) at xi, yi, zi, t
    //! \param[in] xi X-coordinate at which to evaluate the analytic solution
    //! \param[in] yi Y-coordinate at which to evaluate the analytic solution
    //! \param[in] zi Z-coordinate at which to evaluate the analytic solution
    //! \param[in] t Physical time at which to evaluate the analytic solution
    //! \return Vector of analytic solution at given spatial location and time
    std::vector< tk::real >
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return Problem::analyticSolution( m_system, m_ncomp, xi, yi, zi, t ); }

    //! Return analytic solution for conserved variables
    //! \param[in] xi X-coordinate at which to evaluate the analytic solution
    //! \param[in] yi Y-coordinate at which to evaluate the analytic solution
    //! \param[in] zi Z-coordinate at which to evaluate the analytic solution
    //! \param[in] t Physical time at which to evaluate the analytic solution
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    solution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return Problem::initialize( m_system, m_ncomp, xi, yi, zi, t ); }

    //! Return time history field output evaluated at time history points
    //! \param[in] h History point data
    std::vector< std::vector< tk::real > >
    histOutput( const std::vector< HistData >& h,
                const std::vector< std::size_t >&,
                const tk::UnsMesh::Coords&,
                const tk::Fields& ) const
    {
      std::vector< std::vector< tk::real > > Up(h.size()); //punt for now
      return Up;
    }

  private:
    const Physics m_physics;            //!< Physics policy
    const Problem m_problem;            //!< Problem policy
    const ncomp_t m_system;             //!< Equation system index
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    const ncomp_t m_offset;             //!< Offset this PDE operates from
    //! BC configuration
    BCStateFn m_bc;

    //! Evaluate physical flux function for this PDE system
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ugp Numerical solution at the Gauss point at which to
    //!   evaluate the flux
    //! \param[in] v Prescribed velocity evaluated at the Gauss point at which
    //!   to evaluate the flux
    //! \return Flux vectors for all components in this PDE system
    //! \note The function signature must follow tk::FluxFn
    static tk::FluxFn::result_type
    flux( ncomp_t,
          ncomp_t ncomp,
          const std::vector< tk::real >& ugp,
          const std::vector< std::array< tk::real, 3 > >& v )
    {
      Assert( ugp.size() == ncomp, "Size mismatch" );
      Assert( v.size() == ncomp, "Size mismatch" );

      std::vector< std::array< tk::real, 3 > > fl( ugp.size() );

      for (ncomp_t c=0; c<ncomp; ++c)
        fl[c] = {{ v[c][0] * ugp[c], v[c][1] * ugp[c], v[c][2] * ugp[c] }};

      return fl;
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at extrapolation boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    extrapolate( ncomp_t, ncomp_t, const std::vector< tk::real >& ul,
                 tk::real, tk::real, tk::real, tk::real,
                 const std::array< tk::real, 3 >& )
    {
      return {{ ul, ul }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at extrapolation boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    inlet( ncomp_t, ncomp_t, const std::vector< tk::real >& ul,
           tk::real, tk::real, tk::real, tk::real,
           const std::array< tk::real, 3 >& )
    {
      auto ur = ul;
      std::fill( begin(ur), end(ur), 0.0 );
      return {{ ul, std::move(ur) }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at outlet boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    outlet( ncomp_t, ncomp_t, const std::vector< tk::real >& ul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& )
    {
      return {{ ul, ul }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at Dirichlet boundaries
    //! \param[in] system Equation system index
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ul Left (domain-internal) state
    //! \param[in] x X-coordinate at which to compute the states
    //! \param[in] y Y-coordinate at which to compute the states
    //! \param[in] z Z-coordinate at which to compute the states
    //! \param[in] t Physical time
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    dirichlet( ncomp_t system, ncomp_t ncomp, const std::vector< tk::real >& ul,
               tk::real x, tk::real y, tk::real z, tk::real t,
               const std::array< tk::real, 3 >& )
    {
      return {{ ul, Problem::initialize( system, ncomp, x, y, z, t ) }};
    }
};

} // dg::
} // inciter::

#endif // DGTransport_h
