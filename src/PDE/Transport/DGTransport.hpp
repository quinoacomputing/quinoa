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
#include "PrefIndicator.hpp"
#include "EoS/EOS.hpp"
#include "FunctionPrototypes.hpp"
#include "ConfigureTransport.hpp"

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
    explicit Transport() :
      m_physics( Physics() ),
      m_problem( Problem() ),
      m_ncomp( g_inputdeck.get< tag::ncomp >() )
    {
      // associate boundary condition configurations with state functions, the
      // order in which the state functions listed matters, see ctr::bc::Keys
      brigand::for_each< ctr::bclist::Keys >( ConfigBC( m_bc,
        // BC State functions
        { dirichlet
        , invalidBC  // Symmetry BC not implemented
        , inlet
        , outlet
        , invalidBC  // Characteristic BC not implemented
        , extrapolate
        , invalidBC },      // No slip wall BC not implemented
        // BC Gradient functions
        { noOpGrad
        , noOpGrad
        , noOpGrad
        , noOpGrad
        , noOpGrad
        , noOpGrad
        , noOpGrad }
        ) );
      m_problem.errchk( m_ncomp );
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
        numEqDof.push_back(g_inputdeck.get< tag::ndof >());
      }
    }

    //! Find how 'stiff equations', which we currently
    //! have none for Transport
    //! \return number of stiff equations
    std::size_t nstiffeq() const
    { return 0; }

    //! Find how 'nonstiff equations', which we currently
    //! don't use for Transport
    //! \return number of non-stiff equations
    std::size_t nnonstiffeq() const
    { return 0; }

    //! Locate the stiff equations. Unused for transport.
    //! \param[out] stiffEqIdx list
    void setStiffEqIdx( std::vector< std::size_t >& stiffEqIdx ) const
    {
      stiffEqIdx.resize(0);
    }

    //! Locate the nonstiff equations. Unused for transport.
    //! \param[out] nonStiffEqIdx list
    void setNonStiffEqIdx( std::vector< std::size_t >& nonStiffEqIdx ) const
    {
      nonStiffEqIdx.resize(0);
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
      const std::unordered_map< std::size_t, std::set< std::size_t > >&,
      tk::Fields& unk,
      tk::real t,
      const std::size_t nielem ) const
    {
      tk::initialize( m_ncomp, m_mat_blk, L, inpoel, coord,
                      Problem::initialize, unk, t, nielem );
    }

    //! Compute density constraint for a given material
    // //! \param[in] nelem Number of elements
    // //! \param[in] unk Array of unknowns
    //! \param[out] densityConstr Density Constraint: rho/(rho0*det(g))
    void computeDensityConstr( std::size_t /*nelem*/,
                               tk::Fields& /*unk*/,
                               std::vector< tk::real >& densityConstr) const
    {
      densityConstr.resize(0);
    }

    //! Compute the left hand side mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const {
      const auto ndof = g_inputdeck.get< tag::ndof >();
      tk::mass( m_ncomp, ndof, geoElem, l );
    }

    //! Update the interface cells to first order dofs
    //! \details This function resets the high-order terms in interface cells,
    //!   and is currently not used in transport.
    void updateInterfaceCells( tk::Fields&,
      std::size_t,
      std::vector< std::size_t >&,
      std::vector< std::size_t >& ) const {}

    //! Update the primitives for this PDE system
    //! \details This function computes and stores the dofs for primitive
    //!   quantities, which are currently unused for transport.
    void updatePrimitives( const tk::Fields&,
                           const tk::Fields&,
                           const tk::Fields&,
                           tk::Fields&,
                           std::size_t,
                           std::vector< std::size_t >& ) const {}

    //! Clean up the state of trace materials for this PDE system
    //! \details This function cleans up the state of materials present in trace
    //!   quantities in each cell. This is currently unused for transport.
    void cleanTraceMaterial( tk::real,
                             const tk::Fields&,
                             tk::Fields&,
                             tk::Fields&,
                             std::size_t ) const {}

    //! Reconstruct second-order solution from first-order
//    //! \param[in] t Physical time
//    //! \param[in] geoFace Face geometry array
//    //! \param[in] geoElem Element geometry array
//    //! \param[in] fd Face connectivity and boundary conditions object
//    //! \param[in] esup Elements-surrounding-nodes connectivity
//    //! \param[in] inpoel Element-node connectivity
//    //! \param[in] coord Array of nodal coordinates
//    //! \param[in,out] U Solution vector at recent time step
//    //! \param[in,out] P Primitive vector at recent time step
    void reconstruct( tk::real,
                      const tk::Fields&,
                      const tk::Fields&,
                      const inciter::FaceData&,
                      const std::map< std::size_t, std::vector< std::size_t > >&,
                      const std::vector< std::size_t >&,
                      const tk::UnsMesh::Coords&,
                      tk::Fields&,
                      tk::Fields&,
                      const bool,
                      const std::vector< std::size_t >& ) const
    {
      // do reconstruction only if P0P1
      if (g_inputdeck.get< tag::rdof >() == 4 &&
        g_inputdeck.get< tag::ndof >() == 1)
        Throw("P0P1 not supported for Transport.");
    }

    //! Limit second-order solution
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] esup Elements surrounding points
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] ndofel Vector of local number of degrees of freedome
//    //! \param[in] gid Local->global node id map
//    //! \param[in] bid Local chare-boundary node ids (value) associated to
//    //!   global node ids (key)
//    //! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
//    //!   variables
    //! \param[in,out] U Solution vector at recent time step
    void limit( [[maybe_unused]] tk::real t,
                [[maybe_unused]] const bool pref,
                [[maybe_unused]] const tk::Fields& geoFace,
                const tk::Fields&,
                const inciter::FaceData& fd,
                const std::map< std::size_t, std::vector< std::size_t > >& esup,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& ndofel,
                const std::vector< std::size_t >&,
                const std::unordered_map< std::size_t, std::size_t >&,
                const std::vector< std::vector<tk::real> >&,
                const std::vector< std::vector<tk::real> >&,
                const std::vector< std::vector<tk::real> >&,
                tk::Fields& U,
                tk::Fields&,
                std::vector< std::size_t >& ) const
    {
      const auto limiter = g_inputdeck.get< tag::limiter >();

      if (limiter == ctr::LimiterType::WENOP1)
        WENO_P1( fd.Esuel(), U );
      else if (limiter == ctr::LimiterType::SUPERBEEP1)
        Superbee_P1( fd.Esuel(), inpoel, ndofel, coord, U );
      else if (limiter == ctr::LimiterType::VERTEXBASEDP1)
        VertexBasedTransport_P1( esup, inpoel, ndofel, fd.Esuel().size()/4,
          coord, U );
    }

    //! Update the conservative variable solution for this PDE system
    //! \details This function computes the updated dofs for conservative
    //!   quantities based on the limited solution and is currently not used in
    //!   transport.
    void CPL( const tk::Fields&,
              const tk::Fields&,
              const std::vector< std::size_t >&,
              const tk::UnsMesh::Coords&,
              tk::Fields&,
              std::size_t ) const {}

    //! Return cell-average deformation gradient tensor (no-op for transport)
    //! \details This function is a no-op in transport.
    std::array< std::vector< tk::real >, 9 > cellAvgDeformGrad(
      const tk::Fields&,
      std::size_t ) const
    {
      return {};
    }

    //! Reset the high order solution for p-adaptive scheme
    //! \details This function reset the high order coefficient for p-adaptive
    //!   solution polynomials and is currently not used in transport.
    void resetAdapSol( const inciter::FaceData&,
                       tk::Fields&,
                       tk::Fields&,
                       const std::vector< std::size_t >& ) const {}

    //! Compute right hand side
    //! \param[in] t Physical time
    //! \param[in] pref Indicator for p-adaptive algorithm
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] U Solution vector at recent time step
    //! \param[in] P Primitive vector at recent time step
    //! \param[in] ndofel Vector of local number of degrees of freedom
    //! \param[in] dt Delta time
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              const bool pref,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const std::vector< std::unordered_set< std::size_t > >&,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              const tk::Fields& P,
              const std::vector< std::size_t >& ndofel,
              const tk::real dt,
              tk::Fields& R ) const
    {
      const auto ndof = g_inputdeck.get< tag::ndof >();
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto intsharp = g_inputdeck.get< tag::transport,
        tag::intsharp >();

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
      std::vector< std::size_t > solidx(1, 0);
      tk::surfInt( pref, m_ncomp, m_mat_blk, t, ndof, rdof,
                   inpoel, solidx, coord, fd, geoFace, geoElem, Upwind::flux,
                   Problem::prescribedVelocity, U, P, ndofel, dt, R, vriem,
                   riemannLoc, riemannDeriv, intsharp );

      if(ndof > 1)
        // compute volume integrals
        tk::volInt( m_ncomp, t, m_mat_blk, ndof, rdof,
                    fd.Esuel().size()/4, inpoel, coord, geoElem, flux,
                    Problem::prescribedVelocity, U, P, ndofel, R, intsharp );

      // compute boundary surface flux integrals
      for (const auto& b : m_bc)
        tk::bndSurfInt( pref, m_ncomp, m_mat_blk, ndof, rdof,
          std::get<0>(b), fd, geoFace, geoElem, inpoel, coord, t, Upwind::flux,
          Problem::prescribedVelocity, std::get<1>(b), U, P, ndofel, R, vriem,
          riemannLoc, riemannDeriv, intsharp );
    }

    //! Evaluate the adaptive indicator and mark the ndof for each element
    //! \param[in] nunk Number of unknowns
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] unk Array of unknowns
    //! \param[in] prim Array of primitive quantities
    //! \param[in] indicator p-refinement indicator type
    //! \param[in] ndof Number of degrees of freedom in the solution
    //! \param[in] ndofmax Max number of degrees of freedom for p-refinement
    //! \param[in] tolref Tolerance for p-refinement
    //! \param[in,out] ndofel Vector of local number of degrees of freedome
    void eval_ndof( std::size_t nunk,
                    [[maybe_unused]] const tk::UnsMesh::Coords& coord,
                    [[maybe_unused]] const std::vector< std::size_t >& inpoel,
                    const inciter::FaceData& fd,
                    const tk::Fields& unk,
                    [[maybe_unused]] const tk::Fields& prim,
                    inciter::ctr::PrefIndicatorType indicator,
                    std::size_t ndof,
                    std::size_t ndofmax,
                    tk::real tolref,
                    std::vector< std::size_t >& ndofel ) const
    {
      const auto& esuel = fd.Esuel();

      if(indicator == inciter::ctr::PrefIndicatorType::SPECTRAL_DECAY)
        spectral_decay( 1, nunk, esuel, unk, ndof, ndofmax, tolref, ndofel );
      else
        Throw( "No such adaptive indicator type" );
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

    //! Compute stiff terms for a single element, not implemented here
    // //! \param[in] e Element number
    // //! \param[in] geoElem Element geometry array
    // //! \param[in] inpoel Element-node connectivity
    // //! \param[in] coord Array of nodal coordinates
    // //! \param[in] U Solution vector at recent time step
    // //! \param[in] P Primitive vector at recent time step
    // //! \param[in] ndofel Vector of local number of degrees of freedom
    // //! \param[in,out] R Right-hand side vector computed
    void stiff_rhs( std::size_t /*e*/,
                    const tk::Fields& /*geoElem*/,
                    const std::vector< std::size_t >& /*inpoel*/,
                    const tk::UnsMesh::Coords& /*coord*/,
                    const tk::Fields& /*U*/,
                    const tk::Fields& /*P*/,
                    const std::vector< std::size_t >& /*ndofel*/,
                    tk::Fields& /*R*/ ) const
    {}

    //! Return a map that associates user-specified strings to functions
    //! \return Map that associates user-specified strings to functions that
    //!  compute relevant quantities to be output to file
    std::map< std::string, tk::GetVarFn > OutVarFn() const {
      std::map< std::string, tk::GetVarFn > OutFnMap;
      OutFnMap["material_indicator"] = transport::matIndicatorOutVar;

      return OutFnMap;
    }

    //! Return analytic field names to be output to file
    //! \return Vector of strings labelling analytic fields output in file
    std::vector< std::string > analyticFieldNames() const {
      std::vector< std::string > n;
      auto depvar = g_inputdeck.get< tag::depvar >()[0];
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

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const {
      std::vector< std::string > n;
      const auto& depvar =
      g_inputdeck.get< tag::depvar >().at(0);
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
    { return Problem::analyticSolution( m_ncomp, m_mat_blk, xi, yi,
                                        zi, t ); }

    //! Return analytic solution for conserved variables
    //! \param[in] xi X-coordinate at which to evaluate the analytic solution
    //! \param[in] yi Y-coordinate at which to evaluate the analytic solution
    //! \param[in] zi Z-coordinate at which to evaluate the analytic solution
    //! \param[in] t Physical time at which to evaluate the analytic solution
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    solution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return Problem::initialize( m_ncomp, m_mat_blk, xi, yi, zi, t ); }

    //! Return time history field output evaluated at time history points
    //! \param[in] h History point data
    std::vector< std::vector< tk::real > >
    histOutput( const std::vector< HistData >& h,
                const std::vector< std::size_t >&,
                const tk::UnsMesh::Coords&,
                const tk::Fields&,
                const tk::Fields& ) const
    {
      std::vector< std::vector< tk::real > > Up(h.size()); //punt for now
      return Up;
    }

    //! Return cell-averaged total component mass per unit volume for an element
    //! \param[in] e Element id for which total energy is required
    //! \param[in] unk Vector of conserved quantities
    //! \return Cell-averaged total component mass per unit volume for given
    //!   element. Since transport does not have an associated total energy,
    //!   return total mass.
    tk::real sp_totalenergy(std::size_t e, const tk::Fields& unk) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();

      tk::real sp_m(0.0);
      for (std::size_t c=0; c<m_ncomp; ++c) {
        sp_m += unk(e,c*rdof);
      }
      return sp_m;
    }

  private:
    const Physics m_physics;            //!< Physics policy
    const Problem m_problem;            //!< Problem policy
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    //! BC configuration
    BCStateFn m_bc;
    //! \brief EOS material block - This PDE does not require an EOS block,
    //! thus this variable has not been intialized.
    std::vector< EOS > m_mat_blk;

    //! Evaluate physical flux function for this PDE system
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ugp Numerical solution at the Gauss point at which to
    //!   evaluate the flux
    //! \param[in] v Prescribed velocity evaluated at the Gauss point at which
    //!   to evaluate the flux
    //! \return Flux vectors for all components in this PDE system
    //! \note The function signature must follow tk::FluxFn
    static tk::FluxFn::result_type
    flux( ncomp_t ncomp,
          const std::vector< EOS >&,
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
    extrapolate( ncomp_t, const std::vector< EOS >&,
                 const std::vector< tk::real >& ul, tk::real, tk::real,
                 tk::real, tk::real, const std::array< tk::real, 3 >& )
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
    inlet( ncomp_t, const std::vector< EOS >&,
           const std::vector< tk::real >& ul, tk::real, tk::real, tk::real,
           tk::real, const std::array< tk::real, 3 >& )
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
    outlet( ncomp_t, const std::vector< EOS >&,
            const std::vector< tk::real >& ul, tk::real, tk::real, tk::real,
            tk::real, const std::array< tk::real, 3 >& )
    {
      return {{ ul, ul }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at Dirichlet boundaries
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
    dirichlet( ncomp_t ncomp, 
               const std::vector< EOS >& mat_blk,
               const std::vector< tk::real >& ul, tk::real x, tk::real y,
               tk::real z, tk::real t, const std::array< tk::real, 3 >& )
    {
      return {{ ul, Problem::initialize( ncomp, mat_blk, x, y, z, t ) }};
    }

  //----------------------------------------------------------------------------
  // Boundary Gradient functions
  //----------------------------------------------------------------------------

  //! \brief Boundary gradient function copying the left gradient to the right
  //!   gradient at a face
  //! \param[in] dul Left (domain-internal) state
  //! \return Left and right states for all scalar components in this PDE
  //!   system
  //! \note The function signature must follow tk::StateFn.
  static tk::StateFn::result_type
  noOpGrad( ncomp_t,
            const std::vector< EOS >&,
            const std::vector< tk::real >& dul,
            tk::real, tk::real, tk::real, tk::real,
            const std::array< tk::real, 3 >& )
  {
    return {{ dul, dul }};
  }
};

} // dg::
} // inciter::

#endif // DGTransport_h
