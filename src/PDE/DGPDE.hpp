// *****************************************************************************
/*!
  \file      src/PDE/DGPDE.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Partial differential equation base for discontinuous Galerkin PDEs
  \details   This file defines a generic partial differential equation (PDE)
    class for PDEs that use discontinuous Galerkin spatial discretization.
    The class uses runtime polymorphism without client-side inheritance:
    inheritance is confined to the internals of the class, invisible to
    client-code. The class exclusively deals with ownership enabling client-side
    value semantics. Credit goes to Sean Parent at Adobe:
    https://github.com/sean-parent/sean-parent.github.com/wiki/
    Papers-and-Presentations.
*/
// *****************************************************************************
#ifndef DGPDE_h
#define DGPDE_h

#include <array>
#include <string>
#include <vector>
#include <memory>
#include <unordered_set>
#include <functional>

#include "Types.hpp"
#include "Fields.hpp"
#include "FaceData.hpp"
#include "UnsMesh.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "FunctionPrototypes.hpp"
#include "History.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

using ncomp_t = tk::ncomp_t;
using BCStateFn =
  std::vector< std::pair< std::vector< std::size_t >, tk::StateFn > >;

//! Extract BC configuration ignoring if BC not specified
//! \note A more preferable way of catching errors such as this function
//!   hides is during parsing, so that we don't even get here if BCs are
//!   not correctly specified. For now we simply ignore if BCs are not
//!   specified by allowing empty BC vectors from the user input.
struct ConfigBC {
  BCStateFn& state;    //!< BC state config: sidesets + statefn
  const std::vector< tk::StateFn >& fn;    //!< BC state functions
  std::size_t c;       //!< Counts BC types configured
  //! Constructor
  ConfigBC( BCStateFn& s,
            const std::vector< tk::StateFn >& f ) :
    state(s), fn(f), c(0) {}
  //! Function to call for each BC type
  template< typename U > void operator()( brigand::type_<U> ) {
    std::vector< std::size_t > cfg, v;
    // collect sidesets across all meshes
    for (const auto& ibc : g_inputdeck.get< tag::bc >()) {
      v.insert(v.end(), ibc.get< U >().begin(), ibc.get< U >().end());
    }
    if (v.size() > 0) cfg = v;
    Assert( fn.size() > c, "StateFn missing for BC type" );
    state.push_back( { cfg, fn[c++] } );
  }
};

//! State function for invalid/un-configured boundary conditions
[[noreturn]] tk::StateFn::result_type
invalidBC( ncomp_t, const std::vector< EOS >&,
           const std::vector< tk::real >&, tk::real, tk::real, tk::real,
           tk::real, const std::array< tk::real, 3> & );

//! \brief Partial differential equation base for discontinuous Galerkin PDEs
//! \details This class uses runtime polymorphism without client-side
//!   inheritance: inheritance is confined to the internals of the this class,
//!   invisible to client-code. The class exclusively deals with ownership
//!   enabling client-side value semantics. Credit goes to Sean Parent at Adobe:
//!   https://github.com/sean-parent/sean-parent.github.com/wiki/
//!   Papers-and-Presentations. For example client code that models a DGPDE,
//!   see inciter::CompFlow.
class DGPDE {

  private:
    using ncomp_t = tk::ncomp_t;

  public:
    //! Default constructor taking no arguments for Charm++
    explicit DGPDE() = default;

    //! \brief Constructor taking an object modeling Concept.
    //! \details The object of class T comes pre-constructed.
    //! \param[in] x Instantiated object of type T given by the template
    //!   argument.
    template< typename T > explicit DGPDE( T x ) :
      self( std::make_unique< Model<T> >( std::move(x) ) ) {}

    //! \brief Constructor taking a function pointer to a constructor of an
    //!   object modeling Concept.
    //! \details Passing std::function allows late execution of the constructor,
    //!   i.e., as late as inside this class' constructor, and thus usage from
    //!   a factory. Note that there are at least two different ways of using
    //!   this constructor:
    //!   - Bind T's constructor arguments and place it in std::function<T()>
    //!   and passing no arguments as args.... This case then instantiates the
    //!   model via its constructor and stores it in here.
    //!   - Bind a single placeholder argument to T's constructor and pass it in
    //!   as host's args..., which then forwards it to model's constructor. This
    //!   allows late binding, i.e., binding the argument only here.
    //! \see See also the wrapper tk::recordModel() which does the former and
    //!   tk::recordModelLate() which does the latter, both defined in
    //!   src/Base/Factory.h.
    //! \param[in] x Function pointer to a constructor of an object modeling
    //!    Concept.
    //! \param[in] args Zero or more constructor arguments
    template< typename T, typename...Args >
    explicit DGPDE( std::function<T(Args...)> x, Args&&... args ) :
      self( std::make_unique< Model<T> >(
              std::move( x( std::forward<Args>(args)... ) ) ) ) {}

    //! Public interface to find number of primitive quantities for the diff eq
    std::size_t nprim() const
    { return self->nprim(); }

    //! Public interface to find number of materials for the diff eq
    std::size_t nmat() const
    { return self->nmat(); }

    //! Public interface to find Dofs for each equation in pde system
    void numEquationDofs(std::vector< std::size_t >& numEqDof) const
    { return self->numEquationDofs(numEqDof); }

    //! Public interface to determine elements that lie inside the IC box
    void IcBoxElems( const tk::Fields& geoElem,
      std::size_t nielem,
      std::vector< std::unordered_set< std::size_t > >& inbox ) const
    { self->IcBoxElems( geoElem, nielem, inbox ); }

    //! Public interface to setting the initial conditions for the diff eq
    void initialize(
      const tk::Fields& L,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const std::vector< std::unordered_set< std::size_t > >& inbox,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&
        elemblkid,
      tk::Fields& unk,
      tk::real t,
      const std::size_t nielem ) const
    { self->initialize( L, inpoel, coord, inbox, elemblkid, unk, t, nielem ); }

    //! Public interface to computing the left-hand side matrix for the diff eq
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const
    { self->lhs( geoElem, l ); }

    //! Public interface to updating the interface cells for the diff eq
    void updateInterfaceCells( tk::Fields& unk,
                               std::size_t nielem,
                               std::vector< std::size_t >& ndofel ) const
    { self->updateInterfaceCells( unk, nielem, ndofel ); }

    //! Public interface to updating the primitives for the diff eq
    void updatePrimitives( const tk::Fields& unk,
                           const tk::Fields& L,
                           const tk::Fields& geoElem,
                           tk::Fields& prim,
                           std::size_t nielem ) const
    { self->updatePrimitives( unk, L, geoElem, prim, nielem ); }

    //! Public interface to cleaning up trace materials for the diff eq
    void cleanTraceMaterial( tk::real t,
                             const tk::Fields& geoElem,
                             tk::Fields& unk,
                             tk::Fields& prim,
                             std::size_t nielem ) const
    { self->cleanTraceMaterial( t, geoElem, unk, prim, nielem ); }

    //! Public interface to reconstructing the second-order solution
    void reconstruct( tk::real t,
                      const tk::Fields& geoFace,
                      const tk::Fields& geoElem,
                      const inciter::FaceData& fd,
                      const std::map< std::size_t, std::vector< std::size_t > >&
                        esup,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      tk::Fields& U,
                      tk::Fields& P ) const
    {
      self->reconstruct( t, geoFace, geoElem, fd, esup, inpoel, coord, U, P );
    }

    //! Public interface to limiting the second-order solution
    void limit( tk::real t,
                const tk::Fields& geoFace,
                const tk::Fields& geoElem,
                const inciter::FaceData& fd,
                const std::map< std::size_t, std::vector< std::size_t > >& esup,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& ndofel,
                const std::vector< std::size_t >& gid,
                const std::unordered_map< std::size_t, std::size_t >& bid,
                const std::vector< std::vector<tk::real> >& uNodalExtrm,
                const std::vector< std::vector<tk::real> >& pNodalExtrm,
                const std::vector< std::vector<tk::real> >& mtInv,
                tk::Fields& U,
                tk::Fields& P,
                std::vector< std::size_t >& shockmarker ) const
    {
      self->limit( t, geoFace, geoElem, fd, esup, inpoel, coord, ndofel, gid,
                   bid, uNodalExtrm, pNodalExtrm, mtInv, U, P, shockmarker );
    }

    //! Public interface to update the conservative variable solution
    void CPL( const tk::Fields& prim,
              const tk::Fields& geoElem,
              const std::vector< std::size_t >& inpoel,
              const tk::UnsMesh::Coords& coord,
              tk::Fields& unk,
              std::size_t nielem ) const
    {
      self->CPL( prim, geoElem, inpoel, coord, unk, nielem );
    }

    //! Public interface to getting the cell-averaged deformation gradients
    std::array< std::vector< tk::real >, 9 > cellAvgDeformGrad(
      const tk::Fields& U,
      std::size_t nielem ) const
    {
      return self->cellAvgDeformGrad( U, nielem );
    }

    //! Public interface to computing the P1 right-hand side vector
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const std::vector< std::unordered_set< std::size_t > >& boxelems,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              const tk::Fields& P,
              const std::vector< std::size_t >& ndofel,
              const tk::real dt,
              tk::Fields& R ) const
    {
      self->rhs( t, geoFace, geoElem, fd, inpoel, boxelems, coord, U, P,
                 ndofel, dt, R );
    }

    //! Evaluate the adaptive indicator and mark the ndof for each element
    void eval_ndof( std::size_t nunk,
                    const tk::UnsMesh::Coords& coord,
                    const std::vector< std::size_t >& inpoel,
                    const inciter::FaceData& fd,
                    const tk::Fields& unk,
                    const tk::Fields& prim,
                    inciter::ctr::PrefIndicatorType indicator,
                    std::size_t ndof,
                    std::size_t ndofmax,
                    tk::real tolref,
                    std::vector< std::size_t >& ndofel ) const
    {
      self->eval_ndof( nunk, coord, inpoel, fd, unk, prim, indicator, ndof,
        ndofmax, tolref, ndofel );
    }

    //! Public interface for computing the minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const inciter::FaceData& fd,
                 const tk::Fields& geoFace,
                 const tk::Fields& geoElem,
                 const std::vector< std::size_t >& ndofel,
                 const tk::Fields& U,
                 const tk::Fields& P,
                 const std::size_t nielem ) const
    { return self->dt( coord, inpoel, fd, geoFace, geoElem, ndofel, U,
                       P, nielem ); }

    //! Public interface to returning maps of output var functions
    std::map< std::string, tk::GetVarFn > OutVarFn() const
    { return self->OutVarFn(); }

    //! Public interface to returning analytic field output labels
    std::vector< std::string > analyticFieldNames() const
    { return self->analyticFieldNames(); }

    //! Public interface to returning time history field output labels
    std::vector< std::string > histNames() const { return self->histNames(); }

    //! Public interface to returning variable names
    std::vector< std::string > names() const { return self->names(); }

    //! Public interface to returning surface field output
    std::vector< std::vector< tk::real > >
    surfOutput( const std::map< int, std::vector< std::size_t > >& bnd,
                tk::Fields& U ) const
    { return self->surfOutput( bnd, U ); }

    //! Public interface to return point history output
    std::vector< std::vector< tk::real > >
    histOutput( const std::vector< HistData >& h,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const tk::Fields& U,
      const tk::Fields& P ) const
    { return self->histOutput( h, inpoel, coord, U, P ); }

    //! Public interface to returning analytic solution
    tk::InitializeFn::result_type
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return self->analyticSolution( xi, yi, zi, t ); }

    //! Public interface to returning the analytic solution for conserved vars
    tk::InitializeFn::result_type
    solution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return self->solution( xi, yi, zi, t ); }

    //! Public interface to returning the specific total energy
    tk::real
    sp_totalenergy( std::size_t e, const tk::Fields& unk ) const
    { return self->sp_totalenergy( e, unk ); }

    //! Copy assignment
    DGPDE& operator=( const DGPDE& x )
    { DGPDE tmp(x); *this = std::move(tmp); return *this; }
    //! Copy constructor
    DGPDE( const DGPDE& x ) : self( x.self->copy() ) {}
    //! Move assignment
    DGPDE& operator=( DGPDE&& ) noexcept = default;
    //! Move constructor
    DGPDE( DGPDE&& ) noexcept = default;

  private:
    //! \brief Concept is a pure virtual base class specifying the requirements
    //!   of polymorphic objects deriving from it
    struct Concept {
      Concept() = default;
      Concept( const Concept& ) = default;
      virtual ~Concept() = default;
      virtual Concept* copy() const = 0;
      virtual std::size_t nprim() const = 0;
      virtual std::size_t nmat() const = 0;
      virtual void numEquationDofs(std::vector< std::size_t >&) const = 0;
      virtual void IcBoxElems( const tk::Fields&,
        std::size_t,
        std::vector< std::unordered_set< std::size_t > >& ) const = 0;
      virtual void initialize(
        const tk::Fields&,
        const std::vector< std::size_t >&,
        const tk::UnsMesh::Coords&,
        const std::vector< std::unordered_set< std::size_t > >&,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&,
        tk::Fields&,
        tk::real,
        const std::size_t nielem ) const = 0;
      virtual void lhs( const tk::Fields&, tk::Fields& ) const = 0;
      virtual void updateInterfaceCells( tk::Fields&,
                                         std::size_t,
                                         std::vector< std::size_t >& ) const = 0;
      virtual void updatePrimitives( const tk::Fields&,
                                     const tk::Fields&,
                                     const tk::Fields&,
                                     tk::Fields&,
                                     std::size_t ) const = 0;
      virtual void cleanTraceMaterial( tk::real,
                                       const tk::Fields&,
                                       tk::Fields&,
                                       tk::Fields&,
                                       std::size_t ) const = 0;
      virtual void reconstruct( tk::real,
                                const tk::Fields&,
                                const tk::Fields&,
                                const inciter::FaceData&,
                                const std::map< std::size_t,
                                  std::vector< std::size_t > >&,
                                const std::vector< std::size_t >&,
                                const tk::UnsMesh::Coords&,
                                tk::Fields&,
                                tk::Fields& ) const = 0;
      virtual void limit( tk::real,
                          const tk::Fields&,
                          const tk::Fields&,
                          const inciter::FaceData&,
                          const std::map< std::size_t,
                            std::vector< std::size_t > >&,
                          const std::vector< std::size_t >&,
                          const tk::UnsMesh::Coords&,
                          const std::vector< std::size_t >&,
                          const std::vector< std::size_t >&,
                          const std::unordered_map< std::size_t, std::size_t >&,
                          const std::vector< std::vector<tk::real> >&,
                          const std::vector< std::vector<tk::real> >&,
                          const std::vector< std::vector<tk::real> >&,
                          tk::Fields&,
                          tk::Fields&,
                          std::vector< std::size_t >& ) const = 0;
      virtual void CPL( const tk::Fields&,
                        const tk::Fields&,
                        const std::vector< std::size_t >&,
                        const tk::UnsMesh::Coords&,
                        tk::Fields&,
                        std::size_t ) const = 0;
      virtual std::array< std::vector< tk::real >, 9 > cellAvgDeformGrad(
        const tk::Fields&,
        std::size_t ) const = 0;
      virtual void rhs( tk::real,
                        const tk::Fields&,
                        const tk::Fields&,
                        const inciter::FaceData&,
                        const std::vector< std::size_t >&,
                        const std::vector< std::unordered_set< std::size_t > >&,
                        const tk::UnsMesh::Coords&,
                        const tk::Fields&,
                        const tk::Fields&,
                        const std::vector< std::size_t >&,
                        const tk::real,
                        tk::Fields& ) const = 0;
      virtual void eval_ndof( std::size_t,
                              const tk::UnsMesh::Coords&,
                              const std::vector< std::size_t >&,
                              const inciter::FaceData&,
                              const tk::Fields&,
                              const tk::Fields&,
                              inciter::ctr::PrefIndicatorType,
                              std::size_t,
                              std::size_t,
                              tk::real,
                              std::vector< std::size_t >& ) const = 0;
      virtual tk::real dt( const std::array< std::vector< tk::real >, 3 >&,
                           const std::vector< std::size_t >&,
                           const inciter::FaceData&,
                           const tk::Fields&,
                           const tk::Fields&,
                           const std::vector< std::size_t >&,
                           const tk::Fields&,
                           const tk::Fields&,
                           const std::size_t ) const = 0;
      virtual std::map< std::string, tk::GetVarFn > OutVarFn() const = 0;
      virtual std::vector< std::string > analyticFieldNames() const = 0;
      virtual std::vector< std::string > histNames() const = 0;
      virtual std::vector< std::string > names() const = 0;
      virtual std::vector< std::vector< tk::real > > surfOutput(
        const std::map< int, std::vector< std::size_t > >&,
        tk::Fields& ) const = 0;
      virtual std::vector< std::vector< tk::real > > histOutput(
        const std::vector< HistData >&,
        const std::vector< std::size_t >&,
        const tk::UnsMesh::Coords&,
        const tk::Fields&,
        const tk::Fields& ) const = 0;
      virtual tk::InitializeFn::result_type analyticSolution(
        tk::real xi, tk::real yi, tk::real zi, tk::real t ) const = 0;
      virtual tk::InitializeFn::result_type solution(
        tk::real xi, tk::real yi, tk::real zi, tk::real t ) const = 0;
      virtual tk::real sp_totalenergy(
        std::size_t, const tk::Fields& ) const = 0;
    };

    //! \brief Model models the Concept above by deriving from it and overriding
    //!   the virtual functions required by Concept
    template< typename T >
    struct Model : Concept {
      explicit Model( T x ) : data( std::move(x) ) {}
      Concept* copy() const override { return new Model( *this ); }
      std::size_t nprim() const override
      { return data.nprim(); }
      std::size_t nmat() const override
      { return data.nmat(); }
      void numEquationDofs(std::vector< std::size_t >& numEqDof) const override
      { data.numEquationDofs(numEqDof); }
      void IcBoxElems( const tk::Fields& geoElem,
        std::size_t nielem,
        std::vector< std::unordered_set< std::size_t > >& inbox )
      const override { data.IcBoxElems( geoElem, nielem, inbox ); }
      void initialize(
        const tk::Fields& L,
        const std::vector< std::size_t >& inpoel,
        const tk::UnsMesh::Coords& coord,
        const std::vector< std::unordered_set< std::size_t > >& inbox,
        const std::unordered_map< std::size_t, std::set< std::size_t > >&
          elemblkid,
        tk::Fields& unk,
        tk::real t,
        const std::size_t nielem )
      const override { data.initialize( L, inpoel, coord, inbox, elemblkid, unk,
        t, nielem ); }
      void lhs( const tk::Fields& geoElem, tk::Fields& l ) const override
      { data.lhs( geoElem, l ); }
      void updateInterfaceCells( tk::Fields& unk,
                                 std::size_t nielem,
                                 std::vector< std::size_t >& ndofel )
      const override { data.updateInterfaceCells( unk, nielem, ndofel ); }
      void updatePrimitives( const tk::Fields& unk,
                             const tk::Fields& L,
                             const tk::Fields& geoElem,
                             tk::Fields& prim,
                             std::size_t nielem )
      const override { data.updatePrimitives( unk, L, geoElem, prim, nielem ); }
      void cleanTraceMaterial( tk::real t,
                               const tk::Fields& geoElem,
                               tk::Fields& unk,
                               tk::Fields& prim,
                               std::size_t nielem )
      const override { data.cleanTraceMaterial( t, geoElem, unk, prim, nielem ); }
      void reconstruct( tk::real t,
                        const tk::Fields& geoFace,
                        const tk::Fields& geoElem,
                        const inciter::FaceData& fd,
                        const std::map< std::size_t,
                          std::vector< std::size_t > >& esup,
                        const std::vector< std::size_t >& inpoel,
                        const tk::UnsMesh::Coords& coord,
                        tk::Fields& U,
                        tk::Fields& P ) const override
      {
        data.reconstruct( t, geoFace, geoElem, fd, esup, inpoel, coord, U, P );
      }
      void limit( tk::real t,
                  const tk::Fields& geoFace,
                  const tk::Fields& geoElem,
                  const inciter::FaceData& fd,
                  const std::map< std::size_t, std::vector< std::size_t > >&
                    esup,
                  const std::vector< std::size_t >& inpoel,
                  const tk::UnsMesh::Coords& coord,
                  const std::vector< std::size_t >& ndofel,
                  const std::vector< std::size_t >& gid,
                  const std::unordered_map< std::size_t, std::size_t >& bid,
                  const std::vector< std::vector<tk::real> >& uNodalExtrm,
                  const std::vector< std::vector<tk::real> >& pNodalExtrm,
                  const std::vector< std::vector<tk::real> >& mtInv,
                  tk::Fields& U,
                  tk::Fields& P,
                  std::vector< std::size_t >& shockmarker ) const override
      {
        data.limit( t, geoFace, geoElem, fd, esup, inpoel, coord, ndofel, gid,
                    bid, uNodalExtrm, pNodalExtrm, mtInv, U, P, shockmarker );
      }
      void CPL( const tk::Fields& prim,
                const tk::Fields& geoElem,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                tk::Fields& unk,
                std::size_t nielem ) const override
      {
        data.CPL( prim, geoElem, inpoel, coord, unk, nielem );
      }
      std::array< std::vector< tk::real >, 9 > cellAvgDeformGrad(
        const tk::Fields& U,
        std::size_t nielem ) const override
      {
        return data.cellAvgDeformGrad( U, nielem );
      }
      void rhs(
        tk::real t,
        const tk::Fields& geoFace,
        const tk::Fields& geoElem,
        const inciter::FaceData& fd,
        const std::vector< std::size_t >& inpoel,
        const std::vector< std::unordered_set< std::size_t > >& boxelems,
        const tk::UnsMesh::Coords& coord,
        const tk::Fields& U,
        const tk::Fields& P,
        const std::vector< std::size_t >& ndofel,
        const tk::real dt,
        tk::Fields& R ) const override
      {
        data.rhs( t, geoFace, geoElem, fd, inpoel, boxelems, coord, U, P,
                  ndofel, dt, R );
      }
      void eval_ndof( std::size_t nunk,
                      const tk::UnsMesh::Coords& coord,
                      const std::vector< std::size_t >& inpoel,
                      const inciter::FaceData& fd,
                      const tk::Fields& unk,
                      const tk::Fields& prim,
                      inciter::ctr::PrefIndicatorType indicator,
                      std::size_t ndof,
                      std::size_t ndofmax,
                      tk::real tolref,
                      std::vector< std::size_t >& ndofel ) const override
      { data.eval_ndof( nunk, coord, inpoel, fd, unk, prim, indicator, ndof,
                        ndofmax, tolref, ndofel ); }
      tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                   const std::vector< std::size_t >& inpoel,
                   const inciter::FaceData& fd,
                   const tk::Fields& geoFace,
                   const tk::Fields& geoElem,
                   const std::vector< std::size_t >& ndofel,
                   const tk::Fields& U,
                   const tk::Fields& P,
                   const std::size_t nielem ) const override
      { return data.dt( coord, inpoel, fd, geoFace, geoElem, ndofel,
                        U, P, nielem ); }
      std::map< std::string, tk::GetVarFn > OutVarFn() const override
      { return data.OutVarFn(); }
      std::vector< std::string > analyticFieldNames() const override
      { return data.analyticFieldNames(); }
      std::vector< std::string > histNames() const override
      { return data.histNames(); }
      std::vector< std::string > names() const override
      { return data.names(); }
      std::vector< std::vector< tk::real > > surfOutput(
        const std::map< int, std::vector< std::size_t > >& bnd,
        tk::Fields& U ) const override
      { return data.surfOutput( bnd, U ); }
      std::vector< std::vector< tk::real > > histOutput(
        const std::vector< HistData >& h,
        const std::vector< std::size_t >& inpoel,
        const tk::UnsMesh::Coords& coord,
        const tk::Fields& U,
        const tk::Fields& P ) const override
      { return data.histOutput( h, inpoel, coord, U, P ); }
      tk::InitializeFn::result_type
      analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t )
       const override { return data.analyticSolution( xi, yi, zi, t ); }
      tk::InitializeFn::result_type
      solution( tk::real xi, tk::real yi, tk::real zi, tk::real t )
       const override { return data.solution( xi, yi, zi, t ); }
      tk::real sp_totalenergy( std::size_t e, const tk::Fields& unk )
       const override { return data.sp_totalenergy( e, unk ); }
      T data;
    };

    std::unique_ptr< Concept > self;    //!< Base pointer used polymorphically
};

} // inciter::

#endif // DGPDE_h
