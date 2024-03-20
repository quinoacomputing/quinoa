// *****************************************************************************
/*!
  \file      src/PDE/PDEStack.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Stack of differential equations
  \details   This file declares class PDEStack, which implements various
    functionality related to registering and instantiating partial differential
    equation types. Registration and instantiation use a partial differential
    equation factory, which is a std::map (an associative container),
    associating unique partial differential equation keys to their constructor
    calls. For more details, see the in-code documentation of the constructor.
*/
// *****************************************************************************
#ifndef PDEStack_h
#define PDEStack_h

#include <map>
#include <set>
#include <string>
#include <vector>

#include "NoWarning/back.hpp"
#include "NoWarning/front.hpp"

#include "Tags.hpp"
#include "Types.hpp"
#include "Exception.hpp"
#include "Factory.hpp"
#include "CGPDE.hpp"
#include "DGPDE.hpp"
#include "FVPDE.hpp"
#include "PDEFactory.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"

namespace inciter {

extern ctr::New2InputDeck g_newinputdeck;

//! \brief Partial differential equations stack
class PDEStack {

  private:
    using ncomp_t = tk::ncomp_t;

  public:
    //! Constructor: register partial differential equations into factory
    explicit PDEStack();

    //! Instantiate selected PDEs using continuous Galerkin discretization
    std::vector< CGPDE > selectedCG() const;

    //! Instantiate selected PDEs using discontinuous Galerkin discretization
    std::vector< DGPDE > selectedDG() const;

    //! Instantiate selected PDEs using finite volume discretization
    std::vector< FVPDE > selectedFV() const;

    //! Constant accessor to CGPDE factory
    //! \return Constant reference to the CGPDE factory
    const CGFactory& cgfactory() const { return m_cgfactory; }

    //! Constant accessor to DGPDE factory
    //! \return Constant reference to the DGPDE factory
    const DGFactory& dgfactory() const { return m_dgfactory; }

    //! Constant accessor to FVPDE factory
    //! \return Constant reference to the FVPDE factory
    const FVFactory& fvfactory() const { return m_fvfactory; }

    //! Return info on selected partial differential equations
    std::vector< std::vector< std::pair< std::string, std::string > > > info()
    const;

    //! Return number of unique equation types registered into the CG factory
    //! \return The number of unique equation types registered into the CG
    //!   factory the factory
    std::size_t cgntypes() const { return m_cgEqTypes.size(); }

    //! Return number of unique equation types registered into the DG factory
    //! \return The number of unique equation types registered into the DG
    //!   factory the factory
    std::size_t dgntypes() const { return m_dgEqTypes.size(); }

    //! Return number of unique equation types registered into the FV factory
    //! \return The number of unique equation types registered into the FV
    //!   factory the factory
    std::size_t fvntypes() const { return m_fvEqTypes.size(); }

  private:
    //! \brief Instantiate a partial differential equation
    //! \details The template argument, EqTag, is used to find the given
    //!   partial differential equation in the input deck, the hierarchical data
    //!   filled during control file parsing, containing user input. The
    //!   template argument Factory specifies which factory we search in. The
    //!   template argument PDE specifies the "base" PDE type that the
    //!   instantiated "child" PDE class object is used polymorphically with.
    //! \param[in] f Factory in which to search PDE in
    //! \param[in] eq The unique partial differential equation key whose object
    //!   to instantiate.
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   partial differential equations by type.
    template< class EqTag, class Factory, class PDE >
    PDE createPDE( const Factory& f,
                   ctr::PDEType eq,
                   std::map< ctr::PDEType, ncomp_t >& cnt ) const
    {
      auto c = ++cnt[ eq ];   // count eqs
      --c;                    // used to index vectors starting with 0
      const auto& nc = g_newinputdeck.get< newtag::ncomp >();
      if ( nc ) {
        // re-create key and search for it
        ctr::PDEKey key{{ eq,
          g_newinputdeck.get< newtag::physics >(),
          g_newinputdeck.get< EqTag, newtag::problem >() }};
        const auto it = f.find( key );
        Assert( it != end( f ),
                "Can't find PDE with key('" +
                  ctr::PDE().name( key.get< tag::pde >() ) + "', '" +
                  ctr::Physics().name( key.get< tag::physics >() ) + "', '" +
                  ctr::Problem().name( key.get< tag::problem >() ) +
                  + "') in factory" );
        // Associate equation system index (value) to all variable offsets
        for (ncomp_t i=0; i<nc; ++i) g_newinputdeck.get<newtag::sys>()[i] = c;
        // instantiate and return PDE object
        return it->second();
      } else Throw ( "Can't create PDE with zero components" );
    }

    //! Wrapper of createPDE specialized for registering CG PDEs
    //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   partial differential equations by type.
    //! \details The sole reason for this function is to simplify client-code
    //!   calling createPDE specialized to CG PDEs
    template< class EqTag >
    CGPDE createCG( ctr::PDEType t, std::map< ctr::PDEType, ncomp_t >& cnt )
    const {
      return createPDE< EqTag, CGFactory, CGPDE >( m_cgfactory, t, cnt );
    }

    //! Wrapper of createPDE specialized for registering DG PDEs
    //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   partial differential equations by type.
    //! \details The sole reason for this function is to simplify client-code
    //!   calling createPDE specialized to DG PDEs
    template< class EqTag >
    DGPDE createDG( ctr::PDEType t, std::map< ctr::PDEType, ncomp_t >& cnt )
    const {
      return createPDE< EqTag, DGFactory, DGPDE >( m_dgfactory, t, cnt );
    }

    //! Wrapper of createPDE specialized for registering FV PDEs
    //! \param[in] t Enum selecting PDE type, Control/Inciter/Options/PDE.h
    //! \param[in,out] cnt Counter, a std::map, that counts all instantiated
    //!   partial differential equations by type.
    //! \details The sole reason for this function is to simplify client-code
    //!   calling createPDE specialized to FV PDEs
    template< class EqTag >
    FVPDE createFV( ctr::PDEType t, std::map< ctr::PDEType, ncomp_t >& cnt )
    const {
      return createPDE< EqTag, FVFactory, FVPDE >( m_fvfactory, t, cnt );
    }

    //! PDE factory for continuous Galerkin discretization
    CGFactory m_cgfactory;
    //! PDE factory for discontinuous Galerkin discretization
    DGFactory m_dgfactory;
    //! PDE factory for finite volume discretization
    FVFactory m_fvfactory;
    //! Counters for equation types registered into the CG factory
    std::set< ctr::PDEType > m_cgEqTypes;
    //! Counters for equation types registered into the DG factory
    std::set< ctr::PDEType > m_dgEqTypes;
    //! Counters for equation types registered into the FV factory
    std::set< ctr::PDEType > m_fvEqTypes;
};

} // inciter::

#endif // PDEStack_h
