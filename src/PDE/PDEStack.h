// *****************************************************************************
/*!
  \file      src/PDE/PDEStack.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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

#include "NoWarning/back.h"
#include "NoWarning/front.h"

#include "Tags.h"
#include "Exception.h"
#include "Factory.h"
#include "CGPDE.h"
#include "DGPDE.h"
#include "PDEFactory.h"
#include "SystemComponents.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief Partial differential equations stack
class PDEStack {

  private:
    using ncomp_t = tk::ctr::ncomp_type;

  public:
    //! Constructor: register partial differential equations into factory
    explicit PDEStack();

    //! Instantiate selected PDEs using continuous Galerkin discretization
    std::vector< CGPDE > selectedCG() const;

    //! Instantiate selected PDEs using discontinuous Galerkin discretization
    std::vector< DGPDE > selectedDG() const;

    //! Constant accessor to CGPDE factory
    //! \return Constant reference to the CGPDE factory
    const CGFactory& factory() const { return m_cgfactory; }

    //! Return info on selected partial differential equations
    std::vector< std::vector< std::pair< std::string, std::string > > > info()
    const;

    //! \brief Return number of unique equation types registered
    //! \return The number of unique equation types registered in the factory
    std::size_t ntypes() const { return m_eqTypes.size(); }

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
      Assert( c < (g_inputdeck.get< tag::component, EqTag >().size()),
              "The number of scalar components is unspecified for the PDE to "
              "be instantiated. This is most likely a grammar error in the "
              "parser. The parser should not allow the user to select a PDE "
              "without configuring the number of scalar components the "
              "equation consists of. See inciter::deck::check_eq." );
      if ( g_inputdeck.get< tag::component, EqTag >()[c] ) {
        // re-create key and search for it
        ctr::PDEKey key{ eq,
          g_inputdeck.get< tag::param, EqTag, tag::physics >()[c],
          g_inputdeck.get< tag::param, EqTag, tag::problem >()[c] };
        const auto it = f.find( key );
        Assert( it != end( f ),
                "Can't find PDE with key('" +
                  ctr::PDE().name( key.get< tag::pde >() ) + "', '" +
                  ctr::Physics().name( key.get< tag::physics >() ) + "', '" +
                  ctr::Problem().name( key.get< tag::problem >() ) +
                  + "') in factory" );
        return it->second( c );    // instantiate and return PDE object
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

    //! PDE factory for continuous Galerkin discretization
    CGFactory m_cgfactory;
    //! PDE factory for discontinuous Galerkin discretization
    DGFactory m_dgfactory;
    //! Counters for equation types
    std::set< ctr::PDEType > m_eqTypes;
};

} // inciter::

#endif // PDEStack_h
