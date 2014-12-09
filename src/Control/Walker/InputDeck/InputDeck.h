//******************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 08:29:28 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Walker's input deck
  \details   Walker's input deck
*/
//******************************************************************************
#ifndef WalkerInputDeck_h
#define WalkerInputDeck_h

#include <limits>

#include <Control.h>
#include <Walker/CmdLine/CmdLine.h>
#include <Walker/Components.h>

namespace walker {
namespace ctr {

//! InputDeck : Control< specialized to Walker >, see Types.h,
class InputDeck :
  public tk::Control< // tag           type
                      tag::title,      std::string,
                      tag::selected,   selects,
                      tag::discr,      discretization,
                      tag::component,  ncomps,
                      tag::interval,   intervals,
                      tag::cmd,        CmdLine,
                      tag::param,      parameters,
                      tag::stat,       std::vector< tk::ctr::Product >,
                      tag::pdf,        std::vector< tk::ctr::Probability >,
                      tag::error,      std::vector< std::string > > {

  private:
    //! Function object for extracting dependent variable vectors from
    //! components
    struct depvar {
      const InputDeck* const m_host;
      std::vector< std::vector< char > >& m_vars;
      depvar( const InputDeck* const host,
              std::vector< std::vector< char > >& vars ) :
        m_host( host ), m_vars( vars ) {}
      template< typename U > void operator()( U ) {
        m_vars.push_back( m_host->get< tag::param, U, tag::depvar >() );
      }
    };

  public:
    //! Constructor: set all defaults
    InputDeck() {
      // Default discretization parameters
      set< tag::discr, tag::npar >( 1 );
      set< tag::discr, tag::nstep >( std::numeric_limits< uint64_t >::max() );
      set< tag::discr, tag::term >( 1.0 );
      set< tag::discr, tag::dt >( 0.5 );
      set< tag::discr, tag::precision >( std::cout.precision() );
      // Default intervals
      set< tag::interval, tag::tty >( 1 );
      set< tag::interval, tag::stat >( 1 );
      set< tag::interval, tag::pdf >( 1 );
      // Default requested statistics
      set< tag::stat >( std::vector< tk::ctr::Product >() );
    }

    //! Extract data on whether to plot ordinary moments of requested statistics
    std::vector< bool > plotOrdinary() const {
      std::vector< bool > plot;
      for (const auto& product : get< tag::stat >()) {
        if (ordinary( product )) {
          plot.emplace_back( false );
          for (const auto& term : product) if (term.plot) plot.back() = true;
        }
      }
      return plot;
    }

    //! Extract moment names of requested statistics
    std::vector< std::string > momentNames( std::function<
      bool ( const std::vector< tk::ctr::Term >& ) > momentType ) const
    {
      std::vector< std::string > names;
      for (const auto& product : get< tag::stat >()) {
        if (momentType( product )) {
          names.emplace_back( std::string() );
          for (const auto& term : product)
            names.back() += tk::ctr::FieldVar( term.var, term.field );
        }
      }
      return names;
    }

    // Count number requested PDFs with given number of sample space dimensions
    template< std::size_t d >
    std::size_t npdf() const {
      std::size_t n = 0;
      for (const auto& bs : get< tag::discr, tag::binsize >())
        if (bs.size() == d) ++n;
      return n;
    }

    //! Extract vector of vector of dependent variables from components
    std::vector< std::vector< char > > depvars() const {
      std::vector< std::vector< char > > vars;
      boost::mpl::for_each< ncomps::tags >( depvar( this, vars ) );
      return vars;
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      tk::Control< tag::title,      std::string,
                   tag::selected,   selects,
                   tag::discr,      discretization,
                   tag::component,  ncomps,
                   tag::interval,   intervals,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tag::stat,       std::vector< tk::ctr::Product >,
                   tag::pdf,        std::vector< tk::ctr::Probability >,
                   tag::error,      std::vector< std::string > >::pup(p);
    }
    friend void operator|( PUP::er& p, InputDeck& c ) { c.pup(p); }
};

} // ctr::
} // Walker::

#endif // WalkerInputDeck_h
