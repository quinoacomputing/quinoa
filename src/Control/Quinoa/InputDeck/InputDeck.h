//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Sat 27 Sep 2014 09:58:49 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's input deck
  \details   Quinoa's input deck
*/
//******************************************************************************
#ifndef QuinoaInputDeck_h
#define QuinoaInputDeck_h

#include <limits>

#include <Control.h>
#include <Quinoa/CmdLine/CmdLine.h>
#include <Quinoa/Components.h>

namespace quinoa {
namespace ctr {

//! InputDeck : Control< specialized to Quinoa >, see Types.h,
class InputDeck :
  public tk::Control< // tag           type
                      tag::title,      std::string,
                      tag::selected,   selects,
                      tag::discr,      discretization,
                      tag::component,  ncomps< unsigned int >,
                      tag::interval,   intervals,
                      tag::cmd,        CmdLine,
                      tag::param,      parameters,
                      tag::stat,       std::vector< Product >,
                      tag::pdf,        std::vector< Probability >,
                      tk::tag::error,  std::vector< std::string > > {

  public:
    //! Constructor: set all defaults
    InputDeck() {
      // Default discretization parameters
      set< tag::discr, tag::npar >( 1 );
      set< tag::discr, tag::nstep >( std::numeric_limits< uint64_t >::max() );
      set< tag::discr, tag::term >( 1.0 );
      set< tag::discr, tag::dt >( 0.5 );
      set< tag::discr, tag::dt >( 0.5 );
      // Default intervals
      set< tag::interval, tag::tty >( 1 );
      set< tag::interval, tag::dump >( 1 );
      set< tag::interval, tag::stat >( 1 );
      set< tag::interval, tag::pdf >( 1 );
      set< tag::interval, tag::glob >( 1 );
      // Default beta mass model parameters
      set< tag::param, tag::beta, tag::atwood >( 0.5 );
      // Default gamma mix model parameters
      set< tag::param, tag::gamma, tag::c1 >( 0.5 );
      set< tag::param, tag::gamma, tag::c2 >( 0.73 );
      set< tag::param, tag::gamma, tag::c3 >( 5.0 );
      set< tag::param, tag::gamma, tag::c4 >( 0.25 );
      // Default simplified Langevin hydro model parameters
      set< tag::param, tag::slm, tag::c0 >( 2.1 );
      // Default generalized Langevin hydro model parameters
      set< tag::param, tag::slm, tag::c0 >( 2.1 );
      // Default requested statistics
      set< tag::stat >( std::vector< Product >() );
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
      bool ( const std::vector< ctr::Term >& ) > momentType ) const
    {
      std::vector< std::string > names;
      for (const auto& product : get< tag::stat >()) {
        if (momentType( product )) {
          names.emplace_back( std::string() );
          for (const auto& term : product)
            names.back() += ctr::FieldVar( term.var, term.field );
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

    // PDF information
    struct PDFInfo {
      const std::string& name;                  //!< identifier
      const std::vector< tk::real >& exts;      //!< extents
    };

    //! Find PDF information given the sample space dimension and its index
    //! \param[in]  idx  Index of the PDF with given sample space dimension
    template< std::size_t d >
    PDFInfo pdf( long int idx ) const {
      const auto& binsizes = get< tag::discr, tag::binsize >();
      const auto& names = get< tag::cmd, tag::io, tag::pdfnames >();
      const auto& exts = get< tag::discr, tag::extent >();
      Assert( binsizes.size() == names.size(),
              "Number of binsizes vector and the number of PDF names must "
              "equal in InputDeck::pdfname()." );
      long int n = -1;
      long int i = 0;
      for (const auto& bs : binsizes) {
        if (bs.size() == d) ++n;
        if (n == idx) return { names[i], exts[i] };
        ++i;
      }
      Throw( "Cannot find PDF name." );
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      tk::Control< tag::title,      std::string,
                   tag::selected,   selects,
                   tag::discr,      discretization,
                   tag::component,  ncomps< unsigned int >,
                   tag::interval,   intervals,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tag::stat,       std::vector< Product >,
                   tag::pdf,        std::vector< Probability >,
                   tk::tag::error,  std::vector< std::string > >::pup(p);
    }
    friend void operator|( PUP::er& p, InputDeck& c ) { c.pup(p); }
};

} // ctr::
} // quinoa::

#endif // QuinoaInputDeck_h
