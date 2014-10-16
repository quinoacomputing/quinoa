//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Fri 10 Oct 2014 02:10:12 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Statistics
  \details   Computing ordinary and central moments
*/
//******************************************************************************
#ifndef Statistics_h
#define Statistics_h

#include <Types.h>
#include <Quinoa/InputDeck/InputDeck.h>
#include <ParticleProperties.h>
#include <UniPDF.h>
#include <BiPDF.h>
#include <TriPDF.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck;

//! Statistics estimator
class Statistics {

  public:
    //! Constructor
    explicit Statistics( const ParProps& particles );

    //! Accumulate (i.e., only do the sum for) ordinary moments
    void accumulateOrd();

    //! Accumulate (i.e., only do the sum for) central moments
    void accumulateCen( const std::vector< tk::real >& ord );

    //! Accumulate (i.e., only do the sum for) ordinary PDFs
    void accumulateOrdPDF();

    //! Accumulate (i.e., only do the sum for) central PDFs
    void accumulateCenPDF();

    //! Ordinary moments accessor
    const std::vector< tk::real >& ord() const noexcept { return m_ordinary; }

    //! Central moments accessor
    const std::vector< tk::real >& ctr() const noexcept { return m_central; }

    //! Ordinary univariate PDFs accessor
    const std::vector< UniPDF >& oupdf() const noexcept { return m_ordupdf; }

    //! Ordinary bivariate PDFs accessor
    const std::vector< BiPDF >& obpdf() const noexcept { return m_ordbpdf; }

    //! Ordinary trivariate PDFs accessor
    const std::vector< TriPDF >& otpdf() const noexcept { return m_ordtpdf; }

    //! Central univariate PDFs accessor
    const std::vector< UniPDF >& cupdf() const noexcept { return m_cenupdf; }

    //! Central bivariate PDFs accessor
    const std::vector< BiPDF >& cbpdf() const noexcept { return m_cenbpdf; }

    //! Central trivariate PDFs accessor
    const std::vector< TriPDF >& ctpdf() const noexcept { return m_centpdf; }

  private:
    // Get offset type out of InputDeck for computing map< depvar, offset >
    using ncomps = ctr::InputDeck::nT< tag::component >;

    // Map associating a dependend variable to its offset in the particle data.
    // We use a case-insensitive character comparison functor for the offset
    // map, since the keys (dependent variables) in the offset map are only used
    // to indicate the equation's dependent variable, however, queries (using
    // find) are fired up for both ordinary and central moments (which are upper
    // and lower case) for which the offset (for the same depvar) should be the
    // same.
    using OffsetMap =
      std::map< char, ncomps::ncomp, ctr::CaseInsensitiveCharLess >;

    //! Function object for creating a run-time std::map< depvar, offset >
    struct depvar {
      OffsetMap& m_map;
      explicit depvar( OffsetMap& map ) : m_map(map) {}
      template< typename U > void operator()( U ) {
        // Loop through the vector of dependent variables for a differential
        // equation type given by type U and insert an entry to the offset map
        // for each equation, i.e., handle all potentially multiple equations of
        // the same type.
        ncomps::ncomp c = 0;
        for (const auto& v : g_inputdeck.get< tag::param, U, tag::depvar >())
          m_map[ v ] =
            g_inputdeck.get< tag::component >().template offset< U >( c++ );
      }
    };

    //! Setup ordinary moments
    void setupOrdinary( const OffsetMap& offset );

    //! Setup central moments
    void setupCentral( const OffsetMap& offset );

    //! Setup PDFs
    void setupPDF( const OffsetMap& offset );

    //! Return mean for fluctuation
    std::size_t mean(const ctr::Term& term) const;

    const ParProps& m_particles;              //!< Particle properties

    // Statistics

    //! Instantaneous variable pointers for computing ordinary moments
    std::vector< std::vector< const tk::real* > > m_instOrd;
    std::vector< tk::real > m_ordinary;        //!< Ordinary moments
    std::vector< ctr::FieldVar > m_ordFieldVar;//!< Ordinary moment field names
    std::size_t m_nord;                        //!< Number of ordinary moments

    //! Instantaneous variable pointers for computing central moments
    std::vector< std::vector< const tk::real* > > m_instCen;
    std::vector< tk::real > m_central;         //!< Central moments
    //! Ordinary moments about which to compute central moments
    std::vector< std::vector< const tk::real* > > m_ctr;
    std::size_t m_ncen;                        //!< Number of central moments

    // Univariate probability density functions

    //! Instantaneous variable pointers for computing ordinary univariate PDFs
    std::vector< std::vector< const tk::real* > > m_instOrdUniPDF;
    std::vector< UniPDF > m_ordupdf;           //!< Ordinary univariate PDFs

    //! Instantaneous variable pointers for computing central univariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenUniPDF;
    std::vector< UniPDF > m_cenupdf;           //!< Central univariate PDFs
    //! Ordinary moments about which to compute central univariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrUniPDF;

    // Bivariate probability density functions

    //! Instantaneous variable pointers for computing ordinary bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instOrdBiPDF;
    std::vector< BiPDF > m_ordbpdf;           //!< Ordinary bivariate PDFs

    //! Instantaneous variable pointers for computing central bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenBiPDF;
    std::vector< BiPDF > m_cenbpdf;           //!< Central bivariate PDFs
    //! Ordinary moments about which to compute central bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrBiPDF;

    // Trivariate probability density functions

    //! Instantaneous variable pointers for computing ordinary trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instOrdTriPDF;
    std::vector< TriPDF > m_ordtpdf;          //!< Ordinary trivariate PDFs

    //! Instantaneous variable pointers for computing central trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenTriPDF;
    std::vector< TriPDF > m_centpdf;          //!< Central trivariate PDFs
    //! Ordinary moments about which to compute central trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrTriPDF;
};

} // quinoa::

#endif // Statistics_h

