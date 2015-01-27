//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Wed 21 Jan 2015 03:40:09 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Statistics
  \details   Computing ordinary and central moments
*/
//******************************************************************************
#ifndef Statistics_h
#define Statistics_h

#include <Types.h>
#include <ParticleProperties.h>
#include <UniPDF.h>
#include <BiPDF.h>
#include <TriPDF.h>
#include <Components.h>

namespace tk {

//! Statistics estimator
class Statistics {

  public:
    //! Constructor
    explicit Statistics( const tk::ParProps& particles,
                         const ctr::OffsetMap& offset,
                         const std::vector< ctr::Product >& stat,
                         const std::vector< ctr::Probability >& pdf,
                         const std::vector< std::vector< tk::real > >& binsize );

    //! Accumulate (i.e., only do the sum for) ordinary moments
    void accumulateOrd();

    //! Accumulate (i.e., only do the sum for) central moments
    void accumulateCen( const std::vector< tk::real >& ord );

    //! Accumulate (i.e., only do the sum for) ordinary PDFs
    void accumulateOrdPDF();

    //! Accumulate (i.e., only do the sum for) central PDFs
    void accumulateCenPDF( const std::vector< tk::real >& ord );

    //! Ordinary moments accessor
    const std::vector< tk::real >& ord() const noexcept { return m_ordinary; }

    //! Central moments accessor
    const std::vector< tk::real >& ctr() const noexcept { return m_central; }

    //! Ordinary univariate PDFs accessor
    const std::vector< tk::UniPDF >& oupdf() const noexcept { return m_ordupdf; }

    //! Ordinary bivariate PDFs accessor
    const std::vector< tk::BiPDF >& obpdf() const noexcept { return m_ordbpdf; }

    //! Ordinary trivariate PDFs accessor
    const std::vector< tk::TriPDF >& otpdf() const noexcept { return m_ordtpdf; }

    //! Central univariate PDFs accessor
    const std::vector< tk::UniPDF >& cupdf() const noexcept { return m_cenupdf; }

    //! Central bivariate PDFs accessor
    const std::vector< tk::BiPDF >& cbpdf() const noexcept { return m_cenbpdf; }

    //! Central trivariate PDFs accessor
    const std::vector< tk::TriPDF >& ctpdf() const noexcept { return m_centpdf; }

  private:
    //! Setup ordinary moments
    void setupOrdinary( const ctr::OffsetMap& offset,
                        const std::vector< ctr::Product >& stat );

    //! Setup central moments
    void setupCentral( const ctr::OffsetMap& offset,
                       const std::vector< ctr::Product >& stat );

    //! Setup PDFs
    void setupPDF( const ctr::OffsetMap& offset,
                   const std::vector< ctr::Probability >& pdf,
                   const std::vector< std::vector< tk::real > >& binsize  );

    //! Return mean for fluctuation
    std::size_t mean(const tk::ctr::Term& term) const;

    //! Particle properties
    const tk::ParProps& m_particles;

    //! Map used to lookup ordinary moments
    std::map< tk::ctr::Product, const tk::real* > m_ordLookup;
    //! Map used to lookup central moments
    std::map< tk::ctr::Product, const tk::real* > m_cenLookup;

    // Statistics

    //! Instantaneous variable pointers for computing ordinary moments
    std::vector< std::vector< const tk::real* > > m_instOrd;
    //! Ordinary moments
    std::vector< tk::real > m_ordinary;
    //! Ordinary moment Terms, used to find means for fluctuations
    std::vector< tk::ctr::Term > m_ordTerm;
    //! Number of ordinary moments
    std::size_t m_nord;

    //! Instantaneous variable pointers for computing central moments
    std::vector< std::vector< const tk::real* > > m_instCen;
    //! Central moments
    std::vector< tk::real > m_central;
    //! Ordinary moments about which to compute central moments
    std::vector< std::vector< const tk::real* > > m_ctr;
    //! Number of central moments
    std::size_t m_ncen;

    // Univariate probability density functions

    //! Instantaneous variable pointers for computing ordinary univariate PDFs
    std::vector< std::vector< const tk::real* > > m_instOrdUniPDF;
    //! Ordinary univariate PDFs
    std::vector< tk::UniPDF > m_ordupdf;

    //! Instantaneous variable pointers for computing central univariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenUniPDF;
    //! Central univariate PDFs
    std::vector< tk::UniPDF > m_cenupdf;
    //! Ordinary moments about which to compute central univariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrUniPDF;

    // Bivariate probability density functions

    //! Instantaneous variable pointers for computing ordinary bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instOrdBiPDF;
    //! Ordinary bivariate PDFs
    std::vector< tk::BiPDF > m_ordbpdf;

    //! Instantaneous variable pointers for computing central bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenBiPDF;
    //! Central bivariate PDFs
    std::vector< tk::BiPDF > m_cenbpdf;
    //! Ordinary moments about which to compute central bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrBiPDF;

    // Trivariate probability density functions

    //! Instantaneous variable pointers for computing ordinary trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instOrdTriPDF;
    //! Ordinary trivariate PDFs
    std::vector< tk::TriPDF > m_ordtpdf;

    //! Instantaneous variable pointers for computing central trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenTriPDF;
    //! Central trivariate PDFs
    std::vector< tk::TriPDF > m_centpdf;
    //! Ordinary moments about which to compute central trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrTriPDF;
};

} // tk::

#endif // Statistics_h

