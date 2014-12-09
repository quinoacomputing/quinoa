//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 11:21:42 AM MST
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

    const tk::ParProps& m_particles;              //!< Particle properties

    // Statistics

    //! Instantaneous variable pointers for computing ordinary moments
    std::vector< std::vector< const tk::real* > > m_instOrd;
    std::vector< tk::real > m_ordinary;        //!< Ordinary moments
    std::vector< tk::ctr::FieldVar > m_ordFieldVar;//!< Ordinary moment field names
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
    std::vector< tk::UniPDF > m_ordupdf;           //!< Ordinary univariate PDFs

    //! Instantaneous variable pointers for computing central univariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenUniPDF;
    std::vector< tk::UniPDF > m_cenupdf;           //!< Central univariate PDFs
    //! Ordinary moments about which to compute central univariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrUniPDF;

    // Bivariate probability density functions

    //! Instantaneous variable pointers for computing ordinary bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instOrdBiPDF;
    std::vector< tk::BiPDF > m_ordbpdf;           //!< Ordinary bivariate PDFs

    //! Instantaneous variable pointers for computing central bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenBiPDF;
    std::vector< tk::BiPDF > m_cenbpdf;           //!< Central bivariate PDFs
    //! Ordinary moments about which to compute central bivariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrBiPDF;

    // Trivariate probability density functions

    //! Instantaneous variable pointers for computing ordinary trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instOrdTriPDF;
    std::vector< tk::TriPDF > m_ordtpdf;          //!< Ordinary trivariate PDFs

    //! Instantaneous variable pointers for computing central trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_instCenTriPDF;
    std::vector< tk::TriPDF > m_centpdf;          //!< Central trivariate PDFs
    //! Ordinary moments about which to compute central trivariate PDFs
    std::vector< std::vector< const tk::real* > > m_ctrTriPDF;
};

} // tk::

#endif // Statistics_h

