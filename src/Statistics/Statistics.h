// *****************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Statistics class declaration
  \details   This file implements a statistics class that can be used to
    estimate statistics from an ensemble. Supported at this time are ordinary
    and central statistical moments of arbitrary-length products and arbitrary
    number of 1D, 2D, and 3D probability density functions (PDF) with sample
    spaces of ordinary and/or central sample space variables.

    _Definitions and nomenclature:_

    - Upper-case letters denote a full random variable, e.g., X

    - Lower-case letters denote a fluctuation about the mean, i.e., x = X - <X>

    - Letters can be augmented by a field ID, i.e., X2 is the full variable of
      the second component of the vector X, while x1 = X1 - <X1> is the
      fluctuation about the mean of the first component of vector X.

    - If the field ID is unspecified, it defaults to the first field, i.e.,
      X = X1, x = x1, etc.

    - Statistical moments of arbitrary-length products can be computed.

      Examples:
      - <X> - mean,
      - <xx> - variance,
      - <xxx> - third central moment,
      - <xy> - covariance of X and Y,
      - <x1y2> - covariance of the first component of vector X and the second
        component of vector Y

    - In general, arbitrary-length products can be estimated that make up a
      statistical moment, using any number and combinations of upper and
      lower-case letters and their field IDs < [A-Za-z][1-9] ... >.

    - A statistical moment is ordinary if and only if all of its terms are
      ordinary. A central moment has at least one term that is central, i.e., a
      fluctuation about its mean.

      - Examples of ordinary moments: <X>, <XX>, <XYZ>, etc.
      - Examples of central moments: <x1x2>, <Xy>, <XYz>, etc.

    - Estimation of the PDFs can be done using either ordinary or central sample
      space variables.

      Examples:
      - p(X) denotes the univariate PDF of the full variable X,
      - f(x1,x2) denotes the bivariate joint PDF of the fluctuations of the
        variables x1 and x2 about their respective means,
      - g(X,y,Z2) denotes the trivariate joint PDF of variables X, y = Y - <Y>,
        and Z2.
*/
// *****************************************************************************
#ifndef Statistics_h
#define Statistics_h

#include <vector>
#include <cstddef>

#include "Types.h"
#include "StatCtr.h"
#include "Particles.h"
#include "SystemComponents.h"
#include "UniPDF.h"
#include "BiPDF.h"
#include "TriPDF.h"

namespace tk {

//! Statistics estimator
class Statistics {

  public:
    //! Constructor
    explicit Statistics( const tk::Particles& particles,
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
    /** @name Setup functions, called from the constructor */
    ///@{
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
    ///@}

    //! Return mean for fluctuation
    std::size_t mean(const tk::ctr::Term& term) const;

    //! Particle properties
    const tk::Particles& m_particles;

    /** @name Data for statistical moment estimation */
    ///@{
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
    ///@}

    /** @name Data for univariate probability density function estimation */
    ///@{
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
    ///@}

    /** @name Data for bivariate probability density function estimation */
    ///@{
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
    ///@}

    /** @name Data for trivariate probability density function estimation */
    ///@{
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
    ///@}
};

} // tk::

#endif // Statistics_h

