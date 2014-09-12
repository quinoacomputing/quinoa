//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Fri 12 Sep 2014 06:42:15 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Statistics
  \details   Computing ordinary and central moments
*/
//******************************************************************************
#ifndef Statistics_h
#define Statistics_h

#include <Types.h>
#include <Quinoa/InputDeck/InputDeck.h>
#include <ParticleProperties.h>
#include <PDF.h>

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

    //! Accumulate (i.e., only do the sum for) PDFs
    void accumulatePDF();

    //! Ordinary moments accessor
    const std::vector< tk::real >& ord() const noexcept { return m_ordinary; }

    //! Central moments accessor
    const std::vector< tk::real >& ctr() const noexcept { return m_central; }

    //! PDF accessor
    const std::vector< PDF >& pdf() const noexcept { return m_pdf; }

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
    int mean(const ctr::Term& term) const;

    const ParProps& m_particles;              //!< Particle properties

    //! Instantaneous variable pointers for computing ordinary moments
    std::vector< std::vector< const tk::real* > > m_instOrd;
    std::vector< tk::real > m_ordinary;        //!< Ordinary moments
    std::vector< ctr::FieldVar > m_ordFieldVar;//!< Ordinary moment field names
    int m_nord;                                //!< Number of ordinary moments

    //! Instantaneous variable pointers for computing central moments
    std::vector< std::vector< const tk::real* > > m_instCen;
    std::vector< tk::real > m_central;         //!< Central moments
    //! Ordinary moments about which to compute central moments
    std::vector< std::vector< const tk::real* > > m_center;
    int m_ncen;                                //!< Number of central moments

    //! Instantaneous variable pointers for computing PDFs
    std::vector< std::vector< const tk::real* > > m_instPDF;
    std::vector< PDF > m_pdf;                  //!< PDFs
};

} // quinoa::

#endif // Statistics_h

