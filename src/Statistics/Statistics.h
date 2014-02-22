//******************************************************************************
/*!
  \file      src/Statistics/Statistics.h
  \author    J. Bakosi
  \date      Sat 22 Feb 2014 07:11:29 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics
  \details   Statistics
*/
//******************************************************************************
#ifndef Statistics_h
#define Statistics_h

#include <Types.h>
#include <Base.h>
#include <Distribution.h>
#include <Quinoa/Types.h>

namespace quinoa {

//! Statistics estimator
class Statistics {

  public:
    //! Constructor
    explicit Statistics( const Base& base, const ParProps& particles );

    //! Destructor
    virtual ~Statistics() = default;

    //! Accumulate statistics
    void accumulate();

    //! Number of ordinary moments accessor
    int nord() const noexcept { return m_nord; }

    //! Number of central moments accessor
    int ncen() const noexcept { return m_ncen; }

    //! Ordinary moments accessor
    const tk::real* ordinary() const noexcept { return m_ordinary.get(); }

    //! Central moments accessor
    const tk::real* central() const noexcept { return m_central.get(); }

    //! Find out whether ordinary moment is to be plotted
    bool plotOrdinary( int m ) const;

    //! Return the name of ordinary moment
    const std::string& nameOrdinary( int m ) const;

    //! Return the name of central moment
    const std::string& nameCentral( int m ) const;

  private:
    //! Don't permit copy constructor
    Statistics(const Statistics&) = delete;
    //! Don't permit copy assigment
    Statistics& operator=(const Statistics&) = delete;
    //! Don't permit move constructor
    Statistics(Statistics&&) = delete;
    //! Don't permit move assigment
    Statistics& operator=(Statistics&&) = delete;

    // Get types out of InputDeck needed for computing map< depvar, offset >
    using ncomps = ctr::InputDeck::nT< tag::component >;
    using ncomp = ncomps::ncomp;
    using tags = ncomps::tags;

    // Case-insensitive character comparison for offset map
    // The reason we need a case-insensitive character less-than operator is
    // that the keys (dependent variables) in the offset map are only used to
    // indicated the equation, however, queries (using find) are fired up for
    // both ordinary and central moments (which are upper and lower case) for
    // which the offset (for the same depvar) should be the same.
    struct CaseInsensitiveCharLessThan {
      bool operator() ( char lhs, char rhs ) const {
        return tolower(lhs) < tolower(rhs);
      }
    };

    using OffsetMap = std::map< char, ncomp, CaseInsensitiveCharLessThan >;

    //! Function object for creating a run-time std::map< depvar, offset >
    struct depvar {
      const Base& m_base;
      OffsetMap& m_map;
      depvar( const Base& base, OffsetMap& map ) : m_base(base), m_map(map) {}
      template< typename U > void operator()( U ) {
        m_map[ m_base.control.get< tag::param, U, tag::depvar >() ] =
          m_base.control.get< tag::component >().template offset< U >();
      }
    };

    //! Setup ordinary moments
    void setupOrdinary( const Base& base, const OffsetMap& offset );

    //! Setup central moments
    void setupCentral( const Base& base, const OffsetMap& offset );

    //! Estimate ordinary moments
    void estimateOrdinary();

    //! Estimate central moments
    void estimateCentral();

    //! Find out whether product only contains ordinary moment terms
    bool ordinary( const std::vector< ctr::Term >& product ) const;

    //! Return mean for fluctuation
    int mean(const ctr::Term& term) const;

    //! Convert string to upper case
    std::string toUpper( const std::string& s ) const;

    //! Return true if string is all lower case
    bool isLower( const std::string&s ) const;

    const ParProps& m_particles;              //!< Particle properties
    const uint64_t m_nthreads;                //!< Number of threads
    const uint64_t m_npar;                    //!< Number of particles

    //! Instantaneous variable pointers for computing ordinary moments
    std::vector< std::vector< const tk::real* > > m_instOrd;
    std::unique_ptr< tk::real[] > m_ordinary; //!< Ordinary moments
    std::vector< bool > m_plotOrdinary;       //!< Whether to plot ord moments
    //! Ordinary moment field names
    std::vector< ctr::FieldVar > m_ordFieldVar;
    std::vector< std::string > m_nameOrdinary;//!< Ordinary moment names
    int m_nord;                               //!< Number of ordinary moments

    //! Instantaneous variable pointers for computing central moments
    std::vector< std::vector< const tk::real* > > m_instCen;
    std::unique_ptr< tk::real[] > m_central;  //!< Central moments
    //! Ordinary moments about which to compute central moments
    std::vector< std::vector< const tk::real* > > m_center;
    std::vector< std::string > m_nameCentral; //!< Central moment names
    int m_ncen;                               //!< Number of central moments
};

} // quinoa::

#endif // Statistics_h
