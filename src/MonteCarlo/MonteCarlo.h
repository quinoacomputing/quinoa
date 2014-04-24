//******************************************************************************
/*!
  \file      src/MonteCarlo/MonteCarlo.h
  \author    J. Bakosi
  \date      Thu Apr 24 11:12:44 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Monte Carlo
  \details   Monte Carlo
*/
//******************************************************************************
#ifndef MonteCarlo_h
#define MonteCarlo_h

#include <boost/mpl/at.hpp>

#include <Base.h>
#include <Statistics.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>
#include <Factory.h>

namespace quinoa {

//! MonteCarlo

class MonteCarlo {

  public:
    //! Constructor
    explicit MonteCarlo( const Base& b ) :
      m_base( b ),
      m_npar( b.control.get< tag::incpar, tag::npar >() ),
      m_term( b.control.get< tag::incpar, tag::term >() ),
      m_totalTime( b.timer.create("Total solution") ),
      m_particles( m_npar, b.control.get< tag::component >().nprop() ),
      m_statistics( b, m_particles ),
      m_glob( b.control.get< tag::cmd, tag::io, tag::glob >() ),
      m_stat( b.control.get< tag::cmd, tag::io, tag::stat >(), m_statistics )
    {}

    //! Destructor
    virtual ~MonteCarlo() = default;

    //! Run
    virtual void run() = 0;

    //! Constant accessor to base
    //! \return Pointer to base
    const Base& base() const noexcept { return m_base; }

    //! Constant accessor to control object
    //! \return Control object
    const ctr::InputDeck& control() const noexcept { return m_base.control; }

    //! Accessor to particle properties pointer
    //! \return Particle properties array
    const ParProps& particles() const noexcept { return m_particles; }

  protected:

    //! Constant accessor to print object
    //! \return Print object
    const QuinoaPrint& print() const noexcept { return m_base.print; }

    //! Constant accessor to timer object pointer
    //! \return Pointer to timer object
    const tk::Timer& timer() const noexcept { return m_base.timer; }

    //! Accessor to statistics estimator
    //! \return Pointer to statistics estimator
    Statistics& statistics() noexcept { return m_statistics; }

    //! Accessor to glob file writer
    //! \return Pointer to glob file writer
    GlobWriter& globWriter() noexcept { return m_glob; }

    //! Accessor to statistics file writer
    //! \return Pointer to statistics file writer
    TxtStatWriter& statWriter() noexcept { return m_stat; }

    //! Accessor to max run time
    //! \return Max run time
    const tk::real& term() const noexcept { return m_term; }

    //! Standard header at start of time stepping that children can override
    virtual void header() const;

    //! Standard one-liner report that children can override
    virtual void report( uint64_t it,
                         uint64_t nstep,
                         tk::real t,
                         tk::real dt,
                         bool wroteJpdf,
                         bool wroteGlob,
                         bool wrotePlot );

    //! Echo information on RNGs that children can override
    virtual void echoRNGs() const;

    //! Echo information on IO filenames that children can override
    virtual void echoIO() const;

    //! Echo information on increment parameters that children can override
    virtual void echoIncpar() const;

    //! Echo information on intervals that children can override
    virtual void echoIntervals() const;

    //! Echo information on statistics that children can override
    virtual void echoStatistics() const;

    //! Function object for register an SDE into an SDE factory - repeatedly
    //! called by mpl::cartesian_product sweeping all combinations of the SDE
    //! policies
    template< class Host,
              template<class,class> class SDE,
              class ncomp,
              class SDEType >
    struct registerSDE {

      Host* const m_host;
      const SDEType m_type;

      registerSDE( Host* const host, SDEType type ) :
        m_host( host ), m_type( type ) {}

      template< typename U > void operator()( U ) {
        namespace mpl = boost::mpl;

        // Get Initialization policy: 1st type of mpl::vector U
        using InitPolicy = typename mpl::at< U, mpl::int_<0> >::type;
        // Get coefficients policy: 2nd type of mpl::vector U
        using CoeffPolicy = typename mpl::at< U, mpl::int_<1> >::type;

        // Build SDE key
        ctr::SDEKey key;
        key.get< tag::sde >() = m_type;
        key.get< tag::initpolicy >() = InitPolicy().type();
        key.get< tag::coeffpolicy >() = CoeffPolicy().type();

        // Register SDE (with policies given by mpl::vector U) into SDE factory
        const auto& comp = m_host->control().template get< tag::component >();
        tk::record< SDE< InitPolicy, CoeffPolicy > >
                  ( m_host->factory(), key,
                    m_host->base(),
                    std::cref( m_host->particles() ),
                    comp.template offset< ncomp >(),
                    comp.template get< ncomp >() );
      }
    };

    const Base& m_base;                             //!< Essentials
    const uint64_t m_npar;                          //!< Number of particles
    const tk::real m_term;                          //!< Maximum run time
    const tk::TimerId m_totalTime;                  //!< Timer for total run
    const ParProps m_particles;                     //!< Particle properties

  private:
    //! Don't permit copy constructor
    MonteCarlo(const MonteCarlo&) = delete;
    //! Don't permit copy assigment
    MonteCarlo& operator=(const MonteCarlo&) = delete;
    //! Don't permit move constructor
    MonteCarlo(MonteCarlo&&) = delete;
    //! Don't permit move assigment
    MonteCarlo& operator=(MonteCarlo&&) = delete;

    Statistics m_statistics;                        //!< Statistics estimator
    GlobWriter m_glob;                              //!< Glob file writer
    TxtStatWriter m_stat;                           //!< Statistics file writer
};

} // quinoa::

#endif // MonteCarlo_h
