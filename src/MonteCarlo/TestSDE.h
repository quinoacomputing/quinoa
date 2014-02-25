//******************************************************************************
/*!
  \file      src/MonteCarlo/TestSDE.h
  \author    J. Bakosi
  \date      Mon 24 Feb 2014 08:46:53 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SDE testbed
  \details   SDE testbed
*/
//******************************************************************************
#ifndef TestSDE_h
#define TestSDE_h

#include <MonteCarlo.h>
#include <Timer.h>

namespace quinoa {

//! TestSDE : MonteCarlo
class TestSDE : public MonteCarlo {

  public:
    //! Constructor
    explicit TestSDE( const Base& base );

    //! Destructor
    ~TestSDE() override = default;

    //! Run
    void run() override;

    //! Factor accessor
    ctr::SDEFactory& factory() noexcept { return m_SDEFactory; }

  private:
    //! Don't permit copy constructor
    TestSDE(const TestSDE&) = delete;
    //! Don't permit copy assigment
    TestSDE& operator=(const TestSDE&) = delete;
    //! Don't permit move constructor
    TestSDE(TestSDE&&) = delete;
    //! Don't permit move assigment
    TestSDE& operator=(TestSDE&&) = delete;

    //! Instantiate an SDE
    template< class SDE, class SDEType >
    void instantiateSDE( const SDEType& sdetype ) {
      if ( control().get< tag::component, SDE >() ) {
        ctr::SDEKey key{ sdetype,
                         control().get< tag::param, SDE, tag::initpolicy >(),
                         control().get< tag::param, SDE, tag::coeffpolicy >() };
        m_sde.push_back( tk::instantiate( m_SDEFactory, key ) );
      }
    }

    //! Advance
    void advance( tk::real dt );

    //! Output joint scalar PDF
    void outJpdf( tk::real t );

    //! Initialize factories
    void initFactories( const QuinoaPrint& print );

    //! Echo information on SDE test bed
    void echo();

    //! Factory
    ctr::SDEFactory m_SDEFactory;                       //!< SDE factory

    //! Pointers to SDEs
    std::vector< std::unique_ptr< Model > > m_sde;
};

} // quinoa::

#endif // TestSDE_h
