//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.C
  \author    J. Bakosi
  \date      Wed 27 Nov 2013 07:59:43 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************

#include <iostream>

extern "C" {
  #include <unif01.h>
  #include <bbattery.h>
  #include <swrite.h>
}

#include <SmallCrush.h>
#include <MKLRNGWrappers.h>
#include <TestU01.h>

namespace rngtest {

extern std::unique_ptr< tk::RNG > g_rng;

} // rngtest::

using rngtest::SmallCrush;

SmallCrush::SmallCrush(const Base& base) : TestU01Suite(base)
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate RNGs
  ctr::RNGType r = m_base.control.get< ctr::selected, ctr::rng >()[0];
  if (r != ctr::RNGType::NO_RNG) {
    // Save name of random number generator
    quinoa::ctr::RNG rng;
    m_rngname = rng.name( r );
    // Instantiate random number generator
    g_rng = std::unique_ptr< tk::RNG >( m_base.rng[r]() );
  }

  // Create TestU01 external generator
  char* const rngname = const_cast<char*>( m_rngname.c_str() );
  m_gen =
    Gen01Ptr( unif01_CreateExternGen01( rngname,
                                        MKLRNGUniform ) );

  // Type specializations
  struct BirthdaySpacings_info {
    static const char* name() { return "Marsaglia's BirthdaySpacings"; }
  };
  using BirthdaySpacingsTest = TestU01< sres_Poisson,
                                        sres_CreatePoisson,
                                        sres_DeletePoisson,
                                        BirthdaySpacings,
                                        BirthdaySpacings_info >;
  struct Collision_info {
    static const char* name() { return "Knuth's Collision"; }
  };
  using CollisionTest = TestU01< sknuth_Res2,
                                 sknuth_CreateRes2,
                                 sknuth_DeleteRes2,
                                 Collision,
                                 Collision_info >;
  struct Gap_info {
    static const char* name() { return "Knuth's Gap"; }
  };
  using GapTest = TestU01< sres_Chi2,
                           sres_CreateChi2,
                           sres_DeleteChi2,
                           Gap,
                           Gap_info >;

  // Add statistical tests to battery
  add< BirthdaySpacingsTest >( m_tests, m_gen );
  add< CollisionTest >( m_tests, m_gen );
  add< GapTest >( m_tests, m_gen );
}

void
SmallCrush::run()
//******************************************************************************
//  Run SmallCrush battery
//! \author  J. Bakosi
//******************************************************************************
{
  swrite_Basic = FALSE; // no putput from TestU01

  //bbattery_SmallCrush( m_gen.get() );

  for (size_t i=0; i<m_tests.size(); ++i) {
    m_tests[i].pval = m_tests[i].ptr->run();
    std::cout << m_tests[i].ptr->name() << " pval = " << m_tests[i].pval << std::endl;
  }

//   for (auto& test : m_tests) {
//     test.pval = test.ptr->run();
//     std::cout << test.ptr->name() << " pval = " << test.pval << std::endl;
//   }
}
