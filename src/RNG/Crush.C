//******************************************************************************
/*!
  \file      src/RNG/Crush.C
  \author    J. Bakosi
  \date      Fri 06 Dec 2013 01:35:28 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Crush battery
  \details   Crush battery
*/
//******************************************************************************

#include <Crush.h>
#include <TestU01.h>

using rngtest::Crush;

Crush::Crush(const Base& base) :
  TestU01Suite( base,
    tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::CRUSH ) )
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  setupRNGs( *this );
}

void
Crush::addTests( const quinoa::ctr::RNGType& rng, const Gen01Ptr& gen )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
  // Marsaglia Serial Over, t = 2
  add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Serial Over t=2"} ),
       SerialOver, 1L, 500L * MILLION, 0, 4096L, 2 );

  // Marsaglia Serial Over, t = 4
  add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Serial Over t=4"} ),
       SerialOver, 1L, 300L * MILLION, 0, 64L, 4 );

  // Marsaglia Collision Over, t = 2, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=2 r=0"} ),
       CollisionOver, 10L, 10L * MILLION, 0, 1024L * 1024, 2 );

  // Marsaglia Collision Over, t = 2, r = 10
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=2 r=10"} ),
       CollisionOver, 10L, 10L * MILLION, 10, 1024L * 1024, 2 );

  // Marsaglia Collision Over, t = 4, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=4 r=0"} ),
       CollisionOver, 10L, 10L * MILLION, 0, 1024L, 4 );

  // Marsaglia Collision Over, t = 4, r = 20
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=4 r=20"} ),
       CollisionOver, 10L, 10L * MILLION, 20, 1024L, 4 );

  // Marsaglia Collision Over, t = 8, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=8 r=0"} ),
       CollisionOver, 10L, 10L * MILLION, 0, 32L, 8 );

  // Marsaglia Collision Over, t = 8, r = 25
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=8 r=25"} ),
       CollisionOver, 10L, 10L * MILLION, 25, 32L, 8 );

  // Marsaglia Collision Over, t = 20, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=20 r=0"} ),
       CollisionOver, 10L, 10L * MILLION, 0, 4L, 20 );

  // Marsaglia Collision Over, t = 20, r = 28
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=20 r=28"} ),
       CollisionOver, 10L, 10L * MILLION, 28, 4L, 20 );

  #ifdef USE_LONG_LONG

  // Marsaglia Birthday Spacings, t = 2, r = 0
  #if LONG_MAX <= 2147483647L
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
       BirthdaySpacings, 10L, 10L * MILLION, 0, 1073741824L, 2, 1 );
  #else // LONG_MAX
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
       BirthdaySpacings, 5L, 20L * MILLION, 0, 2L*1073741824L, 2, 1 );
  #endif // LONG_MAX

  // Marsaglia Birthday Spacings, t = 3, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=3 r=0"} ),
       BirthdaySpacings, 5L, 20L * MILLION, 0, 2097152L, 3, 1 );

  // Marsaglia Birthday Spacings, t = 4, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=0"} ),
       BirthdaySpacings, 5L, 20L * MILLION, 0, 65536L, 4, 1 );

  // Marsaglia Birthday Spacings, t = 7, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=7 r=0"} ),
       BirthdaySpacings, 3L, 20L * MILLION, 0, 512L, 7, 1 );

  // Marsaglia Birthday Spacings, t = 7, r = 7
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=7 r=7"} ),
       BirthdaySpacings, 3L, 20L * MILLION, 7, 512L, 7, 1 );

  // Marsaglia Birthday Spacings, t = 8, r = 14
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=8 r=14"} ),
       BirthdaySpacings, 3L, 20L * MILLION, 14, 256L, 8, 1 );

  // Marsaglia Birthday Spacings, t = 8, r = 22
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=8 r=22"} ),
       BirthdaySpacings, 3L, 20L * MILLION, 22, 256L, 8, 1 );

  #else // USE_LONG_LONG

  // Marsaglia Birthday Spacings, t = 2, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
       BirthdaySpacings, 200L, 4L * MILLION / 10, 0, 67108864L, 2, 1 );

  // Marsaglia Birthday Spacings, t = 3, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=3 r=0"} ),
       BirthdaySpacings, 100L, 4L * MILLION / 10, 0, 131072L, 3, 1 );

  // Marsaglia Birthday Spacings, t = 4, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=0"} ),
       BirthdaySpacings, 200L, 4L * MILLION / 10, 0, 1024L * 8, 4, 1 );

  // Marsaglia Birthday Spacings, t = 13, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=0"} ),
       BirthdaySpacings, 100L, 4L * MILLION / 10, 0, 16L, 13, 1 );

  // Marsaglia Birthday Spacings, t = 13, r = 10
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=10"} ),
       BirthdaySpacings, 100L, 4L * MILLION / 10, 10, 16L, 13, 1 );

  // Marsaglia Birthday Spacings, t = 13, r = 20
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=20"} ),
       BirthdaySpacings, 100L, 4L * MILLION / 10, 20, 16L, 13, 1 );

  // Marsaglia Birthday Spacings, t = 13, r = 26
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=26"} ),
       BirthdaySpacings, 100L, 4L * MILLION / 10, 26, 16L, 13, 1 );

  #endif // USE_LONG_LONG



}
