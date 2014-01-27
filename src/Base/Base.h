//******************************************************************************
/*!
  \file      src/Base/Base.h
  \author    J. Bakosi
  \date      Mon 27 Jan 2014 04:28:17 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Collection of essentials
  \details   Collection of essentials
*/
//******************************************************************************
#ifndef Base_h
#define Base_h

#include <QuinoaPrint.h>
#include <RNGTestPrint.h>
#include <Paradigm.h>
#include <Timer.h>
#include <Driver.h>
#include <LayoutPolicy.h>

namespace quinoa {

//! Base: collection of essentials
struct Base {
  QuinoaPrint& print;
  tk::Paradigm& paradigm;
  ctr::InputDeck& control;
  tk::Timer& timer;
  tk::RNGFactory& rng;

  //! Initializer constructor
  Base( QuinoaPrint& Print,
        tk::Paradigm& Paradigm,
        ctr::InputDeck& Control,
        tk::Timer& Timer,
        tk::RNGFactory& Rng ) : print( Print ),
                                paradigm( Paradigm ),
                                control( Control ),
                                timer( Timer ),
                                rng( Rng ) {}
};

//! Select data layout policy for particle properties
#if   defined LAYOUT_PAR_EQ_COMP
using ParProps = ParticleProperties< ParEqComp >;
#elif defined LAYOUT_PAR_COMP_EQ
using ParProps = ParticleProperties< ParCompEq >;
#elif defined LAYOUT_EQ_COMP_PAR
using ParProps = ParticleProperties< EqCompPar >;
#elif defined LAYOUT_EQ_PAR_COMP
using ParProps = ParticleProperties< EqParComp >;
#elif defined LAYOUT_COMP_EQ_PAR
using ParProps = ParticleProperties< CompEqPar >;
#elif defined LAYOUT_COMP_PAR_EQ
using ParProps = ParticleProperties< CompParEq >;
#endif

} // quinoa::

namespace rngtest {

//! Base: collection of essentials
struct Base {
  RNGTestPrint& print;
  tk::Paradigm& paradigm;
  ctr::InputDeck& control;
  tk::Timer& timer;
  tk::RNGFactory& rng;

  //! Initializer constructor
  Base( RNGTestPrint& Print,
        tk::Paradigm& Paradigm,
        ctr::InputDeck& Control,
        tk::Timer& Timer,
        tk::RNGFactory& Rng ) : print( Print ),
                                paradigm( Paradigm ),
                                control( Control ),
                                timer( Timer ),
                                rng( Rng ) {}
};

} // rngtest::

#endif // Base_h
