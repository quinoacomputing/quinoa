// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/New2InputDeck.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's new input deck definition
  \details   This file defines the heterogeneous struct that is used for storing
     the data from user input during the control file parsing of the
     computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************
#ifndef New2InputDeck_h
#define New2InputDeck_h

#include <getopt.h>
#include "SimpTaggedTuple.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/Components.hpp"
#include "Inciter/Options/PDE.hpp"
#include "Inciter/Options/Scheme.hpp"
#include "Options/PartitioningAlgorithm.hpp"

namespace newtag {
  constexpr int ic = 0;

  //namespace ic {
    constexpr int materialid = 0;
    constexpr int pressure = 1;
    constexpr int temperature = 2;
    constexpr int velocity = 3;
    constexpr int box = 4;

    //namespace box {
      constexpr int volume = 0;
      constexpr int mass = 1;
      constexpr int density = 2;
      constexpr int energy = 3;
      constexpr int energy_content = 4;
      constexpr int initiate = 5;
      constexpr int point = 6;
      constexpr int init_time = 7;
      constexpr int front_width = 8;
      constexpr int xmin = 9;
      constexpr int xmax = 10;
      constexpr int ymin = 11;
      constexpr int ymax = 12;
      constexpr int zmin = 13;
      constexpr int zmax = 14;
      constexpr int orientation = 15;
    //}
  //}

  constexpr int dt = 1;
} // newtag::

namespace inciter {

namespace ctr {

class NwDeck {

  public:

  std::tuple<

    std::string,  // title

    // time stepping options
    uint64_t,     // nstep
    tk::real,     // term
    tk::real,     // t0
    tk::real,     // dt
    tk::real,     // cfl
    uint32_t,     // ttyi

    // ic block
    std::tuple<
      std::size_t, //materialid
      tk::real, //pressure
      tk::real, //temperature
      tk::real, //velocity
      std::vector< //box
        std::tuple<
          std::size_t, //materialid
          tk::real, //volume
          tk::real, //mass
          tk::real, //density
          tk::real, //velocity
          tk::real, //pressure
          tk::real, //energy
          tk::real, //energy_content
          tk::real, //temperature
          tk::real, //xmin
          tk::real, //xmax
          tk::real, //ymin
          tk::real, //ymax
          tk::real, //zmin
          tk::real, //zmax
          std::vector< tk::real >, //orientation
          InitiateType, //initiate
          std::vector< tk::real >, //point
          tk::real, //init_time
          tk::real //front_width
        >
      >
    >,

    tk::real // dt
  > myDeck;

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | myDeck;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i NwDeck object reference
  friend void operator|( PUP::er& p, NwDeck& i ) { i.pup(p); }
  //@}

};

} // ctr::
} // inciter::

#endif // New2InputDeck_h
