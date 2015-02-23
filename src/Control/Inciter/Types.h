//******************************************************************************
/*!
  \file      src/Control/Incitier/Types.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:52:44 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Types for Incitier's parsers
  \details   Types for Incitier's parsers. This file defines the components of the
    tagged tuple that stores heteroegeneous objects in a hierarchical way. These
    components are therefore part of the grammar stack that is filled during
    parsing (both command-line argument parsing and control file parsing).
*/
//******************************************************************************
#ifndef IncitierTypes_h
#define IncitierTypes_h

#include <Tags.h>
#include <Types.h>
#include <RNGParam.h>
#include <Options/InitPolicy.h>
#include <Options/CoeffPolicy.h>
#include <Options/PDFFile.h>
#include <Options/PDFPolicy.h>
#include <Options/PDFCentering.h>
#include <Options/TxtFloatFormat.h>
#include <Options/RNG.h>
#include <PUPUtil.h>

namespace inciter {
namespace ctr {

//! Discretization parameters storage
using discretization = tk::tuple::tagged_tuple<
  tag::nstep,     kw::nstep::info::expect::type, //!< Number of time steps
  tag::term,      kw::term::info::expect::type,  //!< Time to terminate
  tag::dt,        kw::dt::info::expect::type     //!< Size of time step
>;

//! Output intervals storage
using intervals = tk::tuple::tagged_tuple<
  tag::tty,  kw::ttyi::info::expect::type,  //!< TTY output interval
  tag::dump, uint32_t,  //!< Dump output interval
  tag::glob, uint32_t   //!< Glob output interval
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,     kw::control::info::expect::type,  //!< Control filename
  tag::input,       std::string,  //!< Input filename
  tag::output,      std::string   //!< Output filename
>;

//! PEGTL location type to use throughout Incitier's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // inciter::

#endif // IncitierTypes_h
