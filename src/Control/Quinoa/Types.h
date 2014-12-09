//******************************************************************************
/*!
  \file      src/Control/Quinoa/Types.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 09:34:51 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Types for Quinoa's parsers
  \details   Types for Quinoa's parsers
*/
//******************************************************************************
#ifndef QuinoaTypes_h
#define QuinoaTypes_h

#include <Tags.h>
#include <Types.h>
#include <ControlTypes.h>
#include <Options/InitPolicy.h>
#include <Options/CoeffPolicy.h>
#include <Options/PDFFile.h>
#include <Options/PDFPolicy.h>
#include <Options/PDFCentering.h>
#include <Options/TxtFloatFormat.h>
#include <Options/RNG.h>
#include <Quinoa/Options/MonteCarlo.h>
#include <Quinoa/Options/Position.h>
#include <Quinoa/Options/Mass.h>
#include <Quinoa/Options/Hydro.h>
#include <Quinoa/Options/Energy.h>
#include <Quinoa/Options/Mix.h>
#include <Quinoa/Options/Frequency.h>
#include <Quinoa/Options/MixRate.h>
#include <PUPUtil.h>

namespace quinoa {
namespace ctr {

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::montecarlo,   MonteCarloType,   //!< Physics
  tag::position,     PositionType,     //!< Position model
  tag::mass,         MassType,         //!< Mass model
  tag::hydro,        HydroType,        //!< Hydrodynamics model
  tag::energy,       EnergyType,       //!< Internal energy model
  tag::mix,          MixType,          //!< Material mix model
  tag::frequency,    FrequencyType,    //!< Turbulence frequency model
  tag::mixrate,      MixRateType,      //!< Material mix rate model
  tag::rng,          std::vector< tk::ctr::RNGType >, //!< RNGs
  tag::pdffiletype,  tk::ctr::PDFFileType,      //!< PDF output file type
  tag::pdfpolicy,    tk::ctr::PDFPolicyType,    //!< PDF output file policy
  tag::pdfctr,       tk::ctr::PDFCenteringType, //!< PDF output file centering
  tag::float_format, tk::ctr::TxtFloatFormatType//!< Text floating-point format
>;

//! Discretization parameters storage
using discretization = tk::tuple::tagged_tuple<
  tag::npar,      uint64_t,  //!< Total number of particles
  tag::nstep,     uint64_t,  //!< Number of time steps to take
  tag::term,      tk::real,  //!< Time to terminate time stepping
  tag::dt,        tk::real,  //!< Size of time step
  tag::binsize,   std::vector< std::vector< tk::real > >, //!< PDF binsizes
  tag::extent,    std::vector< std::vector< tk::real > >, //!< PDF extents
  tag::precision, std::streamsize  //!< Precision in digits
>;

//! Output intervals storage
using intervals = tk::tuple::tagged_tuple<
  tag::tty,  uint32_t,  //!< TTY output interval
  tag::dump, uint32_t,  //!< Dump output interval
  tag::stat, uint32_t,  //!< Statistics output interval
  tag::pdf,  uint32_t,  //!< PDF output interval
  tag::glob, uint32_t   //!< Glob output interval
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,     std::string,  //!< Control filename
  tag::input,       std::string,  //!< Input filename
  tag::output,      std::string,  //!< Output filename
  tag::pdf,         std::string,  //!< PDF filename
  tag::glob,        std::string,  //!< Glob filename
  tag::stat,        std::string,  //!< Statistics filename
  tag::pdfnames,    std::vector< std::string >  //!< PDF identifiers
>;

//! Position parameters storage
using PositionParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Mass parameters storage
using MassParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Hydro parameters storage
using HydroParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Mix parameters storage
using MixParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Frequency parameters storage
using FrequencyParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Simplified Langevin hydro model parameters storage
using SLMParameters = tk::tuple::tagged_tuple<
  tag::c0, tk::real
>;

//! Generalized Langevin hydro model parameters storage
using GLMParameters = tk::tuple::tagged_tuple<
  tag::c0, tk::real
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  #ifdef HAS_MKL
  tag::rngmkl,       tk::ctr::RNGMKLParameters,   //!< MKL RNG parameters
  #endif
  tag::rngsse,       tk::ctr::RNGSSEParameters,   //!< RNGSSE RNG parameters
  tag::position,     PositionParameters,
  tag::mass,         MassParameters,
  tag::hydro,        HydroParameters,
  tag::mix,          MixParameters,
  tag::frequency,    FrequencyParameters,
  tag::slm,          SLMParameters,
  tag::glm,          GLMParameters
>;

//! PEGTL location type to use throughout Quinoa's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // quinoa::

#endif // QuinoaTypes_h
