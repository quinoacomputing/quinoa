// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/NewInputDeck.hpp
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
#ifndef NewInputDeck_h
#define NewInputDeck_h

#include <limits>
#include <iomanip>
#include <iostream>

#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/Components.hpp"
#include "Inciter/Options/PDE.hpp"
#include "Inciter/Options/Scheme.hpp"
#include "Options/PartitioningAlgorithm.hpp"

namespace inciter {

namespace ctr {

//! Struct for storing info
struct info_t {
  std::string keyword;
  std::string shortDescription;
  std::string longDescription;
  std::string typeDescription;

  // constructor
  info_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    keyword(kw),
    shortDescription(sDescr),
    longDescription(lDescr),
    typeDescription(tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | keyword;
    p | shortDescription;
    p | longDescription;
    p | typeDescription;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i info_t object reference
  friend void operator|( PUP::er& p, info_t& i ) { i.pup(p); }
  //@}
};

//! Struct for storing entries in the input deck
template< typename T >
struct entry_t {
  info_t info;
  T data;

  // constructor
  entry_t< T >(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | data;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i entry_t object reference
  friend void operator|( PUP::er& p, entry_t& i ) { i.pup(p); }
  //@}
};

//! Struct for storing options in the input deck
struct option_t {
  info_t info;
  std::string data;

  // constructor
  option_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr,
    std::string d="" ) :
    info(kw, sDescr, lDescr, tDescr),
    data(d)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | data;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i option_t object reference
  friend void operator|( PUP::er& p, option_t& i ) { i.pup(p); }
  //@}
};

//! Struct for storing vector entries
template< typename T >
struct vecentry_t {
  info_t info;
  std::vector< T > data;

  // constructor
  vecentry_t< T >(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr,
    std::vector< T > d={} ) :
    info(kw, sDescr, lDescr, tDescr),
    data(d)
    {}

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | info;
      p | data;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i vecentry_t object reference
    friend void operator|( PUP::er& p, vecentry_t& i ) { i.pup(p); }
    //@}
};

// -----------------------------------------------------------------------------
//! Structs for blocks in inputdeck
// -----------------------------------------------------------------------------

//! Transport PDE block
// -----------------------------------------------------------------------------
struct transport_t {
  info_t info;

  // entries inside the block
  entry_t< std::size_t > ncomp {"ncomp",
    "Set number of scalar components for a system of transport equations",
    R"(This keyword is used to specify the number of scalar
    components of transport (linear advection) equations.)", "uint"};

  entry_t< ProblemType > problem {"problem",
    "Specify problem configuration for partial differential equation solver",
    R"(This keyword is used to specify the problem configuration for the
    partial differential equation solver.)", "string"};

  // constructor
  transport_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | ncomp;
    p | problem;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i transport_t object reference
  friend void operator|( PUP::er& p, transport_t& i ) { i.pup(p); }
  //@}
};

//! CompFlow PDE block
// -----------------------------------------------------------------------------
struct compflow_t {
  info_t info;

  // entries inside the block
  entry_t< ProblemType > problem {"problem",
    "Specify problem configuration for partial differential equation solver",
    R"(This keyword is used to specify the problem configuration for the
    partial differential equation solver.)", "string"};

  // constructor
  compflow_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | problem;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i compflow_t object reference
  friend void operator|( PUP::er& p, compflow_t& i ) { i.pup(p); }
  //@}
};

//! MultiMat PDE block
// -----------------------------------------------------------------------------
struct multimat_t {
  info_t info;

  // entries inside the block
  entry_t< std::size_t > nmat {"nmat",
    "Set number of materials for the multi-material system",
    R"(This keyword is used to specify the number of materials for
    multi-material flow, see also the keyword 'multimat'.)", "uint"};

  entry_t< uint64_t > prelax {"prelax",
    "Toggle multi-material finite pressure relaxation",
    R"(This keyword is used to turn finite pressure relaxation between
    multiple materials on/off. It is used only for the multi-material solver,
    and has no effect when used for the other PDE types.)", "uint 0/1"};

  entry_t< tk::real > prelax_timescale {"prelax_timescale",
    "Time-scale for multi-material finite pressure relaxation",
    R"(This keyword is used to specify the time-scale at which finite pressure
    relaxation between multiple materials occurs. The default value of 0.25
    corresponds to a relaxation time that is 4 times the time required for a
    sound wave to pass through a computational element. It is used only for
    multimat, and has no effect for the other PDE types.)", "real"};

  entry_t< int > intsharp {"intsharp",
    "Toggle multi-material interface sharpening",
    R"(This keyword is used to turn interface sharpening on/off. It uses the
    multi-material THINC interface reconstruction.
    Ref. Pandare A. K., Waltz J., & Bakosi J. (2021) Multi-Material
    Hydrodynamics with Algebraic Sharp Interface Capturing. Computers &
    Fluids, doi: https://doi.org/10.1016/j.compfluid.2020.104804. It is used
    for the multi-material and the transport solver, and has no effect when
    used for the other PDE types.)", "uint 0/1"};

  entry_t< tk::real > intsharp_param {"intsharp_param",
    "Parameter for multi-material interface sharpening",
    R"(This keyword is used to specify the parameter for the interface
    sharpening. This parameter affects how many cells the material interfaces
    span, after the use of sharpening. It is used for multimat and transport,
    and has no effect for the other PDE types.)", "real" };

  entry_t< ProblemType > problem {"problem",
    "Specify problem configuration for partial differential equation solver",
    R"(This keyword is used to specify the problem configuration for the
    partial differential equation solver.)", "string"};

  // constructor
  multimat_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | nmat;
    p | prelax;
    p | prelax_timescale;
    p | intsharp;
    p | intsharp_param;
    p | problem;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i multimat_t object reference
  friend void operator|( PUP::er& p, multimat_t& i ) { i.pup(p); }
  //@}
};

//! Material block
// -----------------------------------------------------------------------------
struct material_t {
  info_t info{"material",
    "Start configuration block for material (eos) properties",
    R"(This keyword is used to introduce a material block, used to
    specify material properties.)", "block-title"};

  // entries inside the block
  vecentry_t< uint64_t > id {"id", "ID",
    R"(This keyword is used to specify an ID, a positive integer. Inside the
    material block, it is used to specify a block consisting of IDs of
    materials of that EOS type)", "uint(s)"};

  entry_t< MaterialType > eos {"eos", "Select equation of state (type)",
    R"(This keyword is used to select an equation of state for a material.)",
    "string"};

  vecentry_t< tk::real > gamma {"gamma", "ratio of specific heats",
    R"(This keyword is used to specify the material property, ratio of
    specific heats.)", "real"};

  vecentry_t< tk::real > pstiff {"pstiff", "EoS stiffness parameter",
    R"(This keyword is used to specify the material property, stiffness
    parameter in the stiffened gas equation of state.)", "real"};

  vecentry_t< tk::real > w_gru {"w_gru", "Grueneisen coefficient",
    R"(This keyword is used to specify the material property, Gruneisen
    coefficient for the Jones-Wilkins-Lee equation of state.)", "real"};

  vecentry_t< tk::real > A_jwl {"A_jwl", "JWL EoS A parameter",
    R"(This keyword is used to specify the material property A (units: Pa) for
    the Jones-Wilkins-Lee equation of state.)", "real"};

  vecentry_t< tk::real > B_jwl {"B_jwl", "JWL EoS B parameter",
    R"(This keyword is used to specify the material property B (units: Pa) for
    the Jones-Wilkins-Lee equation of state.)", "real"};

  vecentry_t< tk::real > C_jwl {"C_jwl", "JWL EoS C parameter",
    R"(This keyword is used to specify the material property C (units: Pa) for
    the Jones-Wilkins-Lee equation of state.)", "real"};

  vecentry_t< tk::real > R1_jwl {"R1_jwl", "JWL EoS R1 parameter",
    R"(This keyword is used to specify the material property R1 for the
    Jones-Wilkins-Lee equation of state.)", "real"};

  vecentry_t< tk::real > R2_jwl {"R2_jwl", "JWL EoS R2 parameter",
    R"(This keyword is used to specify the material property R2 for the
    Jones-Wilkins-Lee equation of state.)", "real"};

  vecentry_t< tk::real > rho0_jwl {"rho0_jwl", "JWL EoS rho0 parameter",
    R"(This keyword is used to specify the material property rho0, which is the
    density of initial state (units: kg/m3) for the Jones-Wilkins-Lee
    equation of state.)", "real"};

  vecentry_t< tk::real > de_jwl {"de_jwl", "JWL EoS de parameter",
    R"(This keyword is used to specify the material property de, which is the
    heat of detonation for products; and for reactants, it is chosen such that
    the ambient internal energy (e0) is 0 (units: J/kg). Used for the
    Jones-Wilkins-Lee equation of state.)", "real"};

  vecentry_t< tk::real > rhor_jwl {"rhor_jwl", "JWL EoS rhor parameter",
    R"(This keyword is used to specify the material property rhor, which is the
    density of reference state (units: kg/m3) for the Jones-Wilkins-Lee
    equation of state.)", "real"};

  vecentry_t< tk::real > Tr_jwl {"Tr_jwl", "JWL EoS Tr parameter",
    R"(This keyword is used to specify the material property Tr, which is the
    temperature of reference state (units: K) for the Jones-Wilkins-Lee
    equation of state.)", "real"};

  vecentry_t< tk::real > Pr_jwl {"Pr_jwl", "JWL EoS er parameter",
    R"(This keyword is used to specify the material property Pr, which is the
    pressure at the reference state (units: Pa) for the Jones-Wilkins-Lee
    equation of state. It is used to calculate the reference temperature for
    the EoS.)", "real"};

  vecentry_t< tk::real > mu {"mu", "shear modulus/dynamic viscosity",
    R"(This keyword is used to specify the material property, shear modulus
    for solids, or dynamic viscosity for fluids.)", "real"};

  vecentry_t< tk::real > cv {"cv", "specific heat at constant volume",
    R"(This keyword is used to specify the material property, specific heat at
    constant volume.)", "real"};

  vecentry_t< tk::real > k {"k", "heat conductivity",
    R"(This keyword is used to specify the material property, heat
    conductivity.)", "real"};

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | id;
    p | eos;
    p | gamma;
    p | pstiff;
    p | w_gru;
    p | A_jwl;
    p | B_jwl;
    p | C_jwl;
    p | R1_jwl;
    p | R2_jwl;
    p | rho0_jwl;
    p | de_jwl;
    p | rhor_jwl;
    p | Tr_jwl;
    p | Pr_jwl;
    p | mu;
    p | cv;
    p | k;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i material_t object reference
  friend void operator|( PUP::er& p, material_t& i ) { i.pup(p); }
  //@}
};

//! Material index table
// -----------------------------------------------------------------------------
struct matidxmap_t {
  info_t info;

  std::vector< std::size_t > eosidx;
  std::vector< std::size_t > matidx;
  std::vector< std::size_t > solidx;

  // constructor
  matidxmap_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | eosidx;
    p | matidx;
    p | solidx;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i matidxmap_t object reference
  friend void operator|( PUP::er& p, matidxmap_t& i ) { i.pup(p); }
  //@}
};

//! Field output block
// -----------------------------------------------------------------------------
struct fieldoutput_t {
  info_t info;

  // entries inside the block
  entry_t< uint32_t > iter_interval {"interval",
    "Set interval (in units of iteration count)",
    R"(This keyword is used to specify an interval in units of iteration count
    (i.e., number of time steps). This must be used within a relevant
    block.)", "uint"};

  entry_t< tk::real > time_interval {"time_interval",
    "Set interval (in units of physics time)",
    R"(This keyword is used to specify an interval in units of physics time.
    This must be used within a relevant block.)", "real"};

  vecentry_t< tk::real > time_range {"time_range",
    "Configure physics time range for output (in units of physics time)",
    R"(This keyword is used to configure field-, or history-output, specifying
    a start time, a stop time, and an output frequency in physics time units.
    Example: 'time_range 0.2 0.3 0.001 end', which specifies that from t=0.2 to
    t=0.3 output should happen at physics time units of dt=0.001. This must be
    used within a relevant block.)", "3 reals"};

  entry_t< bool > refined {"refined", "Toggle refined field output on/off",
    R"(This keyword can be used to turn on/off refined field output, which
    refines the mesh and evaluates the solution on the refined mesh for saving
    the solution.)", "bool"};

  entry_t< tk::ctr::FieldFileType > filetype {"filetype", "Select output file type",
    R"(This keyword is used to specify the output file type of
    mesh-based field output in a field_output block.)", "string" };

  vecentry_t< std::string > elemvar {"elemvar",
    "Specify list of elem-centered variables for output",
    R"(This keyword is used to specify elem-centered variables for output to
    file. It is used in field_output blocks.)", "string(s)"};

  vecentry_t< std::string > nodevar {"nodevar",
    "Specify list of node-centered variables for output",
    R"(This keyword is used to specify node-centered variables for output to
    file. It is used in field_output blocks.)", "string(s)"};

  vecentry_t< uint64_t > sideset {"sideset",
    "Specify list of side sets for field output",
    R"(This keyword is used to specify side sets on which field output is
    desired.)", "uint(s)"};

  // constructor
  fieldoutput_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | iter_interval;
    p | time_interval;
    p | time_range;
    p | refined;
    p | filetype;
    p | elemvar;
    p | nodevar;
    p | sideset;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i fieldoutput_t object reference
  friend void operator|( PUP::er& p, fieldoutput_t& i ) { i.pup(p); }
  //@}
};

//! Diagnostics output block
// -----------------------------------------------------------------------------
struct diagnostics_t {
  info_t info;

  // entries inside the block
  entry_t< uint32_t > iter_interval {"interval",
    "Set interval (in units of iteration count)",
    R"(This keyword is used to specify an interval in units of iteration count
    (i.e., number of time steps). This must be used within a relevant
    block.)", "uint"};

  entry_t< tk::ctr::ErrorType > error {"error", "Select an error norm",
    R"(This keyword is used to select the estimation of an error norm.
    The keyword is used in a diagnostics block.)", "string"};

  entry_t< tk::ctr::TxtFloatFormatType > format {"format",
    "Specify the ASCII floating-point output format",
    R"(This keyword is used to select the
    floating-point output format for ASCII floating-point number output.
    Valid options are 'default', 'fixed', and 'scientific'. For more info on
    these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)", "string"};

  entry_t< uint32_t > precision {"precision",
    "Precision in digits for ASCII floating-point output",
    R"(This keyword is used to select
    the precision in digits for ASCII floating-point real number output.
    Example: "precision 10", which selects ten digits for floating-point
    output, e.g., 3.141592654. The number of digits must be larger than zero
    and lower than the maximum representable digits for the given
    floating-point type. For more info on setting the precision in C++, see
    http://en.cppreference.com/w/cpp/io/manip/setprecision, and
    http://en.cppreference.com/w/cpp/types/numeric_limits/digits10)", "uint"};

  // constructor
  diagnostics_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | iter_interval;
    p | error;
    p | format;
    p | precision;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i diagnostics_t object reference
  friend void operator|( PUP::er& p, diagnostics_t& i ) { i.pup(p); }
  //@}
};

//! History point block
// -----------------------------------------------------------------------------
struct point_t {
  info_t info{"point",
    "Start configuration block for history point",
    R"(This keyword is used to introduce a point block, used to
    specify probes for history output.)", "block-title"};

  entry_t< std::string > id {"id", "Specify point identifier",
    R"(This keyword is used to specify an identifier for a history-point which
    is used in the corresponding history file name)", "string"};

  vecentry_t< tk::real > coord {"coord", "Specify point coordinates",
    R"(This keyword is used to specify coordinates of the history-point.)",
    "3 reals"};

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | id;
    p | coord;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i point_t object reference
  friend void operator|( PUP::er& p, point_t& i ) { i.pup(p); }
  //@}
};

//! History output block
// -----------------------------------------------------------------------------
struct historyoutput_t {
  info_t info;

  // entries inside the block
  entry_t< uint32_t > iter_interval {"interval",
    "Set interval (in units of iteration count)",
    R"(This keyword is used to specify an interval in units of iteration count
    (i.e., number of time steps). This must be used within a relevant
    block.)", "uint"};

  entry_t< tk::real > time_interval {"time_interval",
    "Set interval (in units of physics time)",
    R"(This keyword is used to specify an interval in units of physics time.
    This must be used within a relevant block.)", "real"};

  vecentry_t< tk::real > time_range {"time_range",
    "Configure physics time range for output (in units of physics time)",
    R"(This keyword is used to configure field-, or history-output, specifying
    a start time, a stop time, and an output frequency in physics time units.
    Example: 'time_range 0.2 0.3 0.001 end', which specifies that from t=0.2 to
    t=0.3 output should happen at physics time units of dt=0.001. This must be
    used within a relevant block.)", "3 reals"};

  entry_t< uint32_t > precision {"precision",
    "Precision in digits for ASCII floating-point output",
    R"(This keyword is used to select
    the precision in digits for ASCII floating-point real number output.
    Example: "precision 10", which selects ten digits for floating-point
    output, e.g., 3.141592654. The number of digits must be larger than zero
    and lower than the maximum representable digits for the given
    floating-point type. For more info on setting the precision in C++, see
    http://en.cppreference.com/w/cpp/io/manip/setprecision, and
    http://en.cppreference.com/w/cpp/types/numeric_limits/digits10)", "uint"};

  std::vector< point_t > point;

  // constructor
  historyoutput_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | iter_interval;
    p | time_interval;
    p | time_range;
    p | precision;
    p | point;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i historyoutput_t object reference
  friend void operator|( PUP::er& p, historyoutput_t& i ) { i.pup(p); }
  //@}
};

//! ALE move sideset block
// -----------------------------------------------------------------------------
struct move_t {
  info_t info {"move",
    "Start configuration block configuring surface movement",
    R"(This keyword is used to introduce a move block, used to
    configure surface movement for ALE simulations.)", "block-title"};

  // entries inside the block
  entry_t< uint64_t > sideset {"sideset", "Sidesets which are to be moved",
    R"(Specify the sidesets that are to be moved in ALE using user defined
    functions.)", "uint"};

  entry_t< tk::ctr::UserTableType > fntype {"fntype",
    "Select how a user-defined function is interpreted",
    R"(This keyword is used to select how a user-defined function should be
    interpreted.)", "string"};

  vecentry_t< tk::real > fn {"fn", "Specify a discrete user-defined function",
    R"(This keyword is used to specify a user-defined function block with
    discrete points, listed inside the fn block. Used in ale mesh motion and
    time-dependent BC specification)", "reals" };

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | sideset;
    p | fntype;
    p | fn;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i move_t object reference
  friend void operator|( PUP::er& p, move_t& i ) { i.pup(p); }
  //@}
};

//! ALE block
// -----------------------------------------------------------------------------
struct ale_t {
  info_t info;

  // entries inside the block
  entry_t< bool > ale {"ale", "Boolean indicating ALE",
    R"(AUTO-GENERATED boolean indicating if ALE mesh-motion is activated)",
    "bool"};

  entry_t< MeshVelocitySmootherType > smoother {"smoother",
    "Select mesh velocity smoother",
    R"(This keyword is used to select a mesh velocity smoother option, used for
    Arbitrary-Lagrangian-Eulerian (ALE) mesh motion. Valid options are
    'laplace', 'helmholtz', and 'none')", "string"};

  entry_t< MeshVelocityType > mesh_velocity {"mesh_velocity",
    "Select mesh velocity",
    R"(This keyword is used to select a mesh velocity option, used for
    Arbitrary-Lagrangian-Eulerian (ALE) mesh motion. Valid options are
    'sine', 'fluid', and 'user_defined".)", "string"};

  vecentry_t< std::size_t > mesh_motion {"mesh_motion",
    "List of dimension indices that are allowed to move in ALE calculations",
    R"(This keyword is used to specify a list of integers (0, 1, or 2) whose
    coordinate directions corresponding to x, y, or z are allowed to move with
    the mesh velocity in ALE calculations. Useful for 1D/2D problems.)",
    "uint(s)"};

  vecentry_t< tk::real > meshforce {"meshforce",
    "Set ALE mesh force model parameter(s)",
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a mesh force model for ALE. The length of the vector must
    exactly be 4.)", "4 reals"};

  std::vector< move_t > move;

  entry_t< tk::real > dvcfl {"dvcfl",
    "Set the volume-change Courant-Friedrichs-Lewy (CFL) coefficient",
    R"(This keyword is used to specify the volume-change (dV/dt) CFL coefficient
    for variable-time-step-size simulations due to volume change in time in
    arbitrary-Lagrangian-Eulerian (ALE) calculations. Setting 'dvcfl' only has
    effect in ALE calculations and used together with 'cfl'. See also J. Waltz,
    N.R. Morgan, T.R. Canfield, M.R.J. Charest, L.D. Risinger, J.G. Wohlbier, A
    three-dimensional finite element arbitrary Lagrangian–Eulerian method for
    shock hydrodynamics on unstructured grids, Computers & Fluids, 92: 172-187,
    2014.)", "real"};

  entry_t< tk::real > vortmult {"vortmult",
    "Configure vorticity multiplier for ALE mesh velocity",
    R"(This keyword is used to configure the multiplier for the vorticity term
    in the mesh velocity smoother (mesh_velocity=fluid) or for the potential
    gradient for the Helmholtz mesh velocity (mesh_velocity=helmholtz) for ALE
    mesh motion. For 'fluid' this is coefficient c2 in Eq.(36) of Waltz,
    Morgan, Canfield, Charest, Risinger, Wohlbier, A three-dimensional finite
    element arbitrary Lagrangian–Eulerian method for shock hydrodynamics on
    unstructured grids, Computers & Fluids, 2014, and for 'helmholtz', this is
    coefficient a1 in Eq.(23) of Bakosi, Waltz, Morgan, Improved ALE mesh
    velocities for complex flows, International Journal for Numerical Methods
    in Fluids, 2017. )", "real"};

  entry_t< std::size_t > maxit {"maxit",
    "Set the max number of iterations for the ALE mesh velocity linear solve",
    R"(This keyword is used to specify the maximum number of linear solver
    iterations taken to converge the mesh velocity linear solve in
    arbitrary-Lagrangian-Eulerian (ALE) calculations. See also J. Waltz,
    N.R. Morgan, T.R. Canfield, M.R.J. Charest, L.D. Risinger, J.G. Wohlbier, A
    three-dimensional finite element arbitrary Lagrangian–Eulerian method for
    shock hydrodynamics on unstructured grids, Computers & Fluids, 92: 172-187,
    2014.)", "uint"};

  entry_t< tk::real > tolerance {"tolerance",
    "Set the tolerance for the ALE mesh velocity linear solve",
    R"(This keyword is used to specify the tolerance to converge the mesh
    velocity linear solve for in
    arbitrary-Lagrangian-Eulerian (ALE) calculations. See also J. Waltz,
    N.R. Morgan, T.R. Canfield, M.R.J. Charest, L.D. Risinger, J.G. Wohlbier, A
    three-dimensional finite element arbitrary Lagrangian–Eulerian method for
    shock hydrodynamics on unstructured grids, Computers & Fluids, 92: 172-187,
    2014.)", "real"};

  // constructor
  ale_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | ale;
    p | smoother;
    p | mesh_velocity;
    p | mesh_motion;
    p | meshforce;
    p | move;
    p | dvcfl;
    p | vortmult;
    p | maxit;
    p | tolerance;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i ale_t object reference
  friend void operator|( PUP::er& p, ale_t& i ) { i.pup(p); }
  //@}
};

//! AMR coords block
// -----------------------------------------------------------------------------
struct coords_t {
  info_t info;

  // entries inside the block
  entry_t< tk::real > xminus {"xminus",
    "Configure initial refinement for coordinates lower than an x-normal plane",
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are less than the x coordinate of a plane perpendicular to
    coordinate x in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'x- <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie less than (-)
    the real number given. Example: 'x- 0.5' refines all edges whose end-point
    coordinates are less than 0.5.)", "real"};

  entry_t< tk::real > xplus {"xplus",
    "Configure initial refinement for coordinates larger than an x-normal plane",
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are larger than the x coordinate of a plane perpendicular
    to coordinate x in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'x+ <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie larger than
    (+) the real number given. Example: 'x+ 0.5' refines all edges whose
    end-point coordinates are larger than 0.5.)", "real"};

  entry_t< tk::real > yminus {"yminus",
    "Configure initial refinement for coordinates lower than an y-normal plane",
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are less than the y coordinate of a plane perpendicular to
    coordinate y in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'y- <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie less than (-)
    the real number given. Example: 'y- 0.5' refines all edges whose end-point
    coordinates are less than 0.5.)", "real"};

  entry_t< tk::real > yplus {"yplus",
    "Configure initial refinement for coordinates larger than an y-normal plane",
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are larger than the y coordinate of a plane perpendicular
    to coordinate y in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'y+ <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie larger than
    (+) the real number given. Example: 'y+ 0.5' refines all edges whose
    end-point coordinates are larger than 0.5.)", "real"};

  entry_t< tk::real > zminus {"zminus",
    "Configure initial refinement for coordinates lower than an z-normal plane",
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are less than the z coordinate of a plane perpendicular to
    coordinate z in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'z- <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie less than (-)
    the real number given. Example: 'z- 0.5' refines all edges whose end-point
    coordinates are less than 0.5.)", "real"};

  entry_t< tk::real > zplus {"zplus",
    "Configure initial refinement for coordinates larger than an z-normal plane",
    R"(This keyword can be used to configure a mesh refinement volume for edges
    whose end-points are larger than the z coordinate of a plane perpendicular
    to coordinate z in 3D space. The keyword must be used in a coords ... end
    block within an amr ... end block with syntax 'z+ <real>'. All edges of the
    input mesh will be tagged for refinement whose end-points lie larger than
    (+) the real number given. Example: 'z+ 0.5' refines all edges whose
    end-point coordinates are larger than 0.5.)", "real"};

  // constructor
  coords_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | xminus;
    p | xplus;
    p | yminus;
    p | yplus;
    p | zminus;
    p | zplus;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i coords_t object reference
  friend void operator|( PUP::er& p, coords_t& i ) { i.pup(p); }
  //@}
};

//! AMR block
// -----------------------------------------------------------------------------
struct amr_t {
  info_t info;

  // entries inside the block
  entry_t< bool > amr {"amr", "Boolean indicating AMR",
    R"(AUTO-GENERATED boolean indicating if Adaptive Mesh Refinement is
    activated)", "bool"};

  entry_t< bool > t0ref {"t0ref", "Enable mesh refinement at t<0",
    R"(This keyword is used to enable initial mesh refinement, which can be
    configured to perform multiple levels of mesh refinement based on various
    refinement criteria and configuration settings.)", "bool"};

  entry_t< bool > dtref {"dtref", "Enable mesh refinement at t>0",
    R"(This keyword is used to enable soution-adaptive mesh refinement during
    time stepping.)", "bool"};

  entry_t< bool > dtref_uniform {"dtref_uniform",
    "Enable mesh refinement at t>0 but only perform uniform refinement",
    R"(This keyword is used to force uniform-only soution-adaptive mesh
    refinement during time stepping.)", "bool"};

  entry_t< std::size_t > dtfreq {"dtfreq",
    "Set mesh refinement frequency during time stepping",
    R"(This keyword is used to configure the frequency of mesh refinement
    during time stepping.)", "uint"};

  entry_t< std::size_t > maxlevels {"maxlevels",
    "Set maximum allowed mesh refinement levels",
    R"(This keyword is used to configure the maximum allowed mesh refinement
    levels.)", "uint"};

  vecentry_t< AMRInitialType > initial {"initial",
    "Configure initial mesh refinement (before time stepping)",
    R"(This keyword is used to add to a list of initial mesh refinement types
    that happens before t = 0. Allowed options are 'uniform',
    'uniform_derefine', 'initial_conditions', 'coords', 'edgelist')", "string"};

  coords_t coords {"coords",
    "Configure initial refinement using coordinate planes",
    R"(This keyword can be used to configure entire volumes on a given side of
    a plane in 3D space. The keyword introduces an coords block within
    an amr block and must contain the either or multiple of the
    following keywords: x- <real>, x+ <real>, y- <real>, y+ <real>, z- <real>,
    z+ <real>. All edges of the input mesh will be tagged for refinement whose
    end-points lie less than (-) or larger than (+) the real number given.
    Example: 'x- 0.5' refines all edges whose end-point coordinates are less
    than 0.5. Multiple specifications are understood by combining with a
    logical AND. That is: 'x- 0.5 y+ 0.3' refines all edges whose end-point x
    coordinates are less than 0.5 AND y coordinates are larger than 0.3.)",
    "block-title"};

  vecentry_t< std::size_t > edgelist {"edgelist",
    "Configure edge-node pairs for initial refinement",
    R"(This keyword can be used to configure a list of edges that are explicitly
    tagged for initial refinement during setup in inciter. The keyword
    introduces an edgelist block within an amr block and must
    contain a list of integer pairs, i.e., the number of ids must be even,
    denoting the end-points of the nodes (=edge) which should be tagged for
    refinement.)", "even number of uints"};

  entry_t< AMRErrorType > error {"error",
    "Configure the error type for solution-adaptive mesh refinement",
    R"(This keyword is used to select the algorithm used to estimate the error
    for solution-adaptive mesh refinement. Available options are 'jump' and
    'hessian')", "string"};

  vecentry_t< char > refvar {"refvar",
    "Configure dependent variables used for adaptive mesh refinement",
    R"(This keyword is used to configured a list of dependent variables that
    trigger adaptive mesh refinement based on estimating their numerical error.
    These refinement variables are used for both initial (i.e., before time
    stepping) mesh refinement as well as during time stepping. Only previously
    (i.e., earlier in the input file) selected dependent variables can be
    configured as refinement variables. Dependent variables are required to be
    defined in all equation system configuration blocks, e.g., transport ...
    end, by using the 'depvar' keyword. Example: transport depvar c end amr
    refvar c end end. Selecting a particular scalar component in a system is
    done by appending the equation number to the refvar: Example: transport
    depvar q ncomp 3 end amr refvar q1 q2 end end, which configures two
    refinement variables: the first and third scalar component of the previously
    configured transport equation system.)", "char"};

  entry_t< tk::real > tol_refine {"tol_refine", "Configure refine tolerance",
    R"(This keyword is used to set the tolerance used to tag an edge for
    refinement if the relative error exceeds this value.)", "real"};

  entry_t< tk::real > tol_derefine {"tol_derefine",
    "Configure derefine tolerance",
    R"(This keyword is used to set the tolerance used to tag an edge for
    derefinement if the relative error is below this value.)", "real"};

  // constructor
  amr_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | amr;
    p | t0ref;
    p | dtref;
    p | dtref_uniform;
    p | dtfreq;
    p | maxlevels;
    p | initial;
    p | coords;
    p | edgelist;
    p | error;
    p | refvar;
    p | tol_refine;
    p | tol_derefine;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i amr_t object reference
  friend void operator|( PUP::er& p, amr_t& i ) { i.pup(p); }
  //@}
};

//! p-refinement block
// -----------------------------------------------------------------------------
struct pref_t {
  info_t info;

  // entries inside the block
  entry_t< bool > pref {"pref", "Indicator for p-refinement",
    R"(This AUTO-GENERATED boolean indicated if p-adaptation is triggered.)",
    "bool"};

  entry_t< PrefIndicatorType > indicator {"indicator",
    "Configure the specific adaptive indicator for p-adaptive DG scheme",
    R"(This keyword can be used to configure a specific type of adaptive
    indicator for p-adaptive refinement  of the DG scheme. The keyword must
    be used in a pref block. Available options are 'pref_spectral_decay' and
    'pref_non_conformity'.)", "string"};

  entry_t< std::size_t > ndofmax {"ndofmax",
    "Configure the maximum number of degree of freedom for p-adaptive DG",
    R"(This keyword can be used to configure a maximum number of degree of
    freedom for p-adaptive refinement  of the DG scheme. The keyword must
    be used in a pref block.)", "uint either 4 or 10"};

  entry_t< tk::real > tolref {"tolref",
    "Configure the tolerance for p-refinement for p-adaptive DG",
    R"(This keyword can be used to configure a tolerance for p-adaptive
    refinement  for the DG scheme. The keyword must be used in a pref
    block. All elements with a refinement indicator larger than this tolerance
    will be p-refined.)", "real between 0 and 1"};

  // constructor
  pref_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | pref;
    p | indicator;
    p | ndofmax;
    p | tolref;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i pref_t object reference
  friend void operator|( PUP::er& p, pref_t& i ) { i.pup(p); }
  //@}
};

//! Time-dependent Boundary conditions (bc) block
// -----------------------------------------------------------------------------
struct bctimedep_t {
  info_t info;

  // entries inside the block
  vecentry_t< uint64_t > sideset {"sideset",
    "Specify list of side sets for time-dependent BCs",
    R"(This keyword is used to specify side sets on time-dependent BCs are
    required.)", "uint(s)"};

  vecentry_t< tk::real > fn {"fn", "Specify a discrete user-defined function",
    R"(This keyword is used to specify a user-defined function block with
    discrete points, listed inside the fn block. Used in ale mesh motion and
    time-dependent BC specification)", "reals" };

  // constructor
  bctimedep_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | sideset;
    p | fn;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i bctimedep_t object reference
  friend void operator|( PUP::er& p, bctimedep_t& i ) { i.pup(p); }
  //@}
};

//! Boundary conditions (bc) block
// -----------------------------------------------------------------------------
struct bc_t {
  info_t info {"bc",
    "Start configuration block for boundary conditions",
    R"(This keyword is used to introduce the bc block, used for
    boundary conditions.)", "block-title"};

  // entries inside the block
  vecentry_t< std::size_t > mesh {"mesh",
    "List meshes on which the following BCs apply",
    R"(This keyword is used to list multiple meshes on which the boundary
    conditions listed in this particular bc-block apply.)", "uint(s)"};

  vecentry_t< uint64_t > dirichlet {"dirichlet",
    "List sidesets with Dirichlet boundary conditions",
    R"(This keyword is used to list Dirichlet sidesets.
    This keyword is used to list multiple sidesets on
    which a prescribed Dirichlet BC is then applied. Such prescribed BCs
    at each point in space and time are evaluated using a built-in function,
    e.g., using the method of manufactured solutions.)", "uint(s)"};

  vecentry_t< uint64_t > symmetry {"symmetry",
    "List sidesets with symmetry boundary conditions",
    R"(This keyword is used to list (multiple) symmetry BC sidesets.)",
    "uint(s)"};

  vecentry_t< uint64_t > inlet {"inlet",
    "List sidesets with inlet boundary conditions",
    R"(This keyword is used to list (multiple) inlet BC sidesets.)",
    "uint(s)"};

  vecentry_t< uint64_t > outlet {"outlet",
    "List sidesets with outlet boundary conditions",
    R"(This keyword is used to list (multiple) outlet BC sidesets.)",
    "uint(s)"};

  vecentry_t< uint64_t > farfield {"farfield",
    "List sidesets with farfield boundary conditions",
    R"(This keyword is used to list (multiple) farfield BC sidesets.
    Keywords allowed in a bc_farfield block are 'density', 'velocity',
    'pressure')", "uint(s)"};

  vecentry_t< uint64_t > extrapolate {"extrapolate",
    "List sidesets with Extrapolation boundary conditions",
    R"(This keyword is used to list (multiple) extrapolate BC sidesets.)",
    "uint(s)"};

  vecentry_t< uint64_t > stag {"stag",
    "List sidesets with stagnation boundary conditions",
    R"(This keyword is used to list (multiple) stagnation BC sidesets.)",
    "uint(s)"};

  vecentry_t< uint64_t > skip {"skip",
    "List sidesets with skip boundary conditions",
    R"(This keyword is used to list (multiple) skip BC sidesets. If a
    mesh point falls into a skip region, configured by a point and a radius,
    any application of boundary conditions on those points will be skipped.
    'point' and 'radius' are required to be specified within this block.)",
    "uint(s)"};

  vecentry_t< uint64_t > sponge {"sponge",
    "List sidesets with a sponge boundary",
    R"(This keyword is used to list (multiple) sponge sidesets, used to
    specify the configuration for applying sponge parameters on boundaries.
    Keywords allowed in a sponge block: 'pressure' and 'velocity'.)",
    "uint(s)"};

  bctimedep_t timedep {"timedep",
    "Start configuration block describing time dependent boundary conditions",
    R"(This keyword is used to introduce a bc_timedep block, used to
    specify the configuration of time dependent boundary conditions for a
    partial differential equation. A discrete function in time t in the form of
    a table with 6 columns (t, pressure(t), density(t), vx(t), vy(t), vz(t)) is
    expected inside a fn ... end block, specified within the bc_timedep
    block. Multiple such bc_timedep blocks can be specified for different
    time dependent BCs on different groups of side sets.)", "block-title"};

  vecentry_t< tk::real > point {"point", "Specify a point",
    R"(This keyword is used to specify a point, used, e.g., in specifying a
    point in 3D space for setting a stagnation (velocity vector = 0).)",
    "3 reals"};

  entry_t< tk::real > radius {"radius", "Specify a radius",
    R"(This keyword is used to specify a radius, used, e.g., in specifying a
    point in 3D space for setting a stagnation (velocity vector = 0).)",
    "real"};

  vecentry_t< tk::real > velocity {"velocity", "Specify velocity",
    R"(This keyword is used to configure a velocity vector, used for, e.g.,
    boundary or initial conditions.)", "3 reals"};

  entry_t< tk::real > pressure {"pressure", "Specify pressure",
    R"(This keyword is used to configure a pressure, used for, e.g., boundary
    or initial conditions or as a keyword that selects pressure in some other
    context-specific way.)", "real"};

  entry_t< tk::real > density {"density", "Specify density",
    R"(This keyword is used to configure a density, used for, e.g., boundary
    or initial conditions.)", "real"};

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | mesh;
    p | dirichlet;
    p | symmetry;
    p | inlet;
    p | outlet;
    p | farfield;
    p | extrapolate;
    p | stag;
    p | skip;
    p | sponge;
    p | timedep;
    p | point;
    p | radius;
    p | velocity;
    p | pressure;
    p | density;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i bc_t object reference
  friend void operator|( PUP::er& p, bc_t& i ) { i.pup(p); }
  //@}
};

//! Box initialization block
// -----------------------------------------------------------------------------
struct box_t {
  info_t info {"box", "Introduce a box block used to assign initial conditions",
    R"(This keyword is used to introduce a box block used to assign
    initial conditions within a box given by spatial coordinates.)",
    "block-title"};

  // entries inside the block
  entry_t< std::size_t > materialid {"materialid", "Specify material id",
    R"(This keyword is used to configure the material id within a box as a
    part of the initialization.)", "uint"};

  entry_t< tk::real > volume {"volume", "Specify volume",
    R"(This keyword is used to configure the volume of a mesh block.)",
    "real"};

  entry_t< tk::real > mass {"mass", "Specify mass",
    R"(This keyword is used to configure the mass within a box.)", "real"};

  entry_t< tk::real > density {"density", "Specify density",
    R"(This keyword is used to configure a density, used for, e.g., boundary
    or initial conditions.)", "real"};

  vecentry_t< tk::real > velocity {"velocity", "Specify velocity",
    R"(This keyword is used to configure a velocity vector, used for, e.g.,
    boundary or initial conditions or as a keyword that selects velocity in some
    other context-specific way, e.g., 'velocity' as opposed to 'position'.)",
    "3 reals"};

  entry_t< tk::real > pressure {"pressure", "Specify pressure",
    R"(This keyword is used to configure a pressure, used for, e.g., boundary
    or initial conditions or as a keyword that selects pressure in some other
    context-specific way.)", "real"};

  entry_t< tk::real > energy {"energy", "Specify energy per unit mass",
    R"(This keyword is used to configure energy per unit mass, used for, e.g.,
    boundary or initial conditions.)", "real"};

  entry_t< tk::real > energy_content {"energy_content",
    "Specify energy per unit volume",
    R"(This keyword is used to configure energy per unit volume, used for
    initial conditions.)", "real"};

  entry_t< tk::real > temperature {"temperature", "Specify temperature",
    R"(This keyword is used to configure temperature, used for, e.g.,
    boundary or initial conditions.)" , "real"};

  entry_t< tk::real > xmin {"xmin", "Minimum x coordinate",
    R"(This keyword used to configure a minimum x coordinate, e.g., to specify
    a box.)", "real"};

  entry_t< tk::real > xmax {"xmax", "Maximum x coordinate",
    R"(This keyword used to configure a maximum x coordinate, e.g., to specify
    a box.)", "real"};

  entry_t< tk::real > ymin {"ymin", "Minimum y coordinate",
    R"(This keyword used to configure a minimum y coordinate, e.g., to specify
    a box.)", "real"};

  entry_t< tk::real > ymax {"ymax", "Maximum y coordinate",
    R"(This keyword used to configure a maximum y coordinate, e.g., to specify
    a box.)", "real"};

  entry_t< tk::real > zmin {"zmin", "Minimum z coordinate",
    R"(This keyword used to configure a minimum z coordinate, e.g., to specify
    a box.)", "real"};

  entry_t< tk::real > zmax {"zmax", "Maximum z coordinate",
    R"(This keyword used to configure a maximum z coordinate, e.g., to specify
    a box.)", "real"};

  vecentry_t< tk::real > orientation {"orientation", "Configure orientation",
    R"(Configure orientation of an IC box for rotation about centroid of box.)",
    "3 reals"};

  entry_t< InitiateType > initiate {"initiate", "Initiation type",
    R"(This keyword is used to select an initiation type to configure how
    values are assigned, e.g., for a box initial condition. This can be used to
    specify, how the values are assigned to mesh nodes within a box. Uses:
    (1) impulse: assign the full values at t=0 for all points in a box,
    (2) linear: use a linear function in time and space, configured with an
    initiation point in space, a constant velocity of the growing spherical
    front in time (and space) linearly, and width of the front and assigns
    values to mesh points falling within the growing spherical shell inside a
    configured box.)", "string"};

  vecentry_t< tk::real > point {"point", "Specify a point",
    R"(This keyword is used to specify a point, used, e.g., in specifying a
    point in 3D space for the 'initiate linear' seed point.)", "3 reals"};

  entry_t< tk::real > init_time {"init_time","Specify the initialization time",
    R"(This keyword is used to specify the time at which the propagating front
    is initialized for a mesh block or box IC, with 'initiate linear' type.
    Delays in initializing separate mesh blocks or boxes can be achieved using
    different initialization times.)", "real"};

  entry_t< tk::real > front_width {"front_width", "Specify a front_width",
    R"(This keyword is used to specify the width of the propagating front for
    a mesh block or box IC, with 'initiate linear' type. The suggested value of
    the front width is about 4-5 times the mesh size inside the mesh block
    or box.)", "real"};

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | materialid;
    p | volume;
    p | mass;
    p | density;
    p | velocity;
    p | pressure;
    p | energy;
    p | energy_content;
    p | temperature;
    p | xmin;
    p | xmax;
    p | ymin;
    p | ymax;
    p | zmin;
    p | zmax;
    p | orientation;
    p | initiate;
    p | point;
    p | init_time;
    p | front_width;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i box_t object reference
  friend void operator|( PUP::er& p, box_t& i ) { i.pup(p); }
  //@}
};

//! Mesh-block initialization block
// -----------------------------------------------------------------------------
struct block_t {
  info_t info {"meshblock",
    "Introduce a meshblock block used to assign initial conditions",
    R"(This keyword is used to introduce a meshblock block used to
    assign initial conditions within a mesh block specified in the mesh file.)",
    "block-title"};

  // entries inside the block
  entry_t< uint64_t > blockid {"blockid", "Specify mesh block id",
    R"(This keyword is used to configure the mesh block id as a part
    of the initialization. It is strongly recommended to use contiguous block
    ids in mesh file starting from 1.)", "uint"};

  entry_t< std::size_t > materialid {"materialid", "Specify material id",
    R"(This keyword is used to configure the material id within a box as a
    part of the initialization.)", "uint"};

  entry_t< tk::real > volume {"volume", "Specify volume",
    R"(This keyword is used to configure the volume of a mesh block.)",
    "real"};

  entry_t< tk::real > mass {"mass", "Specify mass",
    R"(This keyword is used to configure the mass within a box.)", "real"};

  entry_t< tk::real > density {"density", "Specify density",
    R"(This keyword is used to configure a density, used for, e.g., boundary
    or initial conditions.)", "real"};

  vecentry_t< tk::real > velocity {"velocity", "Specify velocity",
    R"(This keyword is used to configure a velocity vector, used for, e.g.,
    boundary or initial conditions or as a keyword that selects velocity in some
    other context-specific way, e.g., 'velocity' as opposed to 'position'.)",
    "3 reals"};

  entry_t< tk::real > pressure {"pressure", "Specify pressure",
    R"(This keyword is used to configure a pressure, used for, e.g., boundary
    or initial conditions or as a keyword that selects pressure in some other
    context-specific way.)", "real"};

  entry_t< tk::real > energy {"energy", "Specify energy per unit mass",
    R"(This keyword is used to configure energy per unit mass, used for, e.g.,
    boundary or initial conditions.)", "real"};

  entry_t< tk::real > energy_content {"energy_content",
    "Specify energy per unit volume",
    R"(This keyword is used to configure energy per unit volume, used for
    initial conditions.)", "real"};

  entry_t< tk::real > temperature {"temperature", "Specify temperature",
    R"(This keyword is used to configure temperature, used for, e.g.,
    boundary or initial conditions.)" , "real"};

  entry_t< InitiateType > initiate {"initiate", "Initiation type",
    R"(This keyword is used to select an initiation type to configure how
    values are assigned, e.g., for a box initial condition. This can be used to
    specify, how the values are assigned to mesh nodes within a box. Uses:
    (1) impulse: assign the full values at t=0 for all points in a box,
    (2) linear: use a linear function in time and space, configured with an
    initiation point in space, a constant velocity of the growing spherical
    front in time (and space) linearly, and width of the front and assigns
    values to mesh points falling within the growing spherical shell inside a
    configured box.)", "string"};

  vecentry_t< tk::real > point {"point", "Specify a point",
    R"(This keyword is used to specify a point, used, e.g., in specifying a
    point in 3D space for setting a stagnation (velocity vector = 0), or
    for the 'initiate linear' seed point.)", "3 reals"};

  entry_t< tk::real > init_time {"init_time","Specify the initialization time",
    R"(This keyword is used to specify the time at which the propagating front
    is initialized for a mesh block or box IC, with 'initiate linear' type.
    Delays in initializing separate mesh blocks or boxes can be achieved using
    different initialization times.)", "real"};

  entry_t< tk::real > front_width {"front_width", "Specify a front_width",
    R"(This keyword is used to specify the width of the propagating front for
    a mesh block or box IC, with 'initiate linear' type. The suggested value of
    the front width is about 4-5 times the mesh size inside the mesh block
    or box.)", "real"};

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | blockid;
    p | materialid;
    p | volume;
    p | mass;
    p | density;
    p | velocity;
    p | pressure;
    p | energy;
    p | energy_content;
    p | temperature;
    p | initiate;
    p | point;
    p | init_time;
    p | front_width;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i block_t object reference
  friend void operator|( PUP::er& p, block_t& i ) { i.pup(p); }
  //@}
};

//! Initial conditions (ic) block
// -----------------------------------------------------------------------------
struct ic_t {
  info_t info;

  // entries inside the block
  entry_t< std::size_t > materialid {"materialid", "Specify material id",
    R"(This keyword is used to configure the material id within a box or in the
    background as a part of the initialization.)", "uint"};

  entry_t< tk::real > pressure {"pressure", "Specify pressure",
    R"(This keyword is used to configure a pressure, used for, e.g., boundary
    or initial conditions or as a keyword that selects pressure in some other
    context-specific way.)", "real"};

  entry_t< tk::real > temperature {"temperature", "Specify temperature",
    R"(This keyword is used to configure temperature, used for, e.g.,
    boundary or initial conditions.)" , "real"};

  vecentry_t< tk::real > velocity {"velocity", "Specify velocity",
    R"(This keyword is used to configure a velocity vector, used for, e.g.,
    boundary or initial conditions or as a keyword that selects velocity in some
    other context-specific way, e.g., 'velocity' as opposed to 'position'.)",
    "3 reals"};

  std::vector< box_t > box;

  std::vector< block_t > meshblock;

  // constructor
  ic_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | materialid;
    p | pressure;
    p | temperature;
    p | velocity;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i ic_t object reference
  friend void operator|( PUP::er& p, ic_t& i ) { i.pup(p); }
  //@}
};

//! Mesh (for overset) block
// -----------------------------------------------------------------------------
struct mesh_t {
  info_t info {"mesh",
    "Start configuration block assigning a mesh to a solver",
    R"(This keyword is used to introduce a mesh block, used to
    assign and configure a mesh to a solver.)", "block-title"};

  // entries inside the block
  entry_t< std::string > filename {"filename", "Set filename",
    R"(Set filename, e.g., mesh filename for solver coupling.)", "string"};

  vecentry_t< tk::real > location {"location", "Configure location",
    R"(Configure location of a mesh relative to another.)", "3 reals"};

  vecentry_t< tk::real > orientation {"orientation", "Configure orientation",
    R"(Configure orientation of a mesh relative to another.)",
    "3 reals"};

  vecentry_t< tk::real > velocity {"velocity", "Specify overset mesh velocity",
    R"(This keyword is used to configure a velocity vector used for
    moving the overset mesh.)",
    "3 reals"};

  /** @name Charm++ pack/unpack serializer member functions */
  ///@{
  //! \brief Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er &p ) {
    p | info;
    p | filename;
    p | location;
    p | orientation;
    p | velocity;
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] i mesh_t object reference
  friend void operator|( PUP::er& p, mesh_t& i ) { i.pup(p); }
  //@}
};

//! The new input deck struct for inciter
//! \details The stack is a tagged tuple, a hierarchical heterogeneous data
//!    structure where all parsed information is stored.
struct NewInputDeck {

  public:
    //! \brief Inciter input deck keywords

    entry_t< std::string > inciter {"inciter",
      "Start configuration block for inciter",
      R"(This keyword is used to select inciter. Inciter, is a continuum-realm
      shock hydrodynamics tool, solving a system of PDEs.)", "string"};

    entry_t< std::string > title {"title", "Title", R"(The title may be
      specified in the input file. It is optional.)", "string"};

    // -------------------------------------------------------------------------
    // time stepping options
    // -------------------------------------------------------------------------

    entry_t< uint64_t > nstep {"nstep", "Set number of time steps to take",
      R"(This keyword is used to specify the number of time steps to take in a
      simulation. The number of time steps are used in conjunction with the
      maximmum time specified by keyword 'term': the simulation stops whichever
      is reached first. Both 'nstep' and 'term' can be left unspecified, in
      which case their default values are used. See also 'term'.)", "uint"};

    entry_t< tk::real > term {"term", "Set maximum physical time to simulate",
      R"(This keyword is used to specify the termination time in a simulation.
      The termination time and number of time steps, specified by 'nstep', are
      used in conjunction to determine when to stop a simulation: whichever is
      reached first. Both 'nstep' and 'term' can be left unspecified, in which
      case their default values are used. See also 'nstep'.)", "real"};

    entry_t< tk::real > t0 {"t0", "Set starting non-dimensional time",
      R"(This keyword is used to specify the starting time in a simulation.)",
      "real"};

    entry_t< tk::real > dt {"dt", "Select constant time step size",
      R"(This keyword is used to specify the time step size that used as a
      constant during simulation. Setting 'cfl' and 'dt' are mutually
      exclusive. If both 'cfl' and 'dt' are set, 'dt' wins.)", "real"};

    entry_t< tk::real > cfl {"cfl",
    "Set the Courant-Friedrichs-Lewy (CFL) coefficient",
    R"(This keyword is used to specify the CFL coefficient for
    variable-time-step-size simulations. Setting 'cfl' and 'dt' are mutually
    exclusive. If both 'cfl' and 'dt' are set, 'dt' wins.)", "real"};

    entry_t< uint32_t > ttyi {"ttyi", "Set screen output interval",
      R"(This keyword is used to specify the interval in time steps for screen
      output during a simulation.)", "uint"};

    // -------------------------------------------------------------------------
    // steady-state solver options
    // -------------------------------------------------------------------------

    entry_t< bool > steady_state {"steady_state", "March to steady state",
      R"(This keyword is used indicate that local time stepping should be used
      to march towards a stationary solution.)", "bool"};

    entry_t< tk::real > residual {"residual",
      "Set the convergence criterion for the residual to reach",
      R"(This keyword is used to specify a convergence criterion for the local
      time stepping marching to steady state, below which the simulation is
      considered converged.)", "real"};

    entry_t< uint32_t > rescomp {"rescomp",
      "Equation system component index for convergence",
      R"(This keyword is used to specify a single integer that is used to denote
      the equation component index in the complete system of equations
      configured, to use for the convergence criterion for local
      time stepping marching towards steady state.)", "uint"};

    // -------------------------------------------------------------------------
    // mesh partitioning and reordering/sorting choices
    // -------------------------------------------------------------------------

    entry_t< tk::ctr::PartitioningAlgorithmType > partitioning {"partitioning",
      "Select mesh partitioning algorithm",
      R"(This keyword is used to select a mesh partitioning algorithm. See
      Control/Options/PartitioningAlgorithm.hpp for valid options.)", "string"};

    option_t rcb {"rcb",
      "Select recursive coordinate bisection mesh partitioner",
      R"(This keyword is used to select the recursive coordinate bisection (RCB)
      mesh partitioner. RCB is a geometry-based partitioner used to distribute
      an input mesh among processing elements. See
      Control/Options/PartitioningAlgorithm.hpp for other valid options.)",
      "string"};

    option_t rib {"rib",
      "Select recursive inertial bisection mesh partitioner",
      R"(This keyword is used to select the recursive inertial bisection (RIB)
      mesh partitioner. RIB is a geometry-based partitioner used to distribute
      an input mesh among processing elements. See
      Control/Options/PartitioningAlgorithm.hpp for other valid options.)",
      "string"};

    option_t hsfc {"hsfc",
      "Select Hilbert Space Filling Curve (HSFC) mesh partitioner",
      R"(This keyword is used to select the Hilbert Space Filling Curve (HSFC)
      mesh partitioner. HSFC is a geometry-based partitioner used to distribute
      an input mesh among processing elements. See
      Control/Options/PartitioningAlgorithm.hpp for other valid options.)",
      "string"};

    option_t phg {"phg",
      "Select parallel hypergraph mesh partitioner",
      R"(This keyword is used to select the parallel hypergraph (PHG)
      mesh partitioner. PHG is a graph-based partitioner used to distribute an
      input mesh among processing elements. See
      Control/Options/PartitioningAlgorithm.hpp for other valid options.)",
      "string"};

    option_t mj {"mj",
      "Select multi-jagged (MJ) mesh partitioner",
      R"(This keyword is used to select the multi-jagged (MJ) mesh partitioner.
      MJ is a geometry-based partitioner used to distribute an input mesh among
      processing elements. See
      Control/Options/PartitioningAlgorithm.hpp for other valid options.)",
      "string"};

    entry_t< bool > pelocal_reorder {"pelocal_reorder",
      "PE-local reorder",
      R"(This keyword is used in inciter as a keyword in the inciter...end block
      as "pelocal_reorder true" (or false) to do (or not do) a global
      distributed mesh reordering across all PEs that yields an approximately
      continuous mesh node ID order as mesh partitions are assigned to PEs after
      mesh partitioning. This reordering is optional.)", "bool"};

    entry_t< bool > operator_reorder {"operator_reorder",
      "Operator-access reorder",
      R"(This keyword is used in inciter as a keyword in the inciter...end block
      as "operator_reorder on" (or off) to do (or not do) a local mesh node
      reordering based on the PDE operator access pattern. This reordering is
      optional.)", "bool"};

    // -------------------------------------------------------------------------
    // discretization scheme choices
    // -------------------------------------------------------------------------

    entry_t< SchemeType > scheme {"scheme", "Select discretization scheme",
      R"(This keyword is used to select a spatial discretization scheme,
      necessarily connected to the temporal discretization scheme. See
      Control/Inciter/Options/Scheme.hpp for valid options.)", "string"};

    option_t diagcg {"diagcg",
      "Select continuous Galerkin + Lax Wendroff with a lumped-mass matrix LHS",
      R"(This keyword is used to select the lumped-mass matrix continuous Galerkin
      (CG) finite element spatial discretiztaion used in inciter. CG is combined
      with a Lax-Wendroff scheme for time discretization and flux-corrected
      transport (FCT) for treating discontinuous solutions. This option selects
      the scheme that stores the left-hand side matrix lumped, i.e., only the
      diagonal elements stored and thus does not require a linear solver. See
      Control/Inciter/Options/Scheme.hpp for other valid options.)", "string"};

    option_t alecg {"alecg",
      "Select continuous Galerkin with ALE + Runge-Kutta",
      R"(This keyword is used to select the continuous Galerkin finite element
      scheme in the arbitrary Lagrangian-Eulerian (ALE) reference frame combined
      with Runge-Kutta (RK) time stepping.
      See Control/Inciter/Options/Scheme.hpp for other valid options.)",
      "string"};

    option_t oversetfe {"oversetfe",
      "Select continuous Galerkin finite element with overset meshes + "
      "Runge-Kutta",
      R"(This keyword is used to select the continuous Galerkin finite element
      scheme with Runge-Kutta (RK) time stepping, combined with overset grids.
      See Control/Inciter/Options/Scheme.hpp for other valid options.)",
      "string"};

    option_t dg {"dg",
      "Select 1st-order discontinuous Galerkin discretization + Runge-Kutta",
      R"(This keyword is used to select the first-order accurate discontinuous
      Galerkin, DG(P0), spatial discretiztaion used in Inciter. As this is first
      order accurate, it is intended for testing and debugging purposes only.
      Selecting this spatial discretization also selects the Runge-Kutta scheme
      for time discretization. See Control/Inciter/Options/Scheme.hpp for other
      valid options.)", "string"};

    option_t p0p1 {"p0p1",
      "Select 2nd-order finite volume discretization + Runge-Kutta",
      R"(This keyword is used to select the second-order accurate finite volume,
      P0P1, spatial discretiztaion used in Inciter. This method uses a
      least-squares procedure to reconstruct the second-order solution from the
      first-order one. Selecting this spatial discretization also selects the
      Runge-Kutta scheme for time discretization.
      See Control/Inciter/Options/Scheme.hpp for other valid options.)",
      "string"};

    option_t dgp1 {"dgp1",
      "Select 2nd-order discontinuous Galerkin discretization + Runge-Kutta",
      R"(This keyword is used to select the second-order accurate discontinuous
      Galerkin, DG(P1), spatial discretiztaion used in Inciter. Selecting this
      spatial discretization also selects the Runge-Kutta scheme for time
      discretization. See Control/Inciter/Options/Scheme.hpp for other
      valid options.)", "string"};

    option_t dgp2 {"dgp2",
      "Select 3nd-order discontinuous Galerkin discretization + Runge-Kutta",
      R"(This keyword is used to select the third-order accurate discontinuous
      Galerkin, DG(P2), spatial discretiztaion used in Inciter. Selecting this
      spatial discretization also selects the Runge-Kutta scheme for time
      discretization. See Control/Inciter/Options/Scheme.hpp for other
      valid options.)", "string"};

    option_t pdg {"pdg",
      "Select p-adaptive discontinuous Galerkin discretization + Runge-Kutta",
      R"(This keyword is used to select the polynomial adaptive discontinuous
      Galerkin spatial discretizaion used in Inciter. Selecting this spatial
      discretization also selects the Runge-Kutta scheme for time
      discretization. See Control/Inciter/Options/Scheme.hpp for other valid
      options.)", "string"};

    option_t fv {"fv",
      "Select 2nd-order finite volume discretization + Runge-Kutta",
      R"(This keyword is used to select the second-order accurate finite volume,
      P0P1, spatial discretiztaion used in Inciter. This method uses a
      least-squares procedure to reconstruct the second-order solution from the
      first-order one. See Control/Inciter/Options/Scheme.hpp for other valid
      options.)", "string"};

    entry_t< std::size_t > ndof {"ndof", "Number of evolved solution DOFs",
      R"(The number of solution DOFs that are evolved.)", "uint"};

    entry_t< std::size_t > rdof {"rdof", "Total number of solution DOFs",
      R"(The total number of solution DOFs, including the reconstructed and the
      evolved ones.)", "uint"};

    // -------------------------------------------------------------------------
    // limiter options
    // -------------------------------------------------------------------------

    entry_t< LimiterType > limiter {"limiter", "Select limiter function",
      R"(This keyword is used to select a limiter function, used for
      discontinuous Galerkin (DG) spatial discretization used in inciter. See
      Control/Inciter/Options/Limiter.hpp for valid options.)", "string"};

    option_t nolimiter {"nolimiter", "No limiter used",
      R"(This keyword is used for discontinuous Galerkin (DG) spatial
      discretization without any limiter in inciter. See
      Control/Inciter/Options/Limiter.hpp for other valid options.)", "string"};

    option_t wenop1 {"wenop1",
      "Select the Weighted Essentially Non-Oscillatory (WENO) limiter for DGP1",
      R"(This keyword is used to select the Weighted Essentially Non-Oscillatory
      limiter used for discontinuous Galerkin (DG) P1 spatial discretization
      used in inciter. See Control/Inciter/Options/Limiter.hpp for other valid
      options.)", "string"};

    entry_t< tk::real > cweight {"cweight",
      "Set value for central linear weight used by WENO, cweight",
      R"(This keyword is used to set the central linear weight used for the
      central stencil in the Weighted Essentially Non-Oscillatory (WENO) limiter
      for discontinuous Galerkin (DG) methods.)", "real"};

    option_t superbeep1 {"superbeep1",
      "Select the Superbee limiter for DGP1",
      R"(This keyword is used to select the Superbee limiter used for
      discontinuous Galerkin (DG) P1 spatial discretization used in inciter.
      See Control/Inciter/Options/Limiter.hpp for other valid options.)",
      "string"};

    entry_t< tk::real > shock_detector_coeff {"shock_detector_coeff",
      "Configure the coefficient used in shock indicator",
      R"(This keyword can be used to configure the coefficient used in the
      threshold calculation for the shock indicator.)", "real"};

    option_t vertexbasedp1 {"vertexbasedp1",
      "Select the vertex-based limiter for DGP1",
      R"(This keyword is used to select the vertex-based limiter used for
      discontinuous Galerkin (DG) P1 spatial discretization used in inciter.
      Ref. Kuzmin, D. (2010). A vertex-based hierarchical slope limiter for
      p-adaptive discontinuous Galerkin methods. Journal of computational and
      applied mathematics, 233(12), 3077-3085.
      See Control/Inciter/Options/Limiter.hpp for other valid options.)",
      "string"};

    entry_t< bool > accuracy_test {"accuracy_test", "Toggle accuracy test setup",
      R"(This keyword is used to specify if the current setup is for an
      order-of-accuracy testing, used for discontinuous Galerkin (DG) spatial
      discretization in inciter. This deactivates certain robustness corrections
      which might impact order-of-accuracy. Only intended for simple test
      problems and not for real problems.)", "bool"};

    entry_t< bool > limsol_projection {"limsol_projection",
      "Toggle limited solution projection",
      R"(This keyword is used to specify limited solution projection.
      This is used for discontinuous Galerkin (DG) spatial discretization in
      inciter, for multi-material hydrodynamics. This uses a projection to
      obtain bulk momentum and material energies from the limited primitive
      quantities. This step is essential to obtain closure-law obeying limited
      quantities. See Pandare et al. (2023). On the Design of Stable,
      Consistent, and Conservative High-Order Methods for Multi-Material
      Hydrodynamics. J Comp Phys (490).)", "bool"};

    entry_t< bool > fct {"fct", "Turn flux-corrected transport on/off",
      R"(This keyword can be used to turn on/off flux-corrected transport (FCT).
      Note that FCT is only used in conjunction with continuous Galerkin finite
      element discretization, configured by scheme diagcg and it has no
      effect when the discontinuous Galerkin (DG) scheme is used, configured by
      'scheme dg'. Also note that even if FCT is turned off, it is still
      performed, only its result is not applied.)", "bool"};

    entry_t< bool > fctclip {"fctclip",
      "Turn on clipping flux-corrected transport on/off",
      R"(This keyword can be used to turn on/off the clipping limiter used for
      flux-corrected transport (FCT). The clipping limiter only looks at the
      current low order solution to determine the allowed solution minima and
      maxima, instead of the minimum and maximum of the low order solution and
      the previous solution.)", "bool"};

    entry_t< tk::real > fcteps {"fcteps",
      "A number that is considered small enough for FCT",
      R"(This keyword is used to set the epsilon (a small number) below which
      FCT quantities are considered small enough to be treated as zero. Setting
      this number to be somewhat larger than the machine zero, e.g., 1.0e-15,
      helps ignoring some noise that otherwise could contaminate the
      solution.)", "real"};

    entry_t< tk::real > ctau {"ctau", "Set FCT mass diffusion coefficient",
      R"(This keyword is used to set the mass diffusion coefficient used in
      flux-corrected transport, used for integrating transport equations.)",
      "real"};

    entry_t< bool > sysfct {"sysfct",
      "Turn on system nature of flux-corrected transport",
      R"(This keyword can be used to enable a system-nature for flux-corrected
      transport (FCT). Note that FCT is only used in conjunction with continuous
      Galerkin finite element discretization, configured by scheme diagcg and it
      has no effect when the discontinuous Galerkin (DG) scheme is used,
      configured by 'scheme dg'. Enabling the system-nature for FCT will choose
      the limiter coefficients for a system of equations, e.g., compressible
      flow, in way that accounts for the system-nature of the equations. An
      example is assinging the minimum of the limit coefficient to all variables
      limited in a computational cell, e.g., density, momentum, and specitic
      total energy. This yields better, more monotonic, results.)", "bool"};

    vecentry_t< std::size_t > sysfctvar {"sysfctvar",
      "Specify a list of scalar component indices considered for system FCT",
      R"(This keyword is used to specify a list of integers that are considered
      for computing the system-nature of flux-corrected transport. Example:
      'sysfctvar 0 1 2 3 end', which means ignoring the energy (by not listing 4)
      when computing the coupled limit coefficient for a system of mass, momentum,
      and energy for single-material compressible flow.)", "uint(s)"};

    // -------------------------------------------------------------------------
    // flux options
    // -------------------------------------------------------------------------

    entry_t< FluxType > flux {"flux", "Select flux function",
      R"(This keyword is used to select a flux function, used for
      discontinuous Galerkin (DG) spatial discretization used in inciter. See
      Control/Inciter/Options/Flux.hpp for valid options.)", "string"};

    option_t laxfriedrichs {"laxfriedrichs",
      "Select Lax-Friedrichs flux function",
      R"(This keyword is used to select the Lax-Friedrichs flux function used
      for discontinuous Galerkin (DG) spatial discretization used in inciter.
      See Control/Inciter/Options/Flux.hpp for other valid options.)",
      "string"};

    option_t hllc {"hllc",
      "Select the Harten-Lax-van Leer-Contact (HLLC) flux function",
      R"(This keyword is used to select the Harten-Lax-van Leer-Contact flux
      function used for discontinuous Galerkin (DG) spatial discretization
      used in inciter. See Control/Inciter/Options/Flux.hpp for other valid
      options.)", "string"};

    option_t upwind {"upwind", "Select the upwind flux function",
      R"(This keyword is used to select the upwind flux
      function used for discontinuous Galerkin (DG) spatial discretization
      used in inciter. It is only usable for scalar transport.
      See Control/Inciter/Options/Flux.hpp for other valid options.)",
      "string"};

    option_t ausm {"ausm",
      "Select the Advection Upstream Splitting Method (AUSM) flux function",
      R"(This keyword is used to select the AUSM flux
      function used for discontinuous Galerkin (DG) spatial discretization
      used in inciter. It is only set up for for multi-material hydro, and
      not selectable for anything else.)", "string"};

    option_t hll {"hll",
      "Select the Harten-Lax-vanLeer (HLL) flux function",
      R"(This keyword is used to select the HLL flux
      function used for discontinuous Galerkin (DG) spatial discretization
      used in inciter. It is only set up for for multi-material hydro, and
      not selectable for anything else.)", "string"};

    entry_t< PDEType > pde {"pde",
      "Storage for the PDE type selected",
      R"(This AUTO-GENERATED string holds the PDE type specified by the user.)",
      "ctr::PDEType"};

    // -------------------------------------------------------------------------
    // Transport object
    // -------------------------------------------------------------------------

    transport_t transport {"transport",
      "Start configuration block for an transport equation",
      R"(This keyword is used to introduce a transport block, used to
      specify the configuration for a transport equation type.)",
      "block-title"};

    // -------------------------------------------------------------------------
    // CompFlow object
    // -------------------------------------------------------------------------

    compflow_t compflow {"compflow",
      "Start configuration block for the compressible flow equations",
      R"(This keyword is used to introduce the compflow block, used to
      specify the configuration for a system of partial differential equations,
      governing single material compressible fluid flow.)", "block-title"};

    // -------------------------------------------------------------------------
    // MultiMat object
    // -------------------------------------------------------------------------

    multimat_t multimat {"multimat",
      "Start configuration block for the compressible multi-material equations",
      R"(This keyword is used to introduce the multimat block,
      used to specify the configuration for a system of partial differential
      equations, governing compressible multi-material hydrodynamics assuming
      velocity equilibrium (single velocity).)", "block-title"};

    // Dependent variable name

    entry_t< std::string > depvar {"depvar",
      "Select dependent variable name for PDE.",
      R"(Select dependent variable name for PDE.)", "string"};

    // -------------------------------------------------------------------------
    // physics choices
    // -------------------------------------------------------------------------

    entry_t< PhysicsType > physics {"physics",
      "Specify the physics configuration for a system of PDEs",
      R"(This keyword is used to select the physics configuration for a
      particular PDE system. Valid options depend on the system of PDEs in which
      the keyword is used.)", "string"};

    option_t advection {"advection",
      "Specify the advection physics",
      R"(This keyword is used to select the advection physics for the transport
      PDE system. Only usable for 'transport'.)", "string"};

    option_t advdiff {"advdiff",
      "Specify the advection + diffusion physics",
      R"(This keyword is used to select the advection + diffusion physics
      for transport PDEs. Only usable for 'transport'.)", "string"};

    option_t navierstokes {"navierstokes",
      "Specify the Navier-Stokes (viscous) compressible flow physics",
      R"(This keyword is used to select the Navier-Stokes (viscous) compressible
      flow physics configuration. Currently setup only for 'compflow'.)",
      "string"};

    option_t euler {"euler",
      "Specify the Euler (inviscid) compressible flow physics",
      R"(This keyword is used to select the Euler (inviscid) compressible
      flow physics configuration. Usable for 'compflow' and 'multimat')",
      "string"};

    option_t energy_pill {"energy_pill",
      "Specify the energy pill physics",
      R"(This keyword is used to select an energy pill initialization as physics
      configuration for multiple material compressible flow. Parameters for the
      linearly traveling front are required to be specified when energy_pill is
      selected. See 'linear' for more details. Currently setup only for
      'multimat')", "string"};

    entry_t< std::size_t > ncomp {"ncomp",
      "Stores the total number of equations to solve.",
      R"(This keyword is used to store the number of equations in the PDE
      system.)", "uint"};

    // -------------------------------------------------------------------------
    // material/eos object
    // -------------------------------------------------------------------------

    std::vector< material_t > material;

    option_t stiffenedgas {"stiffenedgas",
      "Select the stiffened gas equation of state",
      R"(This keyword is used to select the stiffened gas equation of state.)",
      "string"};

    option_t jwl {"jwl", "Select the JWL equation of state",
      R"(This keyword is used to select the Jones, Wilkins, Lee equation of
      state.)", "string"};

    option_t smallshearsolid {"smallshearsolid",
      "Select the SMALLSHEARSOLID equation of state",
      R"(This keyword is used to select the small shear strain equation of state
      for solids. This EOS uses a small-shear approximation for the elastic
      contribution, and a stiffened gas EOS for the hydrodynamic contribution of
      the internal energy See Plohr, J. N., & Plohr, B. J. (2005). Linearized
      analysis of Richtmyer–Meshkov flow for elastic materials. Journal of Fluid
      Mechanics, 537, 55-89 for further details.)", "string"};

    matidxmap_t matidxmap {"matidxmap",
    "Material index map for EOS",
    R"(The following AUTO-GENERATED data structure is used to index into the
    correct material vector entry. This is done using the following three maps:
    1. eosidx: This vector provides the eos-index (value) in the
    vector<tag::material> for the given user-spec material id (index).
    2. matidx: This vector provides the material-index (value) inside the
    vector<tag::material>[eosidx] block for the given user-specified
    material id (index).
    3. solidx: This vector provides the solid-index (value) assigned to
    the given user-specified material id (index). It is 0 for fluids.)", "map"};

    // -------------------------------------------------------------------------
    // output object
    // -------------------------------------------------------------------------

    fieldoutput_t field_output {"field_output",
      "Start of field_output input block",
      R"(This keyword is used to start a block in the input file containing the
      list and settings of requested field output.)", "block-title"};

    diagnostics_t diagnostics {"diagnostics",
      "Specify the diagnostics block",
      R"(This keyword is used to introduce the dagnostics block, used to
      configure diagnostics output.)", "block-title"};

    historyoutput_t history_output {"history_output",
      "Start of history_output input block",
      R"(This keyword is used to start a block in the input file containing the
      descriptions and settings of requested history output.)", "string"};

    option_t exodusii {"exodusii", "Select ExodusII output",
      R"(This keyword is used to select the
      ExodusII output file type readable by, e.g., ParaView of either a requested
      probability density function (PDF) within a pdfs ... end block or for
      mesh-based field output in a field_output ... end block. Example:
      "filetype exodusii", which selects ExodusII file output. For more info on
      ExodusII, see http://sourceforge.net/projects/exodusii.)", "string"};

    option_t outvar_density {"outvar_density", "Request density",
      R"(This keyword is used to request the fluid density as an output
         variable.)", "string"};

    option_t outvar_xmomentum {"outvar_xmomentum",
      "Request x-momentum",
      R"(This keyword is used to request the fluid x-momentum as an output
      variable.)", "string"};

    option_t outvar_ymomentum {"outvar_ymomentum",
      "Request y-momentum",
      R"(This keyword is used to request the fluid y-momentum as an output
      variable.)", "string"};

    option_t outvar_zmomentum {"outvar_zmomentum",
      "Request z-momentum",
      R"(This keyword is used to request the fluid z-momentum as an output
      variable.)", "string"};

    option_t outvar_specific_total_energy {
      "outvar_specific_total_energy", "Request specific total energy",
      R"(This keyword is used to request the specific total energy as an output
      variable.)", "string"};

    option_t outvar_volumetric_total_energy {
      "outvar_volumetric_total_energy", "Request total volumetric energy",
      R"(This keyword is used to request the volumetric total energy as an output
      variable.)", "string"};

    option_t outvar_xvelocity {"outvar_xvelocity",
      "Request x-velocity",
      R"(This keyword is used to request the fluid x-velocity as an output
      variable.)", "string"};

    option_t outvar_yvelocity {"outvar_yvelocity",
      "Request y-velocity",
      R"(This keyword is used to request the fluid y-velocity as an output
      variable.)", "string"};

    option_t outvar_zvelocity {"outvar_zvelocity",
      "Request z-velocity",
      R"(This keyword is used to request the fluid z-velocity as an output
      variable.)", "string"};

    option_t outvar_pressure {"outvar_pressure",
     "Request pressure",
      R"(This keyword is used to request the fluid pressure as an output
      variable.)", "string"};

    option_t outvar_material_indicator {"outvar_material_indicator",
      "Request material_indicator",
      R"(This keyword is used to request the material indicator function as an
      output variable.)", "string"};

    option_t outvar_analytic {"outvar_analytic",
      "Request analytic solution",
      R"(This keyword is used to request the analytic solution (if exist) as an
      output variable.)", "string"};

    option_t l2 {"l2", "Select the L2 norm",
      R"(This keyword is used to enable computing the L2 norm.)", "string"};

    option_t linf {"linf", "Select the L_{infinity} norm",
      R"(This keyword is used to enable computing the L-infinity norm.)",
      "string"};

    option_t txt_float_default {"default",
      "Select the default ASCII floating-point output",
      R"(This keyword is used to select the
      'default' floating-point output format for ASCII floating-point real
      number output. For more info on these various formats, see
      http://en.cppreference.com/w/cpp/io/manip/fixed.)", "string"};

    option_t fixed {"fixed",
      "Select the fixed ASCII floating-point output",
      R"(This keyword is used to select the
      'fixed' floating-point output format for ASCII floating-point real
      number output. For more info on these various formats, see
      http://en.cppreference.com/w/cpp/io/manip/fixed.)", "string"};

    option_t scientific {"scientific",
      "Select the scientific ASCII floating-point output",
      R"(This keyword is used to select the
      'scientific' floating-point output format for ASCII floating-point real
      number output. For more info on these various formats, see
      http://en.cppreference.com/w/cpp/io/manip/fixed.)", "string"};

    // -------------------------------------------------------------------------
    // ALE options
    // -------------------------------------------------------------------------

    ale_t ale {"ale", "Start block configuring ALE",
      R"(This keyword is used to introduce the ale block, used to
      configure arbitrary Lagrangian-Eulerian (ALE) mesh movement.)", "block-title"};

    option_t laplace {"laplace",
      "Select the Laplace mesh velocity smoother for ALE",
      R"(This keyword is used to select the 'Laplace' mesh velocity smoother for
      Arbitrary-Lagrangian-Eulerian (ALE) mesh motion.)", "string"};

    option_t helmholtz {"helmholtz",
      "Select the Helmholtz velocity for ALE",
      R"(This keyword is used to select the a velocity, computed from the
      Helmholtz-decomposition as the mesh velocity for
      Arbitrary-Lagrangian-Eulerian (ALE) mesh motion. See J. Bakosi, J. Waltz,
      N. Morgan, Improved ALE mesh velocities for complex flows, Int. J. Numer.
      Meth. Fl., 1-10, 2017, https://doi.org/10.1002/fld.4403.)", "string"};

    option_t none {"none", "Select none option",
      R"(This keyword is used to select the 'none' option from a list of
      configuration options.)", "string"};

    option_t sine {"sine",
      "Prescribe sinusoidal mesh velocity for ALE",
      R"(This keyword is used to prescribe a sinusoidal mesh velocity
      for Arbitrary-Lagrangian-Eulerian (ALE) mesh motion.)", "string"};

    option_t fluid {"fluid", "Select the fluid velocity for ALE",
      R"(This keyword is used to select the 'fluid' velocity as the mesh velocity
      for Arbitrary-Lagrangian-Eulerian (ALE) mesh motion.)", "string"};

    // -------------------------------------------------------------------------
    // h/p adaptation objects
    // -------------------------------------------------------------------------

    amr_t amr {"amr",
      "Start configuration block configuring adaptive mesh refinement",
      R"(This keyword is used to introduce the amr block, used to
      configure adaptive mesh refinement.)", "block-title"};

    pref_t pref {"pref",
      "Start configuration block configuring p-adaptive refinement",
      R"(This keyword is used to introduce the pref block, to
      configure p-adaptive refinement)", "block-title"};

    option_t amr_uniform {"uniform",
      "Select uniform initial mesh refinement",
      R"(This keyword is used to select uniform initial mesh refinement.)",
      "string"};

    option_t amr_uniform_derefine {"uniform_derefine",
      "Select uniform initial mesh de-refinement",
      R"(This keyword is used to select uniform initial mesh de-refinement.)",
      "string"};

    option_t amr_initial_conditions {"initial_conditions",
      "Select initial-conditions-based initial mesh refinement",
      R"(This keyword is used to select initial-conditions-based initial mesh
      refinement.)", "string"};

    option_t amr_jump {"jump",
      "Error estimation based on the solution jump normalized by solution value",
      R"(This keyword is used to select the jump-based error indicator for
      solution-adaptive mesh refinement. The error is estimated by computing the
      magnitude of the jump in the solution value normalized by the solution
      value.)", "string"};

    option_t amr_hessian {"hessian",
      "Error estimation based on the Hessian normalized by solution value",
      R"(This keyword is used to select the Hessian-based error indicator for
      solution-adaptive mesh refinement. The error is estimated by computing the
      Hessian (2nd derivative matrix) of the solution normalized by sum of the
      absolute values of the gradients at edges-end points.)", "string"};

    option_t pref_spectral_decay {"pref_spectral_decay",
      "Select the spectral-decay indicator for p-adaptive DG scheme",
      R"(This keyword is used to select the spectral-decay indicator used for
      p-adaptive discontinuous Galerkin (DG) discretization used in inciter.
      See Control/Inciter/Options/PrefIndicator.hpp for other valid options.)",
      "string"};

    option_t pref_non_conformity {"pref_non_conformity",
      "Select the non-conformity indicator for p-adaptive DG scheme",
      R"(This keyword is used to select the non-conformity indicator used for
      p-adaptive discontinuous Galerkin (DG) discretization used in inciter.
      See Control/Inciter/Options/PrefIndicator.hpp for other valid options.)",
      "string"};

    // -------------------------------------------------------------------------
    // boundary condition options
    // -------------------------------------------------------------------------

    std::vector< bc_t > bc;

    // -------------------------------------------------------------------------
    // IC object
    // -------------------------------------------------------------------------

    ic_t ic {"ic",
      "Introduce an ic block used to configure initial conditions",
      R"(This keyword is used to introduce an ic block used to set initial
      conditions.)", "block-title"};

    option_t impulse {"impulse",
      "Select the impulse initiation type, for a box/meshblock IC",
      R"(This keyword can be used to select the 'impulse' initiation/assignment
      type for box initial conditions. It simply assigns the prescribed values
      to mesh points within a configured box at t=0.)", "string"};

    option_t linear {"linear",
      "Select the linear initiation type, for a box/meshblock IC",
      R"(This keyword is be used to specify the 'linear' initiation parameters
      for a particular box or meshblock, as a part of the 'energy_pill'
      initialization. Linear initiation uses a linear function in time and space,
      configured with an initiation point in space, a constant velocity of the
      growing spherical front in time (and space) linearly, and width of the front
      and assigns values to mesh points falling within the growing spherical shell
      inside a configured box or meshblock. The following keywords are required
      in a linear block: 'init_time', 'front_width', 'velocity')",
      "block-title"};

    // -------------------------------------------------------------------------
    // Overset mesh object
    // -------------------------------------------------------------------------

    std::vector< mesh_t > mesh;

    // -------------------------------------------------------------------------
    // pre-configured problems
    // -------------------------------------------------------------------------

    option_t user_defined {"user_defined",
      "Select user-defined specification for a problem",
      R"(This keyword is used to select the user-defined specification for an
      option. This could be a 'problem' to be solved by a partial differential
      equation, but can also be a 'user-defined' mesh velocity specification for
      ALE mesh motion.)", "string"};

    option_t shear_diff {"shear_diff", 
      "Select the shear + diffusion test problem ",
      R"(This keyword is used to select the shear diffusion test problem. The
      initial and boundary conditions are specified to set up the test problem
      suitable to exercise and test the advection and diffusion terms of the
      scalar transport equation.)", "string" };

    option_t slot_cyl {"slot_cyl",
      "Select Zalesak's slotted cylinder test problem",
      R"(This keyword is used to select Zalesak's slotted cylinder test
      problem. The initial and boundary conditions are specified to set up the
      test problem suitable to exercise and test the advection and diffusion
      terms of the scalar transport equation.)", "string"};

    option_t gauss_hump {"gauss_hump",
      "Select advection of 2D Gaussian hump test problem",
      R"(This keyword is used to select the advection of 2D Gaussian hump test
      problem. The initial and boundary conditions are specified to set up the
      test problem suitable to exercise and test the advection
      terms of the scalar transport equation.)", "string"};

    option_t cyl_advect {"cyl_advect",
      "Select advection of cylinder test problem",
      R"(This keyword is used to select the advection of cylinder test
      problem. The initial and boundary conditions are specified to set up the
      test problem suitable to exercise and test the advection
      terms of the scalar transport equation.)", "string"};

    option_t cyl_vortex {"cyl_vortex",
      "Select deformation of cylinder in a vortex test problem",
      R"(This keyword is used to select the test problem which deforms a cylinder
      in a vortical velocity field. The initial and boundary conditions are
      specified to set up the test problem suitable to exercise and test the
      advection terms of the scalar transport equation.)", "string"};

    option_t vortical_flow {"vortical_flow",
      "Select the vortical flow test problem ",
      R"(This keyword is used to select the vortical flow test problem. The
      purpose of this test problem is to test velocity errors generated by spatial
      operators in the presence of 3D vorticity and in particluar the
      superposition of planar and vortical flows, analogous to voritcity
      stretching. For more details, see Waltz,
      et. al, "Manufactured solutions for the three-dimensional Euler equations
      with relevance to Inertial Confinement Fusion", Journal of Computational
      Physics 267 (2014) 196-209.)", "string"};

    option_t nl_energy_growth {"nl_energy_growth",
      "Select the nonlinear energy growth test problem",
      R"(This keyword is used to select the nonlinear energy growth test problem.
      The purpose of this test problem is to test nonlinear, time dependent energy
      growth and the subsequent development of pressure gradients due to coupling
      between the internal energy and the equation of state. For more details,
      see Waltz, et. al, "Manufactured
      solutions for the three-dimensional Euler equations with relevance to
      Inertial Confinement Fusion", Journal of Computational Physics 267 (2014)
      196-209.)", "string"};

    option_t rayleigh_taylor {"rayleigh_taylor",
      "Select the Rayleigh-Taylor test problem ",
      R"(This keyword is used to select the Rayleigh-Taylor unstable configuration
      test problem. The purpose of this test problem is to assess time dependent
      fluid motion in the presence of Rayleigh-Taylor unstable conditions, i.e.
      opposing density and pressure gradients.
      For more details, see Waltz, et. al, "Manufactured solutions for the
      three-dimensional Euler equations with relevance to Inertial Confinement
      Fusion", Journal of Computational Physics 267 (2014) 196-209.)",
      "string"};

    option_t taylor_green {"taylor_green",
      "Select the Taylor-Green test problem ",
      R"(This keyword is used to select the Taylor-Green vortex test problem. The
      purpose of this problem is to test time accuracy and the correctness of the
      discretization of the viscous term in the Navier-Stokes equation. For more
      details on the flow, see G.I. Taylor, A.E.
      Green, "Mechanism of the Production of Small Eddies from Large Ones", Proc.
      R. Soc. Lond. A 1937 158 499-521; DOI: 10.1098/rspa.1937.0036. Published 3
      February 1937.)", "string"};

    option_t sod_shocktube {"sod_shocktube",
      "Select the Sod shock-tube test problem ",
      R"(This keyword is used to select the Sod shock-tube test problem. The
      purpose of this test problem is to test the correctness of the
      approximate Riemann solver and its shock and interface capturing
      capabilities. For more details, see
      G. A. Sod, "A Survey of Several Finite Difference Methods for Systems of
      Nonlinear Hyperbolic Conservation Laws", J. Comput. Phys., 27 (1978)
      1–31.)", "string"};

    option_t rotated_sod_shocktube {"rotated_sod_shocktube",
      "Select the rotated Sod shock-tube test problem ",
      R"(This keyword is used to select the rotated Sod shock-tube test problem.
      This the same as Sod shocktube but the geometry is rotated about X, Y, Z
      each by 45 degrees (in that order) so that none of the domain boundary align
      with any of the coordinate directions. The purpose of this test problem is
      to test the correctness of the approximate Riemann solver and its shock and
      interface capturing capabilities in an arbitrarily oriented geometry.
      For more details on the Sod
      problem, see G. A. Sod, "A Survey of Several Finite Difference Methods for
      Systems of Nonlinear Hyperbolic Conservation Laws", J. Comput. Phys., 27
      (1978) 1–31.)", "string"};

    option_t shedding_flow {"shedding_flow",
      "Select the Shedding flow test problem ",
      R"(This keyword is used to select the Shedding flow test problem. It
      describe a quasi-2D inviscid flow over a triangular wedge in tetrahedron
      grid. The purpose of this test problem is to test the capability of DG
      scheme for retaining the shape of vortices and also different error
      indicator behavior for this external flow problem when p-adaptive DG scheme
      is applied.)", "string"};

    option_t sedov_blastwave {"sedov_blastwave",
      "Select the Sedov blast-wave test problem ",
      R"(This keyword is used to select the Sedov blast-wave test problem. The
      purpose of this test problem is to test the correctness of the
      approximate Riemann solver and its strong shock and interface capturing
      capabilities.)", "string"};

    option_t interface_advection {"interface_advection",
      "Select the interface advection test problem ",
      R"(This keyword is used to select the interface advection test problem. The
      purpose of this test problem is to test the well-balancedness of the
      multi-material discretization and its interface capturing
      capabilities.)", "string"};

    option_t gauss_hump_compflow {"gauss_hump_compflow",
      "Select advection of 2D Gaussian hump test problem",
      R"(This keyword is used to select the advection of 2D Gaussian hump test
      problem. The initial and boundary conditions are specified to set up the
      test problem suitable to exercise and test the advection terms of the
      Euler equations. The baseline of the density distribution in this testcase
      is 1 instead of 0 in gauss_hump_transport which enables it to be the
      regression testcase for p-adaptive DG scheme.)", "string"};

    option_t waterair_shocktube{" waterair_shocktube",
      "Select the water-air shock-tube test problem ",
      R"(This keyword is used to select the Water-air shock-tube test problem. The
      purpose of this test problem is to test the correctness of the
      multi-material pressure relaxation procedure and its interface capturing
      capabilities. For more details, see
      Chiapolino, A., Saurel, R., & Nkonga, B. (2017). Sharpening diffuse
      interfaces with compressible fluids on unstructured meshes. Journal of
      Computational Physics, 340, 389-417.)", "string"};

    option_t shock_hebubble {"shock_hebubble",
      "Select the shock He-bubble test problem ",
      R"(This keyword is used to select the shock He-bubble test problem. The
      purpose of this test problem is to test the correctness of the
      multi-material algorithm and its shock-interface interaction
      capabilities. For more details, see
      Quirk, J. J., & Karni, S. (1996). On the dynamics of a shock–bubble
      interaction. Journal of Fluid Mechanics, 318, 129-163.)", "string"};

    option_t underwater_ex {"underwater_ex",
      "Select the underwater explosion test problem ",
      R"(This keyword is used to select the underwater explosion test problem. The
      purpose of this test problem is to test the correctness of the
      multi-material algorithm and its interface capturing capabilities in the
      presence of strong shocks and large deformations.
      For more details, see
      Chiapolino, A., Saurel, R., & Nkonga, B. (2017). Sharpening diffuse
      interfaces with compressible fluids on unstructured meshes. Journal of
      Computational Physics, 340, 389-417.)", "string"};

    option_t shockdensity_wave {"shockdensity_wave",
      "Select the shock-density wave test problem ",
      R"(This keyword is used to select the shock-density wave test problem. The
      purpose of this test problem is to assess the accuracy of high order method
      in predicting the interaction of a density wave with a shock front.
      For more details, see Yu, L., Matthias
      I. (2014). Discontinuous Galerkin method for multicomponent chemically
      reacting flows and combustion. Journal of Computational Physics, 270,
      105-137.)", "string"};

    option_t equilinterface_advect {"equilinterface_advect",
      "Select the advection of equilibrium interface problem ",
      R"(This keyword is used to select the advection of equilibrium interface
      problem. This is a manufactured problem with source terms with nonlinear
      solutions near the material interface. Source terms are used to ensure that
      the conservation laws are satisfied by the manufactured solution.)",
      "string"};

    option_t sinewave_packet {"sinewave_packet",
      "Select the advection of sinewave packet problem ",
      R"(This keyword is used to select the advection of sinewave packet
      problem.)", "string"};

    option_t richtmyer_meshkov {"richtmyer_meshkov",
      "Select the Richtmyer-Meshkov instability problem ",
      R"(This keyword is used to select the Richtmyer-Meshkov instability
      problem. In this problem, a shock hits a perturbed material interface.)",
      "string"};

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | inciter;
      p | title;
      p | nstep;
      p | term;
      p | t0;
      p | dt;
      p | cfl;
      p | ttyi;
      p | steady_state;
      p | residual;
      p | rescomp;
      p | partitioning;
      p | pelocal_reorder;
      p | operator_reorder;
      p | scheme;
      p | ndof;
      p | rdof;
      p | limiter;
      p | cweight;
      p | shock_detector_coeff;
      p | accuracy_test;
      p | limsol_projection;
      p | fct;
      p | fctclip;
      p | fcteps;
      p | ctau;
      p | sysfct;
      p | sysfctvar;
      p | flux;
      p | pde;
      p | transport;
      p | compflow;
      p | multimat;
      p | depvar;
      p | physics;
      p | ncomp;
      p | material;
      p | matidxmap;
      p | field_output;
      p | diagnostics;
      p | history_output;
      p | ale;
      p | amr;
      p | pref;
      p | bc;
      p | ic;
      p | mesh;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i NewInputDeck object reference
    friend void operator|( PUP::er& p, NewInputDeck& i ) { i.pup(p); }
    //@}
};

} // ctr::
} // inciter::

#endif // NewInputDeck_h
