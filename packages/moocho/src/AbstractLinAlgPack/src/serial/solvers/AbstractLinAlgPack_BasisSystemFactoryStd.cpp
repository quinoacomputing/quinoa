// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "AbstractLinAlgPack_BasisSystemFactoryStd.hpp"
#include "AbstractLinAlgPack_BasisSystemPermDirectSparse.hpp"
#include "AbstractLinAlgPack_DirectSparseSolverDense.hpp"
#include "Teuchos_Assert.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "OptionsFromStreamPack_StringToIntMap.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"

#ifdef HAVE_MOOCHO_MA28
#include "AbstractLinAlgPack_DirectSparseSolverMA28.hpp"
#include "AbstractLinAlgPack_DirectSparseSolverMA28SetOptions.hpp"
#endif

namespace AbstractLinAlgPack {

BasisSystemFactoryStd::BasisSystemFactoryStd()
  :direct_linear_solver_type_(
#ifdef SPARSE_SOLVER_PACK_USE_MA48
    LA_MA48                        // If we have MA48 use it as a first choice
#else
#  ifdef HAVE_MOOCHO_MA28
    LA_MA28                        // If we have MA28 use it as a second choice
#  else
    LA_DENSE                       // If we don't have any sparse solvers use dense
#  endif
#endif
    )
{}

// Overridden from BasisSystemFactory

void BasisSystemFactoryStd::set_options( const options_ptr_t& options )
{
  options_ = options;
}

const BasisSystemFactoryStd::options_ptr_t&
BasisSystemFactoryStd::get_options() const
{
  return options_;
}

// Overridden from AbstractFactory

BasisSystemFactoryStd::obj_ptr_t
BasisSystemFactoryStd::create() const
{
  namespace mmp = MemMngPack;;

  // Read in the options
  read_options();

  // Create the direct sparse solver
  Teuchos::RCP<DirectSparseSolver>  direct_sparse_solver;
  switch(direct_linear_solver_type_) {
    case LA_DENSE: {
      Teuchos::RCP<DirectSparseSolverDense>
        dss_dense = Teuchos::rcp(new DirectSparseSolverDense());
      direct_sparse_solver = dss_dense;
      break;
    }
    case LA_MA28: {
#ifdef HAVE_MOOCHO_MA28
      Teuchos::RCP<DirectSparseSolverMA28>
        dss_ma28 = Teuchos::rcp(new DirectSparseSolverMA28());
      if(options_.get()) {
        AbstractLinAlgPack::DirectSparseSolverMA28SetOptions
          opt_setter(dss_ma28.get());
        opt_setter.set_options(*options_);
      }
      direct_sparse_solver = dss_ma28;
#else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error
        ,"Error, HAVE_MOOCHO_MA28 is not defined and therefore MA28 is not supported!" );
#endif
      break;
    }
    case LA_MA48: {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error
        ,"Error, MA48 is not supported yet!" );
      break;
    }
    case LA_SUPERLU: {
#ifdef SPARSE_SOLVER_PACK_USE_SUPERLU
      Teuchos::RCP<DirectSparseSolverSuperLU>
        dss_slu = Teuchos::rcp(new DirectSparseSolverSuperLU());
      // ToDo: Set options from stream!
      direct_sparse_solver = dss_slu;
#else
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error
        ,"Error, SPARSE_SOLVER_PACK_USE_SUPERLU is not defined and therefore SuperLU is not supported!" );
#endif
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should not be called?
  }

  // Return the basis system
  return Teuchos::rcp(new BasisSystemPermDirectSparse(direct_sparse_solver));

}

// private

void BasisSystemFactoryStd::read_options() const
{
  namespace	ofsp = OptionsFromStreamPack;
  using		ofsp::OptionsFromStream;
  typedef		OptionsFromStream::options_group_t		options_group_t;
  using		ofsp::StringToIntMap;
  using		ofsp::StringToBool;

  if(!options_.get())
    return;

  const std::string opt_grp_name = "BasisSystemFactoryStd";
  const OptionsFromStream::options_group_t optgrp = options_->options_group( opt_grp_name );
  if( OptionsFromStream::options_group_exists( optgrp ) ) {

    const int num_opts = 1;
    enum EBasisSystemFactorStd {
      DIRECT_LINEAR_SOLVER
    };
    const char* SBasisSystemFactorStd[num_opts]	= {
      "direct_linear_solver"
    };
    StringToIntMap	map( opt_grp_name, num_opts, SBasisSystemFactorStd );

    options_group_t::const_iterator itr = optgrp.begin();
    for( ; itr != optgrp.end(); ++itr ) {
      switch( (EBasisSystemFactorStd)map( ofsp::option_name(itr) ) ) {
        case DIRECT_LINEAR_SOLVER:
        {
          const std::string &linear_solver = ofsp::option_value(itr);
          if( linear_solver == "DENSE" ) {
            direct_linear_solver_type_ = LA_DENSE;
          } else if( linear_solver == "MA28" ) {
#ifdef HAVE_MOOCHO_MA28
            direct_linear_solver_type_ = LA_MA28;
#else
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error
              ,"BasisSystemFactoryStd::read_options(...) : MA28 is not supported,"
              " you must configure with --enable-moocho-ma28!" );
#endif
          } else if( linear_solver == "MA48" ) {
#ifdef SPARSE_SOLVER_PACK_USE_MA48
            direct_linear_solver_type_ = LA_MA48;
#else
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error
              ,"BasisSystemFactoryStd::read_options(...) : MA48 is not supported,"
              " must define SPARSE_SOLVER_PACK_USE_MA48!" );
#endif
          } else if( linear_solver == "SUPERLU" ) {
#ifdef SPARSE_SOLVER_PACK_USE_SUPERLU
            direct_linear_solver_type_ = LA_SUPERLU;
#else
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error
              ,"BasisSystemFactoryStd::read_options(...) : SUPERLU is not supported,"
              " must define SPARSE_SOLVER_PACK_USE_SUPERLU!" );
#endif
          } else {
            TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::invalid_argument
              ,"BasisSystemFactoryStd::read_options(...) : "
              "Error, incorrect value for \"direct_linear_solver\" "
              "Only the options \'DENSE\', \'MA28\' and \'SUPERLU\' are avalible." );
          }
          break;
        }
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true);	// this would be a local programming error only.
      }
    }
  }
  else {
    // Warning, options group was not found!!!
  }
  
}

}  // end namespace AbstractLinAlgPack
