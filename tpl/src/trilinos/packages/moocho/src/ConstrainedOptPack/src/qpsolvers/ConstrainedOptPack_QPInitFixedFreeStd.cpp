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

#include <assert.h>

#include <sstream>

#include "ConstrainedOptPack_QPInitFixedFreeStd.hpp"
#include "ConstrainedOptPack_initialize_Q_R_Q_X.hpp"

namespace ConstrainedOptPack {
namespace QPSchurPack {

QPInitFixedFreeStd::QPInitFixedFreeStd()
  :n_(0)
  ,n_R_(0)
  ,m_(0)
  ,G_(NULL)
  ,A_(NULL)
  ,Ko_(NULL)
  ,constraints_(NULL)
{}

void QPInitFixedFreeStd::initialize(
  const DVectorSlice                   &g
  ,const MatrixSymOp              &G
  ,const MatrixOp                 *A
  ,size_type                          n_R
  ,const size_type                    i_x_free[]
  ,const size_type                    i_x_fixed[]
  ,const EBounds                      bnd_fixed[]
  ,const DVectorSlice                  &b_X
  ,const MatrixSymOpNonsing    &Ko
  ,const DVectorSlice                  &fo
  ,Constraints                        *constraints
  ,std::ostream                       *out
  ,bool                               test_setup
  ,value_type                         warning_tol
  ,value_type                         error_tol
  ,bool                               print_all_warnings
  )
{
  namespace GPMSTP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;

  if(!constraints)
    throw std::invalid_argument( "QPInitFixedFreeStd::initialize(...) : Error, "
      "constraints == NULL is not allowed." );

  // Validate the sizes of the input arguments
  const size_type
    n = constraints->n(),
    n_X = n - n_R;
  if( n_R > n )
    throw std::invalid_argument( "QPInitFixedFreeStd::initialize(...) : Error, "
      "n_R > constraints->n() is not allowed." );
  if(g.dim() !=n)
    throw std::invalid_argument( "QPInitFixedFreeStd::initialize(...) : Error, "
      "g.dim() != constraints->n()." );
  if(G.rows() != n || G.cols() !=  n)
    throw std::invalid_argument( "QPInitFixedFreeStd::initialize(...) : Error, "
      "G.rows() != constraints->n() or G.cols() !=  constraints->n()." );
  size_type
    m = 0;
  if(A) {
    m = A->cols();
    if( A->rows() != n )
      throw std::invalid_argument( "QPInitFixedFreeStd::initialize(...) : Error, "
        "A->rows() != constraints->n()." );
  }
  if(b_X.dim() != n_X)
    throw std::invalid_argument( "QPInitFixedFreeStd::initialize(...) : Error, "
      "b_X.dim() != constraints->n() - n_R." );
  if(Ko.rows() != n_R+m || Ko.cols() !=  n_R+m)
    throw std::invalid_argument( "QPInitFixedFreeStd::initialize(...) : Error, "
      "Ko.rows() != n_R+A->cols() or Ko.cols() !=  n_R+A->cols()." );
  if(fo.dim() != n_R+m)
    throw std::invalid_argument( "QPInitFixedFreeStd::initialize(...) : Error, "
      "fo.dim() != n_R+A->cols()." );

  // Setup x_init, l_x_X_map, i_x_X_map

  const int NOT_SET_YET = -9999; // Can't be FREE, LOWER, UPPER, or EQUALITY
  if(test_setup)
    x_init_.assign( n, (EBounds)NOT_SET_YET );
  else
    x_init_.resize(n);
  l_x_X_map_.assign(n,0);
  i_x_X_map_.assign(n_X,0);

  // Set free portion of x_init
  if( i_x_free == NULL ) {
    for( size_type i = 0; i < n_R; ++i )
      x_init_[i] = FREE;
  }
  else {
    if(test_setup) {
      for( const size_type *i_x_R = i_x_free; i_x_R != i_x_free + n_R; ++i_x_R ) {
        if( *i_x_R < 1 || *i_x_R > n ) {
          std::ostringstream omsg;
          omsg
            << "QPInitFixedFreeStd::initialize(...) : Error, "
            << "i_x_free[" << i_x_R-i_x_free << "] = "
            << (*i_x_R) << " is out of bounds";
          throw std::invalid_argument( omsg.str() );
        }
        if( x_init_(*i_x_R) != NOT_SET_YET ) {
          std::ostringstream omsg;
          omsg
            << "QPInitFixedFreeStd::initialize(...) : Error, "
            << "Duplicate entries for i_x_free[i] = "
            << (*i_x_R);
          throw std::invalid_argument( omsg.str() );
        }
        x_init_(*i_x_R) = FREE;
      }
    }
    else {
      for( const size_type *i_x_R = i_x_free; i_x_R != i_x_free + n_R; ++i_x_R ) {
        x_init_(*i_x_R) = FREE;
      }
    }
  }

  // Setup the fixed portion of x_init and l_x_X_map and i_x_X_map
  {
    const size_type
      *i_x_X = i_x_fixed;
    const EBounds
      *bnd = bnd_fixed;
    i_x_X_map_t::iterator
      i_x_X_map_itr = i_x_X_map_.begin();
    if(test_setup) {
      for( size_type l = 1; l <= n_X; ++l, ++i_x_X, ++bnd, ++i_x_X_map_itr ) {
        if( *i_x_X < 1 || *i_x_X > n ) {
          std::ostringstream omsg;
          omsg
            << "QPInitFixedFreeStd::initialize(...) : Error, "
            << "i_x_fixed[" << i_x_X-i_x_fixed << "] = "
            << (*i_x_X) << " is out of bounds";
          throw std::invalid_argument( omsg.str() );
        }
        if( *bnd == FREE ) {
          std::ostringstream omsg;
          omsg
            << "QPInitFixedFreeStd::initialize(...) : Error, "
            << "bnd_fixed[" << l-1 << "] can not equal FREE";
          throw std::invalid_argument( omsg.str() );
        }
        if( x_init_(*i_x_X) != NOT_SET_YET ) {
          std::ostringstream omsg;
          omsg
            << "QPInitFixedFreeStd::initialize(...) : Error, "
            << "Duplicate entries for i_x_fixed[i] = "
            << (*i_x_X);
          throw std::invalid_argument( omsg.str() );
        }
        x_init_(*i_x_X) = *bnd;
        l_x_X_map_(*i_x_X) = l;
        *i_x_X_map_itr = *i_x_X;
      }
      // Check that x_init was filled up properly
      for( size_type i = 1; i <= n; ++i ) {
        if( x_init_(i) == NOT_SET_YET ) {
          std::ostringstream omsg;
          omsg
            << "QPInitFixedFreeStd::initialize(...) : Error, "
            << "x_init(" << i << ") has not been set by"
               " i_x_free[] or i_x_fixed[].";
          throw std::invalid_argument( omsg.str() );
        }
      }
    }
    else {
      for( size_type l = 1; l <= n_X; ++l, ++i_x_X, ++bnd, ++i_x_X_map_itr ) {
        x_init_(*i_x_X) = *bnd;
        l_x_X_map_(*i_x_X) = l;
        *i_x_X_map_itr = *i_x_X;
      }
    }
  }

  // Setup Q_R and Q_X
  Q_R_row_i_.resize( i_x_free ? n_R : 0 );
  Q_R_col_j_.resize( i_x_free ? n_R : 0 );
  Q_X_row_i_.resize(n_X);
  Q_X_col_j_.resize(n_X);
  initialize_Q_R_Q_X(
    n_R,n_X,i_x_free,i_x_fixed,test_setup
    ,n_R && i_x_free ? &Q_R_row_i_[0] : NULL
    ,n_R && i_x_free ? &Q_R_col_j_[0] : NULL
    ,&Q_R_
    ,n_X ? &Q_X_row_i_[0] : NULL
    ,n_X ? &Q_X_col_j_[0] : NULL
    ,&Q_X_
  );

  // Setup other arguments
  n_				= n;
  n_R_			= n_R;
  m_				= m;
  g_.bind(const_cast<DVectorSlice&>(g));
  G_				= &G;
  A_				= A;
  b_X_.bind(const_cast<DVectorSlice&>(b_X));
  Ko_				= &Ko;
  fo_.bind(const_cast<DVectorSlice&>(fo));
  constraints_ 	= constraints;
}

// Overridden from QP 

size_type QPInitFixedFreeStd::n() const
{
  assert_initialized();
  return n_;
}

size_type QPInitFixedFreeStd::m() const
{
  assert_initialized();
  return m_;
}

const DVectorSlice QPInitFixedFreeStd::g() const
{
  assert_initialized();
  return g_;
}

const MatrixSymOp& QPInitFixedFreeStd::G() const
{
  assert_initialized();
  return *G_;
}

const MatrixOp& QPInitFixedFreeStd::A() const
{
  TEUCHOS_TEST_FOR_EXCEPT( !( A_ ) );	// ToDo: make this throw an exception
  return *A_;
}

size_type QPInitFixedFreeStd::n_R() const
{
  assert_initialized();
  return n_R_;
}

const QP::x_init_t& QPInitFixedFreeStd::x_init() const
{
  assert_initialized();
  return x_init_;
}

const QP::l_x_X_map_t& QPInitFixedFreeStd::l_x_X_map() const
{
  assert_initialized();
  return l_x_X_map_;
}

const QP::i_x_X_map_t& QPInitFixedFreeStd::i_x_X_map() const
{
  assert_initialized();
  return i_x_X_map_;
}

const DVectorSlice QPInitFixedFreeStd::b_X() const
{
  assert_initialized();
  return b_X_;
}

const GenPermMatrixSlice& QPInitFixedFreeStd::Q_R() const
{
  assert_initialized();
  return Q_R_;
}

const GenPermMatrixSlice& QPInitFixedFreeStd::Q_X() const
{
  assert_initialized();
  return Q_X_;
}

const MatrixSymOpNonsing& QPInitFixedFreeStd::Ko() const
{
  assert_initialized();
  return *Ko_;
}

const DVectorSlice QPInitFixedFreeStd::fo() const
{
  assert_initialized();
  return fo_;
}

Constraints& QPInitFixedFreeStd::constraints()
{
  assert_initialized();
  return *constraints_;
}

const Constraints& QPInitFixedFreeStd::constraints() const
{
  assert_initialized();
  return *constraints_;
}

// private member functions

void QPInitFixedFreeStd::assert_initialized() const
{
  if( !n_ )
    throw std::logic_error( "QPInitFixedFreeStd::assert_initialized(), Error "
      "object not initialized\n" );
}

}	// end namespace QPSchurPack
}	// end namespace ConstrainedOptPack 
