/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef BELOS_PLAYA_ADAPTER_HPP
#define BELOS_PLAYA_ADAPTER_HPP


#include "PlayaDefs.hpp"
#include "PlayaAnasaziAdapter.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"

namespace Belos
{


using Playa::Vector;
using Teuchos::RCP;
using Teuchos::Array;
using Anasazi::SimpleMV;

/** */
template<> 
class MultiVecTraits<double, SimpleMV>
{
public:
  typedef SimpleMV _MV;
  typedef Teuchos::ScalarTraits<double> SCT;
  typedef Anasazi::MultiVecTraits<double, _MV> AMVT;

  static double one() {static double rtn = SCT::one(); return rtn;}
  static double zero() {static double rtn = SCT::zero(); return rtn;}
  
  /** */
  static RCP<_MV> Clone( const  _MV & mv, const int numvecs )
    {return AMVT::Clone(mv, numvecs);}

  /** */
  static RCP< _MV > CloneCopy( const  _MV & mv )
    {return AMVT::CloneCopy(mv);}

  /** */
  static RCP< _MV > CloneCopy( const  _MV & mv, const std::vector<int>& index )
    {return AMVT::CloneCopy(mv, index);}

  /** */     
  static RCP< _MV > CloneViewNonConst(  _MV & mv, 
    const std::vector<int>& index )
    {return AMVT::CloneViewNonConst(mv, index);}   

  /** */     
  static RCP<_MV> 
  CloneViewNonConst (_MV & mv, 
		     const Teuchos::Range1D& index)
  {
    return AMVT::CloneViewNonConst (mv, index);
  }

  /** */
  static RCP<const _MV > CloneView( const _MV & mv, const std::vector<int>& index )
    {return AMVT::CloneView(mv, index);}

  /** */
  static RCP<const _MV > CloneView( const _MV & mv, const Teuchos::Range1D& index )
    {return AMVT::CloneView(mv, index);}
  
  /** Obtain the vector length of \c mv. */
  static int GetVecLength( const  _MV & mv )
    {return AMVT::GetVecLength(mv);}

  /** Obtain the number of vectors in \c mv */
  static int GetNumberVecs( const  _MV & mv )
    {return AMVT::GetNumberVecs(mv);}


  /*! \brief Update \c mv with \f$ \alpha A B + \beta mv \f$ */
  static void MvTimesMatAddMv( const double alpha, const  _MV & A, 
    const Teuchos::SerialDenseMatrix<int,double>& B, 
    const double beta,  _MV & mv )
    {AMVT::MvTimesMatAddMv(alpha, A, B, beta, mv);}

  
  /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$. */
  static void MvAddMv( const double alpha, const  _MV & A, 
    const double beta,  const  _MV & B,  _MV & mv )
    {AMVT::MvAddMv(alpha, A, beta, B, mv);}

  /*! \brief Compute a dense matrix \c B through 
   * the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
   */
  static void MvTransMv( const double alpha, const  _MV & A, const  _MV & mv, 
    Teuchos::SerialDenseMatrix<int,double>& B )
    {AMVT::MvTransMv(alpha, A, mv, B);}

  
  /**
   * Dot product
  */
  static void MvDot( const  _MV & mv, const  _MV & A, std::vector<double> &b )
    {AMVT::MvDot(mv, A, b);}

  
  /** Scale each element of the vectors in \c *this with \c alpha.
   */
  static void MvScale (  _MV & mv, const double alpha )
    {AMVT::MvScale(mv, alpha);}
      

    
  /** Scale each element of the \c i-th vector in \c *this with \c alpha[i].
   */
  static void MvScale (  _MV & mv, const std::vector<double>& alpha ) 
    {AMVT::MvScale(mv, alpha);}
    
  /** Compute the 2-norm of each individual vector of \c mv. */
  static void MvNorm( const  _MV & mv, 
    std::vector<Teuchos::ScalarTraits<double>::magnitudeType> &normvec,
    NormType type = TwoNorm)
    {
      normvec.resize(mv.size());
      for (int i=0; i<mv.size(); i++)
      {
        if (type==OneNorm) 
          normvec[i] = mv[i].norm1();
        else if (type==TwoNorm) 
          normvec[i] = mv[i].norm2();
        else 
          normvec[i] = mv[i].normInf();
      }
    }

  /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
   */
  static void SetBlock( const  _MV & A, const std::vector<int>& index,  
    _MV & mv )
    {AMVT::SetBlock(A, index, mv);}

  //! Assign (deep copy) A into mv.
  static void Assign( const _MV& A, _MV& mv ) {
    AMVT::Assign (A, mv);
  }

  /*! \brief Replace the vectors in \c mv with random vectors.
   */
  static void MvRandom(  _MV & mv )
    {AMVT::MvRandom(mv);}

  /*! \brief Replace each element of the vectors in \c mv with \c alpha.
   */
  static void MvInit(  _MV & mv, 
    double alpha = Teuchos::ScalarTraits<double>::zero() )
    { AMVT::MvInit(mv, alpha);}
    
  /** Print the \c mv multi-vector to the \c os output stream. */
  static void MvPrint( const  _MV & mv, std::ostream& os )
    { AMVT::MvPrint(mv, os);}

#ifdef HAVE_BELOS_TSQR
  /// \typedef tsqr_adaptor_type
  /// \brief TSQR adapter for the multivector type SimpleMV.
  ///
  /// For now, we provide a "stub" implementation.  It has the right
  /// methods and typedefs, but its constructors and methods all throw
  /// std::logic_error.  If you plan to use TSQR in Belos (e.g.,
  /// through TsqrOrthoManager) with SimpleMV, you must implement a
  /// functional TSQR adapter for SimpleMV.  Please refer to
  /// Epetra::TsqrAdapter (for Epetra_MultiVector) or
  /// Tpetra::TsqrAdaptor (for Tpetra::MultiVector) for examples of
  /// how to implement a TSQR adapter.
  typedef Belos::details::StubTsqrAdapter<SimpleMV> tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
};


/**

*/
template <> 
class OperatorTraits < double, SimpleMV, LinearOperator<double> >
{
public:
  typedef Anasazi::SimpleMV _MV;  
  typedef Anasazi::OperatorTraits<double, SimpleMV, LinearOperator<double> > AOPT;
  /**
  */    
  static void Apply ( 
    const LinearOperator< double >& Op, 
    const  _MV & x,  
    _MV & y,
    ETrans trans=NOTRANS
    )
    {
      if (trans==NOTRANS) AOPT::Apply(Op, x, y);
      else AOPT::Apply(Op.transpose(), x, y);
    }
    
};



}

#endif
