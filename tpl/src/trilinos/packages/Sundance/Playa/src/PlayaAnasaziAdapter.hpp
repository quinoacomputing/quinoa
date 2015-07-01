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

#ifndef ANASAZI_PLAYA_ADAPTER_HPP
#define ANASAZI_PLAYA_ADAPTER_HPP

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziConfigDefs.hpp"


#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaSimpleTransposedOpImpl.hpp"
#include "PlayaLinearCombinationImpl.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#include "PlayaTabs.hpp"

namespace Anasazi 
{
using Playa::Vector;
using Teuchos::RCP;
using Teuchos::Array;

class SimpleMV
{
public:
  SimpleMV() : data_() {}

  SimpleMV(int n) : data_(rcp(new Array<Vector<double> >(n))) {}

  SimpleMV(const Array<Vector<double> >& data) 
    : data_(rcp(new Array<Vector<double> >(data.size())))
    {
      for (int i=0; i<data.size(); i++)
      {
        (*data_)[i] = data[i].copy();
      }
    }

  SimpleMV clone() const
    {
      return SimpleMV(*data_);
    }

  Vector<double>& operator[](int i) {return (*data_)[i];}

  const Vector<double>& operator[](int i) const {return (*data_)[i];}

  int size() const {return data_->size();}

  void resize(int n)
    {
      data_->resize(n);
    }

  void randomize() 
    {
      for (int i=0; i<size(); i++) (*data_)[i].randomize();
    }
  
private:
  RCP<Array<Vector<double> > > data_;
};



inline std::ostream& operator<<(std::ostream& os, const SimpleMV& mv)
{
  os << "MV (size=" << mv.size() << ")" << std::endl;
  for (int i=0; i<mv.size(); i++)
  {
    os << "ptr=" << mv[i].ptr().get() << std::endl;
    os << mv[i] << std::endl;
  }

  return os;
}

/** */
template<>
class MultiVecTraits<double, SimpleMV >
{
public:
  typedef SimpleMV _MV;
  typedef Teuchos::ScalarTraits<double> SCT;

  static double one() {static double rtn = SCT::one(); return rtn;}
  static double zero() {static double rtn = SCT::zero(); return rtn;}

  /** \name Creation methods */
  //@{

  /**
   */
  static RCP<_MV> Clone( const  _MV & mv, const int numvecs )
    { 
      //Out::os() << "Clone(nv == " << numvecs << ")" << endl;
      TEUCHOS_TEST_FOR_EXCEPT(mv.size() <= 0);
      TEUCHOS_TEST_FOR_EXCEPT(numvecs <= 0);

      RCP<_MV> rtn = rcp(new _MV(numvecs));
      for (int i=0; i<numvecs; i++)
      {
        (*rtn)[i] = mv[0].copy();
        (*rtn)[i].setToConstant(zero());
      }
      return rtn;
    }

  /**
   *
   */
  static RCP< _MV > CloneCopy( const  _MV & mv )
    { 
      //Out::os() << "CloneCopy()" << endl;
      int numvecs = mv.size();
      TEUCHOS_TEST_FOR_EXCEPT(numvecs <= 0);

      // create the new multivector
      RCP<_MV> rtn = rcp(new _MV(numvecs));
      for (int i=0; i<numvecs; i++)
      {
        (*rtn)[i] = mv[i].copy();
      }
      return rtn;
    }

  /** 
      
  */
  static RCP< _MV > CloneCopy( const  _MV & mv, const std::vector<int>& index )
    { 
      //Out::os() << "CloneCopy() indexed" << endl;
      int numvecs = index.size();
      TEUCHOS_TEST_FOR_EXCEPT(numvecs <= 0);
      TEUCHOS_TEST_FOR_EXCEPT((int) index.size() > mv.size());

//      TEUCHOS_TEST_FOR_EXCEPT(detectRepeatedIndex(index));

      // create the new multivector
      RCP<  _MV  > rtn = rcp(new _MV(numvecs));

      for (int i=0; i<numvecs; i++)
      {
        (*rtn)[i] = mv[index[i]].copy();
      }
      return rtn;
    }

  /**

  */      
  static RCP< _MV > CloneViewNonConst(  _MV & mv, const std::vector<int>& index )
    {
      int numvecs = index.size();
      //Out::os() << "index.size() = " << numvecs << endl;
      //Out::os() << "input size = " << mv.size() << endl;
      TEUCHOS_TEST_FOR_EXCEPT(numvecs <= 0);
      TEUCHOS_TEST_FOR_EXCEPT((int) index.size() > mv.size());

//      TEUCHOS_TEST_FOR_EXCEPT(detectRepeatedIndex(index));

      // create the new multivector
      RCP<  _MV  > rtn = rcp(new _MV(numvecs));

      for (int i=0; i<numvecs; i++)
      {
        (*rtn)[i] = mv[index[i]]; // shallow copy
      }

      return rtn;
    }

  static RCP<_MV> 
  CloneViewNonConst (_MV & mv, const Teuchos::Range1D& index) 
  {
    typedef std::vector<int>::size_type size_type;
    std::vector<int> ind (index.size ());
    for (size_type i = 0; i < static_cast<size_type> (index.size ()); ++i) {
      ind[i] = index.lbound () + i;
    }
    return CloneViewNonConst (mv, ind);
  }

  /**
   *
   */      
  static RCP<const _MV > CloneView( const _MV & mv, const std::vector<int>& index )
    {
      //Out::os() << "CloneView()" << endl;
      int numvecs = index.size();
      //Out::os() << "index size = " << numvecs << endl;
      //Out::os() << "input size = " << mv.size() << endl;
      TEUCHOS_TEST_FOR_EXCEPT(numvecs <= 0);
      TEUCHOS_TEST_FOR_EXCEPT((int) index.size() > mv.size());

//      TEUCHOS_TEST_FOR_EXCEPT(detectRepeatedIndex(index));

      // create the new multivector
      RCP<  const _MV  > rtn = rcp(new _MV(numvecs));

      for (int i=0; i<numvecs; i++)
      {
        (*(rcp_const_cast<_MV>(rtn)))[i] = mv[index[i]]; // shallow copy
      }
      return rtn;
    }

  static RCP<const _MV> 
  CloneView (const _MV & mv, const Teuchos::Range1D& index) 
  {
    typedef std::vector<int>::size_type size_type;
    std::vector<int> ind (index.size ());
    for (size_type i = 0; i < static_cast<size_type> (index.size ()); ++i) {
      ind[i] = index.lbound () + i;
    }
    return CloneView (mv, ind);
  }

  //@}

  /** \name Attribute methods */
  //@{

  /** Obtain the vector length of \c mv. */
  static int GetVecLength( const  _MV & mv )
    {
      TEUCHOS_TEST_FOR_EXCEPT(mv.size() <= 0);
      return mv[0].space().dim(); 
    }

  /** Obtain the number of vectors in \c mv */
  static int GetNumberVecs( const  _MV & mv )
    {
      //Out::os() << "GetNumVec(" << mv.size() << ")" << endl;
      return mv.size(); 
    }

  //@}

  /** \name Update methods */
  //@{

  /*! \brief Update \c mv with \f$ \alpha A B + \beta mv \f$.
   */
  static void MvTimesMatAddMv( const double alpha, const  _MV & A, 
    const Teuchos::SerialDenseMatrix<int,double>& B, 
    const double beta,  _MV & mv )
    {
//      Out::os() << "MvTimesMatAddMv()" << endl;
      int n = B.numCols();
//      Out::os() << "B.numCols()=" << n << endl;

      TEUCHOS_TEST_FOR_EXCEPT(mv.size() != n);

      for (int j=0; j<mv.size(); j++)
      {
        Vector<double> tmp;
        if (beta==one())
        {
          tmp = mv[j].copy();
        }
        else if (beta==zero())
        {
          tmp = mv[j].copy();
          tmp.setToConstant(zero());
        }
        else
        {
          tmp = beta * mv[j];
        }
        if (alpha != zero())
        {
          for (int i=0; i<A.size(); i++)
          {
            tmp = tmp + alpha*B(i,j)*A[i];
          }
        }
        mv[j].acceptCopyOf(tmp);
      }
    }

  /*! \brief Replace \c mv with \f$\alpha A + \beta B\f$.
   */
  static void MvAddMv( const double alpha, const  _MV & A, 
    const double beta,  const  _MV & B,  _MV & mv )
    { 
      TEUCHOS_TEST_FOR_EXCEPT(A.size() != B.size());
      mv.resize(A.size());
      for (int i=0; i<A.size(); i++)
      {
        if (alpha==zero() && beta != zero()) mv[i]=beta*B[i];
        else if (beta==zero() && alpha != zero()) mv[i]=alpha*A[i];
        else if (alpha!=zero() && beta!=zero())
          mv[i]= alpha*A[i] + beta*B[i] ;
        else
        {
          mv[i].acceptCopyOf(A[i]);
          mv[i].setToConstant(zero());
        }
      }
    }

  /*! \brief Compute a dense matrix \c B through the matrix-matrix multiply \f$ \alpha A^Tmv \f$.
   */
  static void MvTransMv( const double alpha, const  _MV & A, const  _MV & mv, 
    Teuchos::SerialDenseMatrix<int,double>& B )
    { 
      // Create a multivector to hold the result (m by n)
      int m = A.size();
      int n = mv.size();
//      B.shape(m, n);
      //Out::os() << "m=" << m << ", n=" << n << endl;
      for (int i=0; i<m; i++)
      {
        for (int j=0; j<n; j++)
        {
          B(i,j) = alpha * (A[i] * mv[j]);
        }
      }
    
    }

  /**
   * Dot product
  */
  static void MvDot( const  _MV & mv, const  _MV & A, std::vector<double> &b )
    {
      //Out::os() << "MvDot()" << endl;
      TEUCHOS_TEST_FOR_EXCEPT(mv.size() != A.size());
      b.resize(A.size());
      for (int i=0; i<mv.size(); i++) 
        b[i] = mv[i] * A[i];
    }

  /** Scale each element of the vectors in \c *this with \c alpha.
   */
  static void MvScale (  _MV & mv, const double alpha )
    { 
      //Out::os() << "MvScale()" << endl;
      for (int i=0; i<mv.size(); i++) mv[i].scale(alpha);
    }
    
  /** Scale each element of the \c i-th vector in \c *this with \c alpha[i].
   */
  static void MvScale (  _MV & mv, const std::vector<double>& alpha ) 
    { 
      //Out::os() << "MvScale() vector" << endl;
      for (int i=0; i<mv.size(); i++) mv[i].scale(alpha[i]);
    }
    
  //@}

  /** \name Norm method */
  //@{

  /** Compute the 2-norm of each individual vector of \c mv. */
  static void MvNorm( const  _MV & mv, 
    std::vector<Teuchos::ScalarTraits<double>::magnitudeType> &normvec )
    { 
//      Out::os() << "MvNorm()" << endl;
      normvec.resize(mv.size());
      for (int i=0; i<mv.size(); i++) 
      {
        normvec[i] = mv[i].norm2();
        //      Out::os() << "i=" << i << " |v|=" << normvec[i] << endl;
      }
      
    }

  //@}

  /** \name Initialization methods */
  //@{

  /*! \brief Copy the vectors in \c A to a set of vectors in \c mv indicated by the indices given in \c index.
   */
  static void SetBlock( const  _MV & A, const std::vector<int>& index,  _MV & mv )
    { 
      //Out::os() << "SetBlock()" << endl;
      TEUCHOS_TEST_FOR_EXCEPT(A.size() < (int) index.size());
//      mv.resize(index.size());
//      TEUCHOS_TEST_FOR_EXCEPT(detectRepeatedIndex(index));
      for (unsigned int i=0; i<index.size(); i++)
      {
        mv[index[i]].acceptCopyOf(A[i]);
      }
    }

  //! Replace the vectors in \c mv with random vectors.
  static void MvRandom(  _MV & mv )
    { 
      for (int i=0; i<mv.size(); i++) mv[i].randomize(); 
    }

  //! Assign (deep copy) A into mv.
  static void Assign( const _MV& A, _MV& mv ) {
    for (unsigned int i = 0; i < mv.size(); ++i) {
      mv[i].acceptCopyOf (A[i]);
    }
  }

  //! Replace each element of the vectors in \c mv with \c alpha.
  static void MvInit(  _MV & mv, double alpha = Teuchos::ScalarTraits<double>::zero() )
    { 
      //Out::os() << "MvInit()" << endl;
      for (int i=0; i<mv.size(); i++) mv[i].setToConstant(alpha); 
    }

  //@}

  /** \name Print method */
  //@{

  /** Print the \c mv multi-vector to the \c os output stream. */
  static void MvPrint( const  _MV & mv, std::ostream& os )
    { 
      os << mv << std::endl;
    }

  //@}

#ifdef HAVE_ANASAZI_TSQR
  /// \typedef tsqr_adaptor_type
  /// \brief TSQR adapter for the multivector type SimpleMV.
  ///
  /// For now, we provide a "stub" implementation.  It has the right
  /// methods and typedefs, but its constructors and methods all throw
  /// std::logic_error.  If you plan to use TSQR in Anasazi (e.g.,
  /// through TsqrOrthoManager) with SimpleMV, you must implement a
  /// functional TSQR adapter for SimpleMV.  Please refer to
  /// Epetra::TsqrAdapter (for Epetra_MultiVector) or
  /// Tpetra::TsqrAdaptor (for Tpetra::MultiVector) for examples of
  /// how to implement a TSQR adapter.
  typedef Anasazi::details::StubTsqrAdapter<SimpleMV> tsqr_adaptor_type;
#endif // HAVE_ANASAZI_TSQR

  /** */
  static bool detectRepeatedIndex(const std::vector<int>& index)
    {
      std::set<int> s;
      bool rtn = false;

      for (unsigned int i=0; i<index.size(); i++)
      {
        if (s.find(index[i]) != s.end())
        {
          //Out::os() << "detected repeated index " << index[i] << endl;
          rtn = true;
        }
        s.insert(index[i]);
      }
      
      return rtn;
    }

};        


/**

*/
template <>  
class OperatorTraits < double, SimpleMV, LinearOperator<double> >
{
public:
  typedef SimpleMV _MV;  
  /**
  */    
  static void Apply ( 
    const LinearOperator< double >& Op, 
    const  _MV & x,  
    _MV & y )
    {
      //Out::os() << "Apply()" << endl;
      y.resize(x.size());
      for (int i=0; i<x.size(); i++) 
      {
//        y[i] = Op * x[i];
        y[i].acceptCopyOf(Op * x[i]);
//        Out::os() << "i=" << i << " x=" << endl;
//        Out::os() << x[i] << endl;
//        Out::os() << "i=" << i << " y=" << endl;
//        Out::os() << y[i] << endl;
//        TEUCHOS_TEST_FOR_EXCEPT(x[i].norm2() < 1.0e-12);
      }
    }
    
};



} // end of Anasazi namespace 

#endif 
