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

#ifndef PLAYA_LINEARCOMBINATIONDECL_HPP
#define PLAYA_LINEARCOMBINATIONDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaVectorDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"



namespace Playa
{

template <class Scalar> class LinearOperator;

template <class Scalar, int N>
class LCNBase
{
public:
  /** */
  LCNBase(){}

  /** */
  virtual Vector<Scalar> eval() const = 0 ;
  
  /** */
  const Vector<Scalar>& vec(int i) const {return x_[i];}

  /** */
  const Scalar& coeff(int i) const {return a_[i];}

  /** Multiply through by a scalar constant */
  void multiply(const Scalar& beta) 
    {for (int i=0; i<N; i++) a_[i] *= beta;}

  /** */
  int size() const {return N;}

  /** Element-by-element product (Matlab dot-star operator) */
  Vector<Scalar> dotStar(const Vector<Scalar>& other) const 
    {return eval().dotStar(other);}

  /** Element-by-element division (Matlab dot-slash operator) */
  Vector<Scalar> dotSlash(const Vector<Scalar>& other) const 
    {return eval().dotSlash(other);}

  /** Return element-by-element reciprocal as a new vector */
  Vector<Scalar> reciprocal() const {return eval().reciprocal();}

  /** Return element-by-element absolute value as a new vector */
  Vector<Scalar> abs() const {return eval().abs();} 

  /** Return Euclidean norm */
  Scalar norm2() const {return eval().norm2();}

  /** Return 1-norm */
  Scalar norm1() const {return eval().norm1();}

  /** Return infinity norm */
  Scalar normInf() const {return eval().normInf();}

  /** Return max element */
  Scalar max() const {return eval().max();}

  /** Return min element */
  Scalar min() const {return eval().min();}

protected:

  Vector<Scalar> x_[N];
  Scalar a_[N];
};

template <class Scalar, int N>
class LCN : public LCNBase<Scalar, N>
{
public:
  /** */
  LCN(){}
  
  /** */
  void set(int i, const Scalar& a, const Vector<Scalar>& x)
    {
      this->a_[i] = a;
      this->x_[i] = x;
    }

  /** Unary minus */
  LCN<Scalar, N> operator-() const 
    {
      LCN<Scalar, N> rtn=*this; 
      rtn.multiply(-1.0);
      return rtn;
    }

  /** */
  Vector<Scalar> eval() const 
    {
      Vector<Scalar> rtn = this->x_[0].copy();
      if (this->a_[0] != 1.0) rtn.scale(this->a_[0]);
      for (int i=1; i<N; i++)
      {
        rtn += this->a_[i]*this->x_[i];
      }
      return rtn;
    }
};


//template<>
template <class Scalar>
class LCN<Scalar, 1> : public LCNBase<Scalar, 1>
{
public:
  /** */
  LCN(const Scalar& a, const Vector<Scalar>& x) : LCNBase<Scalar, 1>()
    {
      this->x_[0] = x;
      this->a_[0] = a;
    }
    
  /** Unary minus */
  LCN<Scalar, 1> operator-() const 
    {
      LCN<Scalar, 1> rtn=*this; 
      rtn.multiply(-1.0);
      return rtn;
    }

  /** */
  Vector<Scalar> eval() const 
    {
      Vector<Scalar> rtn = this->x_[0].copy();
      if (this->a_[0] != 1.0) rtn.scale(this->a_[0]);
      return rtn;
    }
};


//template<>
template <class Scalar>
class LCN<Scalar, 2> : public LCNBase<Scalar, 2>
{
public:
  /** */
  LCN(const Scalar& a, const Vector<Scalar>& x,
    const Scalar& b, const Vector<Scalar>& y) : LCNBase<Scalar, 2>()
    {
      this->x_[0] = x;
      this->a_[0] = a;
      this->x_[1] = y;
      this->a_[1] = b;
    }

  /** */
  LCN(const LCN<Scalar, 1>& ax,
    const Scalar& b, const Vector<Scalar>& y) : LCNBase<Scalar, 2>()
    {
      this->x_[0] = ax.vec(0);
      this->a_[0] = ax.coeff(0);
      this->x_[1] = y;
      this->a_[1] = b;
    }

  /** */
  LCN(const LCN<Scalar, 1>& ax,
    const LCN<Scalar, 1>& by) : LCNBase<Scalar, 2>()
    {
      this->x_[0] = ax.vec(0);
      this->a_[0] = ax.coeff(0);
      this->x_[1] = by.vec(0);
      this->a_[1] = by.coeff(0);
    }
    

  /** */
  LCN(const Scalar& a, const Vector<Scalar>& x,
    const LCN<Scalar, 1>& by) : LCNBase<Scalar, 2>()
    {
      this->x_[0] = x;
      this->a_[0] = a;
      this->x_[1] = by.vec(0);
      this->a_[1] = by.coeff(0);
    }
    
  /** Unary minus */
  LCN<Scalar, 2> operator-() const 
    {
      LCN<Scalar, 2> rtn=*this; 
      rtn.multiply(-1.0);
      return rtn;
    }

  /** */
  Vector<Scalar> eval() const 
    {
      Vector<Scalar> rtn = this->x_[0].copy();
      if (this->a_[0] != 1.0) rtn.scale(this->a_[0]);
      rtn += this->a_[1]*this->x_[1];
      return rtn;
    }
};


//template<>
template <class Scalar>
class LCN<Scalar, 3> : public LCNBase<Scalar, 3>
{
public:
  /** */
  LCN(const LCN<Scalar, 1>& ax,
    const LCN<Scalar, 2>& bycz)
    {
      this->x_[0] = ax.vec(0);
      this->a_[0] = ax.coeff(0);
      this->x_[1] = bycz.vec(0);
      this->a_[1] = bycz.coeff(0);
      this->x_[2] = bycz.vec(1);
      this->a_[2] = bycz.coeff(1);
    }
    

  /** */
  LCN(const LCN<Scalar, 2>& axby,
    const LCN<Scalar, 1>& cz)
    {
      this->x_[0] = axby.vec(0);
      this->a_[0] = axby.coeff(0);
      this->x_[1] = axby.vec(1);
      this->a_[1] = axby.coeff(1);
      this->x_[2] = cz.vec(0);
      this->a_[2] = cz.coeff(0);
    }

  /** Unary minus */
  LCN<Scalar, 3> operator-() const 
    {
      LCN<Scalar, 3> rtn=*this; 
      rtn.multiply(-1.0);
      return rtn;
    }

  /** */
  Vector<Scalar> eval() const 
    {
      Vector<Scalar> rtn = this->x_[0].copy();
      if (this->a_[0] != 1.0) rtn.scale(this->a_[0]);
      rtn += this->a_[1]*this->x_[1] + this->a_[2]*this->x_[2];
      return rtn;
    }
};





/*======================================================================
 *
 *    Write a LC description
 *
 *======================================================================*/

template <class Scalar, int N> inline
std::ostream& operator<<(std::ostream& os, const LCN<Scalar, N>& lc)
{
  os << "{(size=" << lc.size() << ") ";
  for (int i=0; i<lc.size(); i++)
  {
    if (i != 0) os << ", ";
    os << "(" << lc.coeff(i) << ", " << lc.vec(i).description() << ")";
  }
  os << "}";
  return os;
}

/*======================================================================
 *
 *    operator times vector or LC
 *
 *======================================================================*/

/** */
template <class Scalar> 
Vector<Scalar> operator*(const LinearOperator<Scalar>& A,
  const Vector<Scalar>& x);


/** */
template <class Scalar, int N> 
Vector<Scalar> operator*(const LinearOperator<Scalar>& A,
  const LCN<Scalar, N>& x);
  
/** */
template <class Scalar> 
Vector<Scalar> operator*(const LinearOperator<Scalar>& A,
  const LCN<Scalar, 1>& x);


/*======================================================================
 *
 *    scalar times vector
 *
 *======================================================================*/

/** \relates LCN scalar * vec */
template <class Scalar> 
LCN<Scalar, 1> operator*(const Scalar& alpha, 
  const Vector<Scalar>& x);

/** \relates LCN vec * scalar */
template <class Scalar> 
LCN<Scalar, 1> operator*(const Vector<Scalar>& x, 
  const Scalar& alpha);


/** \relates LCN vec / scalar */
template <class Scalar> 
LCN<Scalar, 1> operator/(const Vector<Scalar>& x, 
  const Scalar& alpha);


/*======================================================================
 *
 *    scalar times LC
 *
 *======================================================================*/


/** */
template <class Scalar, int N> 
LCN<Scalar, N> operator*(const LCN<Scalar, N>& lc, const Scalar& beta);

/** */
template <class Scalar, int N> 
LCN<Scalar, N> operator*(const Scalar& beta, const LCN<Scalar, N>& lc);

/** */
template <class Scalar, int N> 
LCN<Scalar, N> operator/(const LCN<Scalar, N>& lc, const Scalar& beta);


  
/*======================================================================
 *
 *  vector plus/minus vector
 *
 *======================================================================*/


/** */
template <class Scalar> 
LCN<Scalar, 2> operator+(const Vector<Scalar>& x, 
  const Vector<Scalar>& y);

/** */
template <class Scalar> 
LCN<Scalar, 2> operator-(const Vector<Scalar>& x, 
  const Vector<Scalar>& y);

/*======================================================================
 *
 *  LC plus LC
 *
 *======================================================================*/

/** */
template <class Scalar, int N, int M> 
LCN<Scalar, N+M> operator+(const LCN<Scalar, N>& f, const LCN<Scalar, M>& g);


/** */
template <class Scalar, int N> 
LCN<Scalar, N+1> operator+(const LCN<Scalar, N>& f, const Vector<Scalar>& g);


/** */
template <class Scalar, int N> 
LCN<Scalar, N+1> operator+(const Vector<Scalar>& f, const LCN<Scalar, N>& g);


/** */
template <class Scalar> 
LCN<Scalar, 2> operator+(const LCN<Scalar, 1>& lc, const Vector<Scalar>& x);

/** */
template <class Scalar> 
LCN<Scalar, 2> operator+(const Vector<Scalar>& x, const LCN<Scalar, 1>& lc);


/** */
template <class Scalar> 
LCN<Scalar, 2> operator+(const LCN<Scalar, 1>& ax, const LCN<Scalar, 1>& by);


/** */
template <class Scalar> 
LCN<Scalar, 3> operator+(const LCN<Scalar, 1>& ax, 
  const LCN<Scalar, 2>& bycz);

/** */
template <class Scalar> 
LCN<Scalar, 3> operator+(const Vector<Scalar>& x, const LCN<Scalar, 2>& bycz);


/** */
template <class Scalar> 
LCN<Scalar, 3> operator+(const LCN<Scalar, 2>& axby, 
  const LCN<Scalar, 1>& cz);


/** */
template <class Scalar> 
LCN<Scalar, 3> operator+(const LCN<Scalar, 2>& axby, 
  const Vector<Scalar>& z);


/*======================================================================
 *
 *  LC minus LC
 *
 *======================================================================*/


/** */
template <class Scalar, int N, int M> 
LCN<Scalar, N+M> operator-(const LCN<Scalar, N>& f, const LCN<Scalar, M>& g);


/** */
template <class Scalar, int N> 
LCN<Scalar, N+1> operator-(const LCN<Scalar, N>& f, const Vector<Scalar>& g);

/** */
template <class Scalar, int N> 
LCN<Scalar, N+1> operator-(const Vector<Scalar>& f, const LCN<Scalar, N>& g);

/** */
template <class Scalar> 
LCN<Scalar, 2> operator-(const LCN<Scalar, 1>& lc, const Vector<Scalar>& x);

/** */
template <class Scalar> 
LCN<Scalar, 2> operator-(const Vector<Scalar>& x, const LCN<Scalar, 1>& lc);

/** */
template <class Scalar> 
LCN<Scalar, 2> operator-(const LCN<Scalar, 1>& ax, const LCN<Scalar, 1>& by);


/** */
template <class Scalar> 
LCN<Scalar, 3> operator-(const LCN<Scalar, 1>& ax, 
  const LCN<Scalar, 2>& bycz);

/** */
template <class Scalar> 
LCN<Scalar, 3> operator-(const Vector<Scalar>& x, const LCN<Scalar, 2>& bycz);

/** */
template <class Scalar> 
LCN<Scalar, 3> operator-(const LCN<Scalar, 2>& axby, const LCN<Scalar, 1>& cz);

/** */
template <class Scalar> 
LCN<Scalar, 3> operator-(const LCN<Scalar, 2>& axby, const Vector<Scalar>& z);


/*======================================================================
 *
 *  Operations on LC
 *
 *======================================================================*/

/** */
template <class Scalar, int N> 
Scalar norm1(const LCN<Scalar, N>& lc) ;

/** */
template <class Scalar, int N> 
Scalar norm2(const LCN<Scalar, N>& lc) ;

/** */
template <class Scalar, int N> 
Scalar normInf(const LCN<Scalar, N>& lc) ;


/** */
template <class Scalar, int N> 
Vector<Scalar> abs(const LCN<Scalar, N>& lc) ;

/** */
template <class Scalar, int N> 
Scalar min(const LCN<Scalar, N>& lc) ;


/** */
template <class Scalar, int N> 
Scalar max(const LCN<Scalar, N>& lc) ;

/** */
template <class Scalar, int N> 
Vector<Scalar> reciprocal(const LCN<Scalar, N>& lc) ;




}



#endif
