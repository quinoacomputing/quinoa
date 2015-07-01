/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
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

#ifndef SUNDANCE_POINT_H
#define SUNDANCE_POINT_H

#include "SundanceDefs.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Utils.hpp"



namespace Sundance
{
  /**
   * Point represents a spatial point.
   */
  class Point
    {
    public:
      /* empty ctor */
      inline Point();
      inline Point(const double& x);
      inline Point(const double& x, const double& y);
      inline Point(const double& x, const double& y, const double& z);
      inline Point(const Point& other);
      inline Point& operator=(const Point& other);

      inline int dim() const {return dim_;}
      inline double& operator[](int i);
      inline const double& operator[](int i) const ;
      inline void resize(int i);

      /* reflexive arithmetic operators */

      inline Point& operator+=(const Point& p) ;
      inline Point& operator-=(const Point& p) ;
      inline Point& operator*=(const double& a) ;
      inline Point& operator/=(const double& a) ;

      /* unary plus and minus */

      inline Point operator+() const ;
      inline Point operator-() const ;

      /* binary operations with Points */

      inline Point operator+(const Point& p) const ;
      inline Point operator-(const Point& p) const ;
      inline double operator*(const Point& p) const ;

      /* binary operations with doubles */

      inline Point operator*(const double& a) const ;
      inline Point operator/(const double& a) const ;

      inline double distance(const Point& x) const ;

      inline std::string toString() const ;

      static bool unitTest() ;

    protected:
      void boundsCheck(int i) const ;
      int dim_;
      double x_[3];
    };
}


namespace std
{
  ostream& operator<<(std::ostream& os, const Sundance::Point& p);

}

namespace Sundance
{
  inline Point operator*(const double& a, const Point& p);

  inline Point cross(const Point& x, const Point& y);



  inline Point::Point()
    : dim_(0)
    {;}

  inline Point::Point(const double& x)
    : dim_(1)
    {
      x_[0] = x;
    }

  inline Point::Point(const double& x, const double& y)
    : dim_(2)
    {
      x_[0] = x;
      x_[1] = y;
    }

  inline Point::Point(const double& x, const double& y, const double& z)
    : dim_(3)
    {
      x_[0] = x;
      x_[1] = y;
      x_[2] = z;
    }

  inline Point::Point(const Point& other)
    : dim_(other.dim_)
    {
      for (int i=0; i<dim_; i++) x_[i] = other.x_[i];
    }

  Point& Point::operator=(const Point& other)
    {
      if (&other==this) return *this;

      dim_ = other.dim_;
      for (int i=0; i<dim_; i++) x_[i] = other.x_[i];
      return *this;
    }

  double& Point::operator[](int i)
    {
#ifndef NOBOUNDSCHECK
      boundsCheck(i);
#endif
      return x_[i];
    }

  const double& Point::operator[](int i) const
    {
#ifndef NOBOUNDSCHECK
      boundsCheck(i);
#endif
      return x_[i];
    }

  void Point::resize(int i)
    {
#ifndef NOBOUNDSCHECK
      TEUCHOS_TEST_FOR_EXCEPTION(i < 0 || i>3, std::runtime_error,
                         "void Point::resize() invalid dimension");
#endif
      dim_ = i;
    }

  Point& Point::operator+=(const Point& p)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(p.dim() != dim_, std::runtime_error,
                         "Point::operator+=() dimension mismatch "
                         "operands are: " << *this << " and " 
                         << p );

      for (int i=0; i<dim_; i++) x_[i] += p.x_[i];
      return *this;
    }

  Point& Point::operator-=(const Point& p)
    {

      TEUCHOS_TEST_FOR_EXCEPTION(p.dim() != dim_, std::runtime_error,
                         "Point::operator-=() dimension mismatch "
                         "operands are: " << *this << " and " 
                         << p );
      
      for (int i=0; i<dim_; i++) x_[i] -= p.x_[i];
      return *this;
    }

  Point& Point::operator*=(const double& a)
    {
      for (int i=0; i<dim_; i++) x_[i] *= a;
      return *this;
    }

  Point& Point::operator/=(const double& a)
    {
      for (int i=0; i<dim_; i++) x_[i] /= a;
      return *this;
    }

  Point Point::operator-() const
    {
      Point rtn(*this);
      for (int i=0; i<dim_; i++) rtn.x_[i] = -rtn.x_[i];
      return rtn;
    }

  Point Point::operator+() const
    {
      return *this;
    }

  Point Point::operator+(const Point& p) const
    {
      Point rtn(*this);
      rtn += p;
      return rtn;
    }

  Point Point::operator-(const Point& p) const
    {
      Point rtn(*this);
      rtn -= p;
      return rtn;
    }

  double Point::operator*(const Point& p) const
    {
      double rtn = 0.0;

      TEUCHOS_TEST_FOR_EXCEPTION(p.dim() != dim_, std::runtime_error,
                         "Point::operator*() dimension mismatch "
                         "operands are: " << *this << " and " 
                         << p );
      
      for (int i=0; i<dim_; i++) rtn += x_[i]*p.x_[i];
      return rtn;
    }

  Point Point::operator*(const double& a) const
    {
      Point rtn(*this);
      rtn *= a;
      return rtn;
    }

  Point Point::operator/(const double& a) const
    {
      Point rtn(*this);
      rtn /= a;
      return rtn;
    }

  Point operator*(const double& a, const Point& p)
    {
      return p.operator*(a);
    }

  inline Point cross(const Point& a, const Point& b)
  {
    return Point( a[1]*b[2] - b[1]*a[2], 
                  -a[0]*b[2] + a[2]*b[0],
                  a[0]*b[1] - a[1]*b[0]);
  }

  inline std::string Point::toString() const
    {
      std::string rtn = "{";
      for (int i=0; i<dim(); i++)
        {
          rtn += Teuchos::toString(x_[i]);
          if (i<dim()-1) rtn += ", ";
        }
      rtn += "}";
      return rtn;
    }


  inline double Point::distance(const Point& x) const 
  {
    Point dx = *this-x;
    return ::sqrt(dx*dx); 
  }
}

namespace Teuchos
{
  inline std::string toString(const Sundance::Point& x)
    {
      return x.toString();
    }
}

namespace std
{
  ostream& operator<<(std::ostream& os, const Sundance::Point& p);
}


#endif





