// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   example_02.hpp
    \brief  Minimize the Gross-Pitaevskii functional and demonstrate 
            the effect of choice of function space of the Gradient on
            convergence. In this version we implement the correct Sobolev
            inner products and Reisz mapping using a finite difference 
            approximation of the derivative terms.          
                

    \details Minimize the one-dimensional Gross-Pitaevskii (GP) energy 
             functional
             \f[ J[\psi] = \int \frac{1}{2} |\nabla\psi|^2 + V(x)|\psi|^2 
                           +g|\psi|^4 \,\mathrm{d}x \f]
             Subject to the equality constraint that the particle density be
             normalized. 
             \f[ e(\psi) = \int |\psi|^2\,\mathrm{d}x - 1 = 0 \f]
             For simplicity, we will assume the wavefunction \f$\psi\f$ to 
             be real-valued, the potential function \f$ V(x)\geq 0\f$,
             the computational domain is the interval \f$[0,1]\f$, and that
             \f$\psi(0)=\psi(1)=0\f$. We also discretize the problem using
             second-order centered finite differences on a uniform grid. 

             \f[
             \psi''(x_i) \approx = \frac{\psi(x_{i-1})-2\psi(x_i)+\psi(x_{i+1})}{\Delta x^2}
             \f]

             The gradient with respect to the \f$L^2\f$ inner product is actually 
             an element of \f$H^{-1}\f$, so if we search in this direction, we
             are actually looking for a solution \f$\psi\f$ in a larger space than
             we should. If we compute the gradient with respect to the \f$H^1\f$ 
             inner product, by solving a Poisson equation, we search in the right space
             and the optimizer converges much faster than example_01.cpp which  
             does not do this.  
             

    \author Greg von Winckel
    \date   Wed Dec  3 16:40:45 MST 2014
*/

#include <iostream>

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"

#include "numerics/FiniteDifference.hpp"


using namespace ROL;

template <class Real, class Element=Real>
class OptStdVector;  // Optimization space.

template <class Real, class Element=Real>
class OptDualStdVector;  // Dual optimization space.

template <class Real, class Element=Real>
class ConStdVector;  // Constraint space.

template <class Real, class Element=Real>
class ConDualStdVector;  // Dual constraint space.

// Vector space definitions:


// Optimization space.
template <class Real, class Element>
class OptStdVector : public Vector<Real> {

  typedef std::vector<Element>       vector;
  typedef typename vector::size_type uint;

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<OptDualStdVector<Real> >  dual_vec_;

Teuchos::RCP<FiniteDifference<Real> > fd_;


public:

OptStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec, Teuchos::RCP<FiniteDifference<Real> >fd) : 
    std_vec_(std_vec), dual_vec_(Teuchos::null), fd_(fd) {}

void plus( const Vector<Real> &x ) {
  const OptStdVector &ex = Teuchos::dyn_cast<const OptStdVector>(x);
  Teuchos::RCP<const vector> xvalptr = ex.getVector();
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  uint dimension = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}


//! Modify the dot product between primal variables to be \f$(u,v)=\int\limits_0^1 \dot u \dot v\,\mathrm{d}x \f$
Real dot( const Vector<Real> &x ) const {
  Real val = 0;
  const OptStdVector<Real, Element> & ex = Teuchos::dyn_cast<const OptStdVector>(x);
  Teuchos::RCP<const vector> xvalptr = ex.getVector();
   
  Teuchos::RCP<vector> kxvalptr = Teuchos::rcp( new vector(std_vec_->size(), 0.0) );

  fd_->apply(xvalptr,kxvalptr);

  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*kxvalptr)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<Vector<Real> > clone() const {
  return Teuchos::rcp( new OptStdVector( Teuchos::rcp( new vector(std_vec_->size()) ),fd_ ) );
}

Teuchos::RCP<const vector> getVector() const {
  return std_vec_;
}

Teuchos::RCP<vector> getVector()  {
  return std_vec_;
}

Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<vector> e_rcp = Teuchos::rcp( new vector(std_vec_->size(),0.0) );
  Teuchos::RCP<OptStdVector> e = Teuchos::rcp( new OptStdVector( e_rcp, fd_ ) );
  (*e_rcp)[i]= 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}


//! Modify the dual of vector u to be \f$\tilde u = -\ddot u\f$
const Vector<Real> & dual() const {
  Teuchos::RCP<vector> dual_vecp = Teuchos::rcp(new vector(*std_vec_));
  dual_vec_ = Teuchos::rcp( new OptDualStdVector<Real>( dual_vecp, fd_ ) );
  fd_->apply(dual_vecp); 
  return *dual_vec_;
}

}; // class OptStdVector


// Dual optimization space.
template <class Real, class Element>
class OptDualStdVector : public Vector<Real> {

  typedef std::vector<Element>       vector;
  typedef typename vector::size_type uint;

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<OptStdVector<Real> >  dual_vec_;
Teuchos::RCP<FiniteDifference<Real> > fd_;

public:

OptDualStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec, Teuchos::RCP<FiniteDifference<Real> >fd) : 
    std_vec_(std_vec), dual_vec_(Teuchos::null), fd_(fd) {}

void plus( const Vector<Real> &x ) {
  const OptDualStdVector &ex = Teuchos::dyn_cast<const OptDualStdVector>(x);
  Teuchos::RCP<const vector> xvalptr = ex.getVector();
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  uint dimension = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const Vector<Real> &x ) const {
  Real val = 0;
  const OptDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<const OptDualStdVector<Real, Element> >(x);
  Teuchos::RCP<const vector> kxvalptr = ex.getVector();
  Teuchos::RCP<vector> xvalptr = Teuchos::rcp( new vector(std_vec_->size(), 0.0) );
  fd_->solve(kxvalptr,xvalptr);
  
  uint dimension  = std_vec_->size();
  for (unsigned i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<Vector<Real> > clone() const {
  return Teuchos::rcp( new OptDualStdVector( Teuchos::rcp( new std::vector<Element>(std_vec_->size()) ), fd_ ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<std::vector<Element> > getVector()  {
  return std_vec_;
}

Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<vector> e_rcp = Teuchos::rcp( new vector(std_vec_->size(), 0.0 ) );
  Teuchos::RCP<OptDualStdVector> e = Teuchos::rcp( new OptDualStdVector( e_rcp,fd_ ) );
  (*e_rcp)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const Vector<Real> & dual() const {
    Teuchos::RCP<vector> dual_vecp = Teuchos::rcp(new vector(*std_vec_)); 
    dual_vec_ = Teuchos::rcp( new OptStdVector<Real>( dual_vecp, fd_ ) );
    
    fd_->solve(dual_vecp);
    return *dual_vec_;
}

}; // class OptDualStdVector




// Constraint space.
template <class Real, class Element>
class ConStdVector : public Vector<Real> {

  typedef std::vector<Element> vector;
  typedef typename vector::size_type uint;

private:
Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<ConDualStdVector<Real> >  dual_vec_;

public:

ConStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(Teuchos::null) {}

void plus( const Vector<Real> &x ) {
  const ConStdVector &ex = Teuchos::dyn_cast<const ConStdVector>(x);
  Teuchos::RCP<const vector> xvalptr = ex.getVector();
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  uint dimension = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const Vector<Real> &x ) const {
  Real val = 0;
  const ConStdVector<Real, Element> & ex = Teuchos::dyn_cast<const ConStdVector<Real, Element> >(x);
  Teuchos::RCP<const vector> xvalptr = ex.getVector();

  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<Vector<Real> > clone() const {
  return Teuchos::rcp( new ConStdVector( Teuchos::rcp(new vector(std_vec_->size())) ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<std::vector<Element> > getVector() {
  return std_vec_;
}

Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<vector> e_rcp = Teuchos::rcp( new vector(std_vec_->size(),0.0) );
  Teuchos::RCP<ConStdVector> e = Teuchos::rcp( new ConStdVector( e_rcp) );
  (*e_rcp)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const Vector<Real> & dual() const {
  dual_vec_ = Teuchos::rcp( new ConDualStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ) ) );
  return *dual_vec_;
}

}; // class ConStdVector


// Dual constraint space.
template <class Real, class Element>
class ConDualStdVector : public Vector<Real> {

  typedef std::vector<Element>       vector;
  typedef typename vector::size_type uint;

private:

Teuchos::RCP<std::vector<Element> >  std_vec_;
mutable Teuchos::RCP<ConStdVector<Real> >  dual_vec_;

public:

ConDualStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec), dual_vec_(Teuchos::null) {}

void plus( const Vector<Real> &x ) {
  const ConDualStdVector &ex = Teuchos::dyn_cast<const ConDualStdVector>(x);
  Teuchos::RCP<const vector> xvalptr = ex.getVector();
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] += (*xvalptr)[i];
  }
}

void scale( const Real alpha ) {
  uint dimension = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    (*std_vec_)[i] *= alpha;
  }
}

Real dot( const Vector<Real> &x ) const {
  Real val = 0;
  const ConDualStdVector<Real, Element> & ex = Teuchos::dyn_cast<const ConDualStdVector<Real, Element> >(x);
  Teuchos::RCP<const vector> xvalptr = ex.getVector();
  uint dimension  = std_vec_->size();
  for (uint i=0; i<dimension; i++) {
    val += (*std_vec_)[i]*(*xvalptr)[i];
  }
  return val;
}

Real norm() const {
  Real val = 0;
  val = std::sqrt( dot(*this) );
  return val;
}

Teuchos::RCP<Vector<Real> > clone() const {
  return Teuchos::rcp( new ConDualStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size())) ) );
}

Teuchos::RCP<const std::vector<Element> > getVector() const {
  return std_vec_;
}

Teuchos::RCP<std::vector<Element> > getVector() {
  return std_vec_;
}

Teuchos::RCP<Vector<Real> > basis( const int i ) const {
  Teuchos::RCP<vector> e_rcp = Teuchos::rcp( new vector(std_vec_->size(),0.0) );
  Teuchos::RCP<ConDualStdVector> e = Teuchos::rcp( new ConDualStdVector( e_rcp ) );
  (*e_rcp)[i] = 1.0;
  return e;
}

int dimension() const {return static_cast<int>(std_vec_->size());}

const Vector<Real> &  dual() const {
  dual_vec_ = Teuchos::rcp( new ConStdVector<Real>( Teuchos::rcp( new std::vector<Element>(*std_vec_) ) ) );
  return *dual_vec_;
}

}; // class ConDualStdVector

/*** End of declaration of four vector spaces. ***/



/** Objective Function Class */
template<class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real> >
class Objective_GrossPitaevskii : public Objective<Real> {

    typedef std::vector<Real> vector;
    typedef typename vector::size_type uint;

    private:

        /** \var Real g_ appearing before quartic term in GP functional    */ 
        Real g_;    

        /** \var int nx_ Number of interior nodes  */ 
        uint  nx_;     

        /*! \var int nx_ Mesh spacing \f$ \Delta x = \frac{1}{n_x+1} \f$  */ 
        Real dx_;     
        
        /*! \var ptr Vp_ Pointer to potential vector  */ 
        Teuchos::RCP<const std::vector<Real> > Vp_;    

        Teuchos::RCP<FiniteDifference<Real> > fd_;

        //! Apply finite difference operator 
        /*! Compute \f$K\psi\f$, where \f$K\f$ is the finite difference approximation 
            of \f$-D_x^2\f$ */
        void applyK(const Vector<Real> &v, Vector<Real> &kv) {

            using Teuchos::RCP;  using Teuchos::dyn_cast;

            // Pointer to direction vector 
            RCP<const vector> vp = dyn_cast<const XPrim>(v).getVector();

            // Pointer to action of Hessian on direction vector 
            RCP<vector> kvp = dyn_cast<XDual>(kv).getVector();

            Real dx2 = dx_*dx_;

            (*kvp)[0] = (2.0*(*vp)[0]-(*vp)[1])/dx2;
  
            for(uint i=1;i<nx_-1;++i) {
                (*kvp)[i] = (2.0*(*vp)[i]-(*vp)[i-1]-(*vp)[i+1])/dx2;
            } 

            (*kvp)[nx_-1] = (2.0*(*vp)[nx_-1]-(*vp)[nx_-2])/dx2;

        } 

    public: 

        Objective_GrossPitaevskii(const Real &g, const Vector<Real> &V, Teuchos::RCP<FiniteDifference<Real> > fd) : g_(g),  
            Vp_((Teuchos::dyn_cast<const StdVector<Real> >(V)).getVector()), fd_(fd)  {

            nx_ = Vp_->size(); 
            dx_ = (1.0/(1.0+nx_));
        }
           
    //! Evaluate \f$J[\psi]\f$
    /*! \f[ J[\psi]=\frac{1}{2} \int\limits_0^1 |\psi'|^2 + 
            V(x)|\psi|^2+g|\psi|^4\,\mathrm{d}x \f] 
          where the integral is approximated with the trapezoidal rule and
          the derivative is approximated using finite differences */
    Real value( const Vector<Real> &psi, Real &tol ) {

        using Teuchos::RCP;  using Teuchos::rcp;  using Teuchos::dyn_cast;

        // Pointer to opt vector 
        RCP<const vector> psip = dyn_cast<const XPrim>(psi).getVector();

        // Pointer to K applied to opt vector 
        RCP<vector> kpsip = rcp( new vector(nx_, 0.0) );
        XDual kpsi(kpsip,fd_);

        Real J = 0;

        applyK(psi,kpsi);

        for(uint i=0;i<nx_;++i) {
            J += (*psip)[i]*(*kpsip)[i] + (*Vp_)[i]*pow((*psip)[i],2) + g_*pow((*psip)[i],4);
        } 
      
        // Rescale for numerical integration by trapezoidal rule
        J *= 0.5*dx_;

        return J;
    }

    //! Evaluate \f$\nabla J[\psi]\f$
    /*! \f[ \nabla J[\psi] = -\psi'' + V(x)\psi+2g|\psi|^3 \f] */
    void gradient( Vector<Real> &g, const Vector<Real> &psi, Real &tol ) {

        using Teuchos::RCP;  using Teuchos::rcp;  using Teuchos::dyn_cast;

        // Pointer to opt vector 
        RCP<const vector> psip = dyn_cast<const XPrim>(psi).getVector();

        // Pointer to gradient vector 
        RCP<vector> gp = dyn_cast<XDual>(g).getVector();

        // Pointer to K applied to opt vector 
        RCP<vector> kpsip = rcp( new vector(nx_, 0.0) );
        XDual kpsi(kpsip,fd_);

        applyK(psi,kpsi);

        for(uint i=0;i<nx_;++i) {
            (*gp)[i] = ((*kpsip)[i] + (*Vp_)[i]*(*psip)[i] + 2.0*g_*pow((*psip)[i],3))*dx_;
        } 
      
    }



    //! Evaluate \f$\nabla^2 J[\psi] v\f$
    /*! \f[ \nabla^2 J[\psi]v = -v'' + V(x)v+6g|\psi|^2 v \f] */
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol ) {

        using Teuchos::RCP;  using Teuchos::dyn_cast;

        // Pointer to opt vector 
        RCP<const vector> psip = dyn_cast<const XPrim>(psi).getVector();

        // Pointer to direction vector 
        RCP<const vector> vp = dyn_cast<const XPrim>(v).getVector();

        // Pointer to action of Hessian on direction vector 
        RCP<vector> hvp = dyn_cast<XDual>(hv).getVector();

        applyK(v,hv);
 
        for(uint i=0;i<nx_;++i) {
            (*hvp)[i] *= dx_;
            (*hvp)[i] += ( (*Vp_)[i] + 6.0*g_*pow((*psip)[i],2) )*(*vp)[i]*dx_;
        } 

   }

};


/** Constraint class */
template<class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real>, class CPrim=StdVector<Real>, class CDual=StdVector<Real> >
class Normalization_Constraint : public EqualityConstraint<Real> {

    typedef std::vector<Real> vector;
    typedef typename vector::size_type uint;

    private:     
    uint nx_;
    Real dx_;
    Teuchos::RCP<FiniteDifference<Real> > fd_;
    bool exactsolve_; 

    public:
    Normalization_Constraint(int n, Real dx, Teuchos::RCP<FiniteDifference<Real> > fd, bool exactsolve) : 
        nx_(n), dx_(dx), fd_(fd), exactsolve_(exactsolve) {}          

    //! Evaluate \f$c[\psi]\f$
    /*! \f[ c[\psi]= \int\limits_0^1 |\psi|^2\,\mathrm{d}x - 1 \f] 
        where the integral is approximated with the trapezoidal rule and
        the derivative is approximated using finite differences. 
        This constraint is a scalar */
    void value(Vector<Real> &c, const Vector<Real> &psi, Real &tol){

        using Teuchos::RCP;  using Teuchos::dyn_cast;

        // Pointer to constraint vector (only one element)
        RCP<vector> cp = dyn_cast<CPrim>(c).getVector();

        // Pointer to opt vector 
        RCP<const vector> psip = dyn_cast<const XPrim>(psi).getVector();

        (*cp)[0] = -1.0;
        for(uint i=0;i<nx_;++i) {
            (*cp)[0] += dx_*pow((*psip)[i],2);
        } 
    }

    //! Evaluate \f$c'[\psi]v\f$
    /*! \f[ c'[\psi]v= 2 \int\limits_0^1 \psi v\,\mathrm{d}x  \f]
         The action of the Jacobian on a vector produces a scalar */
    void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol){

        using Teuchos::RCP;  using Teuchos::dyn_cast;

        // Pointer to action of Jacobian of constraint on direction vector (yields scalar)
        RCP<vector> jvp = dyn_cast<CPrim>(jv).getVector();

        // Pointer to direction vector     
        RCP<const vector> vp = dyn_cast<const XPrim>(v).getVector();

        // Pointer to opt vector 
        RCP<const vector> psip = dyn_cast<const XPrim>(psi).getVector();
     
        (*jvp)[0] = 0;
        for(uint i=0;i<nx_;++i) {
            (*jvp)[0] += 2.0*dx_*(*psip)[i]*(*vp)[i];
        }
    }

    //! Evaluate \f$(c'[\psi])^\ast v\f$
    /*! \f[ (c'[\psi])^\ast v = 2 \int\limits_0^1 \psi v\,\mathrm{d}x  \f] 
         The action of the Jacobian adjoint on a scalar produces a vector */
    void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &psi, Real &tol){

        using Teuchos::RCP;  using Teuchos::dyn_cast;

        // Pointer to action of adjoint of Jacobian of constraint on direction vector (yields vector)
        RCP<vector> ajvp = dyn_cast<XDual>(ajv).getVector();

        // Pointer to direction vector     
        RCP<const vector> vp = dyn_cast<const CDual>(v).getVector();
 
        // Pointer to opt vector 
        RCP<const vector> psip = dyn_cast<const XPrim>(psi).getVector();

        for(uint i=0;i<nx_;++i) {
            (*ajvp)[i] = 2.0*dx_*(*psip)[i]*(*vp)[0];
        }
    }

    //! Evaluate \f$((c''[\psi])^\ast v)u\f$
    /*! \f[ ((c''[\psi])^\ast v)u = 2 v u   \f] 
         The action of the Hessian adjoint on a on a vector v in a direction u produces a vector of
         the same size as \f$\psi\f$ */
    void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, 
                             const Vector<Real> &psi, Real &tol){


        using Teuchos::RCP;  using Teuchos::dyn_cast;

        // The pointer to action of constraint Hessian in u,v inner product
        RCP<vector> ahuvp = dyn_cast<XDual>(ahuv).getVector();

        // Pointer to direction vector u     
        RCP<const vector> up = dyn_cast<const CDual>(u).getVector();

        // Pointer to direction vector v     
        RCP<const vector> vp = dyn_cast<const XPrim>(v).getVector();

        // Pointer to opt vector 
        RCP<const vector> psip = dyn_cast<const XPrim>(psi).getVector();
  
        for(uint i=0;i<nx_;++i) {
            (*ahuvp)[i] = 2.0*dx_*(*vp)[i]*(*up)[0];        
        }  
    }
     
    /** Solve the system \f[ \begin{\pmatrix} K & c'^\ast(\psi)\\ c'(\psi) & 0 \end{pmatrix}
      * \begin{pmatrix} v_1\\v_2 \end{pmatrix}=\begin{pmatrix} b_1\\b_2\end{pmatrix}\f]
      *  In this example, \f$K\f$ is the finite difference Laplacian the constraint is a 
      * scalar and the Jacobian is a vector and the exact inverse can be computed using the
      * Schur complement method */
    std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2, const Vector<Real> &b1, 
                                           const Vector<Real> &b2, const Vector<Real> &psi, Real &tol) {

        using Teuchos::RCP;  using Teuchos::rcp;  using Teuchos::dyn_cast;

        if(exactsolve_) {
	    RCP<vector> v1p = dyn_cast<XPrim>(v1).getVector();    
	    RCP<vector> v2p = dyn_cast<CDual>(v2).getVector();
	    RCP<const vector> b1p = dyn_cast<const XDual>(b1).getVector();
	    RCP<const vector> b2p = dyn_cast<const CPrim>(b2).getVector();
	    RCP<const vector> psip = dyn_cast<const XPrim>(psi).getVector();
	
	    RCP<vector> jacp = rcp( new vector(nx_, 0.0) );
	    RCP<vector> b1dp = rcp( new vector(nx_, 0.0) );

	    for(uint i=0;i<nx_;++i) {
		(*jacp)[i] = (*psip)[i];
		(*b1dp)[i] = (*b1p)[i];
	    }
	 
	    // The Jacobian of the constraint is \f$c'(\psi)=2dx\psi\f$
	    XDual jac(jacp,fd_);
	    jac.scale(2.0*dx_);

	    // A Dual-in-name-only version of b1, so we can compute the desired inner products involving inv(K) 
	    XDual b1d(b1dp,fd_);
	
	    // \f$ (c'K^{-1}*c'^\ast)^{-1} \f$ 
	    Real d = 1.0/jac.dot(jac);
	    Real p = jac.dot(b1d);

	    (*v2p)[0] = d*(p-(*b2p)[0]);
     
	    v1.set(jac.dual());
	    v1.scale(-(*v2p)[0]);
	    v1.plus(b1d.dual());  

	    return std::vector<Real>(0);
	}     
	else{
	    return EqualityConstraint<Real>::solveAugmentedSystem(v1,v2,b1,b2,psi,tol);
	}
    }
};

