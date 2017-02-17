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


#ifndef ROL_SACADO_EQUALITYCONSTRAINT_SIMOPT
#define ROL_SACADO_EQUALITYCONSTRAINT_SIMOPT

#include "Sacado.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_EqualityConstraint_SimOpt.hpp"

using namespace ROL;

//! \brief ROL interface wrapper for Sacado SimOpt Constraint
template<class Real, template<class> class Constr>
class Sacado_EqualityConstraint_SimOpt : public EqualityConstraint_SimOpt<Real> {
 

   
    protected:     
        Constr<Real> constr_;

        template<class ScalarT>
        void applyJacobian_1AD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                               const Vector<ScalarT> &z, Real &tol);

        template<class ScalarT>
        void applyJacobian_2AD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                               const Vector<ScalarT> &z, Real &tol);

        template<class ScalarT>
        void applyAdjointJacobian_1AD(Vector<ScalarT> &ajv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                                      const Vector<ScalarT> &z, Real &tol);

        template<class ScalarT>
        void applyAdjointJacobian_2AD(Vector<ScalarT> &ajv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                                      const Vector<ScalarT> &z, Real &tol);
         
        template<class ScalarT>
        void applyAdjointHessian_11AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                      const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                      const Vector<ScalarT> &z, Real &tol);         
        template<class ScalarT>
        void applyAdjointHessian_12AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                      const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                      const Vector<ScalarT> &z, Real &tol);         
        template<class ScalarT>
        void applyAdjointHessian_21AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                      const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                      const Vector<ScalarT> &z, Real &tol);         
        template<class ScalarT>
        void applyAdjointHessian_22AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                      const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                      const Vector<ScalarT> &z, Real &tol);         
 
    public:
        Sacado_EqualityConstraint_SimOpt() : constr_(Constr<Real>()) {}
        Sacado_EqualityConstraint_SimOpt(Constr<Real> constr) : constr_(constr) { }

        void value(Vector<Real> &c, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
            constr_.value(c,u,z,tol);
        }

        void applyJacobian_1(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u, 
                             const Vector<Real> &z, Real &tol) {
            this->applyJacobian_1AD(jv,v,u,z,tol);
        } 

        void applyJacobian_2(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u, 
                             const Vector<Real> &z, Real &tol) {
            this->applyJacobian_2AD(jv,v,u,z,tol);
        } 

        void applyAdjointJacobian_1(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u, 
                                    const Vector<Real> &z, Real &tol) {
            this->applyAdjointJacobian_1AD(ajv,v,u,z,tol);
        } 

        void applyAdjointJacobian_2(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u, 
                                    const Vector<Real> &z, Real &tol) {
            this->applyAdjointJacobian_2AD(ajv,v,u,z,tol);
        } 

        void applyAdjointHessian_11(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                    const Vector<Real> &u, const Vector<Real> &z, Real &tol){
            this->applyAdjointHessian_11AD(ahwv,w,v,u,z,tol);
        }
 
        void applyAdjointHessian_12(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                    const Vector<Real> &u, const Vector<Real> &z, Real &tol){
            this->applyAdjointHessian_12AD(ahwv,w,v,u,z,tol);
        }
 
        void applyAdjointHessian_21(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                    const Vector<Real> &u, const Vector<Real> &z, Real &tol){
            this->applyAdjointHessian_21AD(ahwv,w,v,u,z,tol);
        }

        void applyAdjointHessian_22(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                                    const Vector<Real> &u, const Vector<Real> &z, Real &tol){
            this->applyAdjointHessian_22AD(ahwv,w,v,u,z,tol);
        }
};



template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_EqualityConstraint_SimOpt<Real,Constr>::applyJacobian_1AD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                                                                      const Vector<ScalarT> &z, Real &tol) {
    //    v in U (n), jv in C (n)    

    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast; 

    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();
    RCP<const vector> zp = dyn_cast<const SV>(z).getVector();
    RCP<vector> jvp = dyn_cast<SV>(jv).getVector(); 

    int n = up->size();
    int m = zp->size();
    
    RCP<Fadvector> c_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> z_fad_rcp = rcp( new Fadvector );

    c_fad_rcp->reserve(n);
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);

    for(int i=0; i<n; ++i) {
        c_fad_rcp->push_back(0);
        u_fad_rcp->push_back(FadType(n,i,(*up)[i]));
    }

    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back((*zp)[j]); 
    }

    StdVector<FadType> c_fad(c_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
  
    // Evaluate constraint 
    constr_.value(c_fad,u_fad,z_fad,tol);

    for(int i=0; i<n; ++i) {
        (*jvp)[i] = 0; 
        for(int j=0; j<n; ++j) {
            (*jvp)[i] += (*vp)[j]*(*c_fad_rcp)[i].dx(j);
        }
    } 
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_EqualityConstraint_SimOpt<Real,Constr>::applyJacobian_2AD(Vector<ScalarT> &jv, const Vector<ScalarT> &v, const Vector<ScalarT> &u, 
                                                                      const Vector<ScalarT> &z, Real &tol) {
    // v in Z (m), jv in C (n)

    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast;

    RCP<vector> jvp = dyn_cast<SV>(jv).getVector();
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();
    RCP<const vector> zp = dyn_cast<const SV>(z).getVector();

    int n = up->size();
    int m = zp->size(); 
    
    RCP<Fadvector> c_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> z_fad_rcp = rcp( new Fadvector );

    c_fad_rcp->reserve(n);
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);

    for(int i=0; i<n; ++i) {
        c_fad_rcp->push_back(0);
        u_fad_rcp->push_back((*up)[i]); 
    }

    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back(FadType(n,j,(*zp)[j]));
    }

    StdVector<FadType> c_fad(c_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
  
    // Evaluate constraint 
    constr_.value(c_fad,u_fad,z_fad,tol);

    for(int i=0; i<n; ++i) {
        (*jvp)[i] = 0; 
        for(int j=0; j<n; ++j) {
            (*jvp)[i] += (*vp)[j]*(*c_fad_rcp)[i].dx(j);
        }
    } 
 
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_EqualityConstraint_SimOpt<Real,Constr>::applyAdjointJacobian_1AD(Vector<ScalarT> &ajv, const Vector<ScalarT> &v,
                                                                             const Vector<ScalarT> &u, const Vector<ScalarT> &z, Real &tol) {

    // v in C* (n), ajv in U* (n)

    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast;  

    RCP<vector> ajvp = dyn_cast<SV>(ajv).getVector();
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();
    RCP<const vector> zp = dyn_cast<const SV>(z).getVector();
    
    int n = up->size();
    int m = zp->size();
    
    RCP<Fadvector> c_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> z_fad_rcp = rcp( new Fadvector );

    c_fad_rcp->reserve(n);
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);

    for(int i=0; i<n; ++i) {
        c_fad_rcp->push_back(0);
        u_fad_rcp->push_back(FadType(n,i,(*up)[i]));
    }

    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back((*zp)[j]); 
    }

    StdVector<FadType> c_fad(c_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
  
    // Evaluate constraint 
    constr_.value(c_fad,u_fad,z_fad,tol);

    FadType vdotc = 0;

    for(int i=0;i<n;++i) {
        vdotc += (*c_fad_rcp)[i]*(*vp)[i]; 
    } 

    for(int i=0;i<n;++i) {
        (*ajvp)[i] = vdotc.dx(i);
    }

}



template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_EqualityConstraint_SimOpt<Real,Constr>::applyAdjointJacobian_2AD(Vector<ScalarT> &ajv, const Vector<ScalarT> &v, 
                                                                             const Vector<ScalarT> &u, const Vector<ScalarT> &z, 
                                                                             Real &tol) {
    // v in C* (n), ajv in Z* (m)

    typedef Sacado::Fad::DFad<ScalarT> FadType;
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast;  

    RCP<vector> ajvp = dyn_cast<SV>(ajv).getVector();
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();   
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();   
    RCP<const vector> zp = dyn_cast<const SV>(z).getVector();   

    int n = up->size();
    int m = zp->size();
    
    RCP<Fadvector> c_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> z_fad_rcp = rcp( new Fadvector );

    c_fad_rcp->reserve(n);
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);

    for(int i=0; i<n; ++i) {
        c_fad_rcp->push_back(0);
        u_fad_rcp->push_back((*up)[i]); 
    }

    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back(FadType(n,j,(*zp)[j]));
    }

    StdVector<FadType> c_fad(c_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> z_fad(z_fad_rcp);
  
    // Evaluate constraint 
    constr_.value(c_fad,u_fad,z_fad,tol);

    FadType vdotc = 0;

    for(int i=0;i<n;++i) {
        vdotc += (*c_fad_rcp)[i]*(*vp)[i]; 
    } 

    for(int j=0;j<m;++j) {
        (*ajvp)[j] = vdotc.dx(j);
    }
}



template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_EqualityConstraint_SimOpt<Real,Constr>::applyAdjointHessian_11AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                                                             const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                                                             const Vector<ScalarT> &z, Real &tol) {
    // w in C* (n), v in U* (n), ahwv in U* (n)

    typedef Sacado::Fad::DFad<ScalarT> FadType;    
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast; 

    RCP<vector> ahwvp = dyn_cast<SV>(ahwv).getVector();
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();
    RCP<const vector> wp = dyn_cast<const SV>(w).getVector();
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();
    RCP<const vector> zp = dyn_cast<const SV>(z).getVector();

    int n = up->size();
    int m = zp->size();

    RCP<Fadvector> v_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> z_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> jv_fad_rcp = rcp( new Fadvector );

    v_fad_rcp->reserve(n);
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);
    jv_fad_rcp->reserve(n);

    for(int i=0; i<n; ++i) {
        v_fad_rcp->push_back((*vp)[i]);
        u_fad_rcp->push_back(FadType(n,i,(*up)[i]));
        jv_fad_rcp->push_back(0);
    }

    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back((*zp)[j]);    
    }

    StdVector<FadType> v_fad(v_fad_rcp);     
    StdVector<FadType> u_fad(u_fad_rcp);     
    StdVector<FadType> z_fad(z_fad_rcp);     
    StdVector<FadType> jv_fad(jv_fad_rcp);     

    this->applyJacobian_1AD(jv_fad,v_fad,u_fad,z_fad,tol);
 
    FadType wjv_fad = 0;
   
    for(int i=0; i<n; ++i) {
        wjv_fad += (*wp)[i]*(*jv_fad_rcp)[i];
    }

    for(int i=0; i<n; ++i) {
        (*ahwvp)[i] = wjv_fad.dx(i);
    }
}




template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_EqualityConstraint_SimOpt<Real,Constr>::applyAdjointHessian_12AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                                                             const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                                                             const Vector<ScalarT> &z, Real &tol) {
    // w in C* (n), v in U* (n), ahwv in Z* (m)

    typedef Sacado::Fad::DFad<ScalarT> FadType;    
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast;  

    RCP<vector> ahwvp = dyn_cast<SV>(ahwv).getVector();
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();
    RCP<const vector> wp = dyn_cast<const SV>(w).getVector();
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();
    RCP<const vector> zp = dyn_cast<const SV>(z).getVector();

    int n = up->size();
    int m = zp->size();

    RCP<Fadvector> v_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> z_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> jv_fad_rcp = rcp( new Fadvector );

    v_fad_rcp->reserve(n);
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);
    jv_fad_rcp->reserve(n);

    for(int i=0; i<n; ++i) {
        v_fad_rcp->push_back((*vp)[i]);
        u_fad_rcp->push_back((*up)[i]);
        jv_fad_rcp->push_back(0);
    }

    for(int j=0; j<m; ++j) {
        z_fad_rcp->push_back(FadType(m,j,(*zp)[j]));    
    }

    StdVector<FadType> v_fad(v_fad_rcp);     
    StdVector<FadType> u_fad(u_fad_rcp);     
    StdVector<FadType> z_fad(z_fad_rcp);     
    StdVector<FadType> jv_fad(jv_fad_rcp);     

    this->applyJacobian_1AD(jv_fad,v_fad,u_fad,z_fad,tol);
    
    FadType wjv_fad = 0;

    for(int i=0; i<n; ++i) {
        wjv_fad += (*wp)[i]*(*jv_fad_rcp)[i]; 
    }
 
    for(int j=0; j<m; ++j) {
        (*ahwvp)[j] = wjv_fad.dx(j);
    }
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_EqualityConstraint_SimOpt<Real,Constr>::applyAdjointHessian_21AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                                                             const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                                                             const Vector<ScalarT> &z, Real &tol) {
    // w in C* (n), v in Z* (m), ahwv in U* (n)

    typedef Sacado::Fad::DFad<ScalarT> FadType;    
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast; 

    RCP<vector> ahwvp = dyn_cast<SV>(ahwv).getVector();
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();
    RCP<const vector> wp = dyn_cast<const SV>(w).getVector();
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();
    RCP<const vector> zp = dyn_cast<const SV>(z).getVector();

    int n = up->size();
    int m = zp->size();

    RCP<Fadvector> v_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> z_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> jv_fad_rcp = rcp( new Fadvector );

    v_fad_rcp->reserve(m);
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);
    jv_fad_rcp->reserve(n);

    for(int i=0; i<n; ++i) {
        u_fad_rcp->push_back(FadType(1,(*up)[i]));
        jv_fad_rcp->push_back(0);
    }

    for(int j=0; j<m; ++j) {
        v_fad_rcp->push_back((*vp)[j]);
        z_fad_rcp->push_back((*zp)[j]);    
    }

    StdVector<FadType> v_fad(v_fad_rcp);     
    StdVector<FadType> u_fad(u_fad_rcp);     
    StdVector<FadType> z_fad(z_fad_rcp);     
    StdVector<FadType> jv_fad(jv_fad_rcp);     

    this->applyJacobian_2AD(jv_fad,v_fad,u_fad,z_fad,tol);

    FadType wjv_fad = 0;

    for(int i=0; i<n; ++i) {
        wjv_fad += (*wp)[i]*(*jv_fad_rcp)[i];
    }

    for(int i=0; i<n; ++i) {
        (*ahwvp)[i] = wjv_fad.dx(i);
    }
}


template<class Real, template<class> class Constr>
template<class ScalarT>
void Sacado_EqualityConstraint_SimOpt<Real,Constr>::applyAdjointHessian_22AD(Vector<ScalarT> &ahwv, const Vector<ScalarT> &w, 
                                                                             const Vector<ScalarT> &v, const Vector<ScalarT> &u,
                                                                             const Vector<ScalarT> &z, Real &tol) {
    // w in C* (n), v in Z* (m), ahwv in Z* (m)

    typedef Sacado::Fad::DFad<ScalarT> FadType;    
    typedef std::vector<FadType>       Fadvector;
    typedef std::vector<ScalarT>       vector;
    typedef StdVector<ScalarT>         SV;

    using Teuchos::RCP;       using Teuchos::rcp;
    using Teuchos::dyn_cast; 

    RCP<vector> ahwvp = dyn_cast<SV>(ahwv).getVector();
    RCP<const vector> vp = dyn_cast<const SV>(v).getVector();
    RCP<const vector> wp = dyn_cast<const SV>(w).getVector();
    RCP<const vector> up = dyn_cast<const SV>(u).getVector();
    RCP<const vector> zp = dyn_cast<const SV>(z).getVector();
 
    int n = up->size();
    int m = zp->size();

    RCP<Fadvector> v_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> u_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> z_fad_rcp = rcp( new Fadvector );
    RCP<Fadvector> jv_fad_rcp = rcp( new Fadvector );

    v_fad_rcp->reserve(m);
    u_fad_rcp->reserve(n);
    z_fad_rcp->reserve(m);
    jv_fad_rcp->reserve(n);

    for(int i=0; i<n; ++i) {
        u_fad_rcp->push_back((*up)[i]);
        jv_fad_rcp->push_back(0);
    }

    for(int j=0; j<m; ++j) {
        v_fad_rcp->push_back((*vp)[j]);
        z_fad_rcp->push_back(FadType(m,j,(*zp)[j]));    
    }

    StdVector<FadType> v_fad(v_fad_rcp);     
    StdVector<FadType> u_fad(u_fad_rcp);     
    StdVector<FadType> z_fad(z_fad_rcp);     
    StdVector<FadType> jv_fad(jv_fad_rcp);     

    this->applyJacobian_2AD(jv_fad,v_fad,u_fad,z_fad,tol);
 
    FadType wjv_fad = 0;

    for(int i=0; i<n; ++i) {
        wjv_fad += (*wp)[i]*(*jv_fad_rcp)[i];
    }

    for(int j=0; j<m; ++j) {
        (*ahwvp)[j] = wjv_fad.dx(j);

    }
}
#endif
