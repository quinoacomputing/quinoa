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

#ifndef SUNDANCE_UNARYFUNCTOR_H
#define SUNDANCE_UNARYFUNCTOR_H

#include "SundanceDefs.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceFunctorDomain.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RefCountPtr.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace Sundance
{
  using namespace Teuchos;

  
  

  /**
   * 
   */
  class UnaryFunctor
  {
  public:
    /** ctor */
    UnaryFunctor(const std::string& name, 
                 const RCP<FunctorDomain>& domain 
                 = rcp(new UnboundedDomain())) 
      : name_(name), h_(fdStep()), domain_(domain) {;}

    /** */
    virtual ~UnaryFunctor(){;}

    /** */
    const std::string& name() const {return name_;}

    /** */
    virtual void eval0(const double* const x, 
                       int nx, 
                       double* f) const = 0 ;
    
    /** */
    virtual void eval1(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx) const ;
    
    /** */
    virtual void eval2(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx,
                       double* d2f_dxx) const ;
    
    /** */
    virtual void eval3(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx,
                       double* d2f_dxx,
                       double* d3f_dxxx) const ;

    

    /** */
    void evalFDDerivs1(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx) const ;
    /** */
    void evalFDDerivs2(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx,
                       double* d2f_dxx) const ;
    /** */
    void evalFDDerivs3(const double* const x, 
                       int nx, 
                       double* f, 
                       double* df_dx,
                       double* d2f_dxx,
                       double* d3f_dxxx) const ;

    /** */
    bool testDerivs(const double& x, const double& tol) const ;

    /** */
    bool testInvalidValue(const double& xBad) const ;

    /** */
    bool test(int nx, const double& tol) const ;
    
    /** Specify whether we should test for NAN or INFINITE results. */
    static bool& checkResults() {static bool rtn = false; return rtn;}

    static double& fdStep() {static double rtn = 1.0e-3; return rtn;}

    const RCP<FunctorDomain>& domain() const 
    {return domain_;}
  private:
    std::string name_;

    double h_;

    RCP<FunctorDomain> domain_;
  };
}

/** */
#define SUNDANCE_UNARY_FUNCTOR(opName, functorName, description, domain, \
                               funcDefinition, firstDerivDefinition,    \
                               secondDerivDefinition)                   \
  class functorName : public Sundance::UnaryFunctor                \
  {                                                                     \
  public:                                                               \
    /** ctor for description functor */                                 \
    functorName() : Sundance::UnaryFunctor(#opName, rcp(new domain)) {;} \
      /** virtual dtor */                                               \
      virtual ~functorName(){;}                                         \
      /** Evaluate function at an array of values */                    \
      void eval0(const double* const x, int nx, double* f) const ;      \
      /** Evaluate function and first derivative at an array of values */ \
      void eval1(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df) const ;                                    \
      /** Evaluate function and first two derivatives at an array of values */ \
      void eval2(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df,                                            \
                 double* d2f_dxx) const ;                               \
  };                                                                    \
  inline void functorName::eval0(const double* const x, int nx, double* f) const \
  {                                                                     \
    if (checkResults())                                                 \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
          {                                                             \
            f[i] = funcDefinition;                                      \
            TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i]), std::runtime_error, "Non-normal floating point result detected in evaluation of unary functor " << name() << " at argument " << x[i]); \
          }                                                             \
     }                                                                  \
   else                                                                 \
     {                                                                  \
       for (int i=0; i<nx; i++) f[i] = funcDefinition;                  \
     }                                                                  \
}                                                                       \
  inline void functorName::eval1(const double* const x,                 \
                               int nx,                                  \
                               double* f,                               \
                               double* df) const                        \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
             df[i] = firstDerivDefinition;                              \
             TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i]), \
                                std::runtime_error,                           \
                                "Non-normal floating point result detected in " \
                                "evaluation of unary functor "          \
                                << name() << " at argument " << x[i]);  \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
             df[i] = firstDerivDefinition;                              \
           }                                                            \
      }                                                                 \
}                                                                       \
  inline void functorName::eval2(const double* const x,                 \
                               int nx,                                  \
                               double* f,                               \
                               double* df,                              \
                               double* d2f) const                       \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
              df[i] = firstDerivDefinition;                             \
              d2f[i] = secondDerivDefinition;                           \
              TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i])||Teuchos::ScalarTraits<double>::isnaninf(d2f[i]), \
                                 std::runtime_error,                          \
                                 "Non-normal floating point result detected in " \
                                 "evaluation of unary functor "         \
                                 << name() << " at argument " << x[i] ); \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
               df[i] = firstDerivDefinition;                            \
               d2f[i] = secondDerivDefinition;                          \
           }                                                            \
      }                                                                 \
}



/** */
#define SUNDANCE_UNARY_FUNCTOR3(opName, functorName, description, domain, \
                               funcDefinition, firstDerivDefinition,    \
                                secondDerivDefinition, thirdDerivDefinition) \
  class functorName : public Sundance::UnaryFunctor                \
  {                                                                     \
  public:                                                               \
    /** ctor for description functor */                                 \
    functorName() : Sundance::UnaryFunctor(#opName, rcp(new domain)) {;} \
      /** virtual dtor */                                               \
      virtual ~functorName(){;}                                         \
      /** Evaluate function at an array of values */                    \
      void eval0(const double* const x, int nx, double* f) const ;      \
      /** Evaluate function and first derivative at an array of values */ \
      void eval1(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df) const ;                                    \
      /** Evaluate function and first two derivatives at an array of values */ \
      void eval2(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df,                                            \
                 double* d2f_dxx) const ;                               \
      /** Evaluate function and first thress derivatives at an array of values */ \
      void eval3(const double* const x,                                 \
                 int nx,                                                \
                 double* f,                                             \
                 double* df,                                            \
                 double* d2f_dxx,                                       \
                 double* d3f_dxxx) const ;                                        \
  };                                                                    \
  inline void functorName::eval0(const double* const x, int nx, double* f) const \
  {                                                                     \
    if (checkResults())                                                 \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
          {                                                             \
            f[i] = funcDefinition;                                      \
            TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i]), std::runtime_error, "Non-normal floating point result detected in evaluation of unary functor " << name() << " at argument " << x[i]); \
          }                                                             \
     }                                                                  \
   else                                                                 \
     {                                                                  \
       for (int i=0; i<nx; i++) f[i] = funcDefinition;                  \
     }                                                                  \
}                                                                       \
  inline void functorName::eval1(const double* const x,                 \
                               int nx,                                  \
                               double* f,                               \
                               double* df) const                        \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
             df[i] = firstDerivDefinition;                              \
             TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i]), \
                                std::runtime_error,                           \
                                "Non-normal floating point result detected in " \
                                "evaluation of unary functor "          \
                                << name() << " at argument " << x[i]);  \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
             df[i] = firstDerivDefinition;                              \
           }                                                            \
      }                                                                 \
}                                                                       \
  inline void functorName::eval2(const double* const x,                 \
                               int nx,                                  \
                               double* f,                               \
                               double* df,                              \
                               double* d2f) const                       \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
              df[i] = firstDerivDefinition;                             \
              d2f[i] = secondDerivDefinition;                           \
              TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i])||Teuchos::ScalarTraits<double>::isnaninf(d2f[i]), \
                                 std::runtime_error,                          \
                                 "Non-normal floating point result detected in " \
                                 "evaluation of unary functor "         \
                                 << name() << " at argument " << x[i] ); \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
               df[i] = firstDerivDefinition;                            \
               d2f[i] = secondDerivDefinition;                          \
           }                                                            \
      }                                                                 \
}                                                                     \
  inline void functorName::eval3(const double* const x,                 \
                                 int nx,                                \
                                 double* f,                             \
                                 double* df,                            \
                                 double* d2f,                           \
                                 double* d3f) const                     \
{                                                                       \
  if (checkResults())                                                   \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
              df[i] = firstDerivDefinition;                             \
              d2f[i] = secondDerivDefinition;                           \
              d3f[i] = thirdDerivDefinition;                           \
              TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<double>::isnaninf(f[i])||Teuchos::ScalarTraits<double>::isnaninf(df[i])||Teuchos::ScalarTraits<double>::isnaninf(d2f[i])||Teuchos::ScalarTraits<double>::isnaninf(d3f[i]), \
                                 std::runtime_error,                          \
                                 "Non-normal floating point result detected in " \
                                 "evaluation of unary functor "         \
                                 << name() << " at argument " << x[i] ); \
           }                                                            \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (int i=0; i<nx; i++)                                        \
           {                                                            \
             f[i] = funcDefinition;                                     \
               df[i] = firstDerivDefinition;                            \
               d2f[i] = secondDerivDefinition;                          \
               d3f[i] = thirdDerivDefinition;                          \
           }                                                            \
      }                                                                 \
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  

#endif
