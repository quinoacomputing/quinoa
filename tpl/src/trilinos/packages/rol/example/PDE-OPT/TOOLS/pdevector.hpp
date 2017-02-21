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

#ifndef ROL_PDEOPT_PDEVECTOR_HPP
#define ROL_PDEOPT_PDEVECTOR_HPP

#include "ROL_TpetraMultiVector.hpp"
#include "ROL_StdVector.hpp"

#include "assembler.hpp"
#include "solver.hpp"

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    Teuchos::RCP<Teuchos::Time> PDEVectorSimRieszConstruct    = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Construction Time");
    Teuchos::RCP<Teuchos::Time> PDEVectorSimRieszApply        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Application Time");
    Teuchos::RCP<Teuchos::Time> PDEVectorSimRieszSolve        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Sim Riesz Solver Solution Time");
    Teuchos::RCP<Teuchos::Time> PDEVectorOptRieszConstruct    = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Construction Time");
    Teuchos::RCP<Teuchos::Time> PDEVectorOptRieszApply        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Application Time");
    Teuchos::RCP<Teuchos::Time> PDEVectorOptRieszSolve        = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: PDE Vector Opt Riesz Solver Solution Time");
  }
}
#endif


template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_PrimalSimVector;

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_DualSimVector;

template <class Real, class LO, class GO, class Node>
class PDE_PrimalSimVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    Teuchos::RCP<Tpetra::CrsMatrix<> > RieszMap_;
    Teuchos::RCP<Tpetra::MultiVector<> > lumpedRiesz_;
    Teuchos::RCP<Solver<Real> > solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable Teuchos::RCP<PDE_DualSimVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),1));
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void applyRiesz(const Teuchos::RCP<Tpetra::MultiVector<> > &out,
                    const Teuchos::RCP<const Tpetra::MultiVector<> > &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszApply);
      #endif
      if ( useRiesz_ ) {
        if (useLumpedRiesz_) {
          out->elementWiseMultiply(static_cast<Real>(1), *(lumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          RieszMap_->apply(*in,*out);
        }
      }
      else {
        out->scale(static_cast<Real>(1),*in);
      }
    }

  public:
    virtual ~PDE_PrimalSimVector() {}

    PDE_PrimalSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(Teuchos::null),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_PrimalSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<Assembler<Real> > &assembler,
                        Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Sim").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Sim").get("Lump Riesz Map", false);
      assembler->assemblePDERieszMap1(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != Teuchos::null);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
        }
        else {
          solver_ = Teuchos::rcp(new Solver<Real>(parlist.sublist("Solver")));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_PrimalSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<Tpetra::CrsMatrix<> > &RieszMap,
                        const Teuchos::RCP<Solver<Real> > &solver,
                        const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isDualInitialized_(false) {
      if (RieszMap_ != Teuchos::null) {
        useLumpedRiesz_ = (lumpedRiesz_ != Teuchos::null);
        useRiesz_ = (solver_ != Teuchos::null) || useLumpedRiesz_;
      }
      else {
        useLumpedRiesz_ = false;
        useRiesz_ = false;
      }
    }

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      applyRiesz(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return Teuchos::rcp(new PDE_PrimalSimVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             RieszMap_, solver_, lumpedRiesz_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = Teuchos::rcp(new PDE_DualSimVector<Real,LO,GO,Node>(
                    Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                    RieszMap_, solver_, lumpedRiesz_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyRiesz(dual_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *dual_vec_;
    }
}; // class PDE_PrimalSimVector

template <class Real, class LO, class GO, class Node>
class PDE_DualSimVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    Teuchos::RCP<Tpetra::CrsMatrix<Real> > RieszMap_;
    Teuchos::RCP<Tpetra::MultiVector<> > lumpedRiesz_;
    Teuchos::RCP<Tpetra::MultiVector<> > recipLumpedRiesz_;
    Teuchos::RCP<Solver<Real> > solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable Teuchos::RCP<PDE_PrimalSimVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = Teuchos::rcp(new Tpetra::Vector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap()),1);
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void invertLumpedRiesz(void) {
      recipLumpedRiesz_ = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),1));
      recipLumpedRiesz_->reciprocal(*lumpedRiesz_);
    }

    void applyRiesz(const Teuchos::RCP<Tpetra::MultiVector<> > &out,
                    const Teuchos::RCP<const Tpetra::MultiVector<> > &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszSolve);
      #endif
      if ( useRiesz_ ) {
        if (useLumpedRiesz_) {
          out->elementWiseMultiply(static_cast<Real>(1), *(recipLumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          solver_->solve(out,in,false);
        }
      }
      else {
        out->scale(static_cast<Real>(1),*in);
      }
    }

  public:
    virtual ~PDE_DualSimVector() {}

    PDE_DualSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<PDE<Real> > &pde,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(Teuchos::null),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_DualSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<PDE<Real> > &pde,
                      const Teuchos::RCP<Assembler<Real> > &assembler,
                      Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorSimRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Sim").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Sim").get("Lump Riesz Map", false);
      assembler->assemblePDERieszMap1(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != Teuchos::null);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
          invertLumpedRiesz();
        }
        else {
          solver_ = Teuchos::rcp(new Solver<Real>(parlist.sublist("Solver")));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_DualSimVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<Tpetra::CrsMatrix<> > &RieszMap,
                      const Teuchos::RCP<Solver<Real> > &solver,
                      const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isDualInitialized_(false) {
      if (RieszMap_ != Teuchos::null) {
        useLumpedRiesz_ = (lumpedRiesz_ != Teuchos::null);
        useRiesz_ = (solver_ != Teuchos::null) || useLumpedRiesz_;
        if (useLumpedRiesz_) {
          invertLumpedRiesz();
        }
      }
      else {
        useLumpedRiesz_ = false;
        useRiesz_ = false;
      }
    }

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > &ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      applyRiesz(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return Teuchos::rcp(new PDE_DualSimVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             RieszMap_, solver_, lumpedRiesz_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = Teuchos::rcp(new PDE_PrimalSimVector<Real,LO,GO,Node>(
                      Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                      RieszMap_, solver_, lumpedRiesz_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyRiesz(primal_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *primal_vec_;
    }
}; // class PDE_DualSimVector

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_PrimalOptVector;

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_DualOptVector;

template <class Real, class LO, class GO, class Node>
class PDE_PrimalOptVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    Teuchos::RCP<Tpetra::CrsMatrix<> > RieszMap_;
    Teuchos::RCP<Tpetra::MultiVector<> > lumpedRiesz_;
    Teuchos::RCP<Solver<Real> > solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable Teuchos::RCP<PDE_DualOptVector<Real> > dual_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),1));
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void applyRiesz(const Teuchos::RCP<Tpetra::MultiVector<> > &out,
                    const Teuchos::RCP<const Tpetra::MultiVector<> > &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszApply);
      #endif
      if ( useRiesz_ ) {
        if (useLumpedRiesz_) {
          out->elementWiseMultiply(static_cast<Real>(1), *(lumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          RieszMap_->apply(*in,*out);
        }
      }
      else {
        out->scale(static_cast<Real>(1),*in);
      }
    }

  public:
    virtual ~PDE_PrimalOptVector() {}

    PDE_PrimalOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(Teuchos::null),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_PrimalOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<PDE<Real> > &pde,
                        const Teuchos::RCP<Assembler<Real> > &assembler,
                        Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Opt").get("Lump Riesz Map", false);
      assembler->assemblePDERieszMap2(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != Teuchos::null);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
        }
        else {
          solver_ = Teuchos::rcp(new Solver<Real>(parlist.sublist("Solver")));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_PrimalOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                        const Teuchos::RCP<Tpetra::CrsMatrix<> > &RieszMap,
                        const Teuchos::RCP<Solver<Real> > &solver,
                        const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isDualInitialized_(false) {
      if (RieszMap_ != Teuchos::null) {
        useLumpedRiesz_ = (lumpedRiesz_ != Teuchos::null);
        useRiesz_ = (solver_ != Teuchos::null) || useLumpedRiesz_;
      }
      else {
        useLumpedRiesz_ = false;
        useRiesz_ = false;
      }
    }

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      applyRiesz(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();
      return Teuchos::rcp(new PDE_PrimalOptVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             RieszMap_, solver_, lumpedRiesz_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        dual_vec_ = Teuchos::rcp(new PDE_DualOptVector<Real,LO,GO,Node>(
                    Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                    RieszMap_, solver_, lumpedRiesz_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyRiesz(dual_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *dual_vec_;
    }
}; // class PDE_PrimalOptVector

template <class Real, class LO, class GO, class Node>
class PDE_DualOptVector : public ROL::TpetraMultiVector<Real,LO,GO,Node> {
  private:
    Teuchos::RCP<Tpetra::CrsMatrix<Real> > RieszMap_;
    Teuchos::RCP<Tpetra::MultiVector<> > lumpedRiesz_;
    Teuchos::RCP<Tpetra::MultiVector<> > recipLumpedRiesz_;
    Teuchos::RCP<Solver<Real> > solver_;

    bool useRiesz_;
    bool useLumpedRiesz_;

    mutable Teuchos::RCP<PDE_PrimalOptVector<Real> > primal_vec_;
    mutable bool isDualInitialized_;

    void lumpRiesz(void) {
      lumpedRiesz_ = Teuchos::rcp(new Tpetra::Vector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap()),1);
      Tpetra::MultiVector<Real,LO,GO,Node> ones(ROL::TpetraMultiVector<Real>::getMap(),1);
      ones.putScalar(static_cast<Real>(1));
      RieszMap_->apply(ones, *lumpedRiesz_);
    }

    void invertLumpedRiesz(void) {
      recipLumpedRiesz_ = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),1));
      recipLumpedRiesz_->reciprocal(*lumpedRiesz_);
    }

    void applyRiesz(const Teuchos::RCP<Tpetra::MultiVector<> > &out,
                    const Teuchos::RCP<const Tpetra::MultiVector<> > &in) const {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszSolve);
      #endif
      if ( useRiesz_ ) {
        if (useLumpedRiesz_) {
          out->elementWiseMultiply(static_cast<Real>(1), *(recipLumpedRiesz_->getVector(0)), *in, static_cast<Real>(0));
        }
        else {
          solver_->solve(out,in,false);
        }
      }
      else {
        out->scale(static_cast<Real>(1),*in);
      }
    }

  public:
    virtual ~PDE_DualOptVector() {}

    PDE_DualOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<PDE<Real> > &pde,
                      const Teuchos::RCP<Assembler<Real> > &assembler)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), solver_(Teuchos::null),
        useRiesz_(false), useLumpedRiesz_(false), isDualInitialized_(false) {}

    PDE_DualOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<PDE<Real> > &pde,
                      const Teuchos::RCP<Assembler<Real> > &assembler,
                      Teuchos::ParameterList &parlist)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec),
        isDualInitialized_(false) {
      #ifdef ROL_TIMERS
        Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::PDEVectorOptRieszConstruct);
      #endif
      useRiesz_       = parlist.sublist("Vector").sublist("Opt").get("Use Riesz Map", false);
      useLumpedRiesz_ = parlist.sublist("Vector").sublist("Opt").get("Lump Riesz Map", false);
      assembler->assemblePDERieszMap2(RieszMap_, pde);
      useRiesz_ = useRiesz_ && (RieszMap_ != Teuchos::null);
      if (useRiesz_) {
        if (useLumpedRiesz_) {
          lumpRiesz();
          invertLumpedRiesz();
        }
        else {
          solver_ = Teuchos::rcp(new Solver<Real>(parlist.sublist("Solver")));
          solver_->setA(RieszMap_);
        }
      }
    }

    PDE_DualOptVector(const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &tpetra_vec,
                      const Teuchos::RCP<Tpetra::CrsMatrix<> > &RieszMap,
                      const Teuchos::RCP<Solver<Real> > &solver,
                      const Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > &lumpedRiesz)
      : ROL::TpetraMultiVector<Real,LO,GO,Node>(tpetra_vec), RieszMap_(RieszMap),
        lumpedRiesz_(lumpedRiesz), solver_(solver), isDualInitialized_(false) {
      if (RieszMap_ != Teuchos::null) {
        useLumpedRiesz_ = (lumpedRiesz_ != Teuchos::null);
        useRiesz_ = (solver_ != Teuchos::null) || useLumpedRiesz_;
        if (useLumpedRiesz_) {
          invertLumpedRiesz();
        }
      }
      else {
        useLumpedRiesz_ = false;
        useRiesz_ = false;
      }
    }

    Real dot( const ROL::Vector<Real> &x ) const {
      TEUCHOS_TEST_FOR_EXCEPTION( (ROL::TpetraMultiVector<Real,LO,GO,Node>::dimension() != x.dimension()),
                                  std::invalid_argument,
                                  "Error: Vectors must have the same dimension." );
      const Teuchos::RCP<const Tpetra::MultiVector<Real,LO,GO,Node> > &ex
        = Teuchos::dyn_cast<const ROL::TpetraMultiVector<Real,LO,GO,Node> >(x).getVector();
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real>::getVector());
      size_t n = ey.getNumVectors();
      // Scale x with 1/scale_vec_
      Teuchos::RCP<Tpetra::MultiVector<Real,LO,GO,Node> > wex
        = Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(), n));
      applyRiesz(wex,ex);
      // Perform Euclidean dot between *this and scaled x for each vector
      Teuchos::Array<Real> val(n,0);
      ey.dot( *wex, val.view(0,n) );
      // Combine dots for each vector to get a scalar
      Real xy(0);
      for (size_t i = 0; i < n; ++i) {
        xy += val[i];
      }
      return xy;
    }

    Teuchos::RCP<ROL::Vector<Real> > clone() const {
      const Tpetra::MultiVector<Real,LO,GO,Node> &ey
        = *(ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      size_t n = ey.getNumVectors();  
      return Teuchos::rcp(new PDE_DualOptVector<Real,LO,GO,Node>(
             Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
             RieszMap_, solver_, lumpedRiesz_));
    }

    const ROL::Vector<Real> & dual() const {
      if ( !isDualInitialized_ ) {
        // Create new memory for dual vector
        size_t n = ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector()->getNumVectors();
        primal_vec_ = Teuchos::rcp(new PDE_PrimalOptVector<Real,LO,GO,Node>(
                      Teuchos::rcp(new Tpetra::MultiVector<Real,LO,GO,Node>(ROL::TpetraMultiVector<Real>::getMap(),n)),
                      RieszMap_, solver_, lumpedRiesz_));
        isDualInitialized_ = true;
      }
      // Scale *this with scale_vec_ and place in dual vector
      applyRiesz(primal_vec_->getVector(),ROL::TpetraMultiVector<Real,LO,GO,Node>::getVector());
      return *primal_vec_;
    }
}; // class PDE_DualOptVector

template <class Real,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class PDE_OptVector : public ROL::Vector<Real> {
private:
  Teuchos::RCP<ROL::TpetraMultiVector<Real,LO,GO,Node> > vec1_;
  Teuchos::RCP<ROL::StdVector<Real> >                    vec2_;
  mutable Teuchos::RCP<ROL::TpetraMultiVector<Real,LO,GO,Node> > dual_vec1_;
  mutable Teuchos::RCP<ROL::StdVector<Real> >                    dual_vec2_;
  mutable Teuchos::RCP<PDE_OptVector<Real,LO,GO,Node> >          dual_vec_;
  mutable bool isDualInitialized_;

public:
  PDE_OptVector(const Teuchos::RCP<ROL::TpetraMultiVector<Real,LO,GO,Node> > &vec1,
                const Teuchos::RCP<ROL::StdVector<Real> >                    &vec2 ) 
    : vec1_(vec1), vec2_(vec2), isDualInitialized_(false) {

    dual_vec1_ = Teuchos::rcp_dynamic_cast<ROL::TpetraMultiVector<Real,LO,GO,Node> >(vec1_->dual().clone());
    dual_vec2_ = Teuchos::rcp_dynamic_cast<ROL::StdVector<Real> >(vec2_->dual().clone());
  }

  PDE_OptVector(const Teuchos::RCP<ROL::TpetraMultiVector<Real,LO,GO,Node> > &vec)
    : vec1_(vec), vec2_(Teuchos::null), dual_vec2_(Teuchos::null), isDualInitialized_(false) {
    dual_vec1_ = Teuchos::rcp_dynamic_cast<ROL::TpetraMultiVector<Real,LO,GO,Node> >(vec1_->dual().clone());
  }

  PDE_OptVector(const Teuchos::RCP<ROL::StdVector<Real> > &vec)
    : vec1_(Teuchos::null), vec2_(vec), dual_vec1_(Teuchos::null), isDualInitialized_(false) {
    dual_vec2_ = Teuchos::rcp_dynamic_cast<ROL::StdVector<Real> >(vec2_->dual().clone());
  }

  void set( const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    if ( vec1_ != Teuchos::null ) {
      vec1_->set(*(xs.getField()));
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->set(*(xs.getParameter()));
    }
  }

  void plus( const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    if ( vec1_ != Teuchos::null ) {
      vec1_->plus(*(xs.getField()));
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->plus(*(xs.getParameter()));
    }
  }

  void scale( const Real alpha ) {
    if ( vec1_ != Teuchos::null ) {
      vec1_->scale(alpha);
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->scale(alpha);
    }
  }

  void axpy( const Real alpha, const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    if ( vec1_ != Teuchos::null ) {
      vec1_->axpy(alpha,*(xs.getField()));
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->axpy(alpha,*(xs.getParameter()));
    }
  }

  Real dot( const ROL::Vector<Real> &x ) const {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    Real val(0);
    if ( vec1_ != Teuchos::null ) {
      val += vec1_->dot(*(xs.getField()));
    }
    if ( vec2_ != Teuchos::null ) {
      val += vec2_->dot(*(xs.getParameter()));
    }
    return val;
  }

  Real norm() const {
    Real val(0);
    if ( vec1_ != Teuchos::null ) {
      Real norm1 = vec1_->norm();
      val += norm1*norm1;
    }
    if ( vec2_ != Teuchos::null ) {
      Real norm2 = vec2_->norm();
      val += norm2*norm2;
    }
    return std::sqrt(val);
  } 

  Teuchos::RCP<ROL::Vector<Real> > clone(void) const {
    if ( vec2_ == Teuchos::null ) {
      return Teuchos::rcp(new PDE_OptVector<Real,LO,GO,Node>(
             Teuchos::rcp_dynamic_cast<ROL::TpetraMultiVector<Real,LO,GO,Node> >(vec1_->clone())));
    }
    if ( vec1_ == Teuchos::null ) {
      return Teuchos::rcp(new PDE_OptVector<Real,LO,GO,Node>(
             Teuchos::rcp_dynamic_cast<ROL::StdVector<Real> >(vec2_->clone())));
    }
    return Teuchos::rcp(new PDE_OptVector<Real,LO,GO,Node>(
           Teuchos::rcp_dynamic_cast<ROL::TpetraMultiVector<Real,LO,GO,Node> >(vec1_->clone()),
           Teuchos::rcp_dynamic_cast<ROL::StdVector<Real> >(vec2_->clone())));
  }

  const ROL::Vector<Real> & dual(void) const {
    if ( !isDualInitialized_ ) {
      if ( vec1_ == Teuchos::null ) {
        dual_vec_ = Teuchos::rcp(new PDE_OptVector<Real>(dual_vec2_));
      }
      else if ( vec2_ == Teuchos::null ) {
        dual_vec_ = Teuchos::rcp(new PDE_OptVector<Real>(dual_vec1_));
      }
      else {
        dual_vec_ = Teuchos::rcp(new PDE_OptVector<Real>(dual_vec1_,dual_vec2_));
      }
      isDualInitialized_ = true;
    }
    if ( vec1_ != Teuchos::null ) {
      dual_vec1_->set(vec1_->dual());
    }
    if ( vec2_ != Teuchos::null ) {
      dual_vec2_->set(vec2_->dual());
    }
    return *dual_vec_;
  }

  Teuchos::RCP<ROL::Vector<Real> > basis( const int i )  const {
    Teuchos::RCP<ROL::Vector<Real> > e;
    if ( vec1_ != Teuchos::null && vec2_ != Teuchos::null ) {
      int n1 = vec1_->dimension();
      Teuchos::RCP<ROL::Vector<Real> > e1, e2;
      if ( i < n1 ) {
        e1 = vec1_->basis(i);
        e2 = vec2_->clone(); e2->zero();
      }
      else {
        e1 = vec1_->clone(); e1->zero();
        e2 = vec2_->basis(i-n1);
      }
      e = Teuchos::rcp(new PDE_OptVector(
        Teuchos::rcp_dynamic_cast<ROL::TpetraMultiVector<Real> >(e1),
        Teuchos::rcp_dynamic_cast<ROL::StdVector<Real> >(e2)));
    }
    if ( vec1_ != Teuchos::null && vec2_ == Teuchos::null ) {
      int n1 = vec1_->dimension();
      Teuchos::RCP<ROL::Vector<Real> > e1;
      if ( i < n1 ) {
        e1 = vec1_->basis(i);
      }
      else {
        e1->zero();
      }
      e = Teuchos::rcp(new PDE_OptVector(
        Teuchos::rcp_dynamic_cast<ROL::TpetraMultiVector<Real> >(e1)));
    }
    if ( vec1_ == Teuchos::null && vec2_ != Teuchos::null ) {
      int n2 = vec2_->dimension();
      Teuchos::RCP<ROL::Vector<Real> > e2;
      if ( i < n2 ) {
        e2 = vec2_->basis(i);
      }
      else {
        e2->zero();
      }
      e = Teuchos::rcp(new PDE_OptVector(
        Teuchos::rcp_dynamic_cast<ROL::StdVector<Real> >(e2)));
    }
    return e;
  }

  void applyUnary( const ROL::Elementwise::UnaryFunction<Real> &f ) {
    if ( vec1_ != Teuchos::null ) {
      vec1_->applyUnary(f);
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->applyUnary(f);
    }
  }

  void applyBinary( const ROL::Elementwise::BinaryFunction<Real> &f, const ROL::Vector<Real> &x ) {
    const PDE_OptVector<Real> &xs = Teuchos::dyn_cast<const PDE_OptVector<Real> >(x);
    if ( vec1_ != Teuchos::null ) {
      vec1_->applyBinary(f,*xs.getField());
    }
    if ( vec2_ != Teuchos::null ) {
      vec2_->applyBinary(f,*xs.getParameter());
    }
  }

  Real reduce( const ROL::Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    if ( vec1_ != Teuchos::null ) {
      r.reduce(vec1_->reduce(r),result);
    }
    if ( vec2_ != Teuchos::null ) {
      r.reduce(vec2_->reduce(r),result);
    }
    return result;
  }

  int dimension() const {
    int dim(0);
    if ( vec1_ != Teuchos::null ) {
      dim += vec1_->dimension();
    }
    if ( vec2_ != Teuchos::null ) {
      dim += vec2_->dimension();
    }
    return dim;
  }

  Teuchos::RCP<const ROL::TpetraMultiVector<Real> > getField(void) const { 
    return vec1_;
  }

  Teuchos::RCP<const ROL::StdVector<Real> > getParameter(void) const { 
    return vec2_; 
  }

  Teuchos::RCP<ROL::TpetraMultiVector<Real> > getField(void) { 
    return vec1_;
  }

  Teuchos::RCP<ROL::StdVector<Real> > getParameter(void) { 
    return vec2_; 
  }

  void setField(const ROL::Vector<Real>& vec) { 
    vec1_->set(vec);
  }
  
  void setParameter(const ROL::Vector<Real>& vec) { 
    vec2_->set(vec); 
  }
}; // class PDE_OptVector

#endif
