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


#ifndef PLAYA_VECTORTESTER_HPP
#define PLAYA_VECTORTESTER_HPP

#include "PlayaVectorDecl.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "Thyra_TestSpecifier.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace Playa;
using namespace PlayaExprTemplates;
using namespace Teuchos;



namespace Playa
{

  /** 
   * Run comparisons between element-wise calculations and Vector member
   * functions.
   */
  template <class Scalar>
  class VectorTester
  {
  public:
    /** \brief Local typedef for promoted scalar magnitude */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** */
    VectorTester(const VectorSpace<Scalar>& space,
                 const TestSpecifier<Scalar>& spec,
                 const Playa::MPIComm& comm = Playa::MPIComm::world());

    /** */
    bool runAllTests() const ;

    /** */
    bool sumTest() const ;

    /** */
    bool dotStarTest() const ;

    /** */
    bool dotSlashTest() const ;

    /** */
    bool scalarMultTest() const ;

    /** */
    bool overloadedUpdateTest() const ;


  private:

    /** */
    void randomizeVec(Vector<Scalar>& x) const ;

    TestSpecifier<Scalar> spec_;

    VectorSpace<Scalar> space_;

    Playa::MPIComm comm_;

  };

  template <class Scalar> 
  inline VectorTester<Scalar>
  ::VectorTester(const VectorSpace<Scalar>& space,
                 const TestSpecifier<Scalar>& spec,
                 const Playa::MPIComm& comm)
    : spec_(spec), space_(space), comm_(comm)
  {;}

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::runAllTests() const
  {
    bool pass = true;

    pass = sumTest() && pass;
    pass = dotStarTest() && pass;
    pass = dotSlashTest() && pass;
    pass = scalarMultTest() && pass;
    pass = overloadedUpdateTest() && pass;

    return pass;
  }

  template <class Scalar> 
  inline void VectorTester<Scalar>
  ::randomizeVec(Vector<Scalar>& x) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;

    /* do the operation elementwise */
    SequentialIterator<Scalar> i;
    for (i=space_.begin(); i != space_.end(); i++)
      {
        x[i] = 2.0*(drand48()-0.5);
      }    
  }

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::sumTest() const 
  {
    if (spec_.doTest())
      {
        std::cerr << "running vector addition test..." << std::endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        x.zero();
        y.zero();
        cout << "x = " << x << std::endl;
        cout << "y = " << y << std::endl;
        randomizeVec(a);
        randomizeVec(b);
        cout << "a = " << a << std::endl;
        cout << "b = " << b << std::endl;

        /* do the operation elementwise */
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = a[i] + b[i];
          }

        /* do the operation with member functions */
        x = a + b ;

        cout << "op   (a+b)=" << std::endl << x << std::endl;
        cout << "loop (a+b)=" << std::endl << y << std::endl;
	
        double err = (x-y).normInf();

        std::cerr << "|sum error|=" << err << std::endl;
        if (err > spec_.errorTol())
          {
            std::cerr << "vector sum test FAILED: tol = " 
                 << spec_.errorTol() << std::endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            std::cerr << "WARNING: vector sum test could not beat tol = " 
                 << spec_.warningTol() << std::endl;
          }
	
      }
    else
      {
        std::cerr << "skipping vector addition test..." << std::endl;
      }
    std::cerr << "vector addition test PASSED: tol = " 
         << spec_.errorTol() << std::endl;
    return true;
  }

  

  

  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::dotStarTest() const 
  {
    if (spec_.doTest())
      {
        std::cerr << "running vector dotStar test..." << std::endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = a.dotStar(b);

        /* do the operation elementwise */
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = a[i] * b[i];
          }

        double err = (x-y).normInf();

        std::cerr << "|dotStar error|=" << err << std::endl;
        if (err > spec_.errorTol())
          {
            std::cerr << "vector dotStar test FAILED: tol = " 
                 << spec_.errorTol() << std::endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            std::cerr << "WARNING: vector dotStar test could not beat tol = " 
                 << spec_.warningTol() << std::endl;
          }
	
      }
    else
      {
        std::cerr << "skipping vector dotStar test..." << std::endl;
      }
    std::cerr << "vector dotStar test PASSED: tol = " 
         << spec_.errorTol() << std::endl;
    return true;
  }


  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::dotSlashTest() const 
  {
    if (spec_.doTest())
      {
        std::cerr << "running vector dotSlash test..." << std::endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = a.dotSlash(b);


        /* do the operation elementwise */
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = a[i] / b[i];
          }
	
        double err = (x-y).normInf();

        std::cerr << "|dotSlash error|=" << err << std::endl;
        if (err > spec_.errorTol())
          {
            std::cerr << "vector dotSlash test FAILED: tol = " 
                 << spec_.errorTol() << std::endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            std::cerr << "WARNING: vector dotSlash test could not beat tol = " 
                 << spec_.warningTol() << std::endl;
          }
	
      }
    else
      {
        std::cerr << "skipping vector dotSlash test..." << std::endl;
      }
    std::cerr << "vector dotSlash test PASSED: tol = " 
         << spec_.errorTol() << std::endl;
    return true;
  }

  
  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::scalarMultTest() const 
  {
    if (spec_.doTest())
      {
        std::cerr << "running vector scalarMult test..." << std::endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);


        /* do the operation with member functions */
        x = 3.14*a;

        /* do the operation elementwise */
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = 3.14 * a[i];
          }

        double err = (x-y).normInf();

        std::cerr << "|scalarMult error|=" << err << std::endl;
        if (err > spec_.errorTol())
          {
            std::cerr << "vector scalarMult test FAILED: tol = " 
                 << spec_.errorTol() << std::endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            std::cerr << "WARNING: vector scalarMult test could not beat tol = " 
                 << spec_.warningTol() << std::endl;
          }
	
      }
    else
      {
        std::cerr << "skipping vector scalarMult test..." << std::endl;
      }
    std::cerr << "vector scalarMult test PASSED: tol = " 
         << spec_.errorTol() << std::endl;
    return true;
  }
 
  template <class Scalar> 
  inline bool VectorTester<Scalar>
  ::overloadedUpdateTest() const 
  {
    if (spec_.doTest())
      {
        std::cerr << "running vector overloadedUpdate test..." << std::endl;

        Vector<Scalar> a = space_.createMember();
        Vector<Scalar> b = space_.createMember();
        Vector<Scalar> x = space_.createMember();
        Vector<Scalar> y = space_.createMember();
        randomizeVec(a);
        randomizeVec(b);


        /* do the operation with member functions */
        x = 3.14*a + 1.4*b;

        /* do the operation elementwise */
        for (SequentialIterator<Scalar> i=space_.begin(); i!=space_.end(); i++)
          {
            y[i] = 3.14*a[i] + 1.4*b[i];
          }

        double err = (x-y).normInf();

        std::cerr << "|overloadedUpdate error|=" << err << std::endl;
        if (err > spec_.errorTol())
          {
            std::cerr << "vector overloadedUpdate test FAILED: tol = " 
                 << spec_.errorTol() << std::endl;
            return false;
          }
        else if (err > spec_.warningTol())
          {
            std::cerr << "WARNING: vector overloadedUpdate test could not beat tol = " 
                 << spec_.warningTol() << std::endl;
          }
	
      }
    else
      {
        std::cerr << "skipping vector overloadedUpdate test..." << std::endl;
      }
    std::cerr << "vector overloadedUpdate test PASSED: tol = " 
         << spec_.errorTol() << std::endl;
    return true;
  }



}
#endif
