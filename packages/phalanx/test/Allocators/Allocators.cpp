// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "mock_object.hpp"

template<typename T>
bool arrayCompare(Teuchos::ArrayRCP<T> a, Teuchos::ArrayRCP<T> b, double tol)
{
  TEUCHOS_TEST_FOR_EXCEPTION(a.size() != b.size(), std::logic_error, "The arrays in checkTolerance() are not of the same size!");

  double error = 0.0;

  for (typename Teuchos::ArrayRCP<T>::Ordinal i = 0; i < a.size(); ++i)
    error += fabs(a[i] - b[i]);

  if (error > tol)
    return true;
  
  return false;
}

template<int size>
class Alignment {
  char a[size];
};


int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  GlobalMPISession mpi_session(&argc, &argv);

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Allocator Testing
    // *********************************************************************
    {
      cout << "Testing  NewAllocator:" << endl;

      cout << "  Testing ctor...";
      NewAllocator na;
      cout << "passed!" << endl;

      cout << "  Testing allocate()...";
      ArrayRCP<double> vec = na.allocate<double>(4);
      TEUCHOS_TEST_FOR_EXCEPTION(vec.size() != 4, std::runtime_error, 
			 "Allocator is broken!");
      cout << "passed!" << endl;

    }

    {
      cout << "\nTesting ContiguousAllocator:" << endl;
      cout << "  Testing ctor...";
      ContiguousAllocator<double> ca;
      cout << "passed!" << endl;

      cout << "  Testing addRequiredChunk()...";
      const int size = 10;
      const int num_bytes = size * sizeof(double);
      ca.addRequiredChunk(sizeof(double), size);
      ca.addRequiredChunk(sizeof(double), size);
      ca.addRequiredChunk(sizeof(double), size);
      ca.addRequiredChunk(sizeof(double), 2 * size);
      cout << "passed!" << endl;
      
      cout << "  Testing getTotalBytes()...";
      const int total_bytes = ca.getTotalBytes();
      TEUCHOS_TEST_FOR_EXCEPTION(total_bytes != 5 * num_bytes, std::logic_error,
			 "addRequiredBytes() is failing!");
      cout << "passed!" << endl;
            
      cout << "  Testing setup()...";
      ca.setup();
      cout << "passed!" << endl;

      cout << "  Testing allocate()...";
      ArrayRCP<double> a = ca.allocate<double>(size);
      ArrayRCP<double> b = ca.allocate<double>(size);
      ArrayRCP<double> c = ca.allocate<double>(size);
      ArrayRCP<double> d = ca.allocate<double>(2 * size);
      TEUCHOS_TEST_FOR_EXCEPTION(a.size() != size, std::runtime_error, 
			 "Allocator is broken!");
      TEUCHOS_TEST_FOR_EXCEPTION(d.size() != 2 * size, std::runtime_error, 
			 "Allocator is broken!");
      cout << "passed!" << endl;

      
      cout << "  Testing array integrity...";
      
      for (ArrayRCP<double>::Ordinal i = 0; i < a.size(); ++i)
	a[i] = 1.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < b.size(); ++i)
	b[i] = 2.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < c.size(); ++i)
	c[i] = 3.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < d.size(); ++i)
	d[i] = 4.0;

//       for (ArrayRCP<double>::Ordinal i = 0; i < a.size(); ++i)
// 	cout << "a[" << i << "] = " << a[i] << endl;
//       for (ArrayRCP<double>::Ordinal i = 0; i < b.size(); ++i)
// 	cout << "b[" << i << "] = " << b[i] << endl;
//       for (ArrayRCP<double>::Ordinal i = 0; i < c.size(); ++i)
// 	cout << "c[" << i << "] = " << c[i] << endl;
//       for (ArrayRCP<double>::Ordinal i = 0; i < d.size(); ++i)
// 	cout << "d[" << i << "] = " << d[i] << endl;

      ArrayRCP<double> a_ref = arcp<double>(size);
      ArrayRCP<double> b_ref = arcp<double>(size);
      ArrayRCP<double> c_ref = arcp<double>(size);
      ArrayRCP<double> d_ref = arcp<double>(2 * size);

      for (ArrayRCP<double>::Ordinal i = 0; i < a.size(); ++i)
	a_ref[i] = 1.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < b.size(); ++i)
	b_ref[i] = 2.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < c.size(); ++i)
	c_ref[i] = 3.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < d.size(); ++i)
	d_ref[i] = 4.0;

      TEUCHOS_TEST_FOR_EXCEPTION(arrayCompare(a, a_ref, 1.0e-8), std::runtime_error,
			 "Array Integrity failed.  Arrays are different!");
      TEUCHOS_TEST_FOR_EXCEPTION(arrayCompare(b, b_ref, 1.0e-8), std::runtime_error,
			 "Array Integrity failed.  Arrays are different!");
      TEUCHOS_TEST_FOR_EXCEPTION(arrayCompare(c, c_ref, 1.0e-8), std::runtime_error,
			 "Array Integrity failed.  Arrays are different!");
      TEUCHOS_TEST_FOR_EXCEPTION(arrayCompare(d, d_ref, 1.0e-8), std::runtime_error,
			 "Array Integrity failed.  Arrays are different!");
      
      cout << "passed!" << endl;


      cout << "  Testing reset()...";
      ca.reset();
      cout << "passed!" << endl;

    }


    //! This is a mock test to verify allocated objects have their dtor
    //called correctly
    {
      cout << "  Testing that ArrayRCP casted objects are destroyed properly:" 
	   << endl;
      {
	ArrayRCP<MockObject<int> > mo_block_1;
	ArrayRCP<MockObject<int> > mo_block_2;
	
	{
	  const int array_size = 10;
	  const int sizeof_mock_obj = sizeof(MockObject<int>);
	  int total_num_bytes = array_size * sizeof_mock_obj;
	  ArrayRCP<char> raw_array = arcp<char>(total_num_bytes);
	  
	  int length = 5 * sizeof_mock_obj;
	  ArrayRCP<char> view_1 = raw_array.persistingView(0, length);
	  ArrayRCP<char> view_2 = raw_array.persistingView(length,length);
	  
	  MockObject<int> mo(8);
	  
	  mo_block_1 = arcp_reinterpret_cast_nonpod<MockObject<int> >(view_1);
	  
	  mo_block_2 = arcp_reinterpret_cast_nonpod<MockObject<int> >(view_2, mo);
	}
      
      }
      print_num_ctor_dtor_calls<int>(std::cout);
      
      TEUCHOS_TEST_FOR_EXCEPTION(MockObject<int>::m_num_times_empty_ctor_called != 1,
			 std::runtime_error,
			 "Empty ctors are not being called correctly!");
      
      TEUCHOS_TEST_FOR_EXCEPTION(MockObject<int>::m_num_times_sized_ctor_called != 1,
			 std::runtime_error,
			 "Sized ctors are not being called correctly!");
      
      TEUCHOS_TEST_FOR_EXCEPTION(MockObject<int>::m_num_times_copy_ctor_called != 10,
			 std::runtime_error,
			 "Copy ctors are not being called correctly!");
      
      TEUCHOS_TEST_FOR_EXCEPTION(MockObject<int>::m_num_times_dtor_called != 12,
			 std::runtime_error,
			 "dtors are not being called correctly!");
      
      cout << "  Testing that ArrayRCP casted objects are destroyed "
	   << "properly...passed" << endl;
    }


    // Reproduce memory leak due to missing dtors.
    // Fixed by switching to special ArrayRCP ctors.
    {
      ContiguousAllocator<double> ca;
      std::size_t size = 10;
      ca.addRequiredChunk(sizeof(MockObject<int>), size); 
      ca.addRequiredChunk(sizeof(MockObject<int>), size);   
      ca.setup();
      ArrayRCP<MockObject<int> > a = ca.allocate<MockObject<int> >(size);  
      ArrayRCP<MockObject<int> > b = ca.allocate<MockObject<int> >(size);      
      
      MockObject<int> mo(8);
      for (ArrayRCP<MockObject<int> >::iterator i = a.begin();
	   i != a.end(); ++i)
	*i = mo;      
    }

    // Alignment
    {
      cout << "  Testing that Alignment works correctly";

      const std::size_t num_elements_a = 8;
      const std::size_t num_elements_b = 1;
      const std::size_t num_elements_c = 4;
      
      //  Align on size 3
      ContiguousAllocator<Alignment<3> > ca;

      Alignment<2> two;
      Alignment<8> eight;
      ca.addRequiredChunk(sizeof(two), num_elements_a); 
      ca.addRequiredChunk(sizeof(two), num_elements_b); 
      ca.addRequiredChunk(sizeof(eight), num_elements_c);   
      ca.setup();
      ArrayRCP<Alignment<2> > a = ca.allocate<Alignment<2> >(num_elements_a);
      ArrayRCP<Alignment<2> > b = ca.allocate<Alignment<2> >(num_elements_b);
      ArrayRCP<Alignment<8> > c = ca.allocate<Alignment<8> >(num_elements_c);

      TEUCHOS_TEST_FOR_EXCEPTION(a.size() != static_cast<int>(num_elements_a), 
			 std::runtime_error,
			 "Error - mismatch in object sizes - a.size() = "
			 << a.size() << " is not equal to " << num_elements_a);

      TEUCHOS_TEST_FOR_EXCEPTION(b.size() != static_cast<int>(num_elements_b), 
			 std::runtime_error,
			 "Error - mismatch in object sizes - b.size() = "
			 << b.size() << " is not equal to " << num_elements_b);

      TEUCHOS_TEST_FOR_EXCEPTION(c.size() != static_cast<int>(num_elements_c), 
			 std::runtime_error,
			 "Error - mismatch in object sizes - c.size() = "
			 << c.size() << " is not equal to " << num_elements_c);

      cout << "...passed" << endl;
    }


    
    // *********************************************************************
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
