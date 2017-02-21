/*
// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Intrepid2_ArrayTools.hpp"
#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>
#include "Intrepid2_FieldContainer.hpp"
#include "Teuchos_ScalarTraits.hpp"
#ifdef KOKKOS_HAVE_CUDA
#include <cublas.h>
#endif
#include <Teuchos_BLAS.hpp>
#include <Teuchos_as.hpp>

namespace KokkosDenseMat{
  using Kokkos::ALL;
  using Kokkos::pair;
  using Kokkos::parallel_for;
  using Kokkos::subview;	  
	template<typename Scalar, typename DeviceType,typename Layout,int rank>
	struct MultiGemm;
	
template<typename Scalar>
	struct MultiGemm<Scalar,void,void,2>{
	static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, int m, int n, int k, Scalar alpha,
          Scalar* A, Scalar* B,
          Scalar beta, Scalar* C){
	Teuchos::BLAS<int,Scalar>blas;

	blas.GEMM(transA, transB, n, m, k, alpha,
                   B, n,
                   A, k,
                   beta, C, n);
					   
					   
	}
	
};

template<typename Scalar>
	struct MultiGemm<Scalar,void,void,1>{
	static Scalar DDOT(int n, Scalar* A, int incA, Scalar* B, int incB){
	Teuchos::BLAS<int,Scalar>blas;
	Scalar result = blas.DOT(n, A, incA, B, incB);
	return result;				   
	}
	
};

template<typename Scalar>
	struct MultiGemm<Scalar,Kokkos::DefaultExecutionSpace,Kokkos::LayoutLeft,2>{
	static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
          Kokkos::View<Scalar**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> A, Kokkos::View<Scalar**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> B,
          Scalar beta, Kokkos::View<Scalar**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> C){
	Teuchos::BLAS<int,Scalar>blas;
	    const int m = static_cast<int> (C.dimension_0 ()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());

	blas.GEMM(transA, transB, m, n, k, alpha,
                   A.ptr_on_device(), k,
                   B.ptr_on_device(), n,
                   beta, C.ptr_on_device(), n);
					   
					   
	}
	
};

template<typename Scalar>
	struct MultiGemm<Scalar,Kokkos::DefaultExecutionSpace,Kokkos::LayoutRight,2>{
	static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
          Kokkos::View<Scalar**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> A, Kokkos::View<Scalar**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> B,
          Scalar beta, Kokkos::View<Scalar**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> C){
	Teuchos::BLAS<int,Scalar>blas;
	    const int m = static_cast<int> (C.dimension_0 ()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());
   
	blas.GEMM(transA, transB, n, m, k, alpha,
                   B.ptr_on_device(), n,
                   A.ptr_on_device(), k,
                   beta, C.ptr_on_device(), n);
					   
					   
	}
	
};



typedef Kokkos::View<double***, Kokkos::LayoutLeft> left_type;
typedef Kokkos::View<double***, Kokkos::LayoutRight> right_type;
typedef Kokkos::View<double**, Kokkos::LayoutStride> left_subtype;
template<class Scalar>
struct blasOpenMPBatchLeft {
 Teuchos::BLAS<int,Scalar>blas;
 left_type A;
 left_type B;
 left_type C;
  
  

  int msize,nsize,ksize;
  Teuchos::ETransp transA;
  Teuchos::ETransp transB;
  Scalar alpha;
  Scalar beta;
  blasOpenMPBatchLeft (left_type a_, left_type b_, left_type c_, const int m_, const int n_, const int k_,
  Teuchos::ETransp transA_, Teuchos::ETransp transB_, Scalar alpha_, Scalar beta_) :
    A(a_), B(b_), C(c_), msize(m_), nsize(n_), ksize(k_), transA(transA_), transB(transB_),
    alpha(alpha_), beta(beta_)
  {}
 
  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {

	 left_subtype subA=subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subB=subview(B, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subC=subview(C, i, Kokkos::ALL(), Kokkos::ALL());

//#pragma loop(ivdep)

for (unsigned int i = 0u; i != msize; ++i) {
  
    for (unsigned int k = 0u; k != ksize; ++k) {

        const Scalar r = subA(i,k);

        for (unsigned int j = 0u; j != nsize; ++j) {

		subC(i,j)=beta*subC(i,j)+alpha*r*subB(k,j);

        }

    }

} 
                        
  }
};

template<class Scalar>
struct blasOpenMPBatchRight {
 Teuchos::BLAS<int,Scalar>blas;
 right_type A;
 right_type B;
 right_type C;
  
  

  int msize,nsize,ksize;
  Teuchos::ETransp transA;
  Teuchos::ETransp transB;
  Scalar alpha;
  Scalar beta;
  blasOpenMPBatchRight (right_type a_, right_type b_, right_type c_, const int m_, const int n_, const int k_,
  Teuchos::ETransp transA_, Teuchos::ETransp transB_, Scalar alpha_, Scalar beta_) :
    A(a_), B(b_), C(c_), msize(m_), nsize(n_), ksize(k_), transA(transA_), transB(transB_),
    alpha(alpha_), beta(beta_)
  {}
 
  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {

	 left_subtype subA=subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subB=subview(B, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subC=subview(C, i, Kokkos::ALL(), Kokkos::ALL());

//#pragma loop(ivdep)

for (index_type i = 0u; i != static_cast<index_type>(msize); ++i) {
  
    for (index_type k = 0u; k != static_cast<index_type>(ksize); ++k) {

        const Scalar r = subA(i,k);

        for (index_type j = 0u; j != static_cast<index_type>(nsize); ++j) {

		subC(i,j)=beta*subC(i,j)+alpha*r*subB(k,j);

        }

    }

} 
                        
  }
};

/*

template<class Scalar>
struct blasOpenMPLeft {
 Teuchos::BLAS<int,Scalar>blas;
 left_type A;
 left_type B;
 left_type C;
  
  

  int m,n,k;
  Teuchos::ETransp transA;
  Teuchos::ETransp transB;
  Scalar alpha;
  Scalar beta;
  blasOpenMPLeft (left_type a_, left_type b_, left_type c_, const int m_, const int n_, const int k_,
  Teuchos::ETransp transA_, Teuchos::ETransp transB_, Scalar alpha_, Scalar beta_) :
    A(a_), B(b_), C(c_), m(m_), n(n_), k(k_), transA(transA_), transB(transB_),
    alpha(alpha_), beta(beta_)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {

	 left_subtype subA=subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subB=subview(B, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subC=subview(C, i, Kokkos::ALL(), Kokkos::ALL());

 	blas.GEMM(transA, transB, m, n, k, alpha,
                   subA.ptr_on_device(), k,
                   subB.ptr_on_device(), n,
                   beta, subC.ptr_on_device(), n);

       
                   
  }
};*/

template<typename Scalar>
	struct MultiGemm<Scalar,Kokkos::DefaultExecutionSpace,Kokkos::LayoutLeft,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
           Kokkos::View<Scalar***,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> A,  Kokkos::View<Scalar***,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> B,
          Scalar beta, Kokkos::View<Scalar***,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());
   //     printf("m:%d,n:%d,k:%d",m,n,k);
	Kokkos::parallel_for(C.dimension(0),blasOpenMPBatchLeft<Scalar>(A,B,C,m,n,k,transA,transB,alpha,beta));
	}
	
};
template<typename Scalar>
	struct MultiGemm<Scalar,Kokkos::DefaultExecutionSpace,Kokkos::LayoutRight,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
           Kokkos::View<Scalar***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> A,  Kokkos::View<Scalar***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> B,
          Scalar beta, Kokkos::View<Scalar***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());
   //     printf("m:%d,n:%d,k:%d",m,n,k);
	Kokkos::parallel_for(C.dimension(0),blasOpenMPBatchRight<Scalar>(A,B,C,m,n,k,transA,transB,alpha,beta));
	}
	
};

// The Following stuff doesn't compile at all with Cuda. And it can't for example Scalar is not defined because its explicit specialisations.
#if false
	
	struct MultiGemm<double,Kokkos::Cuda,Kokkos::LayoutLeft,2>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
           Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> A,  Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_0()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());
  cublasDgemm(transA, transB, m, n, k, alpha, A.ptr_on_device(), k, B.ptr_on_device(), n, beta, C.ptr_on_device(), n);
	}
	
};
	struct MultiGemm<float,Kokkos::Cuda,Kokkos::LayoutLeft,2>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
           Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> A,  Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_0()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());
  cublasSgemm(transA, transB, m, n, k, alpha, A.ptr_on_device(), k, B.ptr_on_device(), n, beta, C.ptr_on_device(), n);
	}
	
};


	struct MultiGemm<double,Kokkos::Cuda,Kokkos::LayoutRight,2>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
           Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> A,  Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_0()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());

cublasDgemm(transA, transB, n, m, k, alpha, B.ptr_on_device(), n, A.ptr_on_device(), k, beta, C.ptr_on_device(), n);
	}
	
};


	struct MultiGemm<float,Kokkos::Cuda,Kokkos::LayoutRight,2>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
           Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> A,  Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_0()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());

  cublasSgemm(transA, transB, n, m, k, alpha, B.ptr_on_device(), n, A.ptr_on_device(), k, beta, C.ptr_on_device(), n);

	}
	
};

	struct MultiGemm<double,Kokkos::Cuda,Kokkos::LayoutLeft,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
           Kokkos::View<double***,Kokkos::LayoutLeft,Kokkos::Cuda> A,  Kokkos::View<double***,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          double beta, Kokkos::View<double***,Kokkos::LayoutLeft,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ()),
        batchCount=C.dimension_0();
      
cublasDgemmBatched(transA,  transB,
                                  m, n, k,
                                  alpha,
                                  A.ptr_on_device(), k,
                                  B.ptr_on_device(), n,
                                  beta,
                                  C.ptr_on_device(), n, batchCount);

	}
	
};
	
	struct MultiGemm<float,Kokkos::Cuda,Kokkos::LayoutLeft,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
           Kokkos::View<float***,Kokkos::LayoutLeft,Kokkos::Cuda> A,  Kokkos::View<float***,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          float beta, Kokkos::View<float***,Kokkos::LayoutLeft,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ()),
        batchCount=C.dimension_0();
cublasSgemmBatched(transA,  transB,
                                  m, n, k,
                                  alpha,
                                  A.ptr_on_device(), k,
                                  B.ptr_on_device(), n,
                                  beta,
                                  C.ptr_on_device(), n, batchCount);

	}
	
};

	struct MultiGemm<double,Kokkos::Cuda,Kokkos::LayoutRight,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
           Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> A,  Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> B,
          double beta, Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());
#if defined( KOKKOS_HAVE_OPENMP ) 
	Kokkos::parallel_for(C.dimension(0),blasOpenMPBatchRight<Scalar>(A,B,C,m,n,k,transA,transB,alpha,beta));
#endif
	}
	
};




	struct MultiGemm<float,Kokkos::Cuda,Kokkos::LayoutRight,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
           Kokkos::View<float***,Kokkos::LayoutRight,Kokkos::Cuda> A,  Kokkos::View<float***,Kokkos::LayoutRight,Kokkos::Cuda> B,
          float beta, Kokkos::View<float***,Kokkos::LayoutRight,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());
#if defined( KOKKOS_HAVE_OPENMP ) 
	Kokkos::parallel_for(C.dimension(0),blasOpenMPBatchRight<Scalar>(A,B,C,m,n,k,transA,transB,alpha,beta));
#endif
	}
	
};
#endif

}





int main(){

   Kokkos::initialize();
  //initialize viewsto random values

   {
	Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> inputview1("X",5,5);
	Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> inputview2("Y",5,5);
	Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> outputview2("Z",5,5);

  //fill with random numbers
	Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
    Kokkos::fill_random(inputview1,rand_pool64,100);	
    Kokkos::fill_random(inputview2,rand_pool64,100);
    KokkosDenseMat::MultiGemm<double,Kokkos::DefaultExecutionSpace,Kokkos::LayoutLeft,2>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,inputview1,inputview2,1,outputview2);
    /* 
std::cout <<"A:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<inputview1(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<inputview2(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"C:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<outputview2(i,j)<<",";
		
	}
	std::cout <<std::endl;
}*/
}// scope 1

{// scope 2 

	Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> inputview1("X",5,5);
	Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> inputview2("Y",5,5);
	Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> outputview2("Z",5,5);

  //fill with random numbers
	Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
    Kokkos::fill_random(inputview1,rand_pool64,100);	
    Kokkos::fill_random(inputview2,rand_pool64,100);
    KokkosDenseMat::MultiGemm<double,Kokkos::DefaultExecutionSpace,Kokkos::LayoutRight,2>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,inputview1,inputview2,1,outputview2);
/*std::cout <<"A:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<inputview1(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<inputview2(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"C:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<outputview2(i,j)<<",";
		
	}
	std::cout <<std::endl;
}*/
}

//begin scope 4
{
Intrepid2::FieldContainer<double>testcontainerA(5,5);
	for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		testcontainerA(i,j)=Teuchos::ScalarTraits<double>::random();
	}
}
Intrepid2::FieldContainer<double>testcontainerB(5,5);
	for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		testcontainerB(i,j)=Teuchos::ScalarTraits<double>::random();
	}
}
	Intrepid2::FieldContainer<double>testcontainerC(5,5);
	KokkosDenseMat::MultiGemm<double,void,void,2>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,5,5,5,1.0,&testcontainerA[0],&testcontainerB[0],1,&testcontainerC[0]);
	/*
std::cout <<"A:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerA(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerB(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"C:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerC(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<std::endl;	
*/
}//end scope 4

{
	Intrepid2::FieldContainer<double>testcontainerA(2,5,5);
	for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		testcontainerA(i,j,k)=Teuchos::ScalarTraits<double>::random();
	     }
	}
}
Intrepid2::FieldContainer<double>testcontainerB(2,5,5);
	for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		testcontainerB(i,j,k)=Teuchos::ScalarTraits<double>::random();
		}
	}
}
	Intrepid2::FieldContainer<double>testcontainerC(2,5,5);	
for(int i=0;i<2;i++){
	KokkosDenseMat::MultiGemm<double,void,void,2>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,5,5,5,1.0,&testcontainerA[i*5*5],&testcontainerB[i*5*5],1,&testcontainerC[i*5*5]);
}	
/*
std::cout <<"A:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		std::cout <<testcontainerA(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		std::cout <<testcontainerB(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}
	std::cout <<"C:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		std::cout <<testcontainerC(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}	*/
	
}//end scope5


//scope6
{
	Intrepid2::FieldContainer<double>testcontainerA(5);
	for(int i=0;i<5;i++){
		testcontainerA(i)=Teuchos::ScalarTraits<double>::random();
}
Intrepid2::FieldContainer<double>testcontainerB(5);
	for(int i=0;i<5;i++){
		testcontainerB(i)=Teuchos::ScalarTraits<double>::random();
}
/*double result=	KokkosDenseMat::MultiGemm<double,void,void,1>::DDOT(5, &testcontainerA[0], 1, &testcontainerB[0], 1);

std::cout <<"A: "<<std::endl;
	for(int i=0;i<5;i++){
		std::cout <<testcontainerA(i)<<",";
}
std::cout <<std::endl;

std::cout <<"B: "<<std::endl;
	for(int i=0;i<5;i++){
		std::cout <<testcontainerB(i)<<std::endl;
}
std::cout <<std::endl;

std::cout <<"result: "<<result<<std::endl;

*/
}//end scope6
//scope 7
{
Intrepid2::FieldContainer<double>testcontainerA(2,5);
	for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		testcontainerA(i,j)=Teuchos::ScalarTraits<double>::random();
	}
}
Intrepid2::FieldContainer<double>testcontainerB(2,5);
	for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		testcontainerB(i,j)=Teuchos::ScalarTraits<double>::random();
	}
}

Intrepid2::FieldContainer<double>testcontainerC(2);
for(int i=0;i<2;i++){
	testcontainerC(i)=KokkosDenseMat::MultiGemm<double,void,void,1>::DDOT(5, &testcontainerA[i*5], 1, &testcontainerB[i*5], 1);
}
/*
std::cout <<"A:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerA(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerB(i,j)<<std::endl;
		
	}
	std::cout <<std::endl;
}
std::cout <<"C:"<<std::endl;
for(int i=0;i<2;i++){
		std::cout <<testcontainerC(i)<<",";		
	std::cout <<std::endl;
}*/
	
}//end scope 7

{//end scope 8
    Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> inputview1("X",1000,5,6);
	Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> inputview2("Y",1000,6,5);
	Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> outputview2("Z",1000,5,5);

  //fill with random numbers
	Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
    Kokkos::fill_random(inputview1,rand_pool64,100);	
    Kokkos::fill_random(inputview2,rand_pool64,100);
    KokkosDenseMat::MultiGemm<double,Kokkos::DefaultExecutionSpace,Kokkos::LayoutRight,3>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,inputview1,inputview2,1,outputview2);
/*	std::cout <<"A:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<6;k++){
		std::cout <<inputview1(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<6;j++){
		for (int k=0;k<5;k++){
		std::cout <<inputview2(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}
	std::cout <<"C:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		std::cout <<outputview2(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}*/

}

   Kokkos::finalize();
  
  
	
	return 0;
}
