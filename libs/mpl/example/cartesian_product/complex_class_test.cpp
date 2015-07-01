// Copyright George van Venrooij 2008
// http://www.organicvectory.com
//
// Distributed under the Boost Software License, Version 1.0. 
// (See http://www.boost.org/LICENSE_1_0.txt)
//

//  Contains the complex_class test example described in the documentation
//
//  Tested with Visual Studio 2005 (SP1) on Windows XP x64

#include <boost/mpl/cartesian_product.hpp>

#include <iostream>

namespace mpl = boost::mpl;

//  Complex class
template <typename T1, typename T2, typename T3>
class complex_class
{
    public:
        void test()
        {
            std::cout << "complex_class::test()" << std::endl;
        }
};

//  Some dummy types
struct  t11 {}; struct  t12 {}; struct  t13 {};
struct  t21 {}; struct  t22 {}; struct  t23 {};
struct  t31 {}; struct  t32 {}; struct  t33 {};

//  Variations in types
typedef mpl::vector<t11, t12, t13>  T1Variants;
typedef mpl::vector<t21, t22, t23>  T2Variants;
typedef mpl::vector<t31, t32, t33>  T3Variants;

//  Generic test function
template <typename ComplexClass>
void test_complex_class()
{
    //  Create instance
    ComplexClass instance;

    //  Perform tests...
    std::cout << "Generic test: " << typeid(ComplexClass).name() << " ";

    instance.test();
}

//  Specialized test function
template <>
void test_complex_class<complex_class<t11, t22, t33> >()
{
    //  Create instance
    complex_class<t11, t22, t33> instance;

    //  Perform specialized test...
    std::cout << "Specialized test: " << typeid(complex_class<t11, t22, t33>).name() << " ";

    instance.test();
}

//  Functor that determines the type of the complex class and calls a test function for it
struct tester
{
    template <typename Sequence>
    void operator()(Sequence) const
    {
        //  Iterate over the sequence and get the first three elements
        //  Use them to determine the type of the complex class
        typedef typename mpl::begin<Sequence>::type  first;
        typedef typename mpl::next<first>::type      second;
        typedef typename mpl::next<second>::type     third;

        //  Call the test
        test_complex_class
            <   complex_class
                    <   typename mpl::deref<first >::type
                    ,   typename mpl::deref<second>::type
                    ,   typename mpl::deref<third >::type
                    >
            >
            ();
    }
};

int main()
{
    //  Combine the type lists
    typedef mpl::vector<T1Variants, T2Variants, T3Variants>   type_collection;

    //  Print out all generated combinations
    mpl::cartesian_product< type_collection >(tester());

    return 0;
}
