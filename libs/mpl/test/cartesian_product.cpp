// Based on boost::mpl::for_each.cpp test file,
// Copyright Aleksey Gurtovoy 2000-2008
//------------------------------------------------------------
// This file Copyright George van Venrooij 2008
// http://www.organicvectory.com
//
// Distributed under the Boost Software License, Version 1.0. 
// (See http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mpl/cartesian_product.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>

namespace mpl = boost::mpl;

//  Functor that prints out the type name of the type with which it is called
struct type_printer
{
    type_printer(std::ostream& s) : f_stream(&s) {}
    template< typename U > void operator()(U)
    {
        *f_stream << typeid(U).name() << std::endl;
    }

 private:
    std::ostream* f_stream;
};

int main()
{
    //  Define some type lists
    typedef mpl::list<short,long>       int_types;
    typedef mpl::list<float,double>     real_types;
    typedef mpl::list<char, wchar_t>    char_types;

    //  Combine them
    typedef mpl::list<int_types, real_types, char_types>    type_collection;

    //  Print out all generated combinations
    mpl::cartesian_product< type_collection >(type_printer(std::cout));

    return 0;
}
