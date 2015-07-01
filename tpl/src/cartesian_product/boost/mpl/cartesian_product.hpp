#ifndef BOOST_MPL_CARTESIAN_PRODUCT_HPP_INCLUDED
#define BOOST_MPL_CARTESIAN_PRODUCT_HPP_INCLUDED

// Based on boost::mpl::for_each.hpp,
// Copyright Aleksey Gurtovoy 2000-2008
//------------------------------------------------------------
// This file Copyright George van Venrooij 2008
// http://www.organicvectory.com
//
// Distributed under the Boost Software License, Version 1.0. 
// (See http://www.boost.org/LICENSE_1_0.txt)
//
//  Documentation:
//
//  cartesian_product is a run-time algorithm that works on a
//  sequence of sequences.
//
//  While executing, it constructs a mpl::vector containing
//  one element of each of the sequences and it does so until
//  it exhausts all possible combinations.
//
//  See http://en.wikipedia.org/wiki/Cartesian_product for
//  the official specification.
//
//  For example:
//
//      typedef boost::mpl::vector<short, long>     t1;
//      typedef boost::mpl::vector<char, wchar_t>   t2;
//      typedef boost::mpl::vector<float, double>   t3;
//
//      typedef boost::mpl::vector<t1, t2, t3>      tt;
//
//      template <typename Sequence>
//      void F()
//      {
//      }
//
//      cartesian_product<tt>(F);
//
//      
//      F() is called for the following types:
//
//          boost::mpl::vector<short, char, float>
//          boost::mpl::vector<short, char, double>
//
//          boost::mpl::vector<short, wchar_t, float>
//          boost::mpl::vector<short, wchar_t, double>
//
//          boost::mpl::vector<long, char, float>
//          boost::mpl::vector<long, char, double>
//
//          boost::mpl::vector<long, wchar_t, float>
//          boost::mpl::vector<long, wchar_t, double>
//
//
//  Intended use
//  ------------
//
//  This design evolved out of another design where the behavior of software system was specified
//  using various type sequences. A generic implementation needed to be overloaded for specific 
//  combinations of types and the search for a way to do that indicated the need for such an
//  algorithm.  
//
//  Another use for the algorithm was in testing a complex class where specific functionality of
//  that class could be specified in external functors that are then passed as type parameters to
//  the class template.
//  
//  Using cartesian_product is an easy way to test various configurations in a common manner.
//
//  A short piece of code demonstrates this:
//
//  template <typename T1, typename T2, typename T3>
//  class complex_class
//  {
//  };
//
//  //  Test function that creates a complex class from a type sequence and then tests it
//  template <typename Sequence>
//  void
//  test_complex_class()
//  {
//      typedef complex_class
//          <   at<Sequence, 0>::type
//          ,   at<Sequence, 1>::type
//          ,   at<Sequence, 2>::type
//          >   test_class_t;
//
//      test_class_t  my_class;
//
//      //  Perform some tests on my_class
//  }
//
//  //  Variants of the template parameters
//  typedef boost::mpl::vector<t11, t12, t13>   T1Variants;
//  typedef boost::mpl::vector<t21, t22, t23>   T2Variants;
//  typedef boost::mpl::vector<t31, t32, t33>   T3Variants;
//
//  //  Combined sequence
//  typedef boost::mpl::vector<T1Variants, T2Variants, T3Variants>  TT;
//
//  test()
//  {
//      boost::mpl::cartesian_product<TT>(test_complex_class);
//  }

#include <boost/mpl/is_sequence.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/aux_/unwrap.hpp>

#include <boost/utility/value_init.hpp>

namespace boost
{ 
    namespace mpl
    {
        namespace aux
        {
            //  Stop condition for the recursion
            //  If done == true, then the end of the sequence has been reached and nothing needs to be done
            template <bool done = true>
            struct cartesian_product_inner_impl
            {
                template
                    <   typename CurrentType
                    ,   typename LastType
                    ,   typename CurrentSequence
                    ,   typename LastSequence
                    ,   typename ArgumentSequence   //  Sequence under construction
                    ,   typename F
                    >
                static void execute
                    (   CurrentType*
                    ,   LastType*
                    ,   CurrentSequence*
                    ,   LastSequence*
                    ,   ArgumentSequence*
                    ,   F
                    )
                {
                }
            };

            //  Stop condition for the recursion
            //  If done == true, then the end of the sequence has been reached and nothing needs to be done
            template <bool done = true>
            struct cartesian_product_outer_impl
            {
                template
                    <   typename CurrentSequence
                    ,   typename LastSequence
                    ,   typename ArgumentSequence
                    ,   typename F
                    >
                static void execute
                    (   CurrentSequence*
                    ,   LastSequence*
                    ,   ArgumentSequence*
                    ,   F                   f
                    )
                {
                    value_initialized<ArgumentSequence> x;
                    aux::unwrap(f, 0)(boost::get(x));
                }
            };

            //  The specialization takes care of the case where the recursion is NOT yet done
            template<>
            struct cartesian_product_inner_impl<false>
            {
                template
                    <   typename CurrentType        //  Iterator to current type in sequence
                    ,   typename LastType           //  Iterator to last/end type in sequence
                    ,   typename CurrentSequence    //  Iterator to current outer sequence
                    ,   typename LastSequence       //  Iterator to last/end outer seuqence
                    ,   typename ArgumentSequence   //  Sequence under construction
                    ,   typename F                  //  Functor to call
                    >
                static void execute
                    (   CurrentType*
                    ,   LastType*
                    ,   CurrentSequence*
                    ,   LastSequence*
                    ,   ArgumentSequence*
                    ,   F f
                    )
                {
                    //  Retrieve the type from the sequence
                    typedef typename deref<CurrentType>::type item;

                    //  Add it to the argument sequence
                    typedef typename push_back<ArgumentSequence, item>::type    NewSequence;

                    //  Step to next element in the sequence
                    typedef typename mpl::next<CurrentType>::type NextType;

                    //  Start working on the next sequence
                    cartesian_product_outer_impl
                        <   boost::is_same<CurrentSequence, LastSequence>::value
                        >::execute
                            (   (CurrentSequence*)  0
                            ,   (LastSequence*)     0
                            ,   (NewSequence*)      0
                            ,                       f
                            );

                    //  Continue with next element of current sequence
                    cartesian_product_inner_impl
                        <   boost::is_same<NextType, LastType>::value
                        >::execute
                        (   (NextType*)         0
                        ,   (LastType*)         0
                        ,   (CurrentSequence*)  0
                        ,   (LastSequence*)     0
                        ,   (ArgumentSequence*) 0
                        ,                       f
                        );
                }
            };

            //  The specialization takes care of the case where the recursion is NOT yet done
            // 
            //  This is called for every outer sequence
            template<>
            struct cartesian_product_outer_impl<false>
            {
                template
                    <   typename CurrentSequence    //  Iterator to current outer sequence
                    ,   typename LastSequence       //  Iterator to last outer sequence
                    ,   typename ArgumentSequence   //  Sequence under construction
                    ,   typename F                  //  Functor
                    >
                static void execute
                    (   CurrentSequence*
                    ,   LastSequence*
                    ,   ArgumentSequence*
                    ,   F f
                    )
                {
                    //  Retrieve the inner sequence from the current outer sequence
                    typedef typename deref<CurrentSequence>::type InnerSequence;
                
                    //  Check if the inner sequence passed into this call is actually a sequence
                    BOOST_MPL_ASSERT(( is_sequence<InnerSequence> ));

                    //  Step to next element in the sequence
                    typedef typename mpl::next<CurrentSequence>::type NextSequence;
                
                    //  Inner sequence iterators
                    typedef typename begin<InnerSequence>::type FirstType;
                    typedef typename end  <InnerSequence>::type LastType;

                    //  Process inner sequence types
                    cartesian_product_inner_impl
                        <   boost::is_same<FirstType, LastType>::value
                        >::execute
                        (   (FirstType*)        0
                        ,   (LastType*)         0
                        ,   (NextSequence*)     0
                        ,   (LastSequence*)     0
                        ,   (ArgumentSequence*) 0
                        ,                       f
                        );
                }
            };

        }   //  namespace aux

        //  The cartesian_product algorithm
        template
            <   typename SequenceOfSequences
            ,   typename F
            >
        inline
        void cartesian_product
            (   F                       f
            ,   SequenceOfSequences*        = 0
            )
        {
            //  Check if the sequence passed into this call is actually a sequence
            BOOST_MPL_ASSERT(( is_sequence<SequenceOfSequences> ));

            //  Get iterators to start and end of sequence
            typedef typename begin<SequenceOfSequences>::type   FirstSequence;
            typedef typename end  <SequenceOfSequences>::type   LastSequence;

            //  Use recursion to iterate over all outer sequence elements
            //  The recursion stops as soon first == last
            aux::cartesian_product_outer_impl
                <   boost::is_same<FirstSequence, LastSequence>::value 
                >::execute
                (   (FirstSequence*)0
                ,   (LastSequence*) 0
                ,   (vector0<>*)    0
                ,                   f
                );
        }

    }   //  namespace mpl
}   //  namespace boost
#endif  //  BOOST_MPL_CARTESIAN_PRODUCT_HPP_INCLUDED
