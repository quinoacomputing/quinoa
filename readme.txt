DESCRIPTION
-----------

cartesian_product is an extension to the Boost.MPL library and as such
requires a version of the Boost libraries on your system.

Go to http://www.boost.org to download one.

cartesian_product has been tested with Boost 1.36.0 but should be fairly
compatible with older versions of Boost that include the MPL.


INSTALLATION
------------

Just copy the file boost/boost/mpl/cartesian_product.hpp into your boost
installation folder:

	<BOOST_DIR>/boost/mpl/cartesian_product.hpp


Copyright George van Venrooij 2008
Organic Vectory
http://www.organicvectory.com


WHAT IS IT & HOW TO USE
-----------------------
http://tinyurl.com/55wdww


REVISION HISTORY
----------------

2015.07.02  J. Bakosi: Add 'what is it and how to use' section above

2014.02.05  J. Bakosi: Moved

            'template <bool done = true> cartesian_product_outer_impl'

            to immediately after

            'template <bool done = true> struct cartesian_product_inner_impl'

            to avoid compiler errors:

            boost/mpl/cartesian_product.hpp:177:21: error: use of undeclared identifier 'cartesian_product_outer_impl'
                    cartesian_product_outer_impl
                    ^
            boost/mpl/cartesian_product.hpp:179:28: error: no member named 'execute' in the global namespace
                        >::execute

2008.12.05  Renamed algorithm to "cartesian_product" after James Philbin
            correctly identified it as such after posting the algorithm 
            on the boost developers newsgroup.
            
2008.11.21  Initial version named "for_all".
