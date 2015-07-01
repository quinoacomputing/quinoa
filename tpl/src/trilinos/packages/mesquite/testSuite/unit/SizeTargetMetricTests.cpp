#define TARGET_TEST_GROUP "SizeTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "AWSizeNB1.hpp"
#include "AWSizeB1.hpp"
#include "TSizeNB1.hpp"
#include "TSizeB1.hpp"
#include "TTau.hpp"


//                     NAME       !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS( AWSizeNB1,  true, false,  true, false, 0.0 );
TEST_METRIC_WITH_GRAD( AWSizeB1,   true, false,  true,  true, 0.0 );
TEST_METRIC_WITH_HESS( TSizeNB1,   true, false,  true, false, 0.0 );
TEST_METRIC_WITH_HESS( TSizeB1,    true, false,  true,  true, 0.0 );


TEST_NON_QUALITY_METRIC_WITH_HESS( TTau );
