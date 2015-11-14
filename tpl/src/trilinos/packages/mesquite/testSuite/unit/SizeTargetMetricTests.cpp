#define TARGET_TEST_GROUP "SizeTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "Mesquite_AWSizeNB1.hpp"
#include "Mesquite_AWSizeB1.hpp"
#include "Mesquite_TSizeNB1.hpp"
#include "Mesquite_TSizeB1.hpp"
#include "Mesquite_TTau.hpp"


//                     NAME       !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS( AWSizeNB1,  true, false,  true, false, 0.0 );
TEST_METRIC_WITH_GRAD( AWSizeB1,   true, false,  true,  true, 0.0 );
TEST_METRIC_WITH_HESS( TSizeNB1,   true, false,  true, false, 0.0 );
TEST_METRIC_WITH_HESS( TSizeB1,    true, false,  true,  true, 0.0 );


TEST_NON_QUALITY_METRIC_WITH_HESS( TTau );
