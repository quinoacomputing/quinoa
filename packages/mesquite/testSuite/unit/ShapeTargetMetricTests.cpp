#define TARGET_TEST_GROUP "ShapeTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "AWShape2DNB1.hpp"
#include "AWShape2DNB2.hpp"
#include "AWShape2DB1.hpp"
#include "TInverseMeanRatio.hpp"
#include "TShape2DNB2.hpp"
#include "TShape3DB2.hpp"
#include "TShapeB1.hpp"
#include "TShapeNB1.hpp"

//                               NAME       !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_NO_DERIVS_2D( AWShape2DB1,       false, true,  true, true , 0.0 );
TEST_METRIC_NO_DERIVS_2D( AWShape2DNB1,      false, true,  true, false, 0.0 );
TEST_METRIC_WITH_GRAD_2D( AWShape2DNB2,      false, true,  true, false, 0.0 );
TEST_METRIC_WITH_HESS   ( TInverseMeanRatio, false, true,  true, true , 0.0 );
TEST_METRIC_WITH_HESS_2D( TShape2DNB2,       false, true,  true, false, 0.0 );
TEST_METRIC_WITH_HESS_3D( TShape3DB2,        false, true,  true, true , 0.0 );
TEST_METRIC_WITH_HESS   ( TShapeB1,          false, true,  true, true , 0.0 );
TEST_METRIC_WITH_HESS   ( TShapeNB1,         false, true,  true, false, 0.0 );
