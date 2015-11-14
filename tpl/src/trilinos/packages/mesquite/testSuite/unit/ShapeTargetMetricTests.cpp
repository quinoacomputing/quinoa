#define TARGET_TEST_GROUP "ShapeTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "Mesquite_AWShape2DNB1.hpp"
#include "Mesquite_AWShape2DNB2.hpp"
#include "Mesquite_AWShape2DB1.hpp"
#include "Mesquite_TInverseMeanRatio.hpp"
#include "Mesquite_TShape2DNB2.hpp"
#include "Mesquite_TShape3DB2.hpp"
#include "Mesquite_TShapeB1.hpp"
#include "Mesquite_TShapeNB1.hpp"

//                               NAME       !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_NO_DERIVS_2D( AWShape2DB1,       false, true,  true, true , 0.0 );
TEST_METRIC_NO_DERIVS_2D( AWShape2DNB1,      false, true,  true, false, 0.0 );
TEST_METRIC_WITH_GRAD_2D( AWShape2DNB2,      false, true,  true, false, 0.0 );
TEST_METRIC_WITH_HESS   ( TInverseMeanRatio, false, true,  true, true , 0.0 );
TEST_METRIC_WITH_HESS_2D( TShape2DNB2,       false, true,  true, false, 0.0 );
TEST_METRIC_WITH_HESS_3D( TShape3DB2,        false, true,  true, true , 0.0 );
TEST_METRIC_WITH_HESS   ( TShapeB1,          false, true,  true, true , 0.0 );
TEST_METRIC_WITH_HESS   ( TShapeNB1,         false, true,  true, false, 0.0 );
