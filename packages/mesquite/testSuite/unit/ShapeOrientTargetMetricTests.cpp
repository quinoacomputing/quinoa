#define TARGET_TEST_GROUP "ShapeOrientTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "Mesquite_AWShapeOrientNB1.hpp"
#include "Mesquite_TShapeOrientB1.hpp"
#include "Mesquite_TShapeOrientB2.hpp"
#include "Mesquite_TShapeOrientNB1.hpp"
#include "Mesquite_TShapeOrientNB2.hpp"

//                            NAME     !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_GRAD( AWShapeOrientNB1,false, true,false,false, 0.0 );
TEST_METRIC_WITH_HESS( TShapeOrientB1,  false, true,false,true , 0.0 );
TEST_METRIC_WITH_HESS( TShapeOrientB2,  false, true,false,true , 0.0 );
TEST_METRIC_WITH_HESS( TShapeOrientNB1, false, true,false,false, 0.0 );
TEST_METRIC_WITH_HESS( TShapeOrientNB2, false, true,false,false, 0.0 );
