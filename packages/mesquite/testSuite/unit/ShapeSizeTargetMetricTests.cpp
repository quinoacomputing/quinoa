#define TARGET_TEST_GROUP "ShapeSizeTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "AWShapeSizeB1.hpp"
#include "TShapeSize2DB2.hpp"
#include "TShapeSize2DNB1.hpp"
#include "TShapeSize2DNB2.hpp"
#include "TShapeSize3DNB1.hpp"
#include "TShapeSize3DB2.hpp"
#include "TShapeSize3DB4.hpp"
#include "TShapeSizeB1.hpp"
#include "TShapeSizeB3.hpp"
#include "TShapeSizeNB3.hpp"

//                            NAME     !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_GRAD   ( AWShapeSizeB1,  false,false, true,true , 0.0 );
TEST_METRIC_WITH_HESS_2D( TShapeSize2DB2, false,false, true,true , 0.0 );
TEST_METRIC_WITH_HESS_2D( TShapeSize2DNB1,false,false, true,false, 0.0 );
TEST_METRIC_WITH_HESS_2D( TShapeSize2DNB2,false,false, true,false, 0.0 );
TEST_METRIC_WITH_HESS_3D( TShapeSize3DB2, false,false, true,true , 0.0 );
TEST_METRIC_WITH_HESS_3D( TShapeSize3DB4, false,false, true,true , 0.0 );
TEST_METRIC_WITH_HESS_3D( TShapeSize3DNB1,false,false, true,false, 0.0 );
TEST_METRIC_WITH_HESS   ( TShapeSizeB1,   false,false, true,true , 0.0 );
TEST_METRIC_WITH_HESS   ( TShapeSizeB3,   false,false, true,true , 0.0 );
TEST_METRIC_WITH_HESS   ( TShapeSizeNB3,  false,false, true,false, 0.0 );

