#define TARGET_TEST_GROUP "ShapeSizeTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "Mesquite_AWShapeSizeB1.hpp"
#include "Mesquite_TShapeSize2DB2.hpp"
#include "Mesquite_TShapeSize2DNB1.hpp"
#include "Mesquite_TShapeSize2DNB2.hpp"
#include "Mesquite_TShapeSize3DNB1.hpp"
#include "Mesquite_TShapeSize3DB2.hpp"
#include "Mesquite_TShapeSize3DB4.hpp"
#include "Mesquite_TShapeSizeB1.hpp"
#include "Mesquite_TShapeSizeB3.hpp"
#include "Mesquite_TShapeSizeNB3.hpp"

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

