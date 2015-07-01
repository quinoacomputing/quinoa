#define TARGET_TEST_GROUP "UntangleTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "AWUntangleBeta.hpp"
#include "TSizeNB1.hpp"
#include "TShapeSize2DNB1.hpp"
#include "TShapeSize3DNB1.hpp"
#include "TMixed.hpp"
#include "TScale.hpp"
#include "TUntangleBeta.hpp"
#include "TUntangle1.hpp"
#include "TUntangleMu.hpp"

class TUntangleShSz : public TUntangleMu
{
public:
  TShapeSize2DNB1 SS2D;
  TShapeSize3DNB1 SS3D;
  TScale SS2DS; // scale 2D value so that it is sensitive to shape deformation
  TMixed mBase;
  TUntangleShSz() : TUntangleMu(&mBase), SS2DS(10,&SS2D), mBase(&SS2DS,&SS3D) {}
};

class TUntangleSz : public TUntangleMu
{
public:
  TSizeNB1 mBase;
  TUntangleSz() : TUntangleMu(&mBase) {}
};
  

//                               NAME                   !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_GRAD   ( AWUntangleBeta,                 true,  true,  true, false, 0.0 );
TEST_METRIC_WITH_HESS   ( TUntangleBeta,                  true,  true,  true, false, 0.0 );
TEST_METRIC_WITH_HESS   ( TUntangle1,                     true,  true,  true, false, 0.0 );
TEST_NAMED_METRIC_WITH_HESS( TUntangleSz,   TUntangleMu,  true,  false, true, false, 0.0 );
TEST_NAMED_METRIC_WITH_HESS( TUntangleShSz, TUntangleMu, false, false, true, false, 0.0 );
