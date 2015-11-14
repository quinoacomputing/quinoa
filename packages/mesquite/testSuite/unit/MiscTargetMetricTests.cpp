#define TARGET_TEST_GROUP "MiscTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace Mesquite;

#include "Mesquite_TSizeNB1.hpp"
#include "Mesquite_TSizeB1.hpp"
#include "Mesquite_TSquared.hpp"

#include "Mesquite_TOffset.hpp"
#include "Mesquite_TPower2.hpp"
#include "Mesquite_TScale.hpp"
#include "Mesquite_TSum.hpp"

class TOffset_TSizeNB1_2 : public TOffset
{
  public:
  TSizeNB1 mBase;
  TOffset_TSizeNB1_2() : TOffset( 2.0, &mBase ) {}
};

class TPower2_TSizeNB1 : public TPower2
{
  public:
  TSizeNB1 mBase;
  TPower2_TSizeNB1() : TPower2( &mBase ) {}
};

class TScale_TSizeNB1_half : public TScale
{
  public:
  TSizeNB1 mBase;
  TScale_TSizeNB1_half() : TScale( 0.5, &mBase ) {}
};

class TSum_TSize_TSize : public TSum
{
  public:
  TSizeNB1 mu1;
  TSizeB1 mu2;
  TSum_TSize_TSize() : TSum(&mu1,&mu2) {}
};


TEST_NON_QUALITY_METRIC_WITH_HESS( TSquared );

//                           METRIC                NAME    !SHAPE !SIZE !ORIENT BARRIER
TEST_NAMED_METRIC_WITH_HESS( TOffset_TSizeNB1_2,   TOffset, true, false,  true, false, 2.0 );
TEST_NAMED_METRIC_WITH_HESS( TPower2_TSizeNB1,     TPower2, true, false,  true, false, 0.0 );
TEST_NAMED_METRIC_WITH_HESS( TScale_TSizeNB1_half, TScale,  true, false,  true, false, 0.0 );
TEST_NAMED_METRIC_WITH_HESS( TSum_TSize_TSize,     TSum,    true, false,  true,  true, 0.0 );
