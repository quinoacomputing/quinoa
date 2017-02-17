/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/diag/EntityKey.hpp>
#include <stk_search/diag/IdentProc.hpp>
#include <stk_search_util/stk_mesh/PrintBoundingBox.hpp>

using namespace use_case;

using namespace stk_classic::diag;


typedef stk_classic::search::ident::IdentProc<uint64_t,unsigned> IdentProc;
typedef stk_classic::search::box::AxisAlignedBoundingBox<IdentProc,float,3> BoundingBox;
typedef std::vector<std::pair<IdentProc, IdentProc> > IdentProcRelation;

static const int box_count = 100;

void
use_case_4_driver(
  stk_classic::ParallelMachine  comm)
{
  stk_classic::diag::WriterThrowSafe _write_throw_safe(dw());

  dw().m(LOG_SEARCH) << "Use case 4" << stk_classic::diag::push << stk_classic::diag::dendl;

  int parallel_rank = stk_classic::parallel_machine_rank(comm);
//  int parallel_size = stk_classic::parallel_machine_size(comm);

  std::vector<BoundingBox> range_vector;

  for (int i = parallel_rank*box_count; i < (parallel_rank + 1)*box_count; ++i) {
    float box[6] = {i,0.0f,0.0f,i+1.0f,1.0f,1.0f };

    BoundingBox   range;
    range.key.ident = i;
    range.set_box(box);

    range_vector.push_back(range);
  }

  std::vector<BoundingBox> domain_vector;

  for (int i = parallel_rank*box_count; i < (parallel_rank + 1)*box_count; ++i) {
    float box[6] = {i,0.0f,0.0f,i+1.0f,1.0f,1.0f };

    BoundingBox   domain;
    domain.key.ident = i;
    domain.set_box(box);

    domain_vector.push_back(domain);
  }

  dw().m(LOG_SEARCH) << "range " << range_vector << dendl;
  dw().m(LOG_SEARCH) << "domain " << domain_vector << dendl;

  stk_classic::search::FactoryOrder order;
  order.m_communicator = comm;


  dw().m(LOG_SEARCH) << "Search algorithm " << order.m_algorithm << dendl;

  IdentProcRelation relation;

  stk_classic::search::coarse_search(relation, range_vector, domain_vector,order);

  dw().m(LOG_SEARCH) << "relation " << relation << dendl;

  dw().m(LOG_SEARCH) << stk_classic::diag::pop;
}
