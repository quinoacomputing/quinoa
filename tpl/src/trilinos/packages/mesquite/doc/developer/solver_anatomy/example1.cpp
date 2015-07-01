#include "TerminationCriterion.hpp"
#include "FeasibleNewton.hpp"
#include "LPtoPTemplate.hpp"
#include "IdealWeightInverseMeanRatio.hpp"

using namespace Mesquite;

int main()
{
  IdealWeightInverseMeanRatio metric;
  LPtoPTemplate obj_func( 2, &metric );
  
  MsqError err;
  TerminationCriterion inner_term, outer_term;
  inner_term.add_criterion_type_with_double( 
TerminationCriterion::VERTEX_MOVEMENT_ABSOLUTE, 0.01, err );
  inner_term.add_criterion_type_with_int( 
TerminationCriterion::CPU_TIME, 60, err );
  outer_term.add_criterion_type_with_int( 
TerminationCriterion::NUMBER_OF_ITERATES, 1, err );

  FeasibleNewton solver( &obj_func );
  solver.set_inner_termination_criterion( &inner_term );
  solver.set_outer_termination_criterion( &outer_term );
  
  solver.use_global_patch();
  
  return 0;
}

