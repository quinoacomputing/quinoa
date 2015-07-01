# @HEADER
# ************************************************************************
# 
# Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
#                  Copyright (2003) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package moocho subdirectory src built successfully."
	@echo ""

# ------------------------------------------------------------------------
# Install the script that creates a Moocho.opt file automatically.  This
# allows users to easily access this script and to use it to create moocho
# solver options.
# ------------------------------------------------------------------------

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(exec_prefix)/tools/moocho ; \
	$(mkinstalldirs) $(DESTDIR)$(exec_prefix)/tools/moocho/Shared ; \
	$(mkinstalldirs) $(DESTDIR)$(exec_prefix)/tools/moocho/MamaJama ; \
	$(mkinstalldirs) $(DESTDIR)$(exec_prefix)/tools/moocho/IpConfig ; \
	cp $(srcdir)/MoochoPack/configurations/generate-opt-file.pl $(DESTDIR)$(exec_prefix)/tools/moocho/. ; \
	chmod a+x $(exec_prefix)/tools/moocho/generate-opt-file.pl ; \
	cp $(srcdir)/MoochoPack/configurations/Moocho.opt.MoochoSolver $(DESTDIR)$(exec_prefix)/tools/moocho/. ; \
	cp $(srcdir)/MoochoPack/configurations/Shared/Moocho.opt.DecompositionSystemStateStepBuilderStd $(DESTDIR)$(exec_prefix)/tools/moocho/Shared/. ; \
	cp $(srcdir)/MoochoPack/configurations/MamaJama/Moocho.opt.NLPAlgoConfigMamaJama $(DESTDIR)$(exec_prefix)/tools/moocho/MamaJama/. ; \
	cp $(srcdir)/MoochoPack/configurations/IpConfig/Moocho.opt.NLPAlgoConfigIP $(DESTDIR)$(exec_prefix)/tools/moocho/IpConfig/.

# Above, I need to just copy the generate-opt-file.pl file stright over to
# keep its execute permissions!


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#

MOOCHO_UTILITIES_HPP = \
  $(srcdir)/MoochoUtilities/src/BLAS_Cpp_Types.hpp\
  $(srcdir)/MoochoUtilities/src/FortranTypes_CppFortranStrings.hpp\
  $(srcdir)/MoochoUtilities/src/FortranTypes_f_open_file.hpp\
  $(srcdir)/MoochoUtilities/src/InputStreamHelperPack_EatInputComment.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_CommandLineOptionsFromStreamProcessor.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_OptionsFromStreamExceptions.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_OptionsFromStream.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_SetOptionsFromStream.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_SetOptionsFromStreamNode.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_SetOptionsToTargetBase.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_StringToBool.hpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_StringToIntMap.hpp\
  $(srcdir)/MoochoUtilities/src/ProfileHackPack_profile_hack.hpp\
  $(srcdir)/MoochoUtilities/src/ReleaseResource.hpp\
  $(srcdir)/MoochoUtilities/src/ReleaseResource_ref_count_ptr.hpp\
  $(srcdir)/MoochoUtilities/src/RangePack_Range1D.hpp\
  $(srcdir)/MoochoUtilities/src/SerializationPack_Serializable.hpp\
  $(srcdir)/MoochoUtilities/src/SpecialASCII.hpp\
  $(srcdir)/MoochoUtilities/src/StandardAggregationMacros.hpp\
  $(srcdir)/MoochoUtilities/src/StandardCompositionRelationshipsPack.hpp\
  $(srcdir)/MoochoUtilities/src/StopWatchPack_stopwatch.hpp\
  $(srcdir)/MoochoUtilities/src/StrideIterPack_StrideIter.hpp\
  $(srcdir)/MoochoUtilities/src/TestingHelperPack_update_success.hpp

MOOCHO_UTILITIES_CPP = \
  $(srcdir)/MoochoUtilities/src/FortranTypes_CppFortranStrings.cpp\
  $(srcdir)/MoochoUtilities/src/FortranTypes_f_open_file.cpp\
  $(srcdir)/MoochoUtilities/src/InputStreamHelperPack_EatInputComment.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_CommandLineOptionsFromStreamProcessor.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_OptionsFromStream.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_SetOptionsFromStreamNode.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_StringToBool.cpp\
  $(srcdir)/MoochoUtilities/src/OptionsFromStreamPack_StringToIntMap.cpp\
  $(srcdir)/MoochoUtilities/src/ProfileHackPack_profile_hack.cpp\
  $(srcdir)/MoochoUtilities/src/RangePack_Range1D.cpp\
  $(srcdir)/MoochoUtilities/src/StopWatchPack_stopwatch.cpp\
  $(srcdir)/MoochoUtilities/src/TestingHelperPack_update_success.cpp

MOOCHO_UTILITIES_F = \
  $(srcdir)/MoochoUtilities/src/ff_open_file.f\
  $(srcdir)/MoochoUtilities/src/FortranCppStrings.f\
  $(srcdir)/MoochoUtilities/src/FSeconds.f

ITERATION_PACK_HPP = \
  $(srcdir)/IterationPack/src/IterationPack_Algorithm.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmSetOptions.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmState.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmStep.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmTracker.hpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmTrackerComposite.hpp \
  $(srcdir)/IterationPack/src/IterationPack_CastIQMember.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantity.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantityAccess.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantityAccessContiguous.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantityAccessContiguousDecl.hpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantityAccessContiguousDef.hpp \
  $(srcdir)/IterationPack/src/IterationPack_PublicTypes.ud \
  $(srcdir)/IterationPack/src/IterationPack_Types.hpp \
  $(srcdir)/IterationPack/src/IterationPack_cast_iq.hpp \
  $(srcdir)/IterationPack/src/IterationPack_print_algorithm_step.hpp

ITERATION_PACK_CPP = \
  $(srcdir)/IterationPack/src/IterationPack_Algorithm.cpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmSetOptions.cpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmState.cpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmTracker.cpp \
  $(srcdir)/IterationPack/src/IterationPack_AlgorithmTrackerComposite.cpp \
  $(srcdir)/IterationPack/src/IterationPack_CastIQMember.cpp \
  $(srcdir)/IterationPack/src/IterationPack_IterQuantity.cpp \
  $(srcdir)/IterationPack/src/IterationPack_cast_iq.cpp \
  $(srcdir)/IterationPack/src/IterationPack_print_algorithm_step.cpp

RTOP_PACK_BASE_H = \
  $(srcdir)/RTOpPack/src/base/RTOp.h\
  $(srcdir)/RTOpPack/src/base/RTOpToMPI.h\
  $(srcdir)/RTOpPack/src/base/RTOp_MPI_config.h\
  $(srcdir)/RTOpPack/src/base/RTOp_SparseSubVector.h\
  $(srcdir)/RTOpPack/src/base/RTOp_apply_op_mpi.h\
  $(srcdir)/RTOpPack/src/base/RTOp_apply_op_serial.h\
  $(srcdir)/RTOpPack/src/base/RTOp_mpi.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_free_free.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_index_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_null_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_index_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_value_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_values_vtbl.h\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_max_value.h\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_min_value.h\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_sum_value.h\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_sum_values.h\
  $(srcdir)/RTOpPack/src/base/check_nan_inf.h

RTOP_PACK_BASE_C = \
  $(srcdir)/RTOpPack/src/base/RTOp.c\
  $(srcdir)/RTOpPack/src/base/RTOpToMPI.c\
  $(srcdir)/RTOpPack/src/base/RTOp_SparseSubVector.c\
  $(srcdir)/RTOpPack/src/base/RTOp_apply_op_mpi.c\
  $(srcdir)/RTOpPack/src/base/RTOp_apply_op_serial.c\
  $(srcdir)/RTOpPack/src/base/RTOp_mpi.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_free_free.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_index_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_null_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_index_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_value_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_value_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_obj_values_vtbl.c\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_max_value.c\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_min_value.c\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_sum_value.c\
  $(srcdir)/RTOpPack/src/base/RTOp_reduct_sum_values.c

RTOP_PACK_BASE_HPP = \
  $(srcdir)/RTOpPack/src/base/RTOpPack_OldTypes.hpp\
  $(srcdir)/RTOpPack/src/base/RTOpPack_RTOpC.hpp\
  $(srcdir)/RTOpPack/src/base/RTOpPack_RTOpCPostMod.hpp\
  $(srcdir)/RTOpPack/src/base/RTOp_get_reduct_op.hpp\
  $(srcdir)/RTOpPack/src/base/RTOpPack_print_sub_vector.hpp

RTOP_PACK_BASE_CPP = \
  $(srcdir)/RTOpPack/src/base/RTOpPack_RTOpC.cpp\
  $(srcdir)/RTOpPack/src/base/check_nan_inf.cpp\
  $(srcdir)/RTOpPack/src/base/RTOpPack_print_sub_vector.cpp

RTOP_PACK_STD_OPS_H = \
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_combined_nu_comp_err.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_comp_err_with_mu.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_dot_prod.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_find_nan_inf.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_fraction_to_boundary.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_fraction_to_zero_boundary.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_get_ele.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_get_sub_vector.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_log_bound_barrier.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_abs_ele.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_inequ_viol.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_near_feas_step.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_rel_step.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_step.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_norms.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_num_bounded.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_num_nonzeros.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_sum.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_sum_abs.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_add_scalar.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_assign_scalar.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_assign_vectors.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_axpy.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_correct_multipliers.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_divide.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_prod.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_sqrt.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_force_in_bounds.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_inv_of_difference.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_max_abs_vec_scalar.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_max_vec_scalar.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_multiplier_step.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_random_vector.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_scale_vector.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_set_ele.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_set_sub_vector.h\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_sign.h

RTOP_PACK_STD_OPS_C = \
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_combined_nu_comp_err.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_comp_err_with_mu.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_dot_prod.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_find_nan_inf.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_fraction_to_boundary.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_fraction_to_zero_boundary.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_get_ele.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_get_sub_vector.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_log_bound_barrier.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_abs_ele.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_inequ_viol.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_near_feas_step.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_rel_step.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_max_step.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_norms.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_num_bounded.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_num_nonzeros.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_sum.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_ROp_sum_abs.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_add_scalar.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_assign_scalar.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_assign_vectors.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_axpy.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_correct_multipliers.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_divide.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_prod.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_ele_wise_sqrt.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_force_in_bounds.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_inv_of_difference.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_max_abs_vec_scalar.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_max_vec_scalar.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_multiplier_step.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_random_vector.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_scale_vector.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_set_ele.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_set_sub_vector.c\
  $(srcdir)/RTOpPack/src/std_ops/RTOp_TOp_sign.c

RTOP_PACK_STD_OPS_HPP =

# $(srcdir)/RTOpPack/src/std_ops/RTOpPack_initRTOpStdOpsLibCpp.hpp

RTOP_PACK_STD_OPS_CPP =

# $(srcdir)/RTOpPack/src/std_ops/RTOpPack_initRTOpStdOpsLibCpp.cpp

DENSE_LIN_ALG_PACK_HPP = \
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgLAPack.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgLAPack_Types.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgLAPack_Types.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_AssertOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_AssertOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_BLAS_Cpp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_BLAS_Cpp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixAsTriSym.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixAsTriSym.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixAssign.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixAssign.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixClass.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixClass.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixInFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixInFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOpBLAS.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOpBLAS.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOut.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOut.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOutFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOutFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorAssign.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorAssign.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClass.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClass.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassExt.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassExt.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassTmpl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassTmpl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorInFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOut.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOutFunc.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_ElementaryMatVec.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IO.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IOBasic.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IOFormat.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IO_NameLookups.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IVector.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InFormat.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InFormatDecl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InFormatDef.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InvCholUpdate.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LAPACK_C_Decl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LAPACK_Cpp.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LinAlgOpPack.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LinAlgOpPackDecl.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LinAlgOpPackDef.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_MatVecCompare.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_MatlabPack.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_Options.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_OutFormat.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermIn.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermOut.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermVecMat.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PublicTypes.ud\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_Types.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_assert_print_nan_inf.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_delete_row_col.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_extended_value_type.hpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_random_vector.hpp

DENSE_LIN_ALG_PACK_CPP = \
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_BLAS_Cpp.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixClass.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixInFunc.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOpBLAS.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DMatrixOutFunc.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_InvCholUpdate.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_AssertOp.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_IOFormat.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_MatVecCompare.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_MatlabPack.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermIn.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermOut.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_PermVecMat.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorClassTmpl.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorInFunc.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOpBLAS.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_DVectorOutFunc.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_random_vector.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_delete_row_col.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_assert_print_nan_inf.cpp \
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgPack_LAPACK_C_Decl.cpp\
  $(srcdir)/DenseLinAlgPack/src/DenseLinAlgLAPack.cpp

if ENABLE_MOOCHO_MA28
ABSTRACT_LIN_ALG_PACK_MA28_HPP = \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverMA28.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverMA28SetOptions.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28CommonBlockEncap.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28Solver.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28_CppDecl.hpp
ABSTRACT_LIN_ALG_PACK_MA28_CPP = \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverMA28.cpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverMA28SetOptions.cpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28_CppDecl.cpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28CommonBlockEncap.cpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MA28Solver.cpp
else
ABSTRACT_LIN_ALG_PACK_MA28_HPP =
ABSTRACT_LIN_ALG_PACK_MA28_CPP =
endif

ABSTRACT_LIN_ALG_PACK_HPP = \
  $(srcdir)/AbstractLinAlgPack/src/AbstractLinAlgPack_PublicTypes.ud\
  $(srcdir)/AbstractLinAlgPack/src/AbstractLinAlgPack_Types.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_AssertOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_BasisSystem.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_BasisSystemFactory.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_BasisSystemPerm.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_EtaVector.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSlice.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSliceIterator.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSliceOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_InnerProduct.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_InnerProductDot.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_LinAlgOpPack.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixBase.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixNonsing.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpNonsing.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpNonsingAggr.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpSubView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixPermAggr.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymDiag.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymInitDiag.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymNonsing.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymOpNonsing.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymSecant.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MultiVector.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MultiVectorMutable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_Permutation.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_PermutationOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVecIndexLookupClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVecIndexLookupClassDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVecIndexLookupClassDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseElement.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseVectorClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseVectorClassDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseVectorClassDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_Vector.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorApplyOpSerialBase.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorMutable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorMutableSubView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpace.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpaceFactory.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpaceSubSpace.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorStdOps.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSubView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_apply_op_helper.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_compare_element_indexes.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemComposite.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemTester.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemTesterSetOptions.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixComposite.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixOpNonsingTester.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixSymDiagStd.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixSymIdent.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixZero.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_TestMatrixSymSecant.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorAuxiliaryOps.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorMutableBlocked.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceBlocked.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceTester.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceTesterSetOptions.hpp\
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_assert_print_nan_inf.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_BFGS_helpers.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMPartitionOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMPartitionOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixIn.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixOutFunc.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixPartitionViewSubclass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixPartitionedViewClass.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixPartitionedViewClassDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixPartitionedViewClassDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplConvertToSparseCompressedColumn.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplConvertToSparseCompressedColumnDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplConvertToSparseCompressedColumnDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOpDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOpDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOutFunc.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOutFuncDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_COOMatrixTmplOutFuncDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixCOORTmplItfc.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSparseCOORSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSymDiagSparseStd.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSymPosDefCholFactor.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixVectorTemplateOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixVectorTemplateOpDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixWithOpConcreteEncap.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_PermutationSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SortByDescendingAbsValue.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SpVectorOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseCOOElementTransView.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseCOOPtrElement.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseCOOReadMatrix.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseElementOut.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparsePtrElement.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOpDecl.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOpDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOutput.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorOutputDef.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SparseVectorSliceOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_TransSparseCOOElementViewIter.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_VectorDenseEncap.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_dense_Vp_StPtMtV.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_rank_2_chol_update.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_sparse_bounds.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_sparse_bounds_diff.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_GenPermMatrixSliceOp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_LinAlgOpPackHack.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixConvertToSparse.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixConvertToSparseEncap.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixExtractInvCholFactor.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixExtractSparseElements.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixLoadSparseElements.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixNonsingSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpGetGMS.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpGetGMSMutable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpGetGMSTri.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpNonsingSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymAddDelUpdateable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymDenseInitialize.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymDiagSparse.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymNonsingSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpGetGMSSym.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpGetGMSSymMutable.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpNonsingSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MultiVectorMutableCols.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MultiVectorMutableDense.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorMutableDense.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorSpaceFactorySerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorSpaceSerial.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_BasisSystemFactoryStd.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_BasisSystemPermDirectSparse.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolver.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverDense.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverImp.hpp\
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_MatrixScaling_Strategy.hpp\
  $(ABSTRACT_LIN_ALG_PACK_MA28_HPP)

#  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_SuperLUSolver.hpp
#  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverSuperLU.hpp

ABSTRACT_LIN_ALG_PACK_CPP = \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_AssertOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_AbstractLinAlgOpPack.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_apply_op_helper.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_BasisSystem.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSlice.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSliceIterator.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_GenPermMatrixSliceOut.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_InnerProductDot.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixBase.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixNonsing.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpNonsing.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpNonsingAggr.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixOpSubView.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymNonsing.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixSymOpNonsing.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MatrixPermAggr.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MultiVector.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_MultiVectorMutable.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_Vector.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorApplyOpSerialBase.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpace.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpaceSubSpace.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSpaceFactory.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorSubView.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorMutable.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorMutableSubView.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_VectorStdOps.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SparseVectorClass.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorOut.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/interfaces/AbstractLinAlgPack_SpVectorView.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_assert_print_nan_inf.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemComposite.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemTester.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_BasisSystemTesterSetOptions.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixOpNonsingTester.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixSymIdent.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixSymDiagStd.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixComposite.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_MatrixZero.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceTester.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceTesterSetOptions.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorSpaceBlocked.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorMutableBlocked.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_VectorAuxiliaryOps.cpp \
  $(srcdir)/AbstractLinAlgPack/src/abstract/tools/AbstractLinAlgPack_TestMatrixSymSecant.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_LinAlgOpPackHack.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_GenPermMatrixSliceOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorSpaceSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorSpaceFactorySerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_VectorMutableDense.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MultiVectorMutableDense.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MultiVectorMutableCols.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixConvertToSparseEncap.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixOpSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixNonsingSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymOpSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymNonsingSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixSymDiagSparse.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/interfaces/AbstractLinAlgPack_MatrixExtractSparseElements.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_BFGS_helpers.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSymDiagSparseStd.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSymPosDefCholFactor.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_MatrixSparseCOORSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_PermutationSerial.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_SpVectorOp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/AbstractLinAlgPack_rank_2_chol_update.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_BasisSystemPermDirectSparse.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_BasisSystemFactoryStd.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverImp.cpp \
  $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverDense.cpp \
  $(ABSTRACT_LIN_ALG_PACK_MA28_CPP)

# $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_DirectSparseSolverSuperLU.cpp
# $(srcdir)/AbstractLinAlgPack/src/serial/solvers/AbstractLinAlgPack_SuperLUSolver.cpp

ABSTRACT_LIN_ALG_PACK_F = \
  $(srcdir)/AbstractLinAlgPack/src/serial/implementations/dchud.f

# ToDo: Add configure-time logic to only build this if it is missing from $LIBS

NLP_INTERFACE_PACK_HPP = \
  $(srcdir)/NLPInterfacePack/src/NLPInterfacePack_PublicTypes.ud\
  $(srcdir)/NLPInterfacePack/src/NLPInterfacePack_Types.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLP.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPDirect.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPFirstOrder.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPObjGrad.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPSecondOrder.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPVarReductPerm.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPDirectTester.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPDirectTesterSetOptions.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPFirstDerivTester.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPFirstDerivTesterSetOptions.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPTester.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPTesterSetOptions.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_test_basis_system.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_test_nlp_direct.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_test_nlp_first_order.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/tools/NLPInterfacePack_CalcFiniteDiffProd.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/tools/NLPInterfacePack_CalcFiniteDiffProdSetOptions.hpp\
  $(srcdir)/NLPInterfacePack/src/abstract/tools/NLPInterfacePack_NLPBarrier.hpp\
  $(srcdir)/NLPInterfacePack/src/serial/NLPInterfacePack_NLPSerialPreprocess.hpp\
  $(srcdir)/NLPInterfacePack/src/serial/NLPInterfacePack_NLPSerialPreprocessExplJac.hpp

NLP_INTERFACE_PACK_CPP = \
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLP.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPDirect.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPFirstOrder.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPObjGrad.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/interfaces/NLPInterfacePack_NLPSecondOrder.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPDirectTester.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPDirectTesterSetOptions.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPFirstDerivTester.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPFirstDerivTesterSetOptions.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPTester.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_NLPTesterSetOptions.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_test_basis_system.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_test_nlp_direct.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/test/NLPInterfacePack_test_nlp_first_order.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/tools/NLPInterfacePack_CalcFiniteDiffProd.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/tools/NLPInterfacePack_CalcFiniteDiffProdSetOptions.cpp\
  $(srcdir)/NLPInterfacePack/src/abstract/tools/NLPInterfacePack_NLPBarrier.cpp\
  $(srcdir)/NLPInterfacePack/src/serial/NLPInterfacePack_NLPSerialPreprocess.cpp\
  $(srcdir)/NLPInterfacePack/src/serial/NLPInterfacePack_NLPSerialPreprocessExplJac.cpp

CONSTRAINED_OPT_PACK_HPP = \
  $(srcdir)/ConstrainedOptPack/src/ConstrainedOptPack_PublicTypes.ud\
  $(srcdir)/ConstrainedOptPack/src/ConstrainedOptPack_Types.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystem.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemCoordinate.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemOrthogonal.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemTester.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemTesterSetOptions.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemVarReduct.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemVarReductImp.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemVarReductPerm.hpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemVarReductPermStd.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_DirectLineSearchArmQuad_Strategy.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_DirectLineSearchArmQuad_StrategySetOptions.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_DirectLineSearch_Strategy.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncCalc.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncCalc1D.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncCalc1DQuadratic.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncCalcNLE.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncCalcNLP.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLE.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLESqrResid.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLP.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLPL1.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLPModL1.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncPenaltyParam.hpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncPenaltyParams.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixDecompRangeOrthog.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixGenBanded.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixHessianRelaxed.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixHessianSuperBasic.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixHessianSuperBasicInitDiagonal.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixIdentConcat.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixIdentConcatStd.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixKKTFullSpaceRelaxed.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymAddDelBunchKaufman.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymAddDelUpdateableWithOpNonsingular.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymHessianRelaxNonSing.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymIdentitySerial.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymPosDefBandedChol.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymPosDefInvCholFactor.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymPosDefLBFGS.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymPosDefLBFGSSerial.hpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixVarReductImplicit.hpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_ComputeMinMult.hpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_VariableBoundsTester.hpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_VariableBoundsTesterSetOptions.hpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_initialize_Q_R_Q_X.hpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_print_vector_change_stats.hpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_vector_change_stats.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_ConstraintsRelaxedStd.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPInitFixedFreeStd.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPKWIK_Output.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPOPT_CppDecl.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSchur.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSchurInitKKTSystemHessianFixedFree.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSchurInitKKTSystemHessianFull.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSchurInitKKTSystemHessianRelaxed.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSchurInitKKTSystemHessianSuperBasic.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxed.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedLOQO.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedQPKWIK.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedQPOPT.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedQPOPTSOL.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedQPSchurSetOptions.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedTester.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedTesterSetOptions.hpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverStats.hpp

CONSTRAINED_OPT_PACK_CPP = \
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystem.cpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemCoordinate.cpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemOrthogonal.cpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemTester.cpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemTesterSetOptions.cpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemVarReductImp.cpp\
  $(srcdir)/ConstrainedOptPack/src/decompositions/ConstrainedOptPack_DecompositionSystemVarReductPermStd.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_DirectLineSearchArmQuad_Strategy.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_DirectLineSearchArmQuad_StrategySetOptions.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncCalc1DQuadratic.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncCalcNLE.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncCalcNLP.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLESqrResid.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLP.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLPL1.cpp\
  $(srcdir)/ConstrainedOptPack/src/globalization/ConstrainedOptPack_MeritFuncNLPModL1.cpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixDecompRangeOrthog.cpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixIdentConcat.cpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixIdentConcatStd.cpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymAddDelBunchKaufman.cpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymHessianRelaxNonSing.cpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymIdentitySerial.cpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixSymPosDefLBFGS.cpp\
  $(srcdir)/ConstrainedOptPack/src/matrices/ConstrainedOptPack_MatrixVarReductImplicit.cpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_ComputeMinMult.cpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_VariableBoundsTester.cpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_VariableBoundsTesterSetOptions.cpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_initialize_Q_R_Q_X.cpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_print_vector_change_stats.cpp\
  $(srcdir)/ConstrainedOptPack/src/misc/ConstrainedOptPack_vector_change_stats.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_ConstraintsRelaxedStd.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPInitFixedFreeStd.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPKWIK_Output.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSchur.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSchurInitKKTSystemHessianFull.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxed.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedQPKWIK.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedQPSchur.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedQPSchurSetOptions.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedTester.cpp\
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/ConstrainedOptPack_QPSolverRelaxedTesterSetOptions.cpp

# ToDo: Add conditional logic for QPOPT and QPKWIK!

CONSTRAINED_OPT_PACK_F = \
  $(srcdir)/ConstrainedOptPack/src/qpsolvers/QPKWIKNEW.f

#CONSTRAINED_OPT_PACK_F = \
#  $(srcdir)/ConstrainedOptPack/src/qpsolvers/qpoptcpp.f

MOOCHO_PACK_HPP = \
  $(srcdir)/MoochoPack/configurations/MoochoPack_MoochoSolver.hpp \
  $(srcdir)/MoochoPack/configurations/IpConfig/MoochoPack_NLPAlgoConfigIP.hpp \
  $(srcdir)/MoochoPack/configurations/MamaJama/MoochoPack_NLPAlgoConfigMamaJama.hpp \
  $(srcdir)/MoochoPack/configurations/Shared/MoochoPack_DecompositionSystemStateStepBuilderStd.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_IpState.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_Exceptions.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_Types.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgo.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgoClientInterface.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgoConfig.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgoContainer.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgoInterface.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgoState.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPSolverClientInterface.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPSolverClientInterfaceSetOptions.hpp \
  $(srcdir)/MoochoPack/src/MoochoPack_moocho_algo_conversion.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ActSetStats.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ActSetStats_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_BFGSUpdate_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_BFGSUpdate_StrategySetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcDFromYPYZPZ_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcDFromYPY_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcDFromZPZ_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcD_vStep_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcLambdaIndepStd_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcReducedGradLagrangianStd_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergenceIP_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergenceStd_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergenceStd_AddedStepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergenceStd_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergence_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckDecompositionFromPy_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckDecompositionFromRPy_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckDescentQuasiNormalStep_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckSkipBFGSUpdateStd_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckSkipBFGSUpdateStd_StepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_DecompositionSystemHandlerSelectNew_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_DecompositionSystemHandlerStd_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_DecompositionSystemHandlerVarReductPerm_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_DecompositionSystemHandler_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointStd_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointStd_StepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointTailoredApproachCoordinate_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointTailoredApproachOrthogonal_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointTailoredApproach_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointTailoredApproach_StepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_FeasibilityStepReducedStd_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_FeasibilityStepReducedStd_StrategySetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_FeasibilityStep_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_InitFinDiffReducedHessian_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_InitFinDiffReducedHessian_StepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchDirect_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchFailureNewDecompositionSelection_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchFilter_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchFilter_StepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchFullStep_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchNLE_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_DummyUpdate_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_PenaltyParamUpdateGuts_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_PenaltyParamUpdateMultFree_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_PenaltyParamUpdateWithMult_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoAlgorithmStepNames.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoTrackerConsoleStd.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoTrackerStatsStd.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoTrackerSummaryStd.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoTrackerXMLSummary.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_NewDecompositionSelectionStd_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_NewDecompositionSelection_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_PostEvalNewPointBarrier_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_PostProcessBarrierLineSearch_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_PreEvalNewPointBarrier_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_PreProcessBarrierLineSearch_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_QPFailureReinitReducedHessian_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_QuasiNewtonStats.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_QuasiNormalStepStd_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_QuasiRangeSpaceStepStd_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_QuasiRangeSpaceStep_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedGradientStd_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSecantUpdateBFGSFull_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSecantUpdateStd_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSecantUpdate_Strategy.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSerialization_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSerialization_StepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_SetDBoundsStd_AddedStep.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_TangentialStepIP_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_TangentialStepWithInequStd_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_TangentialStepWithInequStd_StepSetOptions.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_TangentialStepWithoutBounds_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_UpdateBarrierParameter_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_UpdateReducedSigma_Step.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_act_set_stats.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_d_bounds_iter_quant.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_qp_solver_stats.hpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_quasi_newton_stats.hpp

MOOCHO_PACK_CPP = \
  $(srcdir)/MoochoPack/configurations/MoochoPack_MoochoSolver.cpp \
  $(srcdir)/MoochoPack/configurations/IpConfig/MoochoPack_NLPAlgoConfigIP.cpp \
  $(srcdir)/MoochoPack/configurations/MamaJama/MoochoPack_NLPAlgoConfigMamaJama.cpp \
  $(srcdir)/MoochoPack/configurations/Shared/MoochoPack_DecompositionSystemStateStepBuilderStd.cpp \
  $(srcdir)/MoochoPack/src/MoochoPack_IpState.cpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgo.cpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgoContainer.cpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPAlgoState.cpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPSolverClientInterface.cpp \
  $(srcdir)/MoochoPack/src/MoochoPack_NLPSolverClientInterfaceSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_BFGSUpdate_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_BFGSUpdate_StrategySetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcDFromYPYZPZ_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcDFromYPY_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcDFromZPZ_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcD_vStep_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CalcReducedGradLagrangianStd_AddedStep.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergenceIP_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergenceStd_AddedStep.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergenceStd_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckConvergence_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckDecompositionFromPy_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckDecompositionFromRPy_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckDescentQuasiNormalStep_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckSkipBFGSUpdateStd_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_CheckSkipBFGSUpdateStd_StepSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_DecompositionSystemHandlerStd_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_DecompositionSystemHandlerVarReductPerm_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointStd_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointStd_StepSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointTailoredApproachCoordinate_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointTailoredApproachOrthogonal_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointTailoredApproach_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_EvalNewPointTailoredApproach_StepSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_FeasibilityStepReducedStd_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_InitFinDiffReducedHessian_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_InitFinDiffReducedHessian_StepSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchDirect_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchFailureNewDecompositionSelection_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchFilter_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchFilter_StepSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchFullStep_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_LineSearchNLE_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_DummyUpdate_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_PenaltyParamUpdateGuts_AddedStep.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_PenaltyParamUpdateMultFree_AddedStep.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoTrackerConsoleStd.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoTrackerStatsStd.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoTrackerSummaryStd.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_MoochoTrackerXMLSummary.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_NewDecompositionSelectionStd_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_PostEvalNewPointBarrier_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_PostProcessBarrierLineSearch_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_PreEvalNewPointBarrier_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_PreProcessBarrierLineSearch_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_QPFailureReinitReducedHessian_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_QuasiNormalStepStd_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_QuasiRangeSpaceStepStd_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedGradientStd_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSecantUpdateBFGSFull_Strategy.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSecantUpdateStd_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSerialization_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_ReducedHessianSerialization_StepSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_SetDBoundsStd_AddedStep.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_TangentialStepIP_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_TangentialStepWithInequStd_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_TangentialStepWithInequStd_StepSetOptions.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_TangentialStepWithoutBounds_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_UpdateBarrierParameter_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_UpdateReducedSigma_Step.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_act_set_stats.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_d_bounds_iter_quant.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_qp_solver_stats.cpp \
  $(srcdir)/MoochoPack/src/std/MoochoPack_quasi_newton_stats.cpp

#
# Specify locations of include files within Trilinos
#

include $(top_builddir)/Makefile.export.moocho

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(MOOCHO_INCLUDES))
else
EXPORT_INCLUDES = $(MOOCHO_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# Moocho library specifications
# --------------------------------------------------------------------

MOOCHO_LIB = libmoocho.a

MOOCHO_H = \
  $(srcdir)/Moocho_ConfigDefs.hpp \
  $(MOOCHO_UTILITIES_HPP) \
  $(ITERATION_PACK_HPP) \
  $(RTOP_PACK_BASE_H) \
  $(RTOP_PACK_BASE_HPP) \
  $(RTOP_PACK_STD_OPS_H) \
  $(RTOP_PACK_STD_OPS_HPP) \
  $(DENSE_LIN_ALG_PACK_HPP) \
  $(ABSTRACT_LIN_ALG_PACK_HPP) \
  $(NLP_INTERFACE_PACK_HPP) \
  $(CONSTRAINED_OPT_PACK_HPP) \
  $(MOOCHO_PACK_HPP)

libmoocho_a_SOURCES = \
  $(MOOCHO_UTILITIES_CPP) \
  $(MOOCHO_UTILITIES_F) \
  $(ITERATION_PACK_CPP) \
  $(RTOP_PACK_BASE_C) \
  $(RTOP_PACK_BASE_CPP) \
  $(RTOP_PACK_STD_OPS_C) \
  $(RTOP_PACK_STD_OPS_CPP) \
  $(DENSE_LIN_ALG_PACK_CPP) \
  $(ABSTRACT_LIN_ALG_PACK_CPP) \
  $(ABSTRACT_LIN_ALG_PACK_F) \
  $(NLP_INTERFACE_PACK_CPP) \
  $(CONSTRAINED_OPT_PACK_CPP) \
  $(CONSTRAINED_OPT_PACK_F) \
  $(MOOCHO_PACK_CPP)

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libmoocho_a_AR = $(ALTERNATE_AR)

else

libmoocho_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

libmoocho_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(MOOCHO_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(MOOCHO_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: MoochoConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(includedir)
nodist_execinclude_HEADERS = Moocho_Config.h

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

EXTRA_DIST = $(srcdir)/MoochoPack/configurations/generate-opt-file.pl \
	$(srcdir)/MoochoPack/configurations/IpConfig/Moocho.opt.NLPAlgoConfigIP \
	$(srcdir)/MoochoPack/configurations/Moocho.opt.MoochoSolver \
	$(srcdir)/MoochoPack/configurations/Shared/Moocho.opt.DecompositionSystemStateStepBuilderStd \
	$(srcdir)/MoochoPack/configurations/MamaJama/Moocho.opt.NLPAlgoConfigMamaJama \
	$(srcdir)/MoochoPack/configurations/sample_option_files/Moocho.sample.opt
