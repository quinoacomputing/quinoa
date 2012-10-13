/*******************************************************************************
!  Copyright(C) 2003-2012 Intel Corporation. All Rights Reserved.
!
!  The source code, information  and  material ("Material") contained herein is
!  owned  by Intel Corporation or its suppliers or licensors, and title to such
!  Material remains  with Intel Corporation  or its suppliers or licensors. The
!  Material  contains proprietary information  of  Intel or  its  suppliers and
!  licensors. The  Material is protected by worldwide copyright laws and treaty
!  provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
!  modified, published, uploaded, posted, transmitted, distributed or disclosed
!  in any way  without Intel's  prior  express written  permission. No  license
!  under  any patent, copyright  or  other intellectual property rights  in the
!  Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
!  implication, inducement,  estoppel or  otherwise.  Any  license  under  such
!  intellectual  property  rights must  be express  and  approved  by  Intel in
!  writing.
!
!  *Third Party trademarks are the property of their respective owners.
!
!  Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
!  this  notice or  any other notice embedded  in Materials by Intel or Intel's
!  suppliers or licensors in any way.
!
!******************************************************************************/

#include <stdlib.h>

void CheckVslError(int num)
{
    switch(num) {
        case VSL_ERROR_CPU_NOT_SUPPORTED: {
            printf("Error: CPU version is not supported. (code %d).\n",num);
            break;
        }
        case VSL_ERROR_FEATURE_NOT_IMPLEMENTED: {
            printf("Error: this feature not implemented yet. (code %d).\n",num);
            break;
        }
        case VSL_ERROR_UNKNOWN: {
            printf("Error: unknown error (code %d).\n",num);
            break;
        }
        case VSL_ERROR_BADARGS: {
            printf("Error: bad arguments (code %d).\n",num);
            break;
        }
        case VSL_ERROR_MEM_FAILURE: {
            printf("Error: memory failure. Memory allocation problem maybe (code %d).\n",num);
            break;
        }
        case VSL_ERROR_NULL_PTR: {
            printf("Error: null pointer (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_INVALID_BRNG_INDEX: {
            printf("Error: invalid BRNG index (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_LEAPFROG_UNSUPPORTED: {
            printf("Error: leapfrog initialization is unsupported (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_SKIPAHEAD_UNSUPPORTED: {
            printf("Error: skipahead initialization is unsupported (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BRNGS_INCOMPATIBLE: {
            printf("Error: BRNGs are not compatible for the operation (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BAD_STREAM: {
            printf("Error: random stream is invalid (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BRNG_TABLE_FULL: {
            printf("Error: table of registered BRNGs is full (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BAD_STREAM_STATE_SIZE: {
            printf("Error: value in StreamStateSize field is bad (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BAD_WORD_SIZE: {
            printf("Error: value in WordSize field is bad (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BAD_NSEEDS: {
            printf("Error: value in NSeeds field is bad (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BAD_NBITS: {
            printf("Error: value in NBits field is bad (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BAD_UPDATE: {
            printf("Error: number of updated entries in buffer is invalid (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_NO_NUMBERS: {
            printf("Error: zero number of updated entries in buffer (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_INVALID_ABSTRACT_STREAM: {
            printf("Error: abstract random stream is invalid (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_NONDETERM_NOT_SUPPORTED: {
            printf("Error: Non-deterministic random number generator is not supported on the CPU which runs the application (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_NONDETERM_NRETRIES_EXCEEDED: {
            printf("Error: Number of retries to generate a random number by using non-deterministic random number generator exceeds threshold (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_FILE_CLOSE: {
            printf("Error: can`t close file (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_FILE_OPEN: {
            printf("Error: can`t open file (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_FILE_WRITE: {
            printf("Error: can`t write to file (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_FILE_READ: {
            printf("Error: can`t read from file (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BAD_FILE_FORMAT: {
            printf("Error: file format is unknown (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_UNSUPPORTED_FILE_VER: {
            printf("Error: unsupported file version (code %d).\n",num);
            break;
        }
        case VSL_RNG_ERROR_BAD_MEM_FORMAT: {
            printf("Error: random stream format is unknown (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_ALLOCATION_FAILURE: {
            printf("Error: memory allocation failure in summary statistics functionality (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_DIMEN: {
            printf("Error: bad dimension value (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_OBSERV_N: {
            printf("Error: bad number of observations (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_STORAGE_NOT_SUPPORTED: {
            printf("Error: storage format is not supported (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_INDC_ADDR: {
            printf("Error: array of indices is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_WEIGHTS: {
            printf("Error: array of weights contains negative values (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MEAN_ADDR: {
            printf("Error: array of means is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_2R_MOM_ADDR: {
            printf("Error: array of 2nd order raw moments is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_3R_MOM_ADDR: {
            printf("Error: array of 3rd order raw moments is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_4R_MOM_ADDR: {
            printf("Error: array of 4th order raw moments is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_2C_MOM_ADDR: {
            printf("Error: array of 2nd order central moments is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_3C_MOM_ADDR: {
            printf("Error: array of 3rd order central moments is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_4C_MOM_ADDR: {
            printf("Error: array of 4th order central moments is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_KURTOSIS_ADDR: {
            printf("Error: array of kurtosis values is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_SKEWNESS_ADDR: {
            printf("Error: array of skewness values is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MIN_ADDR: {
            printf("Error: array of minimum values is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MAX_ADDR: {
            printf("Error: array of maximum values is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_VARIATION_ADDR: {
            printf("Error: array of variation coefficients is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_COV_ADDR: {
            printf("Error: covariance matrix is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_COR_ADDR: {
            printf("Error: correlation matrix is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_QUANT_ORDER_ADDR: {
            printf("Error: array of quantile orders is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_QUANT_ORDER: {
            printf("Error: bad value of quantile order (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_QUANT_ADDR: {
            printf("Error: array of quantiles is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_ORDER_STATS_ADDR: {
            printf("Error: array of order statistics is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_MOMORDER_NOT_SUPPORTED: {
            printf("Error: moment of requested order is not supported (code %d).\n",num);
            break;
        }
        case VSL_SS_NOT_FULL_RANK_MATRIX: {
            printf("Warning: correlation matrix is not of full rank (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_ALL_OBSERVS_OUTLIERS: {
            printf("Error: all observations are outliers (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_ROBUST_COV_ADDR: {
            printf("Error: robust covariance matrix is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_ROBUST_MEAN_ADDR: {
            printf("Error: array of robust means is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_METHOD_NOT_SUPPORTED: {
            printf("Error: requested method is not supported (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_NULL_TASK_DESCRIPTOR: {
            printf("Error: task descriptor is null (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_OBSERV_ADDR: {
            printf("Error: dataset matrix is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_SINGULAR_COV: {
            printf("Error: covariance matrix is singular (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_POOLED_COV_ADDR: {
            printf("Error: pooled covariance matrix is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_POOLED_MEAN_ADDR: {
            printf("Error: array of pooled means is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_GROUP_COV_ADDR: {
            printf("Error: group covariance matrix is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_GROUP_MEAN_ADDR: {
            printf("Error: array of group means is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_GROUP_INDC_ADDR: {
            printf("Error: array of group indices is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_GROUP_INDC: {
            printf("Error: group indices have improper values (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_ADDR: {
            printf("Error: array of parameters for outliers detection algorithm is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_OUTLIERS_PARAMS_N_ADDR: {
            printf("Error: pointer to size of parameter array for outlier detection algorithm is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_OUTLIERS_WEIGHTS_ADDR: {
            printf("Error: output of the outlier detection algorithm is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_ADDR: {
            printf("Error: array of parameters of robust covariance estimation algorithm is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_ROBUST_COV_PARAMS_N_ADDR: {
            printf("Error: pointer to number of parameters of algorithm for robust covariance is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_STORAGE_ADDR: {
            printf("Error: pointer to variable that holds storage format is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_PARTIAL_COV_IDX_ADDR: {
            printf("Error: array that encodes sub-components of random vector for partial covariance algorithm is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_PARTIAL_COV_ADDR: {
            printf("Error: partial covariance matrix is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_PARTIAL_COR_ADDR: {
            printf("Error: partial correlation matrix is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_PARAMS_ADDR: {
            printf("Error: array of parameters for Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_PARAMS_N_ADDR: {
            printf("Error: pointer to number of parameters for Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_BAD_PARAMS_N: {
            printf("Error: bad size of the parameter array of Multiple Imputation method (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_PARAMS: {
            printf("Error: bad parameters of Multiple Imputation method (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_N_ADDR: {
            printf("Error: pointer to number of initial estimates in Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_INIT_ESTIMATES_ADDR: {
            printf("Error: array of initial estimates for Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_SIMUL_VALS_ADDR: {
            printf("Error: array of simulated missing values in Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N_ADDR: {
            printf("Error: pointer to size of the array of simulated missing values in Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_ESTIMATES_N_ADDR: {
            printf("Error: pointer to the number of parameter estimates in Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_ESTIMATES_ADDR: {
            printf("Error: array of parameter estimates in Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_SIMUL_VALS_N: {
            printf("Error: bad size of the array of simulated values in Multiple Imputation method (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_ESTIMATES_N: {
            printf("Error: bad size of array to hold parameter estimates obtained using Multiple Imputation method (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_OUTPUT_PARAMS: {
            printf("Error: array of output parameters in Multiple Imputation method is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_PRIOR_N_ADDR: {
            printf("Error: pointer to the number of prior parameters is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_PRIOR_ADDR: {
            printf("Error: array of prior parameters is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_MI_MISSING_VALS_N: {
            printf("Error: bad number of missing values (code %d).\n",num);
            break;
        }
        case VSL_SS_SEMIDEFINITE_COR: {
            printf("Warning: correlation matrix passed into parametrization function is semidefinite(code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_PARAMTR_COR_ADDR: {
            printf("Error: correlation matrix to be parametrized is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_COR: {
            printf("Error: all eigenvalues of correlation matrix to be parametrized are non-positive (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N_ADDR: {
            printf("Error: pointer to the number of parameters for quantile computation algorithm for streaming data is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_ADDR: {
            printf("Error: array of parameters of quantile computation algorithm for streaming data is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS_N: {
            printf("Error: bad number of parameters of quantile computation algorithm for streaming data (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_STREAM_QUANT_PARAMS: {
            printf("Error: bad parameters of quantile computation algorithm for streaming data (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER_ADDR: {
            printf("Error: array of quantile orders for streaming data is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_STREAM_QUANT_ORDER: {
            printf("Error: bad quantile order for streaming data (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_STREAM_QUANT_ADDR: {
            printf("Error: array of quantiles for streaming data is not defined (code %d).\n",num);
            break;
        }
        case VSL_SS_ERROR_BAD_PARTIAL_COV_IDX: {
            printf("Error: partial covariance indices have improper values (code %d).\n",num);
            break;
        }
    }

    if(num < 0) {
       exit(1);
    }
}
