//******************************************************************************
/*!
  \file      src/Random/VSLException.h
  \author    J. Bakosi
  \date      Sun 27 Jan 2013 12:19:33 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Intel's Vector Statistical Library exception
  \details   Intel's Vector Statistical Library exception
*/
//******************************************************************************
#ifndef VSLException_h
#define VSLException_h

#include <string>
#include <map>

using namespace std;

#include <MKLException.h>

namespace Quinoa {

//! VSL exception types
enum VSLExceptType { VSL_UNIMPLEMENTED=0,
                     VSL_CPU_NOT_SUPPORTED,
                     VSL_FEATURE_NOT_IMPLEMENTED,
                     VSL_UNKNOWN,
                     VSL_BADARGS,
                     VSL_MEM_FAILURE,
                     VSL_NULL_PTR,
                     VSL_INVALID_BRNG_INDEX,
                     VSL_LEAPFROG_UNSUPPORTED,
                     VSL_SKIPAHEAD_UNSUPPORTED,
                     VSL_BRNGS_INCOMPATIBLE,
                     VSL_BAD_STREAM,
                     VSL_BRNG_TABLE_FULL,
                     VSL_BAD_STREAM_STATE_SIZE,
                     VSL_BAD_WORD_SIZE,
                     VSL_BAD_NSEEDS,
                     VSL_BAD_NBITS,
                     VSL_BAD_UPDATE,
                     VSL_NO_NUMBERS,
                     VSL_INVALID_ABSTRACT_STREAM,
                     VSL_NONDETERM_NOT_SUPPORTED,
                     VSL_NONDETERM_NRETRIES_EXCEEDED,
                     VSL_FILE_CLOSE,
                     VSL_FILE_OPEN,
                     VSL_FILE_WRITE,
                     VSL_FILE_READ,
                     VSL_BAD_FILE_FORMAT,
                     VSL_UNSUPPORTED_FILE_VER,
                     VSL_BAD_MEM_FORMAT,
                     VSL_ALLOCATION_FAILURE,
                     VSL_BAD_DIMEN,
                     VSL_BAD_OBSERV_N,
                     VSL_STORAGE_NOT_SUPPORTED,
                     VSL_BAD_INDC_ADDR,
                     VSL_BAD_WEIGHTS,
                     VSL_BAD_MEAN_ADDR,
                     VSL_BAD_2R_MOM_ADDR,
                     VSL_BAD_3R_MOM_ADDR,
                     VSL_BAD_4R_MOM_ADDR,
                     VSL_BAD_2C_MOM_ADDR,
                     VSL_BAD_3C_MOM_ADDR,
                     VSL_BAD_4C_MOM_ADDR,
                     VSL_BAD_KURTOSIS_ADDR,
                     VSL_BAD_SKEWNESS_ADDR,
                     VSL_BAD_MIN_ADDR,
                     VSL_BAD_MAX_ADDR,
                     VSL_BAD_VARIATION_ADDR,
                     VSL_BAD_COV_ADDR,
                     VSL_BAD_COR_ADDR,
                     VSL_BAD_QUANT_ORDER_ADDR,
                     VSL_BAD_QUANT_ORDER,
                     VSL_BAD_QUANT_ADDR,
                     VSL_BAD_ORDER_STATS_ADDR,
                     VSL_MOMORDER_NOT_SUPPORTED,
                     VSL_NOT_FULL_RANK_MATRIX,
                     VSL_ALL_OBSERVS_OUTLIERS,
                     VSL_BAD_ROBUST_COV_ADDR,
                     VSL_BAD_ROBUST_MEAN_ADDR,
                     VSL_METHOD_NOT_SUPPORTED,
                     VSL_NULL_TASK_DESCRIPTOR,
                     VSL_BAD_OBSERV_ADDR,
                     VSL_SINGULAR_COV,
                     VSL_BAD_POOLED_COV_ADDR,
                     VSL_BAD_POOLED_MEAN_ADDR,
                     VSL_BAD_GROUP_COV_ADDR,
                     VSL_BAD_GROUP_MEAN_ADDR,
                     VSL_BAD_GROUP_INDC_ADDR,
                     VSL_BAD_GROUP_INDC,
                     VSL_BAD_OUTLIERS_PARAMS_ADDR,
                     VSL_BAD_OUTLIERS_PARAMS_N_ADDR,
                     VSL_BAD_OUTLIERS_WEIGHTS_ADDR,
                     VSL_BAD_ROBUST_COV_PARAMS_ADDR, 
                     VSL_BAD_ROBUST_COV_PARAMS_N_ADDR,
                     VSL_BAD_STORAGE_ADDR,
                     VSL_BAD_PARTIAL_COV_IDX_ADDR,
                     VSL_BAD_PARTIAL_COV_ADDR,
                     VSL_BAD_PARTIAL_COR_ADDR,
                     VSL_BAD_MI_PARAMS_ADDR,
                     VSL_BAD_MI_PARAMS_N_ADDR,
                     VSL_BAD_MI_BAD_PARAMS_N,
                     VSL_BAD_MI_PARAMS,
                     VSL_BAD_MI_INIT_ESTIMATES_N_ADDR,
                     VSL_BAD_MI_INIT_ESTIMATES_ADDR,
                     VSL_BAD_MI_SIMUL_VALS_ADDR,
                     VSL_BAD_MI_SIMUL_VALS_N_ADDR,
                     VSL_BAD_MI_ESTIMATES_N_ADDR,
                     VSL_BAD_MI_ESTIMATES_ADDR,
                     VSL_BAD_MI_SIMUL_VALS_N,
                     VSL_BAD_MI_ESTIMATES_N,
                     VSL_BAD_MI_OUTPUT_PARAMS,
                     VSL_BAD_MI_PRIOR_N_ADDR,
                     VSL_BAD_MI_PRIOR_ADDR,
                     VSL_BAD_MI_MISSING_VALS_N,
                     VSL_SEMIDEFINITE_COR,
                     VSL_BAD_PARAMTR_COR_ADDR,
                     VSL_BAD_COR,
                     VSL_BAD_STREAM_QUANT_PARAMS_N_ADDR,
                     VSL_BAD_STREAM_QUANT_PARAMS_ADDR,
                     VSL_BAD_STREAM_QUANT_PARAMS_N,
                     VSL_BAD_STREAM_QUANT_PARAMS,
                     VSL_BAD_STREAM_QUANT_ORDER_ADDR,
                     VSL_BAD_STREAM_QUANT_ORDER,
                     VSL_BAD_STREAM_QUANT_ADDR,
                     VSL_BAD_PARTIAL_COV_IDX,
                     NUM_VSL_EXCEPT
};

//! VSL exception error messages
const string VSLMsg[NUM_VSL_EXCEPT] = {
  "VSL exception unimplemented (update VSLException with latest VSL errors)",
  "CPU version is not supported",
  "feature not yet implemented",
  "unknown error",
  "bad arguments",
  "memory allocation problem",
  "null pointer",
  "invalid BRNG index",
  "LeapFrog initialization unsupported",
  "SkipAhead initialization unsupported",
  "BRNGs are not compatible for the operation",
  "Random stream is invalid",
  "table of registered BRNGs is full",
  "value in StreamStateSize field is bad",
  "value in WordSize field is bad",
  "value in NSeeds field is bad",
  "value in NBits field is bad",
  "number of updated entries in buffer is invalid",
  "zero number of updated entries in buffer",
  "abstract random stream is invalid",
  "non-deterministic random number generator is not supported on CPU",
  "number of retries to generate a random number by using non-deterministic "
    "random number generator exceeds threshold",
  "cannot close file",
  "cannot open file",
  "cannot write to file",
  "cannot read from file",
  "file format is unknown",
  "unsupported file version",
  "random stream format is unknown",
  "memory allocation failure in summary statistics functionality",
  "bad dimension value",
  "bad number of observations",
  "storage format is not supported",
  "array of indices is not defined",
  "array of weights contains negative values",
  "array of means is not defined",
  "array of 2nd order raw moments is not defined",
  "array of 3rd order raw moments is not defined",
  "array of 4th order raw moments is not defined",
  "array of 2nd order central moments is not defined",
  "array of 3rd order central moments is not defined",
  "array of 4th order central moments is not defined",
  "array of kurtosis values is not defined",
  "array of skewness values is not defined",
  "array of minimum values is not defined",
  "array of maximum values is not defined",
  "array of variation coefficients is not defined",
  "covariance matrix is not defined",
  "correlation matrix is not defined",
  "array of quantile orders is not defined",
  "bad value of quantile order",
  "array of quantiles is not defined",
  "array of order statistics is not defined",
  "moment of requested order is not supported",
  "correlation matrix is not of full rank",
  "all observations are outliers",
  "robust covariance matrix is not defined",
  "array of robust means is not defined",
  "requested method is not supported",
  "task descriptor is null",
  "dataset matrix is not defined",
  "covariance matrix is singular",
  "pooled covariance matrix is not defined",
  "array of pooled means is not defined",
  "group covariance matrix is not defined",
  "array of group means is not defined",
  "array of group indices is not defined",
  "group indices have improper values",
  "array of parameters for outliers detection algorithm is not defined",
  "pointer to size of parameter array for outlier detection algorithm is not "
    "defined",
  "output of the outlier detection algorithm is not defined",
  "array of parameters of robust covariance estimation algorithm is not "
    "defined",
  "pointer to number of parameters of algorithm for robust covariance is not "
    "defined",
  "pointer to variable that holds storage format is not defined",
  "array that encodes sub-components of random vector for partial covariance "
    "algorithm is not defined",
  "partial covariance matrix is not defined",
  "partial correlation matrix is not defined",
  "array of parameters for Multiple Imputation method is not defined",
  "pointer to number of parameters for Multiple Imputation method is not "
    " defined",
  "bad size of the parameter array of Multiple Imputation method",
  "bad parameters of Multiple Imputation method",
  "pointer to number of initial estimates in Multiple Imputation method is "
    "not defined",
  "array of initial estimates for Multiple Imputation method is not defined",
  "array of simulated missing values in Multiple Imputation method is not "
    "defined",
  "pointer to size of the array of simulated missing values in Multiple "
    "Imputation method is not defined",
  "pointer to the number of parameter estimates in Multiple Imputation method "
    "is not defined",
  "array of parameter estimates in Multiple Imputation method is not defined",
  "bad size of the array of simulated values in Multiple Imputation method",
  "bad size of array to hold parameter estimates obtained using Multiple "
    "Imputation method",
  "array of output parameters in Multiple Imputation method is not defined",
  "pointer to the number of prior parameters is not defined",
  "array of prior parameters is not defined",
  "bad number of missing values",
  "correlation matrix passed into parametrization function is semidefinite",
  "correlation matrix to be parametrized is not defined",
  "all eigenvalues of correlation matrix to be parametrized are non-positive",
  "pointer to the number of parameters for quantile computation algorithm for "
    "streaming data is not defined",
  "array of parameters of quantile computation algorithm for streaming data "
    "is not defined",
  "bad number of parameters of quantile computation algorithm for streaming "
    "data",
  "bad parameters of quantile computation algorithm for streaming data",
  "array of quantile orders for streaming data is not defined",
  "bad quantile order for streaming data",
  "array of quantiles for streaming data is not defined",
  "partial covariance indices have improper values"
};

//! VSLException : RandomException
class VSLException : public MKLException {

  public:
    //! Constructor
    VSLException(ExceptType except,
                 int vslerr,
                 const string& file,
                 const string& func,
                 const unsigned int& line);

    //! Move constructor, necessary for throws, default compiler generated
    VSLException(VSLException&&) = default;

    //! Destructor
    //virtual ~VSLException() {}

    //! Handle VSLException
    virtual ErrCode handleException(Driver* driver);

    //! Get VSLException based on VSLError
    VSLExceptType getException(int vslerr);

  private:
    //! Don't permit copy constructor
    VSLException(const VSLException&) = delete;
    //! Don't permit copy assignment
    VSLException& operator=(const VSLException&) = delete;
    //! Don't permit move assignment
    VSLException& operator=(VSLException&&) = delete;

    //! VSL exception type (VSL_UNIMPLEMENTED, VSL_UNKNOWN, etc.)
    VSLExceptType m_except;

    //! VSLError -> VSLExceptType map
    map<int,VSLExceptType> m_VSLErrMap;
};

} // namespace Quinoa

#endif // VSLException_h
