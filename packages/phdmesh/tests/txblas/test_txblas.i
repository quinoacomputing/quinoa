###
#
# accuracy <#items> <magnitude>
# accuracy 1000000 1e+10
#
# reduce <#items> <magnitude>
# reduce 1000000 1e10
#
# timing_blas1 <Target_Mflops>
# timing_mxv <size> <# bands> <stride-bands> <# Mflops>
#
threadpool 1
timing_blas1 1000
#
timing_mxv    1000 10 1 2000
timing_mxv    2000 10 1 2000
timing_mxv    5000 10 1 2000
timing_mxv   10000 10 1 2000
timing_mxv   20000 10 1 2000
timing_mxv   50000 10 1 2000
timing_mxv  100000 10 1 2000
timing_mxv  200000 10 1 2000
timing_mxv  500000 10 1 2000
timing_mxv 1000000 10 1 2000
###
threadpool 2
timing_blas1 1000
#
timing_mxv    1000 10 1 2000
timing_mxv    2000 10 1 2000
timing_mxv    5000 10 1 2000
timing_mxv   10000 10 1 2000
timing_mxv   20000 10 1 2000
timing_mxv   50000 10 1 2000
timing_mxv  100000 10 1 2000
timing_mxv  200000 10 1 2000
timing_mxv  500000 10 1 2000
timing_mxv 1000000 10 1 2000
###
threadpool 4
timing_blas1 1000
#
timing_mxv    1000 10 1 2000
timing_mxv    2000 10 1 2000
timing_mxv    5000 10 1 2000
timing_mxv   10000 10 1 2000
timing_mxv   20000 10 1 2000
timing_mxv   50000 10 1 2000
timing_mxv  100000 10 1 2000
timing_mxv  200000 10 1 2000
timing_mxv  500000 10 1 2000
timing_mxv 1000000 10 1 2000
###
threadpool 8
timing_blas1 1000
#
timing_mxv    1000 10 1 2000
timing_mxv    2000 10 1 2000
timing_mxv    5000 10 1 2000
timing_mxv   10000 10 1 2000
timing_mxv   20000 10 1 2000
timing_mxv   50000 10 1 2000
timing_mxv  100000 10 1 2000
timing_mxv  200000 10 1 2000
timing_mxv  500000 10 1 2000
timing_mxv 1000000 10 1 2000
###




