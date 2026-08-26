#ifndef EOSDeviceFn_h
#define EOSDeviceFn_h

#if defined(__CUDACC__)
  #define EOS_FN __host__ __device__ inline
#else
  #define EOS_FN inline
#endif

#endif
