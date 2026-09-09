#ifndef EOSDeviceFn_h
#define EOSDeviceFn_h

// Device-capable compilers only: nvcc defines __CUDACC__; hipcc and
// amdclang++ (clang with -x hip) define __HIPCC__ in both the host and the
// device pass. Host-only targets compiled with mpicxx define neither and must
// not see __host__/__device__.
#if defined(__CUDACC__) || defined(__HIPCC__)
  #define EOS_FN __host__ __device__ inline
#else
  #define EOS_FN inline
#endif

// Defined only during the device pass of a device-capable compiler, to gate
// host-only bodies (iostream, Throw, LAPACK) out of device code. nvcc signals
// this with __CUDA_ARCH__; clang does not define __CUDA_ARCH__ when compiling
// HIP and uses __HIP_DEVICE_COMPILE__ instead.
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  #define EOS_DEVICE_PASS 1
#endif

#endif
