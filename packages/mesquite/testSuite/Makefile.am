if ENABLE_IMESH
  IMESH_TEST = imesh
else
  IMESH_TEST =
endif

if ENABLE_IGEOM
  IGEOM_TEST = igeom
else
  IGEOM_TEST =
endif

if WITH_EXODUS
  EXO_TEST = convert 
else
  EXO_TEST = 
endif

if WITH_CPPUNIT
  UNIT_TEST = unit
else
  UNIT_TEST =
endif

if WITH_MPI
  PARALLEL_UNTANGLE_SHAPE = parallel_untangle_shape
  PARALLEL_SMOOTH_LAPLACE = parallel_smooth_laplace
else
  PARALLEL_UNTANGLE_SHAPE = 
  PARALLEL_SMOOTH_LAPLACE = 
endif

SUBDIRS = test_1 \
	  ActiveSetTest \
	  algorithm_test\
	  analytical_grad_3D \
          benchmark_tests \
          paraboloid_domain_test \
          $(EXO_TEST) \
	  headers \
          high_aspect_ratio \
	  higher_order \
	  idft_time \
          $(IGEOM_TEST) \
	  $(IMESH_TEST) \
	  jacobi \
	  laplacian_test\
	  laplacian_polygon_test\
	  $(PARALLEL_SMOOTH_LAPLACE) \
	  $(PARALLEL_UNTANGLE_SHAPE) \
	  nongradient_test\
	  pyramid \
	  simple_hybrid_test\
          slaved \
	  synchronous \
	  transform \
	  tutorial\
	  untangle_test\
	  wrapper_tests\
	  wedge \
	  2d_target \
	  2d_metrics\
	  2d_formulation \
          $(UNIT_TEST)

EXTRA_DIST = CMakeLists.txt

# Override all-recursive so that "make" with no target 
# skips subdirectories in testSuite.
all-recursive:
	@:
