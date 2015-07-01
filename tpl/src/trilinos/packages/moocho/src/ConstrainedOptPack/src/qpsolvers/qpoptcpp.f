*======================================================================
*
* Helper functions to allow QPOPT to be called from C++
*
*======================================================================

*======================================================================
      SUBROUTINE QPOPT_SET_DEFAULTS

*** Reset all of options to thier defaults for QPOPT.

      IMPLICIT NONE

      EXTERNAL QPPRM
      EXTERNAL QPPRMI
      EXTERNAL QPPRMR

      CALL QPPRM( 'Defaults' )

      END

*======================================================================
      SUBROUTINE QPOPT_INT_OPT ( OPTION, VAL )

      IMPLICIT NONE

      INTEGER OPTION, VAL

*** Set an option for QPOPT that can be represented using an INTEGER.
***
***		CHECK_FREQUENCY               = 1
***		EXPAND_FREQUENCY              = 2
***		FEASIBILITY_PHASE_ITER_LIMIT  = 3
***		OPTIMALITY_PHASE_ITER_LIMIT   = 4
***		HESSIAN_ROWS                  = 5
***		ITERATION_LIMIT               = 6
***		MAXIMUM_DEGREES_OF_FREEDOM    = 7
***		PRINT_FILE                    = 8
***		PRINT_LEVEL                   = 9
***		PROBLEM_TYPE                  = 10
***			values of:
***				FP								= 1
***				LP								= 2
***				QP1								= 3
***				QP2								= 4
***				QP3								= 5
***				QP4								= 6
***		SUMMARY_FILE                  = 11
***
*** See the documentation from QPOPT to see a description of
*** these options.
***

      EXTERNAL QPPRM
      EXTERNAL QPPRMI
      EXTERNAL QPPRMR

      IF     ( OPTION .EQ. 1  ) THEN
          CALL QPPRMI( 'Check frequency', VAL )
      ELSEIF ( OPTION .EQ. 2  ) THEN
          CALL QPPRMI( 'Expand frequency', VAL )
      ELSEIF ( OPTION .EQ. 3  ) THEN
          CALL QPPRMI( 'Feasibility Phase Iteration Limit', VAL )
      ELSEIF ( OPTION .EQ. 4  ) THEN
          CALL QPPRMI( 'Optimality Phase Iteration Limit', VAL )
      ELSEIF ( OPTION .EQ. 5  ) THEN
          CALL QPPRMI( 'Hessian rows', VAL )
      ELSEIF ( OPTION .EQ. 6  ) THEN
          CALL QPPRMI( 'Iteration Limit', VAL )
      ELSEIF ( OPTION .EQ. 7  ) THEN
          CALL QPPRMI( 'Maximum degrees of freedom', VAL )
      ELSEIF ( OPTION .EQ. 8  ) THEN
          CALL QPPRMI( 'Print file', VAL )
      ELSEIF ( OPTION .EQ. 9  ) THEN
          CALL QPPRMI( 'Print level', VAL )
      ELSEIF ( OPTION .EQ. 10  ) THEN
          IF     ( VAL .EQ. 1 ) THEN
              CALL QPPRM( 'Problem type      FP' )
          ELSEIF ( VAL .EQ. 2 ) THEN
              CALL QPPRM( 'Problem type      LP' )
          ELSEIF ( VAL .EQ. 3 ) THEN
              CALL QPPRM( 'Problem type      QP1' )
          ELSEIF ( VAL .EQ. 4 ) THEN
              CALL QPPRM( 'Problem type      QP2' )
          ELSEIF ( VAL .EQ. 5 ) THEN
              CALL QPPRM( 'Problem type      QP3' )
          ELSEIF ( VAL .EQ. 6 ) THEN
              CALL QPPRM( 'Problem type      QP4' )
          ENDIF
      ELSEIF ( OPTION .EQ. 11 ) THEN
          CALL QPPRMI( 'Summary file', VAL )
      ENDIF

      END

*======================================================================
      SUBROUTINE QPOPT_LOG_OPT ( OPTION, VAL )

      IMPLICIT NONE

      INTEGER OPTION
      LOGICAL VAL

*** Set an option for QPOPT that can be represented using a LOGICAL.
***
***		WARM_START  = 1
***		LIST        = 2
***		MIN_SUM     = 3
***
*** See the documentation from QPOPT to see a description of
*** these options.
***

      EXTERNAL QPPRM
      EXTERNAL QPPRMI
      EXTERNAL QPPRMR

      IF     ( OPTION .EQ. 1  ) THEN
          IF( VAL ) THEN
              CALL QPPRM( 'Warm start' )
          ELSE
              CALL QPPRM( 'Cold start' )
          ENDIF
      ELSEIF ( OPTION .EQ. 2  ) THEN
          IF( VAL ) THEN
              CALL QPPRM( 'List' )
          ELSE
              CALL QPPRM( 'NoList' )
          ENDIF
      ELSEIF ( OPTION .EQ. 3  ) THEN
          IF( VAL ) THEN
              CALL QPPRM( 'Min sum = Yes' )
          ELSE
              CALL QPPRM( 'Min sum = No' )
          ENDIF
      ENDIF

      END
*======================================================================
      SUBROUTINE QPOPT_REAL_OPT ( OPTION, VAL )

      IMPLICIT NONE

      INTEGER OPTION
      DOUBLE PRECISION VAL

*** Set an option for QPOPT that can be represented using a
*** DOUBLE PRECISION.
***
***		CRASH_TOLERANCE           = 1
***		FEASIBILITY_TOLERANCE     = 2
***		INFINITE_BOUND_SIZE       = 3
***		INFINITE_STEP_SIZE        = 4
***		OPTIMALITY_TOLERANCE      = 5
***		RANK_TOLERANCE            = 6
***
*** See the documentation from QPOPT to see a description of
*** these options.
***

      EXTERNAL QPPRM
      EXTERNAL QPPRMI
      EXTERNAL QPPRMR

      IF     ( OPTION .EQ. 1  ) THEN
          CALL QPPRMR( 'Crash tolerance', VAL )
      ELSEIF ( OPTION .EQ. 2  ) THEN
          CALL QPPRMR( 'Feasibility tolerance', VAL )
      ELSEIF ( OPTION .EQ. 3  ) THEN
          CALL QPPRMR( 'Infinite Bound size', VAL )
      ELSEIF ( OPTION .EQ. 4  ) THEN
          CALL QPPRMR( 'Infinite Step size', VAL )
      ELSEIF ( OPTION .EQ. 5  ) THEN
          CALL QPPRMR( 'Optimality tolerance', VAL )
      ELSEIF ( OPTION .EQ. 6  ) THEN
          CALL QPPRMR( 'Rank tolerance', VAL )
      ENDIF

      END
