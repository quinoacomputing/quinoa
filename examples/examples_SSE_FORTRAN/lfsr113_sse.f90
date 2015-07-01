
PROGRAM LFSR113_Example

   PARAMETER (N=100000000)
   INTEGER i,LFSR113_SSE_generate,sum
   REAL val,rsum

   TYPE LFSR113_state
     INTEGER z(4)
   END TYPE LFSR113_state
   TYPE(LFSR113_state) state

   TYPE LFSR113_SSE_state
     INTEGER z(8)
   END TYPE LFSR113_SSE_state
   TYPE(LFSR113_SSE_state) SSE_state

   CALL LFSR113_init ( state )
   CALL LFSR113_Get_SSE_state (state,SSE_state)
   CALL LFSR113_Print_State (state)

   sum = 0

   DO i=1,N
     sum = sum + LFSR113_SSE_generate (SSE_state)
   END DO

   WRITE(*,1) N
   rsum = sum/4294967296.
   val = LFSR113_SSE_generate (SSE_state)/4294967296.
   IF (rsum<0) THEN
     rsum=rsum+1
   ENDIF
   IF (val<0) THEN
     val=val+1
   ENDIF
   WRITE(*,2) rsum
   WRITE(*,3) val

1  FORMAT(I9," LFSR113_SSE pseudorandom numbers generated.")
2  FORMAT("Fractional part of the total sum of generated numbers: ",F8.6)
3  FORMAT("Next output value: ",F8.6)

END PROGRAM LFSR113_Example

