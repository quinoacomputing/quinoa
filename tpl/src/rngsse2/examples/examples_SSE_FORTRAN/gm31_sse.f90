
PROGRAM GM31_Example

   PARAMETER (N=100000000)
   INTEGER i,GM31_SSE_generate,sum
   REAL val,rsum

   TYPE GM31_state
     INTEGER z(64)
   END TYPE GM31_state
   TYPE(GM31_state) state

   TYPE GM31_sse_state
     INTEGER z(132)
   END TYPE GM31_sse_state
   TYPE(GM31_sse_state) SSE_state

   CALL GM31_init ( state )
   CALL GM31_Get_SSE_state (state,SSE_state)
   CALL GM31_Print_State (state)

   sum = 0

   DO i=1,N
     sum = sum + GM31_SSE_generate (SSE_state)
   END DO

   WRITE(*,1) N
   rsum = sum/4294967296.
   val = GM31_SSE_generate (SSE_state)/4294967296.
   IF (rsum<0) THEN
     rsum=rsum+1
   ENDIF
   IF (val<0) THEN
     val=val+1
   ENDIF
   WRITE(*,2) rsum
   WRITE(*,3) val

1  FORMAT(I9," GM31_SSE pseudorandom numbers generated.")
2  FORMAT("Fractional part of the total sum of generated numbers: ",F8.6)
3  FORMAT("Next output value: ",F8.6)

END PROGRAM GM31_Example

