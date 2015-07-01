
PROGRAM MT19937_Example

   PARAMETER (N=100000000)
   INTEGER i,MT19937_SSE_generate,sum
   REAL val,rsum

   TYPE MT19937_state
     INTEGER z(625)
   END TYPE MT19937_state
   TYPE(MT19937_state) state

   TYPE MT19937_sse_state
     INTEGER z(3770)
   END TYPE MT19937_sse_state
   TYPE(MT19937_sse_state) SSE_state

   CALL MT19937_init ( state )
   CALL MT19937_Get_SSE_state (state,SSE_state)
   CALL MT19937_Print_State (state)

   sum = 0

   DO i=1,N
     sum = sum + MT19937_SSE_generate (SSE_state)
   END DO

   WRITE(*,1) N
   rsum = sum/4294967296.
   val = MT19937_SSE_generate (SSE_state)/4294967296.
   IF (rsum<0) THEN
     rsum=rsum+1
   ENDIF
   IF (val<0) THEN
     val=val+1
   ENDIF
   WRITE(*,2) rsum
   WRITE(*,3) val

1  FORMAT(I9," MT19937_SSE pseudorandom numbers generated.")
2  FORMAT("Fractional part of the total sum of generated numbers: ",F8.6)
3  FORMAT("Next output value: ",F8.6)

END PROGRAM MT19937_Example

