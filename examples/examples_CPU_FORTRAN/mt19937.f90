
PROGRAM MT19937_Example

   PARAMETER (N=100000000)
   INTEGER i,MT19937_generate,sum
   REAL val,rsum

   TYPE MT19937_state
     INTEGER z(625)
   END TYPE MT19937_state
   TYPE(MT19937_state) state

   CALL MT19937_init ( state )
   CALL MT19937_Print_State (state)

   sum = 0

   DO i=1,N
     sum = sum + MT19937_generate (state)
   END DO

   WRITE(*,1) N
   rsum = sum/4294967296.
   val = MT19937_generate (state)/4294967296.
   IF (rsum<0) THEN
     rsum=rsum+1
   ENDIF
   IF (val<0) THEN
     val=val+1
   ENDIF
   WRITE(*,2) rsum
   WRITE(*,3) val

1  FORMAT(I9," MT19937 pseudorandom numbers generated.")
2  FORMAT("Fractional part of the total sum of generated numbers: ",F8.6)
3  FORMAT("Next output value: ",F8.6)

END PROGRAM MT19937_Example

