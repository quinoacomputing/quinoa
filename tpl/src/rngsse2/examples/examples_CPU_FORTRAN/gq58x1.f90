
PROGRAM GQ58X1_Example

   PARAMETER (N=100000000)
   INTEGER i,GQ58X1_generate,sum
   REAL val,rsum

   TYPE GQ58X1_state
     INTEGER z(132)
   END TYPE GQ58X1_state
   TYPE(GQ58X1_state) state

   CALL GQ58X1_init ( state )
   CALL GQ58X1_Print_State (state)

   sum = 0

   DO i=1,N
     sum = sum + GQ58X1_generate (state)
   END DO

   WRITE(*,1) N
   rsum = sum/4294967296.
   val = GQ58X1_generate (state)/4294967296.
   IF (rsum<0) THEN
     rsum=rsum+1
   ENDIF
   IF (val<0) THEN
     val=val+1
   ENDIF
   WRITE(*,2) rsum
   WRITE(*,3) val

1  FORMAT(I9," GQ58X1 pseudorandom numbers generated.")
2  FORMAT("Fractional part of the total sum of generated numbers: ",F8.6)
3  FORMAT("Next output value: ",F8.6)

END PROGRAM GQ58X1_Example

