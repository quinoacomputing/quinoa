
PROGRAM MRG32K3A_Example

   PARAMETER (N=100000000)
   INTEGER i,MRG32K3A_generate,sum
   REAL val,rsum

   TYPE MRG32K3A_state
     INTEGER z(6)
   END TYPE MRG32K3A_state
   TYPE(MRG32K3A_state) state

   CALL MRG32K3A_init ( state )
   CALL MRG32K3A_Print_State (state)

   sum = 0

   DO i=1,N
     sum = sum + MRG32K3A_generate (state)
   END DO

   WRITE(*,1) N
   rsum = sum/4294967296.
   val = MRG32K3A_generate (state)/4294967296.
   IF (rsum<0) THEN
     rsum=rsum+1
   ENDIF
   IF (val<0) THEN
     val=val+1
   ENDIF
   WRITE(*,2) rsum
   WRITE(*,3) val

1  FORMAT(I9," MRG32K3A pseudorandom numbers generated.")
2  FORMAT("Fractional part of the total sum of generated numbers: ",F8.6)
3  FORMAT("Next output value: ",F8.6)

END PROGRAM MRG32K3A_Example

