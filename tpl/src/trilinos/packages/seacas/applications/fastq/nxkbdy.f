C $Id: nxkbdy.f,v 1.2 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: nxkbdy.f,v $
C Revision 1.2  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:12:58  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:12:57  gdsjaar
c Initial revision
c 
C
CC* FILE: [.RENUM]NXKBDY.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NXKBDY (MDIM, NNXK, NPELEM, NXK, KKK, LIST, THREE,
     &   EIGHT, NINE)
C***********************************************************************
C
C  SUBROUTINE NXKBDY = FLAGS ALL SIDES OF ELEMENTS ONLY ONCE BY USE OF
C                      A HASH SCHEME
C
C***********************************************************************
C
C     NXK    = THE NODES PER ELEMENT ARRAY  (CONNECTIVITY)
C               (A NEGATIVE VALUE WILL INDICATE UNIQUENESS OF FOLLOWING
C              SIDE)
C***********************************************************************
C
      LOGICAL THREE, EIGHT, NINE, ITSOK
C
      DIMENSION NXK (NNXK, NPELEM), LIST (MDIM)
C
      DO 100 I = 1, MDIM
         LIST (I)  =  0
  100 CONTINUE
C
      DO 140 K = 1, KKK
         IF ((NXK (3, K) .EQ. 0) .AND. (THREE)) THEN
            ITSOK = .TRUE.
            NEND = 1
         ELSEIF  ((NXK (3, K) .NE. 0) .AND. ((EIGHT) .OR. (NINE))) THEN
            ITSOK = .TRUE.
            NEND = 4
         ELSE
            ITSOK = .FALSE.
         ENDIF
         IF (ITSOK) THEN
            DO 130 N = 1, NEND
C
C  COMPUTE HASH CODE FOR LINE
C
               N2 = N + 1
               IF (N .GE. 4)N2 = 1
               NODE1 = NXK (N, K)
               NODE2 = IABS (NXK (N2, K))
               IF (NODE2 .GT. NODE1) THEN
                  LVAL = NODE1 * 100000 + NODE2
               ELSE
                  LVAL = NODE2 * 100000 + NODE1
               ENDIF
C
C  CALCULATE THE BEGINNING HASH VALUE
C
               HOLD = FLOAT (NODE1 + NODE2) * 3.1830989
               LHASH = INT((HOLD-IFIX (HOLD)) * FLOAT (MDIM) + 1)
               LKEEP = LHASH
C
C  FIRST-TIMERS CLAIM THE NODE
C
  110          CONTINUE
               IF (LIST (LHASH).NE.0) THEN
                  IF (LIST (LHASH) .EQ. LVAL) GOTO 120
                  LHASH = LHASH + 1
                  IF (LHASH .EQ. MDIM) LHASH = 1
                  IF (LHASH.NE.LKEEP) GOTO 110
                  CALL MESAGE ('HASH SCHEME OVERFLOW IN KXNBDY')
                  STOP
               ENDIF
C
               LIST (LHASH) = LVAL
               NXK (N, K) = ISIGN (NXK (N, K), -1)
  120          CONTINUE
  130       CONTINUE
         ENDIF
  140 CONTINUE
      RETURN
      END
