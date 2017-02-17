C $Id: delhol.f,v 1.1 1990/11/30 11:05:55 gdsjaar Exp $
C $Log: delhol.f,v $
C Revision 1.1  1990/11/30 11:05:55  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]DELHOL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DELHOL (IPTR, MXND, LXK, KXL, NXL, LXN, NXH, NUID, NNN,
     &   IAVAIL, NAVAIL, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE DELNOD  =  DELETES ALL LINES, ELEMENTS, ETC. ATTACHED TO
C                        A NODE
C
C***********************************************************************
C
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), NXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), NXH(MXND), NUID(MXND)
      DIMENSION KLIST(20)
C
      LOGICAL ERR, NOROOM
C
      ERR = .FALSE.
      CALL GKXN (MXND, KXL, LXN, IPTR, KS, KLIST, ERR)
      IF (.NOT. ERR) THEN
         ERR = .TRUE.
C
C  DELELE LINES PER ELEMENT; MARK ELEMENT NODES
         DO 110 J = 1, KS
            DO 100 K = 1, 4
               IF (NXH(NXL(1, LXK(K, KLIST(J)))) .EQ. 0)
     &            NXH(NXL(1, LXK(K, KLIST(J)))) = 1
               IF (NXH(NXL(2, LXK(K, KLIST(J)))) .EQ. 0)
     &            NXH(NXL(2, LXK(K, KLIST(J)))) = 1
               IF (KXL(1, LXK(K, KLIST(J))) .EQ. KLIST(J))
     &            KXL(1, LXK(K, KLIST(J))) = 0
               IF (KXL(2, LXK(K, KLIST(J))) .EQ. KLIST(J))
     &            KXL(2, LXK(K, KLIST(J))) = 0
               LXK(K, KLIST(J)) = 0
  100       CONTINUE
  110    CONTINUE
         NXH(IPTR) = -1
C
         DO 120 J = 1, 3
            IF (LXN(J, IPTR) .GT. 0) THEN
C
C  DELETE LINE ATTACHED TO OPPOSITE END NODE
C
               IF (NXL(1, LXN(J, IPTR)) .EQ. IPTR) THEN
                  CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &               NXL(2, LXN(J, IPTR)), LXN(J, IPTR), NNN, ERR,
     &               NOROOM)
               ELSE IF (NXL(2, LXN(J, IPTR)) .EQ. IPTR) THEN
                  CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &               NXL(1, LXN(J, IPTR)), LXN(J, IPTR), NNN, ERR,
     &               NOROOM)
               END IF
               IF (NOROOM) RETURN
C
C  DELETE NODES PER LINE; ELEMENTS PER LINE
C
               NXL(1, LXN(J, IPTR)) = 0
               NXL(2, LXN(J, IPTR)) = 0
               KXL(1, LXN(J, IPTR)) = 0
               KXL(2, LXN(J, IPTR)) = 0
               LXN(J, IPTR) = 0
            END IF
  120    CONTINUE
C
C  FOR LAST LINE, SAVE LINK ON IAVAIL
C
         IF (LXN(4, IPTR) .GT. 0) THEN
            IF (NXL(1, LXN(4, IPTR)) .EQ. IPTR) THEN
               CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &            NXL(2, LXN(4, IPTR)), LXN(4, IPTR), NNN, ERR,
     &            NOROOM)
            ELSE IF (NXL(2, LXN(J, IPTR)) .EQ. IPTR) THEN
               CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &            NXL(1, LXN(4, IPTR)), LXN(4, IPTR), NNN, ERR,
     &            NOROOM)
            END IF
            IF (NOROOM) RETURN
            NXL(1, LXN(4, IPTR)) = 0
            NXL(2, LXN(4, IPTR)) = 0
            KXL(1, LXN(4, IPTR)) = 0
            KXL(2, LXN(4, IPTR)) = 0
         END IF
C
         LXN(4, IPTR) = IAVAIL
         IAVAIL = IPTR
         NAVAIL = NAVAIL + 1
C
         ERR = .FALSE.
      END IF
C
      RETURN
      END
