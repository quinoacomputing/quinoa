C $Id: nps.f,v 1.1 1990/11/30 11:12:52 gdsjaar Exp $
C $Log: nps.f,v $
C Revision 1.1  1990/11/30 11:12:52  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]NPS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS,
     &   ILLIST, LINKL, LINKS, NNPS, ERR)
C***********************************************************************
C
C  SUBROUTINE NPS = GIVES A LIST OF THE NUMBER OF PERIMETER NODES
C                   ON EACH OF A REGION'S SIDES
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES QUAD ELEMENTS
C
C***********************************************************************
C
C  VARIABLES USED:
C     NNPS  = ARRAY OF NUMBER OF NODES PER SIDE
C     CCW   = .TRUE. IF THE SIDE IS ORIENTED CCW
C     KS    = COUNTER OF THE NUMBER OF SIDES
C
C***********************************************************************
C
      DIMENSION NNPS (MNNPS), ISLIST (NS), LINKL (2, ML), LINKS (2, MS)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS * 3), NINT (ML)
C
      LOGICAL ERR, ADDLNK
C
      ERR = .TRUE.
      ADDLNK = .FALSE.
C
      KS = 0
      DO 110 I = 1, NS
         IF (ISLIST (I) .LT. 0) THEN
            KS = KS + 1
            IL = IABS (ISLIST (I))
            CALL LTSORT (ML, LINKL, IL, IPNTR, ADDLNK)
            NNPS (KS) = IABS (NINT (IPNTR)) + 1
         ELSEIF (ISLIST (I) .GT. 0) THEN
            CALL LTSORT (MS, LINKS, ISLIST (I), IPNTR, ADDLNK)
            J1 = IFLINE (IPNTR)
            J2 = J1 + NLPS (IPNTR) - 1
            KS = KS + 1
            NNPS (KS) = 0
            DO 100 J = J1, J2
               IL = ILLIST (J)
               CALL LTSORT (ML, LINKL, IL, IPNTR, ADDLNK)
               NNPS (KS) = NNPS (KS) + IABS (NINT (IPNTR))
  100       CONTINUE
            NNPS (KS) = NNPS (KS) + 1
         ELSE
            RETURN
         ENDIF
  110 CONTINUE
      ERR = .FALSE.
C
      RETURN
C
      END
