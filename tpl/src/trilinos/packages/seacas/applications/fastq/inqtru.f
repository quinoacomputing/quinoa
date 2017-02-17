C $Id: inqtru.f,v 1.1 1990/11/30 11:10:04 gdsjaar Exp $
C $Log: inqtru.f,v $
C Revision 1.1  1990/11/30 11:10:04  gdsjaar
C Initial revision
C
C
      SUBROUTINE INQTRU (PROMPT, IANS)
C***********************************************************************
C
C  SUBROUTINE INQTRU = INPUTS A YES OR NO ANSWER
C
C***********************************************************************
C
      CHARACTER* (*) PROMPT
      CHARACTER*1 RESULT, ANS
      LOGICAL IANS
      DIMENSION ANS (4)
      DATA ANS / 'Y', 'y', 'N', 'n' /
C
  100 CONTINUE
      WRITE (*, 10000)PROMPT
      READ (*, 10010, END = 110, ERR = 120)RESULT
      IF ( (RESULT (1:1) .EQ. ANS (1)) .OR.
     &   (RESULT (1:1) .EQ. ANS (2))) THEN
         IANS = .TRUE.
      ELSEIF ( (RESULT (1:1) .EQ. ANS (3)) .OR.
     &   (RESULT (1:1) .EQ. ANS (4))) THEN
         IANS = .FALSE.
      ELSE
         WRITE (*, 10020)
         GOTO 100
      ENDIF
      RETURN
  110 CONTINUE
      WRITE (*, 10030)
      GOTO 100
  120 CONTINUE
      WRITE (*, 10040)
      GOTO 100
C
10000 FORMAT (' ', A, '? ')
10010 FORMAT (A1)
10020 FORMAT (' RESPONSE MUST BE EITHER YES OR NO  -  TRY AGAIN')
10030 FORMAT (' END OF DATA ENCOUNTERED  -  TRY AGAIN')
10040 FORMAT (' ERROR IN RESPONSE  -  TRY AGAIN')
      END
