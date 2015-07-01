      INTEGER FUNCTION FOR_OPEN_FILE( IUNIT, I_FILE, I_FILE_LEN,
     &	ISTATUS, IFORM, IBLANK, IACCESS, IRECL )

      IMPLICIT NONE

*     *** In Parameters
      INTEGER IUNIT
      INTEGER I_FILE(*), I_FILE_LEN, ISTATUS
      INTEGER IFORM, IBLANK, IACCESS, IRECL
***
*** This function is ment to be called from C++ to open
*** a Fortran file.  The file name is a string stored as an
*** integer array.
***
*** This interface is closly tied with the C++ function
*** in FortranTypes_f_open_file.cpp.
***
*** Arguments:
***	IUNIT		[I] The unit number to give the opened file
***	I_FILE		[I]	Name of the file to open stored as an
***					integer array (maximun length = 100)
***	I_FILE_LEN	[I]	Length the name of I_FILE (<=100)
***	ISTATUS		[I] 0 = 'OLD', 1 = 'NEW', 2 = 'SCRATCH',
***					3 = 'UNKNOWN' (default)
***	IFORM		[I] 0 = 'FORMATTED' (default), 1 = 'UNFORMATTED'
***	IBLANK		[I] 0 = 'NULL' (default), 1 = 'ZERO'
***	IACCESS		[I] 0 = 'SEQUENTIAL' (default), 1 = 'DIRECT'
***	IRECL		[I] > 0 for IACCESS = 1, otherwise who cares
***
*** Return values:
*** 0	    Successfully opened file.
*** < 0.	Filename is not a valid ASCII string and -i will 
***	      be returned where i is the ith character in I_FILE_NAME
***       that could not be converted to a Fortran CHARACTER.
***	> 0		File could not be opened and this is the value of IOSTAT
***       returned from OPEN(...)
***

*     *** External procedures
      INTEGER CONVERT_FROM_C_INT_STR

*     *** Local variables (Ordered as they appear?)
      CHARACTER*100		S_FILE
      INTEGER			RESULT     
      CHARACTER*7		S_STATUS
      CHARACTER*10		S_ACCESS
      CHARACTER*11		S_FORM
      CHARACTER*4		S_BLANK
      INTEGER			IOS

*     *** Executable statements

*     *** Convert parameters to OPEN(...) input

      RESULT = CONVERT_FROM_C_INT_STR( I_FILE, I_FILE_LEN, S_FILE )
      IF( RESULT .NE. 0 ) THEN
          FOR_OPEN_FILE = - RESULT
          RETURN
      ENDIF

*     *** I am going to assume there will not be any errors
*     *** in converting these arguments since I only expect
*     *** The C++ wrapper function to call it and it is not
*     *** for the general user.

      IF		( ISTATUS .EQ. 0 ) THEN
          S_STATUS = 'OLD'
      ELSEIF	( ISTATUS .EQ. 1 ) THEN
          S_STATUS = 'NEW'
      ELSEIF	( ISTATUS .EQ. 2 ) THEN
          S_STATUS = 'SCRATCH'
      ELSEIF	( ISTATUS .EQ. 3 ) THEN
          S_STATUS = 'UNKNOWN'
      ENDIF

      IF		( IFORM .EQ. 0 ) THEN
          S_FORM = 'FORMATTED'
      ELSEIF	( IFORM .EQ. 1 ) THEN
          S_FORM = 'UNFORMATTED'
      ENDIF

      IF		( IBLANK .EQ. 0 ) THEN
          S_BLANK = 'NULL'
      ELSEIF	( IBLANK .EQ. 1 ) THEN
          S_BLANK = 'ZERO'
      ENDIF
   
      IF		( IACCESS .EQ. 0 ) THEN
          S_ACCESS = 'SEQUENTIAL'
      ELSEIF	( IACCESS .EQ. 1 ) THEN
          S_ACCESS = 'DIRECT'
      ENDIF

*     *** Try to open the file
      IF( S_ACCESS .EQ. 'DIRECT' ) THEN
          OPEN( IUNIT, IOSTAT = IOS, FILE = S_FILE, STATUS = S_STATUS
     &        , FORM = S_FORM, BLANK = S_BLANK, ACCESS = 'DIRECT'
     &        , RECL = IRECL ) 
      ELSE
          OPEN( IUNIT, IOSTAT = IOS, FILE = S_FILE, STATUS = S_STATUS
     &        , FORM = S_FORM, BLANK = S_BLANK, ACCESS = 'SEQUENTIAL' ) 
*          OPEN( IUNIT, IOSTAT = IOS, FILE = S_FILE, STATUS = 'UNKNOWN'
*     &        , FORM = 'FORMATTED', BLANK = 'NULL'
*     &        , ACCESS = 'SEQUENTIAL' ) 
      ENDIF

      FOR_OPEN_FILE = IOS

      END

*======================================================================

      SUBROUTINE FOR_CLOSE_FILE( IUNIT, KEEP )

      IMPLICIT NONE

*     *** In Parameters
      INTEGER IUNIT, KEEP
***
*** This function is ment to be called from C++ to close
*** a Fortran file.
***
*** Arguments:
***	IUNIT		[I] The unit number to give the opened file
***	KEEP		[I]	If KEEP == 1 then the file will be kept.

*     *** Local variables
      CHARACTER*7		S_STATUS
      INTEGER			IOS

*     *** Executable statements

      IF ( KEEP .EQ. 1 ) THEN
          S_STATUS = 'KEEP'
      ELSE
          S_STATUS = 'DELETE'
      ENDIF

      CLOSE( IUNIT, IOSTAT = IOS, STATUS = S_STATUS )

      END  

