*** Utility functions for passing strings between C++ and Fortran.
***
*** See FortranTypes_CppFortranStrings.hpp for documentation
***

*** ToDo: 8/9/99: Make sure that the ASCII conversions are legal
*** on the target platform.

      INTEGER FUNCTION CONVERT_FROM_C_INT_STR( I_STRING, STRING_LEN
     & , STRING )

*     *** Formal parameters
*     *** In
      INTEGER I_STRING(*), STRING_LEN
*     *** Out
      CHARACTER*(*) STRING

*     ***
*     *** Convert from an integer representation of a string to
*     *** a fortran CHARACTER string.
*     ***
*     *** Returns 0 if conversion is successful.  Otherwise returns
*     *** the indice of the ith character that was not ASCII converable.
*     ***

*     *** Local variables
      INTEGER I

      STRING = ''
      DO 10 I = 1,STRING_LEN
          STRING(I:I) = CHAR( I_STRING(I) )
10    CONTINUE

      CONVERT_FROM_C_INT_STR = 0

      END

      INTEGER FUNCTION CONVERT_TO_C_INT_STR( STRING, I_STRING
     &, STRING_LEN )

*     *** Formal parameters
*     *** In
      CHARACTER*(*) STRING
*     *** Out
      INTEGER I_STRING(*), STRING_LEN

*     ***
*     *** Convert from a fortran CHARACTER string to an integer
*     *** representation suitable for passing to C++.
*     ***
*     *** Returns 0 if conversion is successful.  Otherwise returns
*     *** the indice of the ith character that was not ASCII converable.
*     ***

*     *** Local variables
      INTEGER I

      STRING_LEN = LEN( STRING )

      DO 10 I = 1,STRING_LEN
          I_STRING(I) = ICHAR( STRING(I:I) )
10    CONTINUE

      CONVERT_TO_C_INT_STR = 0

      END
