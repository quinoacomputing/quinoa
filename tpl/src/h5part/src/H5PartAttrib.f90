
!< \ingroup h5partf_attrib
!! See \ref H5PartWritefileAttribFloat64
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writefileattrib_r8 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    REAL*8, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readfileattrib_r8 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    REAL*8, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritefileAttribFloat32
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writefileattrib_r4 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    REAL*4, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readfileattrib_r4 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    REAL*4, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritefileAttribInt64
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writefileattrib_i8 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    INTEGER*8, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readfileattrib_i8 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    INTEGER*8, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritefileAttribInt32
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writefileattrib_i4 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    INTEGER*4, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readfileattrib_i4 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    INTEGER*4, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritestepAttribFloat64
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writestepattrib_r8 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    REAL*8, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readstepattrib_r8 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    REAL*8, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritestepAttribFloat32
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writestepattrib_r4 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    REAL*4, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readstepattrib_r4 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    REAL*4, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritestepAttribInt64
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writestepattrib_i8 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    INTEGER*8, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readstepattrib_i8 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    INTEGER*8, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION

!< \ingroup h5partf_attrib
!! See \ref H5PartWritestepAttribInt32
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_writestepattrib_i4 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the attribute
    INTEGER*4, INTENT(IN) :: data(*) !< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelem !< the number of elements in the array
END FUNCTION

!< \ingroup h5partf_attrib
!! Read the attribute \c name into the buffer \c data.
!! \return 0 on success or error code
!>
INTEGER*8 FUNCTION h5pt_readstepattrib_i4 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name   !< the name of the attribute
    INTEGER*4, INTENT(OUT) :: data(*) !< buffer to read value into
END FUNCTION
