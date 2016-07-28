! Declaration of subroutines for Fortran Bindings

!> \defgroup h5part_f90_api H5Part F90 API

!> \ingroup h5part_f90_api
!! \defgroup h5partf_open       File Opening and Closing
!<

!> \ingroup h5part_f90_api
!! \defgroup h5partf_model      Setting up the Data Model
!<

!> \ingroup h5part_f90_api
!! \defgroup h5partf_data       Reading and Writing Datasets
!<

!> \ingroup h5part_f90_api
!! \defgroup h5partf_attrib     Reading and Writing Attributes
!<


!!!!!!!! File Opening and Closing !!!!!!!!

!> \ingroup h5partf_open
!! Opens a file for reading. See \ref H5PartOpenFile
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openr ( filename )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for reading
END FUNCTION

!> \ingroup h5partf_open
!! Opens a file for writing in truncate mode. See \ref H5PartOpenFile
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openw ( filename )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for writing
END FUNCTION

!> \ingroup h5partf_open
!! Opens a file for writing in append mode. See \ref H5PartOpenFile
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_opena ( filename )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for appending
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for reading.
!! See \ref H5PartOpenFileParallel
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openr_par ( filename, mpi_communicator )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for reading
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI communicator used by the program
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for writing in truncate mode.
!! See \ref H5PartOpenFileParallel
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openw_par ( filename, mpi_communicator )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for writing
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for writing in append mode.
!! See \ref H5PartOpenFileParallel
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_opena_par ( filename, mpi_communicator )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for appending
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
END FUNCTION

!> \ingroup h5partf_open
!! Opens a file for reading and specifies an HDF5 alignment.
!! See \ref H5PartOpenFileAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openr_align ( filename, align )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for reading
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
END FUNCTION

!> \ingroup h5partf_open
!! Opens a file for writing in truncate mode and specifies an HDF5 alignment.
!! See \ref H5PartOpenFileAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openw_align ( filename, align )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for writing
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
END FUNCTION
 
!> \ingroup h5partf_open
!! Opens a file for writing in append mode and specifies an HDF5 alignment.
!! See \ref H5PartOpenFileAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_opena_align ( filename, align )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for appending
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for reading and specifies an HDF5 alignment.
!! See \ref H5PartOpenFileParallelAlign
!!
!! Flags are specified as a comma separated string that can include:
!!
!! - \c fs_lustre - enable optimizations for the Lustre file system
!! - \c vfd_mpiposix - use the HDF5 MPI-POSIX virtual file driver
!! - \c vfd_mpio_ind - use MPI-IO in indepedent mode
!!
!! See \ref H5PartOpenFileParallelAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openr_par_align ( filename, mpi_communicator, align, flags )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for reading
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
    CHARACTER(LEN=*), INTENT(IN) :: flags       !< additional flags
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for writing in truncate mode and specifies
!! an HDF5 alignment.
!!
!! Flags are specified as a comma separated string that can include:
!!
!! - \c fs_lustre - enable optimizations for the Lustre file system
!! - \c vfd_mpiposix - use the HDF5 MPI-POSIX virtual file driver
!! - \c vfd_mpio_ind - use MPI-IO in indepedent mode
!!
!! See \ref H5PartOpenFileParallelAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_openw_par_align ( filename, mpi_communicator, align, flags )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for writing
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
    CHARACTER(LEN=*), INTENT(IN) :: flags       !< additional flags
END FUNCTION

!> \ingroup h5partf_open
!! Opens a parallel file for writing in append mode and specifies
!! an HDF5 alignment.
!!
!! Flags are specified as a comma separated string that can include:
!!
!! - \c fs_lustre - enable optimizations for the Lustre file system
!! - \c vfd_mpiposix - use the HDF5 MPI-POSIX virtual file driver
!! - \c vfd_mpio_ind - use MPI-IO in indepedent mode
!!
!! See \ref H5PartOpenFileParallelAlign
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_opena_par_align ( filename, mpi_communicator, align, flags )
    CHARACTER(LEN=*), INTENT(IN) :: filename    !< the filename to open for appending
    INTEGER, INTENT(IN) :: mpi_communicator     !< the MPI_Communicator used by the program
    INTEGER*8, INTENT(IN) :: align              !< alignment value in bytes
    CHARACTER(LEN=*), INTENT(IN) :: flags       !< additional flags
END FUNCTION

!> \ingroup h5partf_open
!! Closes a file. See \ref H5PartCloseFile
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_close ( filehandle )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_open
!! See \ref H5PartSetVerbosityLevel
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_set_verbosity_level ( level )
    INTEGER*8, INTENT(IN) :: level      !< the level from 0 (no output) to 5 (most detailed)
END FUNCTION


!!!!!!!! Setting up the Data Model !!!!!!!!

!> \ingroup h5partf_model
!! See \ref H5PartSetNumParticles
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setnpoints ( filehandle, npoints )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: npoints    !< the number of particles on *this* processor
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartSetNumParticlesStrided
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setnpoints_strided ( filehandle, npoints, stride )
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: npoints    !< the number of particles on *this* processor
    INTEGER*8, INTENT(IN) :: stride     !< the stride value (e.g. the number of fields in the particle data array)
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartSetStep
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setstep (filehandle,step)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: step       !< a timestep value >= 1
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetNumSteps
!! \return the number of steps or error code
!<
INTEGER*8 FUNCTION h5pt_getnsteps (filehandle)      
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetNumDatasets
!! \return the number of datasets or error code
!<
INTEGER*8 FUNCTION h5pt_getndatasets (filehandle)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetNumParticles
!! \return the number of particles or error code
!<
INTEGER*8 FUNCTION h5pt_getnpoints (filehandle)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetDatasetName
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getdatasetname (filehandle,index,name)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: index              !< index of dataset to query (starting from 0)
    CHARACTER(LEN=*), INTENT(OUT) :: name       !< buffer to read the dataset name into 
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartSetView
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setview (filehandle,start,end)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: start      !< offset of the first particle in the view
    INTEGER*8, INTENT(IN) :: end        !< offset of the last particle in the view (inclusive)
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartSetViewIndices
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_setview_indices (filehandle,indices,nelem)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: indices(*) !< list of indicies to select in this view
    INTEGER*8, INTENT(IN) :: nelem      !< number of particles in the list
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartResetView
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_resetview (filehandle)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartResetView
!! \return 1 if true, 0 if false, or error code
!<
INTEGER*8 FUNCTION h5pt_hasview (filehandle)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_model
!! See \ref H5PartGetView
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getview (filehandle,start,end)
    INTEGER*8, INTENT(IN) :: filehandle !< the handle returned during file open
    INTEGER*8, INTENT(OUT) :: start     !< buffer to store the offset of the first particle in the view
    INTEGER*8, INTENT(OUT) :: end       !< buffer to store the offset of the last particle in the view (inclusive)
END FUNCTION


!!!!!!!! Reading and Writing Datasets !!!!!!!!

!> \ingroup h5partf_data
!! See \ref H5PartWriteDataFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writedata_r8 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    REAL*8, INTENT(IN) :: data(*)               !< the array of float64 data to write
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartWriteDataFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writedata_r4 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    REAL, INTENT(IN) :: data(*)                 !< the array of float32 data to write
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartWriteDataInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writedata_i8 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    INTEGER*8, INTENT(IN) :: data(*)            !< the array of int64 data to write
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartWriteDataInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writedata_i4 ( filehandle, name, data )
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    INTEGER, INTENT(IN) :: data(*)              !< the array of int32 data to write
END FUNCTION


!> \ingroup h5partf_data
!! See \ref H5PartReadDataFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readdata_r8 (filehandle,name,data)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    REAL*8, INTENT(OUT) :: data(*)              !< array to read float64 data into
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartReadDataFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readdata_r4 (filehandle,name,data)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    REAL, INTENT(OUT) :: data(*)                !< array to read float32 data into
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartReadDataInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readdata_i8 (filehandle,name,data)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    INTEGER*8, INTENT(OUT) :: data(*)           !< array to read int64 data into
END FUNCTION

!> \ingroup h5partf_data
!! See \ref H5PartReadDataInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readdata_i4 (filehandle,name,data)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the dataset
    INTEGER, INTENT(OUT) :: data(*)             !< array to read int32 data into
END FUNCTION


!!!!!!!! Reading and Writing Attributes !!!!!!!!

!> \ingroup h5partf_attrib
!! See \ref H5PartWriteFileAttribString
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writefileattrib_string (filehandle,attrib_name,attrib_value)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the attribute
    CHARACTER(LEN=*), INTENT(IN) :: value       !< the string value to store
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartWriteStepAttribString
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_writestepattrib_string (filehandle,attrib_name,attrib_value)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the attribute
    CHARACTER(LEN=*), INTENT(IN) :: value       !< the string value to store
END FUNCTION

!> \ingroup h5partf_attrib
!! Reads the attribute \c name in the file root ("/")
!! into the string buffer \c value.
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readfileattrib_string (filehandle,attrib_name,attrib_value)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the attribute
    CHARACTER(LEN=*), INTENT(OUT) :: value      !< buffer to read the string value into
END FUNCTION

!> \ingroup h5partf_attrib
!! Reads the attribute \c name in the current step
!! into the string buffer \c value.
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_readstepattrib_string (filehandle,attrib_name,attrib_value)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    CHARACTER(LEN=*), INTENT(IN) :: name        !< the name of the attribute
    CHARACTER(LEN=*), INTENT(OUT) :: value      !< buffer to read the string value into
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartGetNumStepAttribs
!! \return number of attributes or error code
!<
INTEGER*8 FUNCTION h5pt_getnstepattribs (filehandle)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartGetNumFileAttribs
!! \return number of attributes or error code
!<
INTEGER*8 FUNCTION h5pt_getnfileattribs (filehandle)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartGetStepAttribInfo
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getstepattribinfo (filehandle,idx,attrib_name,attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: index              !< index of the attribute to query (starting from 0)
    CHARACTER(LEN=*), INTENT(OUT) :: name       !< buffer to read the attribute name into
    INTEGER*8, INTENT(OUT) :: nelem             !< number of elements in the attribute's array
END FUNCTION

!> \ingroup h5partf_attrib
!! See \ref H5PartGetFileAttribInfo
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5pt_getfileattribinfo (filehandle,idx,attrib_name,attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle         !< the handle returned during file open
    INTEGER*8, INTENT(IN) :: index              !< index of the attribute to query (starting from 0)
    CHARACTER(LEN=*), INTENT(OUT) :: name       !< buffer to read the attribute name into
    INTEGER*8, INTENT(OUT) :: nelem             !< number of elements in the attribute's array
END FUNCTION


