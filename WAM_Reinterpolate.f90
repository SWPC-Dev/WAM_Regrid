PROGRAM WAM_Interpolate

USE netcdf

IMPLICIT NONE

  ! /////////////////////////////////////////////////// !

  INTEGER, PARAMETER :: n_regridded_levels = 50
  REAL(8), PARAMETER :: z0_regridded       = 50.0D0
  REAL(8), PARAMETER :: zM_regridded       = 400.0D0


  ! /////////////////////////////////////////////////// !


  INTEGER :: i, j, k 
  INTEGER :: n_time, n_lat, n_lon, n_pressure_levels

  REAL(8) :: z_regridded(1:n_regridded_levels)
  REAL(8) :: z_regridded_log(1:n_regridded_levels)

  REAL(8), ALLOCATABLE :: A(:,:)
  REAL(8), ALLOCATABLE :: temperature(:,:,:,:)
  REAL(8), ALLOCATABLE :: density(:,:,:,:)
  REAL(8), ALLOCATABLE :: z_pressure(:,:,:,:)
  REAL(8), ALLOCATABLE :: temperature_reordered(:,:,:,:)
  REAL(8), ALLOCATABLE :: density_reordered(:,:,:,:)
  REAL(8), ALLOCATABLE :: z_pressure_reordered(:,:,:,:)
  REAL(8), ALLOCATABLE :: regridded_temperature(:,:,:,:)
  REAL(8), ALLOCATABLE :: regridded_density(:,:,:,:)
  REAL(8), ALLOCATABLE :: regridded_temperature_reordered(:,:,:,:)
  REAL(8), ALLOCATABLE :: regridded_density_reordered(:,:,:,:)
  REAL(8), ALLOCATABLE :: lat(:)
  REAL(8), ALLOCATABLE :: lon(:)
  REAL(8), ALLOCATABLE :: x(:,:,:), y(:,:,:), z(:,:,:)
  INTEGER, ALLOCATABLE :: time(:)
  CHARACTER(200)       :: inputFile, outputFile
  LOGICAL :: success


  CALL Setup( )

  IF( success )THEN


    CALL Read_NetCDF( inputFile )

    lat = lat*90.0/88.5419616699219

!    CALL Convert_to_XYZ( )

    CALL Reorder_Pressure_Arrays( )

    !$OMP PARALLEL DO
    DO k = 1, n_time
      DO j = 1, n_lat
        DO i = 1, n_lon

          CALL Interpolate( z_pressure_reordered(1:n_pressure_levels,i,j,k), z_regridded, temperature_reordered(1:n_pressure_levels,i,j,k), regridded_temperature_reordered(1:n_regridded_levels,i,j,k) )
          CALL Interpolate( z_pressure_reordered(1:n_pressure_levels,i,j,k), z_regridded, density_reordered(1:n_pressure_levels,i,j,k), regridded_density_reordered(1:n_regridded_levels,i,j,k) )


        ENDDO
      ENDDO
    ENDDO
!    !$OMP PARALLEL ENDDO

    CALL Reorder_Regridded_Arrays( )

    CALL Write_NetCDF( outputFile )

    CALL Cleanup( )

  ENDIF



CONTAINS

  SUBROUTINE Setup( )
    ! Local
    INTEGER :: i, nArg, argID
    LOGICAL :: inputGiven, outputGiven, helpNeeded, fileExists
    CHARACTER(200) :: argName
 
     inputGiven  = .FALSE.
     outputGiven = .FALSE.
     helpNeeded  = .FALSE.
     outputFile  = ''
     inputFile   = ''

     success = .FALSE. 

     nArg = command_argument_count( )

     IF( nArg > 0 )THEN
 
        DO argID = 1, nArg
  
          CALL get_command_argument( argID, argName )

          SELECT CASE( TRIM(argName) )

             CASE("-i")
                inputGiven = .TRUE.

             CASE("--input-file")
                inputGiven = .TRUE.

             CASE("-o")
                outputGiven = .TRUE.

             CASE("--output-file")
                outputGiven = .TRUE.

             CASE("-h")

                 PRINT*, '    A tool for regridding WAM output from pressure levels to fixed height grid'
                 PRINT*, '--------------------------------------------------------------'
                 PRINT*, '  Usage : regrid -i input_file_name.nc [optional arguments]                                      '
                 PRINT*, ''
                 PRINT*, ' Required inputs'
                 PRINT*, ''
                 PRINT*, ' -i /full/path/to/wam/netcdf/file '
                 PRINT*, '                OR  '
                 PRINT*, ' --input-file /full/path/to/wam/netcdf/file '
                 PRINT*, ''
                 PRINT*, ' Optional inputs'
                 PRINT*, ''
                 PRINT*, ' -h --help '
                 PRINT*, '  Display this help message.'
                 PRINT*, ''
                 PRINT*, ' -o /full/path/regridded/ouput/file '
                 PRINT*, '                OR  '
                 PRINT*, ' --output-file /full/path/regridded/ouput/file '
                 PRINT*, ''
                 PRINT*, ' --output-directory /full/path/to/output directory '
                 PRINT*, ''
                 PRINT*, ' If the output file is not included, the input file name '
                 PRINT*, ' will be used with "regridded" prepended.'
                 PRINT*, ''
                 PRINT*, '--------------------------------------------------------------'
                 success = .FALSE.
                 RETURN

             CASE("--help")

                 PRINT*, '  regrid '
                 PRINT*, '    A tool for regridding WAM output from pressure levels to fixed height grid'
                 PRINT*, '--------------------------------------------------------------'
                 PRINT*, '  Usage : regrid -i input_file_name.nc [optional arguments]                                      '
                 PRINT*, ''
                 PRINT*, ' Required inputs'
                 PRINT*, ''
                 PRINT*, ' -i /full/path/to/wam/netcdf/file '
                 PRINT*, '                OR  '
                 PRINT*, ' --input-file /full/path/to/wam/netcdf/file '
                 PRINT*, ''
                 PRINT*, ' Optional inputs'
                 PRINT*, ''
                 PRINT*, ' -h --help '
                 PRINT*, '  Display this help message.'
                 PRINT*, ''
                 PRINT*, ' -o /full/path/regridded/ouput/file '
                 PRINT*, '                OR  '
                 PRINT*, ' --output-file /full/path/regridded/ouput/file '
                 PRINT*, ''
                 PRINT*, ' --output-directory /full/path/to/output directory '
                 PRINT*, ''
                 PRINT*, ' If the output file is not included, the input file name '
                 PRINT*, ' will be used with "regridded" prepended.'
                 PRINT*, ''
                 PRINT*, '--------------------------------------------------------------'
                 success = .FALSE.
                 RETURN

             CASE DEFAULT

               IF( inputGiven )THEN

                  inputFile  = TRIM(argName) ! capture list of plasma files
                  inputGiven = .FALSE.

                  success = .TRUE.

               ELSEIF( outputGiven )THEN

                  outputFile  = TRIM(argName) ! capture list of neutral files
                  outputGiven = .FALSE.

               ENDIF

          END SELECT 
        ENDDO

     ENDIF

     IF( TRIM(inputFile) == '' )THEN
       PRINT*, 'No input file given'
       success = .FALSE.
       RETURN
     ENDIF

     INQUIRE( FILE = inputFile, EXIST = fileExists )

     IF( .NOT. fileExists )THEN
       PRINT*, ' File '//TRIM(inputFile)//' does not exist.'
       success = .FALSE.
       RETURN
     ENDIF

     DO i = 1, n_regridded_levels
       z_regridded(i) = z0_regridded + ( zM_regridded - z0_regridded )/REAL( n_regridded_levels-1, 8 ) * REAL(i-1)
       z_regridded_log(i) = log( z_regridded(i) )
     ENDDO


     IF( TRIM(outputfile) == '' )THEN

       outputfile = 'regridded.'//TRIM(inputFile)

     ENDIF
    
     PRINT*, 'Input file  : '//TRIM( inputFile )
     PRINT*, 'Output file : '//TRIM( outputFile )

  END SUBROUTINE Setup
!
  SUBROUTINE Cleanup( )
    ! Clear any allocatable memory, and make sure files are closed
     DEALLOCATE( A, &
                 temperature, &
                 density, &
                 z_pressure, &
                 temperature_reordered, &
                 density_reordered, &
                 z_pressure_reordered, &
                 regridded_temperature, &
                 regridded_density, &
                 regridded_temperature_reordered, &
                 regridded_density_reordered, &
                 lat, &
                 lon, &
                 time )

      DEALLOCATE( x, y, z )
    
  END SUBROUTINE Cleanup
!
  SUBROUTINE Allocate_Arrays( )

  ALLOCATE( A(1:n_pressure_levels, 1:n_regridded_levels), &
            temperature(1:n_lon,1:n_lat,1:n_pressure_levels,1:n_time), &
            density(1:n_lon,1:n_lat,1:n_pressure_levels,1:n_time), &
            z_pressure(1:n_lon,1:n_lat,1:n_pressure_levels,1:n_time), &
            temperature_reordered(1:n_pressure_levels,1:n_lon,1:n_lat,1:n_time), &
            density_reordered(1:n_pressure_levels,1:n_lon,1:n_lat,1:n_time), &
            z_pressure_reordered(1:n_pressure_levels,1:n_lon,1:n_lat,1:n_time), &
            regridded_temperature(1:n_lon,1:n_lat,1:n_regridded_levels,1:n_time), &
            regridded_density(1:n_lon,1:n_lat,1:n_regridded_levels,1:n_time), &
            regridded_temperature_reordered(1:n_regridded_levels,1:n_lon,1:n_lat,1:n_time), &
            regridded_density_reordered(1:n_regridded_levels,1:n_lon,1:n_lat,1:n_time), &
            lat(1:n_lat), &
            lon(1:n_lon), &
            time(1:n_time) )

  ALLOCATE( x(1:n_lon,1:n_lat,1:n_regridded_levels), &
            y(1:n_lon,1:n_lat,1:n_regridded_levels), &
            z(1:n_lon,1:n_lat,1:n_regridded_levels) )

  END SUBROUTINE Allocate_Arrays

  SUBROUTINE Reorder_Pressure_Arrays( )

    ! Local
    INTEGER :: i, j, k 
    REAL(8) :: max_z


    max_z = MAXVAL( z_pressure )
    DO k = 1, n_time
      DO j = 1, n_lat
        DO i = 1, n_lon

          temperature_reordered(1:n_pressure_levels,i,j,k) = temperature(i,j,1:n_pressure_levels,k)
          
          z_pressure_reordered(1:n_pressure_levels,i,j,k)  = z_pressure(i,j,1:n_pressure_levels,k) 
          density_reordered(1:n_pressure_levels,i,j,k)     = density(i,j,1:n_pressure_levels,k)

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE Reorder_Pressure_Arrays

  SUBROUTINE Reorder_Regridded_Arrays( )

    ! Local
    INTEGER :: i, j, k 

    DO k = 1, n_time
      DO j = 1, n_lat
        DO i = 1, n_lon

          regridded_temperature(i,j,1:n_regridded_levels,k) = regridded_temperature_reordered(1:n_regridded_levels,i,j,k) 
          regridded_density(i,j,1:n_regridded_levels,k)     = regridded_density_reordered(1:n_regridded_levels,i,j,k)     

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE Reorder_Regridded_Arrays

  SUBROUTINE Convert_to_XYZ( )
    INTEGER :: i, j, k

    DO k = 1, n_regridded_levels
      DO j = 1, n_lat
        DO i = 1, n_lon
 
          x(i,j,k) = ( z_regridded(k) + 6371.0D0 )*sin( 3.141592653D0*lon(i)/180.0D0 )*cos( 3.141592653D0*lat(i)/180.0D0 )
          y(i,j,k) = ( z_regridded(k) + 6371.0D0 )*cos( 3.141592653D0*lon(i)/180.0D0 )*cos( 3.141592653D0*lat(i)/180.0D0 )
          z(i,j,k) = ( z_regridded(k) + 6371.0D0 )*sin( 3.141592653D0*lat(i)/180.0D0 )

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE Convert_to_XYZ

  SUBROUTINE Read_NetCDF( filename )
    CHARACTER(*), INTENT(in) :: filename
    ! Local 
    CHARACTER(NF90_MAX_NAME) :: nameHolder
    INTEGER :: ncid, dimid, varid
    

      CALL Check( nf90_open( TRIM(filename), NF90_NETCDF4, ncid))


      CALL Check( nf90_inq_dimid( ncid, "time", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, n_time ) )

      CALL Check( nf90_inq_dimid( ncid, "lev", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, n_pressure_levels ) )

      CALL Check( nf90_inq_dimid( ncid, "lat", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, n_lat ) )

      CALL Check( nf90_inq_dimid( ncid, "lon", dimid ) )
      CALL Check( nf90_inquire_dimension( ncid, dimid, nameHolder, n_lon ) )

      CALL Allocate_Arrays( )

      CALL Check( nf90_inq_varid( ncid, "time", varid ) )
      CALL Check( nf90_get_var( ncid, varid, time ) )

      CALL Check( nf90_inq_varid( ncid, "lat", varid ) )
      CALL Check( nf90_get_var( ncid, varid, lat ) )

      CALL Check( nf90_inq_varid( ncid, "lon", varid ) )
      CALL Check( nf90_get_var( ncid, varid, lon ) )

      CALL Check( nf90_inq_varid( ncid, "z", varid ) )
      CALL Check( nf90_get_var( ncid, varid, z_pressure ) )

      CALL Check( nf90_inq_varid( ncid, "den", varid ) )
      CALL Check( nf90_get_var( ncid, varid, density ) )

      CALL Check( nf90_inq_varid( ncid, "t", varid ) )
      CALL Check( nf90_get_var( ncid, varid, temperature ) )

      CALL Check( nf90_close( ncid ) )

  END SUBROUTINE Read_NetCDF
!
  SUBROUTINE Write_NetCDF( filename )
    CHARACTER(*), INTENT(in) :: filename
    ! Local
    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid, z_dimid, time_dimid
    INTEGER :: x_varid, y_varid, z_varid, time_varid
    INTEGER :: temperature_varid, density_varid
    INTEGER :: recStart(1:4), recCount(1:4)


      recStart = (/ 1, 1, 1, 1 /)
      recCount = (/ n_lon, n_lat, n_regridded_levels, n_time /)

      CALL Check( nf90_create( TRIM(filename), NF90_NETCDF4, ncid))

      CALL Check( nf90_def_dim( ncid, "latitude", n_lat, y_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "longitude", n_lon, x_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "z", n_regridded_levels, z_dimid ) ) 
      CALL Check( nf90_def_dim( ncid, "time", NF90_UNLIMITED, time_dimid ) )

      CALL Check( nf90_def_var( ncid, "z", NF90_DOUBLE, z_dimid, z_varid ) )
      CALL Check( nf90_put_att( ncid, z_varid, "long_name", "Radial Distance" ) )
      CALL Check( nf90_put_att( ncid, z_varid, "units", "km" ) )
      CALL Check( nf90_put_att( ncid, z_varid, "units", "km" ) )

      CALL Check( nf90_def_var( ncid, "latitude", NF90_DOUBLE, y_dimid, y_varid ) )
      CALL Check( nf90_put_att( ncid, y_varid, "long_name", "Geographic Latitude" ) )
      CALL Check( nf90_put_att( ncid, y_varid, "units", "degrees north" ) )

      CALL Check( nf90_def_var( ncid, "longitude", NF90_DOUBLE, x_dimid, x_varid ) )
      CALL Check( nf90_put_att( ncid, x_varid, "long_name", "Geographic Longitude" ) )
      CALL Check( nf90_put_att( ncid, x_varid, "units", "degrees east" ) )

      CALL Check( nf90_def_var( ncid, "time", NF90_INT, time_dimid, time_varid ) )
      CALL Check( nf90_put_att( ncid, time_varid, "long_name", "minutes since an unknown time" ) )
      CALL Check( nf90_put_att( ncid, time_varid, "units", "minutes" ) )

      CALL Check( nf90_def_var( ncid, "temperature", NF90_DOUBLE, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , temperature_varid ) )
      CALL Check( nf90_put_att( ncid, temperature_varid, "long_name", "Thermosphere Temperature" ) )
      CALL Check( nf90_put_att( ncid, temperature_varid, "units", "K" ) )

      CALL Check( nf90_def_var( ncid, "density", NF90_DOUBLE, (/ x_dimid, y_dimid, z_dimid, time_dimid /) , density_varid ) )
      CALL Check( nf90_put_att( ncid, density_varid, "long_name", "Volume Averaged Neutral Density" ) )
      CALL Check( nf90_put_att( ncid, density_varid, "units", "kg m^{-3}" ) )

      CALL Check( nf90_enddef(ncid) )

      CALL Check( nf90_put_var( ncid, x_varid, lon ) )
      CALL Check( nf90_put_var( ncid, y_varid, lat ) )
      CALL Check( nf90_put_var( ncid, z_varid, z_regridded ) )
      CALL Check( nf90_put_var( ncid, time_varid, time ) )
      CALL Check( nf90_put_var( ncid, temperature_varid, regridded_temperature, recStart, recCount ) )
      CALL Check( nf90_put_var( ncid, density_varid, regridded_density, recStart, recCount ) )

      CALL Check( nf90_close( ncid ) )



  END SUBROUTINE Write_NetCDF

  SUBROUTINE Interpolate( z_pressure, z_regridded, f_pressure, f_regridded )
    REAL(8), INTENT(in)    :: f_pressure(1:n_pressure_levels)
    REAL(8), INTENT(in)    :: z_pressure(1:n_pressure_levels)
    REAL(8), INTENT(in)    :: z_regridded(1:n_regridded_levels)
    REAL(8), INTENT(inout) :: f_regridded(1:n_regridded_levels)
    ! Local
    REAL(8) :: interp(1:n_pressure_levels)
    INTEGER :: row, col
    REAL(8) :: f_new

      DO row = 1, n_regridded_levels

        interp(1:n_pressure_levels) = Hat_Function( z_pressure, z_regridded(row), n_pressure_levels )
        f_new = 0.0D0

        DO col = 1, n_pressure_levels
     
          f_new = f_new + interp(col)*f_pressure(col)

        ENDDO

        f_regridded(row) = f_new

      ENDDO

  END SUBROUTINE Interpolate

  FUNCTION Hat_Function( z_pressure, z_interp, n_pressure_levels ) RESULT( phi )
    INTEGER :: n_pressure_levels
    REAL(8) :: z_pressure(1:n_pressure_levels)
    REAL(8) :: z_interp
    REAL(8) :: phi(1:n_pressure_levels)
    ! Local
    INTEGER :: i

      i = 1
      IF( z_interp >= z_pressure(i) .AND. z_interp < z_pressure(i+1) )THEN
  
        phi(i) = ( z_interp - z_pressure(i+1) )/( z_pressure(i) - z_pressure(i+1) )
  
      ELSE
  
        phi(i) = 0.0D0
  
      ENDIF

      DO i = 2, n_pressure_levels-1

        IF( z_interp > z_pressure(i-1) .AND. z_interp <= z_pressure(i) )THEN
  
          phi(i) = ( z_interp - z_pressure(i-1) )/( z_pressure(i) - z_pressure(i-1) )
  
        ELSEIF( z_interp >= z_pressure(i) .AND. z_interp < z_pressure(i+1) )THEN
  
          phi(i) = ( z_interp - z_pressure(i+1) )/( z_pressure(i) - z_pressure(i+1) )
  
        ELSE
  
          phi(i) = 0.0D0
  
        ENDIF

      ENDDO

      i = n_pressure_levels
      IF( z_interp > z_pressure(i-1) .AND. z_interp <= z_pressure(i) )THEN

        phi(i) = ( z_interp - z_pressure(i-1) )/( z_pressure(i) - z_pressure(i-1) )

      ELSE

        phi(i) = 0.0D0

      ENDIF


  END FUNCTION Hat_Function

  SUBROUTINE Check(status)
    IMPLICIT NONE
    INTEGER, INTENT (in) :: status
    
     IF(status /= nf90_noerr) THEN 
       PRINT *, trim(nf90_strerror(status))
       STOP "NetCDF Error, Stopped"
     ENDIF
  END SUBROUTINE Check 

END PROGRAM WAM_Interpolate
