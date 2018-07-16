PROGRAM WAM_Interpolate


IMPLICIT NONE

  ! /////////////////////////////////////////////////// !

  INTEGER, PARAMETER :: n_regridded_levels = 60
  REAL(8), PARAMETER :: z0_regridded       = 0.0D0
  REAL(8), PARAMETER :: zM_regridded       = 600.0D0

  ! /////////////////////////////////////////////////// !


  INTEGER :: n_lat, n_lon, n_pressure_levels



CONTAINS

  SUBROUTINE Setup( )

    ! Should 
    ! (1) Read any command line arguments
    ! (2) Load in a netcdf file
    ! (3) Set n_pressure_levels, n_lat, n_lon
    ! (4) Set up the vertical "fixed height grid"
    !
   

  END SUBROUTINE Setup

  SUBROUTINE Cleanup( )
    ! Clear any allocatable memory, and make sure files are closed
    
  END SUBROUTINE Cleanup

  SUBROUTINE Read_NetCDF( )
  END SUBROUTINE Read_NetCDF

  SUBROUTINE Write_NetCDF( )
  END SUBROUTINE Write_NetCDF

  SUBROUTINE Generate_Interpolation_Matrix( z_pressure, z_regridded, interp_matrix )

    REAL(8), INTENT(in)  :: z_pressure(0:n_pressure_levels)
    REAL(8), INTENT(in)  :: z_regridded(0:n_regridded_levels)
    REAL(8), INTENT(out) :: interp_matrix(0:n_regridded_levels,0:n_pressure_levels) !! Or you could do this as a tridiagonal storage...


  END SUBROUTINE Generate_Interpolation_Matrix

  SUBROUTINE Apply_Interpolation_Matrix( f_pressure, f_regridded, interp_matrix )
    REAL(8), INTENT(in)  :: f_pressure(0:n_pressure_levels)
    REAL(8), INTENT(out) :: f_regridded(0:n_pressure_levels)
    REAL(8), INTENT(in)  :: interp_matrix(0:n_regridded_levels,0:n_pressure_levels) !! Or you could do this as a tridiagonal storage...



  END SUBROUTINE Apply_Interpolation_Matrix

  REAL(8) FUNCTION Hat_Function( i, z_pressure, z_interp )

    INTEGER :: i
    REAL(8) :: z_pressure(0:n_pressure_levels)


  END FUNCTION Hat_Function


END PROGRAM WAM_Interpolate
