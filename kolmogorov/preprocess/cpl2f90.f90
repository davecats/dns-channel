!============================================!
!                                            !
!   Test the readbin.f90 module and more     !
!                                            !
!============================================!
!           
! Author: Dr.-Ing. Davide Gatti
! Date  : 28/Jul/2015
! 

PROGRAM cpl2f90

  USE readbin
  USE, intrinsic :: iso_c_binding
  IMPLICIT NONE

  ! Input data
  integer(C_INT), parameter :: nx=384,ny=192,nz=384
  integer(C_INT), parameter :: nmin=200,nmax=200

  ! Internal variables
  integer(C_INT) :: fid,iFIELD,iY,i
  
  CALL allocate_slab(nx,nz)

  DO iFIELD=nmin,nmax
    fid=openfield(iFIELD)
    DO iY=1,ny-1
      DO i=-1,1
        CALL read_plane(iY+i,nx,nz,FIELD(:,:,i)) 
      END DO
    END DO
    CLOSE(fid)
  END DO

  CALL deallocate_slab()

END PROGRAM cpl2f90
