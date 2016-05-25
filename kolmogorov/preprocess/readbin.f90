!============================================!
!                                            !
! Reads preprocess.cpl data with Fortran 08  !
!                                            !
!============================================!
! 
! Author: Dr.-Ing. Davide Gatti
! Date  : 24/May/2016
! 
! This Fortran 2008 module declares a globally
! available ARRAY of STRUCTURE (REAL u,v,w,p),
! called FIELD, large enough to contain three
! wall-parallel planes. This can be useful in 
! case FD should be use to compute derivatives
! in the wall-normal direction
!  

MODULE readbin


  USE, intrinsic :: iso_c_binding
  IMPLICIT NONE


  TYPE :: FIELDTYPE
    real(C_DOUBLE) :: u,v,w,p
  END TYPE FIELDTYPE
  TYPE(FIELDTYPE), allocatable, dimension(:,:,:) :: FIELD


  CONTAINS


  SUBROUTINE allocate_slab(nx,nz)
    integer(C_INT), intent(IN) :: nx,nz
    ALLOCATE(FIELD(0:nz-1,0:nx-1,-1:1))
  END SUBROUTINE allocate_slab


  FUNCTION getfilename(i) RESULT(filename)
    integer(C_INT), intent(IN) :: i
    character(len=40) :: filename
    character(len=30) :: istring
    WRITE(istring,*) i
    filename = "Field.bin."//adjustl(istring) ! XXX add variable length string
  END FUNCTION getfilename


  INTEGER FUNCTION openfield(i) RESULT(fid)
    integer(C_INT), intent(IN) :: i
    fid = 1986 ! XXX add check to see if fid already used
    OPEN(UNIT=fid,FILE=trim(getfilename(i)),access="stream",action="read")
  END FUNCTION openfield


  SUBROUTINE read_plane(iy,nx,nz,planeFIELD)
    integer(C_INT), intent(IN) :: iy,nx,nz
    TYPE(FIELDTYPE), intent(OUT) :: planeFIELD(:,:)
    integer(C_SIZE_T) :: disp
    disp=1+8*4*(iy+1)*nx*nz
    READ(1986,POS=disp) planeFIELD
  END SUBROUTINE read_plane


  SUBROUTINE deallocate_slab()
    DEALLOCATE(FIELD)
  END SUBROUTINE deallocate_slab


END MODULE readbin
