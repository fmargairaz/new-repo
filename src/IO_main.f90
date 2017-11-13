module mod_IO_main

  use decomp_2d
  use decomp_2d_io
  use system_variables, only : u,v,w,p,&
       rhsx,rhsy,rhsz,&
       dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,&
       txx,txy,txz,tyy,tyz,tzz
  use compute_sgs, only : nu_t,cs_opt2

  use parameters_IO

  implicit none
  include 'mpif.h'

  private
  TYPE(DECOMP_INFO),save :: dcp

  public :: io_init,io_finalize
  public :: output_loop

contains !=======================================================================

  subroutine io_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main
    
    return
  end subroutine io_init
  
  subroutine io_finalize()
    implicit none
        
    return
  end subroutine io_finalize
  
  subroutine output_loop(jt,jt_total,me)

    implicit none

    integer,intent(in) :: jt,jt_total,me

    if(nrank==0)then
       write(*,*) 'output at jt=',jt
    endif

    !vtk output
    !call vtk_out(jt,me,decomp) 

    !HDF5 output
    !call add_snap_data(dble(jt*dt),jt)

    !2decomp output 
    call decomp_out(jt,jt_total)

  end subroutine output_loop

  subroutine decomp_out(jt,jt_total)

    implicit none
    integer,intent(in) :: jt,jt_total
    character(512) :: filename

    integer :: ierror, fh
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  
    ! open file for IO system
    write(filename,'(A,I9.9,A)') trim(out_path)//'2decomp_sys_',jt_total,'.dat'
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    
    call decomp_2d_write_var(fh,disp,1,u(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,v(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,w(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,p(:,:,1:dcp%xsz(3)),dcp)

    call decomp_2d_write_var(fh,disp,1,rhsx(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,rhsy(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,rhsz(:,:,1:dcp%xsz(3)),dcp)

    call decomp_2d_write_var(fh,disp,1,dudx(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,dudy(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,dudz(:,:,1:dcp%xsz(3)),dcp)
        
    call decomp_2d_write_var(fh,disp,1,dvdx(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,dvdy(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,dvdz(:,:,1:dcp%xsz(3)),dcp)
        
    call decomp_2d_write_var(fh,disp,1,dwdx(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,dwdy(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,dwdz(:,:,1:dcp%xsz(3)),dcp)
        
    call decomp_2d_write_var(fh,disp,1,txx(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,txy(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,txz(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,tyy(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,tyz(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,tzz(:,:,1:dcp%xsz(3)),dcp)
        
    call decomp_2d_write_var(fh,disp,1,nu_t(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,cs_opt2(:,:,1:dcp%xsz(3)),dcp)
        
    call MPI_FILE_CLOSE(fh,ierror)

  end subroutine decomp_out

end module mod_IO_main
