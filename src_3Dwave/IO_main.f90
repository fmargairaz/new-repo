module mod_IO_main

  !use mod_IO_VTK
  !use mod_IO_HDF5
  use decomp_2d_io
  use physical_system
  use decomp_2d
  
  use parameters_IO

  implicit none
  include 'mpif.h'
    
  private
  public :: io_init,io_finilize,output_loop

contains
  
  subroutine io_init
    !call create_files(vector)
  end subroutine io_init
  
  subroutine io_finilize
    !call close_files
  end subroutine io_finilize

  subroutine output_loop(jt,me,decomp)
    
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: decomp
    integer,intent(in) :: jt,me
    
    !vtk output
    !call vtk_out(jt,me,decomp) 
    
    !HDF5 output
    !call add_snap_data(dble(jt*dt),jt)

    !2decomp output 
    call decomp_out(jt)

  end subroutine output_loop
  
  subroutine decomp_out(jt)
    
    implicit none
    integer,intent(in) :: jt
    character(256) :: filename
    
    write(filename,'(A,I6.6,A)') trim(out_path)//'2decomp_u_',jt,'.dat'
    
    call decomp_2d_write_one(1,u,filename)
    
  end subroutine decomp_out

end module mod_IO_main
