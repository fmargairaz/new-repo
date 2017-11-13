!------------------------------------------------------------------------------
!!  module: checkpoint
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!  - checkpoints
!!
!------------------------------------------------------------------------------

module mod_IO_checkpoint

  use decomp_2d
  use decomp_2d_custom
  use decomp_2d_io

  !! - global variables
  use system_variables, only : u,v,w,rhsx,rhsy,rhsz
  use compute_sgs, only : cs_opt2
  use compute_LagScaleDep, only : F_LM,F_MM,F_QN,F_NN
  use parameters_IO
  
  implicit none
  include 'mpif.h'

  private
  TYPE(DECOMP_INFO),save :: dcp

  public :: checkpoint_io_init,checkpoint_in,checkpoint_out,scheckpoint_out

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: initialise checkpoint IO module
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine checkpoint_io_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main
    dcp = dcp_main

    return
  end subroutine checkpoint_io_init


  subroutine checkpoint_in(jt_total)
    implicit none
    integer,intent(in) :: jt_total
    
    character(128) :: filename
    write(filename,'(A,I9.9)') 'checkpoint_',jt_total
    if(nrank==0)then
       write(*,*) '=============================================================================='
       write(*,'(A,A)') '# reading checkpoint file: ',trim(filename)
       write(*,*) '=============================================================================='
    endif
    call read_checkpoint_file(filename)
    
    !update ghostcell
    call update_ghost(u,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(v,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(w,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(rhsx,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(rhsx,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(rhsx,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_top_ghost(cs_opt2,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    if(Lagran_init==.true.)then 
       call update_top_ghost(F_LM,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_top_ghost(F_MM,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_top_ghost(F_QN,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_top_ghost(F_NN,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    endif
    
    return
  end subroutine checkpoint_in

  subroutine checkpoint_out(jt)
    implicit none
    integer,intent(in) :: jt
    
    character(128) :: filename
    write(filename,'(A,I9.9)') 'checkpoint_',jt
    if(nrank==0)then
       write(*,*) '=============================================================================='
       write(*,'(A,A)') '# writing checkpoint file: ',trim(filename)
       write(*,*) '=============================================================================='
    endif
    call write_checkpoint_file(filename)
    
    return
  end subroutine checkpoint_out

  subroutine scheckpoint_out()
    implicit none
    character(128) :: filename = 'scheckpoint'
    
    call write_checkpoint_file(filename)
    
    return
  end subroutine scheckpoint_out

  !------------------------------------------------------------------------------
  !! subroutine: write checkpoint
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - the checkpoint has ghostcell only at the top of the domain
  !!
  !------------------------------------------------------------------------------
  subroutine write_checkpoint_file(filename)
    
    character(128),intent(IN) :: filename
    
    integer :: ierror,errorcode,fh,newtype
    character(512) :: fileout
    integer (kind=MPI_OFFSET_KIND) :: filesize,disp
    
    integer, dimension(3) :: sizes,subsizes,starts
    
    ! new MPI data type for variable with ghostcell 
    ! at top of the domain
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz+1
    
    subsizes(1) = dcp%xsz(1)
    subsizes(2) = dcp%xsz(2)
    if(dcp%xen(3)==nz)then 
       subsizes(3) = dcp%xsz(3)+1
    else
       subsizes(3) = dcp%xsz(3)
    endif
              
    starts(1) = dcp%xst(1)-1  ! 0-based index
    starts(2) = dcp%xst(2)-1
    starts(3) = dcp%xst(3)-1
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,&
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    
    ! open file for IO checkpoint
    write(fileout,'(A)') trim(out_path)//trim(filename)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileout, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    if(ierror.ne.0)then
       errorcode = 1
       call decomp_2d_abort(errorcode, &
            '[checkpoint error] cannot open checkpoint file')
    end if

    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND

    call write_var_tbg(fh,u,sizes,subsizes,newtype,disp)
    call write_var_tbg(fh,v,sizes,subsizes,newtype,disp)
    call write_var_tbg(fh,w,sizes,subsizes,newtype,disp)
    call write_var_tbg(fh,rhsx,sizes,subsizes,newtype,disp)
    call write_var_tbg(fh,rhsy,sizes,subsizes,newtype,disp)
    call write_var_tbg(fh,rhsz,sizes,subsizes,newtype,disp)
    call write_var_tg(fh,cs_opt2,sizes,subsizes,newtype,disp)
    if(Lagran_init==.true.)then 
       call write_var_tg(fh,F_LM,sizes,subsizes,newtype,disp)
       call write_var_tg(fh,F_MM,sizes,subsizes,newtype,disp)
       call write_var_tg(fh,F_QN,sizes,subsizes,newtype,disp)
       call write_var_tg(fh,F_NN,sizes,subsizes,newtype,disp)
    endif

    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

  end subroutine write_checkpoint_file

  !------------------------------------------------------------------------------
  !! subroutine: read checkpoint
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - the checkpoint has ghostcell only at the top of the domain
  !!
  !------------------------------------------------------------------------------
  subroutine read_checkpoint_file(filename)
    character(128),intent(IN) :: filename

    integer :: ierror,errorcode,fh,newtype
    integer (kind=MPI_OFFSET_KIND) :: disp
    
    integer, dimension(3) :: sizes,subsizes,starts

    ! new MPI data type for variable with ghostcell 
    ! at top of the domain
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz+1
    
    subsizes(1) = dcp%xsz(1)
    subsizes(2) = dcp%xsz(2)
    if(dcp%xen(3)==nz)then 
       subsizes(3) = dcp%xsz(3)+1
    else
       subsizes(3) = dcp%xsz(3)
    endif    
    
    starts(1) = dcp%xst(1)-1  ! 0-based index
    starts(2) = dcp%xst(2)-1
    starts(3) = dcp%xst(3)-1
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,&
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)

    ! read data back in from file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_RDONLY, MPI_INFO_NULL, &
         fh, ierror)

    if(ierror.ne.0)then
       errorcode = 1
       call decomp_2d_abort(errorcode, &
            '[Restart error] cannot open checkpoint file')
    end if

    disp = 0_MPI_OFFSET_KIND

    call read_var_tbg(fh,u,sizes,subsizes,newtype,disp)
    call read_var_tbg(fh,v,sizes,subsizes,newtype,disp)
    call read_var_tbg(fh,w,sizes,subsizes,newtype,disp)
    call read_var_tbg(fh,rhsx,sizes,subsizes,newtype,disp)
    call read_var_tbg(fh,rhsy,sizes,subsizes,newtype,disp)
    call read_var_tbg(fh,rhsz,sizes,subsizes,newtype,disp)
    call read_var_tg(fh,cs_opt2,sizes,subsizes,newtype,disp)
    if(Lagran_init==.true.)then 
       call read_var_tg(fh,F_LM,sizes,subsizes,newtype,disp)
       call read_var_tg(fh,F_MM,sizes,subsizes,newtype,disp)
       call read_var_tg(fh,F_QN,sizes,subsizes,newtype,disp)
       call read_var_tg(fh,F_NN,sizes,subsizes,newtype,disp)
    endif

    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    return
  end subroutine read_checkpoint_file
  
  !------------------------------------------------------------------------------
  !! subroutine: write var with top&bottom ghostcell
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - input variable has top&bottom ghostcell
  !!
  !------------------------------------------------------------------------------
  subroutine write_var_tbg(fh,var,sz,sbsz,newtype,disp)
    implicit none

    integer, intent(IN) :: fh,newtype
    real(mytype), dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1), intent(IN) :: var
    integer, dimension(3),intent(IN) :: sz,sbsz
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,var(1,1,1),sbsz(1)*sbsz(2)*sbsz(3),real_type,MPI_STATUS_IGNORE,ierror)

    ! update displacement for the next write operation
    disp = disp + sz(1)*sz(2)*sz(3)*mytype_bytes

    return
  end subroutine write_var_tbg
  
  !------------------------------------------------------------------------------
  !! subroutine: write var with top ghostcell
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - input variable has only top ghostcell
  !!
  !------------------------------------------------------------------------------
  subroutine write_var_tg(fh,var,sz,sbsz,newtype,disp)
    implicit none

    integer, intent(IN) :: fh,newtype
    real(mytype), dimension(dcp%xsz(1),dcp%xsz(2),1:dcp%xsz(3)+1), intent(IN) :: var
    integer, dimension(3),intent(IN) :: sz,sbsz
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,var(1,1,1),sbsz(1)*sbsz(2)*sbsz(3),real_type,MPI_STATUS_IGNORE,ierror)
    
    ! update displacement for the next write operation
    disp = disp + sz(1)*sz(2)*sz(3)*mytype_bytes

    return
  end subroutine write_var_tg

  !------------------------------------------------------------------------------
  !! subroutine: read var with top&bottom ghostcell
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - input variable has top&bottom ghostcell
  !!
  !------------------------------------------------------------------------------
  subroutine read_var_tbg(fh,var,sz,sbsz,newtype,disp)
    implicit none

    integer, intent(IN) :: fh,newtype
    real(mytype), dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1), intent(IN) :: var
    integer, dimension(3),intent(IN) :: sz,sbsz
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh,var(1,1,1),sbsz(1)*sbsz(2)*sbsz(3),real_type,MPI_STATUS_IGNORE,ierror)
    
    ! update displacement for the next write operation
    disp = disp + sz(1)*sz(2)*sz(3)*mytype_bytes

    return
  end subroutine read_var_tbg

  !------------------------------------------------------------------------------
  !! subroutine: read var with top ghostcell
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - input variable has only top ghostcell
  !!
  !------------------------------------------------------------------------------
  subroutine read_var_tg(fh,var,sz,sbsz,newtype,disp)
    implicit none

    integer, intent(IN) :: fh,newtype
    real(mytype), dimension(dcp%xsz(1),dcp%xsz(2),1:dcp%xsz(3)+1), intent(IN) :: var
    integer, dimension(3),intent(IN) :: sz,sbsz
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh,var(1,1,1),sbsz(1)*sbsz(2)*sbsz(3),real_type,MPI_STATUS_IGNORE,ierror)
    
    ! update displacement for the next write operation
    disp = disp + sz(1)*sz(2)*sz(3)*mytype_bytes

    return
  end subroutine read_var_tg
  
end module mod_IO_checkpoint
