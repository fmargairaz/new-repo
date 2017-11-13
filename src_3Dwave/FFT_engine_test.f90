module FFT_engine

  use parameters_IO
  use physical_system
  
  use decomp_2d
  use bc_mod
  use fft_mod
  use derivatives_mod
  use rhs_mod

  use mod_IO_main

  implicit none
  include 'mpif.h'

  private
  public :: testFFT
  
contains
  
  subroutine testFFT()
    integer :: ierr,me,nproc
    integer, dimension(2) :: dims
    
    integer :: i,j,k
    integer :: jt

    real(rprec) :: time_exec, time_start 
    real(rprec) :: time_x, time_y, time_z, time_out
    
    call MPI_comm_rank(MPI_COMM_WORLD, me, ierr)
    call MPI_comm_size(MPI_COMM_WORLD, nproc, ierr) 

    call read_parameters(me,nproc,ierr)

    !initiate
    call create_system(me,nproc,dims,ierr)
    call fft_init
    if(output)then
       call io_init
    endif
         
    !initial conditions 
    do i=vector%xst(1),vector%xen(1)
       u(i,:,:)=dsin(2.0_rprec*pi*i*dx/lx)
    enddo

    if(output)then
       call output_loop(0,me,vector)
    endif

    time_exec = mpi_wtime()
    time_x = 0.0
    time_y = 0.0
    time_z = 0.0
    time_out = 0.0
        
    !time loop
    do jt=1,nt

       rhs = 0
       
       call MPI_barrier(MPI_COMM_WORLD, ierr)

       !boundary conditions 
       if(vector%xst(3).eq.1)then 
          call set_bc(0,u,'drc',0.0_rprec,1,vector)
       endif
       if(vector%xen(3).eq.nz)then 
          call set_bc(1,u,'neu',0.0_rprec,1,vector)
       endif

       !needed ???
       call MPI_barrier(MPI_COMM_WORLD, ierr)

       !x,y differentiation (spectral) 
       time_start = mpi_wtime()
       call add2rhs_ddx(u,rhs,vector)
       time_x = time_x + (mpi_wtime() - time_start)

       time_start = mpi_wtime()
       if(dims(1).eq.1 .and. dims(2).eq.1)then
          call add2rhs_ddy_notrans(u,rhs,vector)
       else
          call add2rhs_ddy(u,rhs,vector)
       end if
       time_y = time_y + (mpi_wtime() - time_start)
       
       !z differentiation (fdm)
       time_start = mpi_wtime()
       call add2rhs_ddz(u,rhs,vector)
       time_z = time_z + (mpi_wtime() - time_start)

       !advance in time
       unew=u-a*dt*rhs
       
       !update u
       u=unew

       time_start = mpi_wtime()
       !std output
       if(output .and. mod(jt,step_out).eq.0)then
          call output_loop(jt,me,vector) 
       endif
       time_out = time_out + (mpi_wtime() - time_start)
    
    enddo

    time_exec = mpi_wtime() - time_exec

    if(me==0)then
       open(99,file=out_path(1:len_trim(out_path))//'timer')
       write(99,'(I6.6)') p_on_y, p_on_z, nx, ny, nz, nt
       write(99,'(F15.7)') time_x, time_y, time_z, time_out
       write(99,'(F15.7)') time_exec, time_exec - (time_x+time_y+time_z+time_out)
       write(*,*) 'Timers...'
       write(*,*) ' '
       write(*,*) '============================================================'
       write(*,*) 'time x | time y | time z'
       write(*,*) time_x, time_y, time_z
       write(*,*) 'time io | time loop | time rest'
       write(*,*) time_out, time_exec, time_exec - (time_x+time_y+time_z+time_out) 
       write(*,*) '============================================================'
       write(*,*) ' '
    endif

    if(output)then
       call io_finilize
    endif

    call clean_system(me,nproc,ierr)
	
  end subroutine testFFT
  
end module FFT_engine
