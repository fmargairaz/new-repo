!------------------------------------------------------------------------------
!! progam: main
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description: 
!!  - Lagre Eddy Simulation for Atmospheric Boundary Layer
!!  - Hybrid parallelisation using MPI/OpenMP
!!
!! - libraries and modules
!!  - FFTW (version 3.x)
!!  - 2Decomp (version ...)
!!
!------------------------------------------------------------------------------
program main
 
  use LES_engine

  implicit none
  include 'mpif.h'

  integer,parameter :: rprec=kind(1.d0)

  integer :: ierr,me,nproc
  real(rprec) :: time_exec

  call MPI_init(ierr)

  call MPI_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, me, ierr)

  if(me==0)then
     write(*,*) '******************************************************************************'            
     write(*,*) '*                                                                            *'
     write(*,*) '*                                                                            *'
     write(*,*) '*                                                                            *'
     write(*,*) '*                                                                            *'
     write(*,*) '*                                                                            *'
     write(*,*) '*                                                                            *'
     write(*,*) '*             Fabien Magrairaz, University of Utah, 2016                     *'
     write(*,*) '*                                                                            *'
     write(*,*) '******************************************************************************'
     write(*,*) '                                                                              '
  endif

  time_exec = mpi_wtime()
  
  call LES_main_engine()

  time_exec = mpi_wtime()-time_exec

  !call MPI_barrier(MPI_COMM_WORLD, ierr)
  
  if(me==0)then
     write(*,'(A,F10.4,A)') ' wall time of exection:',time_exec,'s' 
  endif

  call MPI_finalize(ierr)

end program main
