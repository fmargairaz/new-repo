module parameters_IO

  implicit none
  include 'mpif.h'

  integer,parameter :: rprec=kind(1.d0)
  real(rprec),parameter :: pi=3.1415926535897932384626433_rprec

  integer :: p_on_y,p_on_z

  integer :: nx,ny,nz
  integer :: nt
  integer :: lhx,lhy,ldx,ldy

  real(rprec) :: dx,dy,dz
  real(rprec) :: dt

  real(rprec) :: a
  real(rprec) :: lx,ly,lz 

  logical :: output
  integer :: step_out
  character(256) :: out_path

  private
  public :: read_parameters
  
  public :: rprec,pi
  public :: p_on_y,p_on_z,nx,ny,nz,nt,lhx,lhy,ldx,ldy
  public :: dx,dy,dz,dt,a,lx,ly,lz
  public :: output,step_out,out_path
  
contains

  subroutine read_parameters(me, nproc,ierr)

    implicit none

    integer, intent(inout) :: ierr, me, nproc

    if(me==0)then
       namelist / numerical_param / p_on_y, p_on_z, nx, ny, nz, nt, dt

       namelist / output_param / output, step_out, out_path

       namelist / physical_param / a, lx, ly, lz

       write(*,*) 'Reading parameters...'
       write(*,*) ' '
       write(*,*) '============================================================'
       read(*,numerical_param)
       write(*,numerical_param)
       write(*,*) '============================================================'
       read(*,output_param)
       write(*,output_param)
       write(*,*) '============================================================'
       read(*,physical_param)
       write(*,physical_param)
       write(*,*) '============================================================'
       write(*,*) ' '

       lhx=nx/2+1
       lhy=ny/2+1
       ldx=2*lhx
       ldy=2*lhy
       
       dx=lx/nx
       dy=ly/ny
       dz=lz/(nz-1)
       
    endif

    call broadcast_parameters(me, nproc, ierr)

    ! needed so all param are bcast
    call MPI_barrier(MPI_COMM_WORLD, ierr)

  end subroutine read_parameters

  subroutine broadcast_parameters(me, nproc,ierr)

    implicit none
    integer, intent (inout) :: ierr, me, nproc
    
    integer,parameter :: bufferlength=1000
    integer :: position
    character :: buffer(bufferlength)

    if(me==0)then
       !! on node 0
       position = 0
       
       call MPI_pack(p_on_y, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(p_on_z, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       
       call MPI_pack(nx, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ny, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(nz, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(nt, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       
       call MPI_pack(lhx, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(lhy, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ldx, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ldy, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       
       call MPI_pack(dx, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(dy, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(dz, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(dt, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)

       
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

    else
       !! on the other nodes
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

       position = 0
       call MPI_unpack(buffer, bufferlength, position, p_on_y, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, p_on_z, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

       call MPI_unpack(buffer, bufferlength, position, nx, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ny, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, nz, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, nt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        
       call MPI_unpack(buffer, bufferlength, position, lhx, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, lhy, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ldx, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ldy, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

       call MPI_unpack(buffer, bufferlength, position, dx, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, dy, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, dz, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, dt, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       
    endif

    if(me==0)then
       !! on node 0
       position = 0
       call MPI_pack(output, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(step_out, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(out_path, 256, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

    else
       !! on the other nodes
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

       position = 0
       call MPI_unpack(buffer, bufferlength, position, output, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, step_out, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, out_path, 256, MPI_CHARACTER, MPI_COMM_WORLD, ierr)
    endif

    
    if(me==0)then
       !! on node 0
       position = 0
       call MPI_pack(a, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(lx, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ly, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(lz, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

    else
       !! on the other nodes
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

       position = 0
       call MPI_unpack(buffer, bufferlength, position, a, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, lx, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ly, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, lz, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    endif

  end subroutine broadcast_parameters

end module parameters_IO
