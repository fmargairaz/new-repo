module physical_system

  use parameters_IO
  use decomp_2d
  
  implicit none

  real(rprec), allocatable, dimension(:,:,:) :: u, unew
  real(rprec), allocatable, dimension(:,:,:) :: rhs 

  TYPE(DECOMP_INFO) :: vector
  TYPE(DECOMP_INFO) :: scalar
    
  private

  public :: create_system, clean_system
  public :: u, unew, rhs
  public :: vector, scalar

contains
  subroutine create_system(me,nproc,dims,ierr)
    
    implicit none
    integer, intent(inout) :: ierr, me, nproc
    integer, dimension(2), intent(out) :: dims

    
    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer :: status, errorcode, ierror

    ! initialise 2d decomp
    call decomp_2d_init(nx,ny,nz,p_on_y,p_on_z)
    
    call decomp_info_init(nx, ny, nz, vector) 
    call decomp_info_init(nx, ny, nz, scalar) 

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, ierror)

    if(me == 0)then
       write(*,*) 'Creating system...'
       write(*,*) ' '
       write(*,*) '============================================================'
       write(*,*) '(nx,ny,nt)',nx,ny,nz
       write(*,*) '2DECOMP on a 2D grid (X-pencil)'
       write(*,*) '-> (p_on_y,p_on_z)',dims
       write(*,*) '============================================================'
       write(*,*) ' '
    end if

    ! allocate vector quantities
    call alloc_x(u, vector, .true.) ! in X-pencil
    call alloc_x(unew, vector,.true.)
    call alloc_x(rhs, vector, .true.)
    
  end subroutine create_system

  subroutine clean_system(me,nproc,ierr)
    
    implicit none
    integer, intent(inout) :: ierr, me, nproc

    deallocate(u)
    deallocate(unew)
    deallocate(rhs)
    
    call decomp_info_finalize(vector)
    call decomp_info_finalize(scalar)
    
    call decomp_2d_finalize

  end subroutine clean_system
  
end module physical_system
