!------------------------------------------------------------------------------
!! module: physical_system
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
!!
!! - description:
!!
!------------------------------------------------------------------------------
module system_variables

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use decomp_2d_fft_2d

  implicit none
  include 'mpif.h'

  ! velocity variables
  real(rprec),allocatable,dimension(:,:,:),public :: u, v, w
  real(rprec),allocatable,dimension(:,:,:),public :: dudx,dudy,dudz
  real(rprec),allocatable,dimension(:,:,:),public :: dvdx,dvdy,dvdz
  real(rprec),allocatable,dimension(:,:,:),public :: dwdx,dwdy,dwdz

  ! pressure variables
  real(rprec),allocatable,dimension(:,:,:),public :: p
  real(rprec),allocatable,dimension(:,:,:),public :: dpdx,dpdy,dpdz
  
  ! right hand side of momentum varibales
  real(rprec),allocatable,dimension(:,:,:),public :: rhsx,rhsy,rhsz
  real(rprec),allocatable,dimension(:,:,:),public :: rhsx_f,rhsy_f,rhsz_f

  ! stress tensor variables
  real(rprec),dimension(:,:,:),allocatable,public :: txx,txy,txz,tyy,tyz,tzz
!  real(rprec),dimension(:,:,:),allocatable,public :: divtx,divty,divtz

  ! boundary condtions vaiables
  real(rprec),allocatable,dimension(:,:),public :: lbc_u,lbc_v,lbc_w
  real(rprec),allocatable,dimension(:,:),public :: ubc_u,ubc_v,ubc_w

  ! bottiom and top wkspace for pressure solver
  real(rprec),allocatable,dimension(:,:),public :: bottom_wk_press,top_wk_press
  
  ! 2decomp info
  TYPE(DECOMP_INFO),public :: dcp,dcpg
  TYPE(DECOMP_INFO),public :: scalar

  
  private

  public :: create_system, clean_system

contains !=======================================================================   
  
  !------------------------------------------------------------------------------
  !! subroutine: initialyse the system
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------

  subroutine create_system(me,nproc,dims,ierr)
    
    implicit none
    integer, intent(in) :: ierr, me, nproc
    integer, dimension(2), intent(out) :: dims

    logical, dimension(3) :: periodic_bc
    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer :: status, errorcode, ierror

    ! to add to param
    integer :: ghost_cell=1

    ! initialise 2d decomp
    periodic_bc(1) = .true.
    periodic_bc(2) = .true.
    periodic_bc(3) = .true.

    call decomp_2d_init(nx,ny,nz,p_on_y,p_on_z,periodic_bc)

    call decomp_info_ghost_init(nx,ny,nz,ghost_cell,p_on_z,dcp,dcpg)

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, dims, dummy_periods, dummy_coords, ierror)
    
   if(me == 0)then
       write(*,*) '## Creating system: summary'
       write(*,*) '=============================================================================='       
       write(*,*) '(nx,ny,nz)',nx,ny,nz
       write(*,*) '2DECOMP on a 2D grid (X-pencil)'
       write(*,*) '(p_on_y,p_on_z)',dims
       write(*,*) '=============================================================================='
       write(*,*) ' '       
    end if

    ! allocate velocity variables
    call alloc_x_ghost(u, dcp)
    call alloc_x_ghost(v, dcp)
    call alloc_x_ghost(w, dcp)

    call alloc_x_ghost(dudx, dcp)
    call alloc_x_ghost(dudy, dcp)
    call alloc_x_ghost(dudz, dcp)
    
    call alloc_x_ghost(dvdx, dcp)
    call alloc_x_ghost(dvdy, dcp)
    call alloc_x_ghost(dvdz, dcp)

    call alloc_x_ghost(dwdx, dcp)
    call alloc_x_ghost(dwdy, dcp)
    call alloc_x_ghost(dwdz, dcp)

    ! allocate pressure variables
    call alloc_x_ghost(p, dcp)

    call alloc_x(dpdx, dcp)
    call alloc_x(dpdy, dcp)
    call alloc_x(dpdz, dcp)

    ! allocate rhs variables
    call alloc_x_ghost(rhsx, dcp)
    call alloc_x_ghost(rhsy, dcp)
    call alloc_x_ghost(rhsz, dcp)

    call alloc_x(rhsx_f, dcp)
    call alloc_x(rhsy_f, dcp)
    call alloc_x(rhsz_f, dcp)
    
    ! allocate stress variables
    call alloc_x_ghost(txx, dcp)
    call alloc_x_ghost(txy, dcp)
    call alloc_x_ghost(txz, dcp)
    call alloc_x_ghost(tyy, dcp)
    call alloc_x_ghost(tyz, dcp)
    call alloc_x_ghost(tzz, dcp)

    ! allocate boundary condition variables
    allocate(lbc_u(dcp%xsz(1),dcp%xsz(2)))
    allocate(lbc_v(dcp%xsz(1),dcp%xsz(2)))
    allocate(lbc_w(dcp%xsz(1),dcp%xsz(2)))
    allocate(ubc_u(dcp%xsz(1),dcp%xsz(2)))
    allocate(ubc_v(dcp%xsz(1),dcp%xsz(2)))
    allocate(ubc_w(dcp%xsz(1),dcp%xsz(2)))

    ! allocate bottom and top wkspace for pressure solver
    allocate(bottom_wk_press(dcp%xsz(1),dcp%xsz(2)))
    allocate(top_wk_press(dcp%xsz(1),dcp%xsz(2)))

    ! safe initialisation of all global variables
    u=0.0_rprec
    v=0.0_rprec
    w=0.0_rprec
    P=0.0_rprec

    dudx=0.0_rprec
    dudy=0.0_rprec
    dudz=0.0_rprec

    dvdx=0.0_rprec
    dvdy=0.0_rprec
    dvdz=0.0_rprec

    dwdx=0.0_rprec
    dwdy=0.0_rprec
    dwdz=0.0_rprec
    dpdx=0.0_rprec

    dpdy=0.0_rprec
    dpdz=0.0_rprec

    rhsx=0.0_rprec
    rhsy=0.0_rprec
    rhsz=0.0_rprec

    rhsx_f=0.0_rprec
    rhsy_f=0.0_rprec
    rhsz_f=0.0_rprec

    txx=0.0_rprec
    txy=0.0_rprec
    txz=0.0_rprec
    tyy=0.0_rprec
    tyz=0.0_rprec
    tzz=0.0_rprec

    bottom_wk_press=0.0_rprec
    top_wk_press=0.0_rprec

    ! set BOGUS to useless nodes
    if(dcp%xst(3)==1)then
       w(:,:,0)=BOGUS
    endif

    call MPI_barrier(MPI_COMM_WORLD, ierror)

    return
  end subroutine create_system

  !------------------------------------------------------------------------------
  !! subroutine: clean the system
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 2/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine clean_system()
    
    implicit none

    integer :: error,ierror

    ! deallocate velocity variables
    deallocate(u,v,w)
    deallocate(dudx,dudy,dudz)
    deallocate(dvdx,dvdy,dvdz)
    deallocate(dwdx,dwdy,dwdz)

    ! deallocate pressure variables
    deallocate(p)
    deallocate(dpdx,dpdy,dpdz)

    ! deallocate rhs variables
    deallocate(rhsx,rhsy,rhsz)
    deallocate(rhsx_f,rhsy_f,rhsz_f)

    ! deallocate stress variables
    deallocate(txx,txy,txz,tyy,tyz,tzz)

    ! deallocate boundary condition variables
    deallocate(lbc_u,lbc_v,lbc_w)
    deallocate(ubc_u,ubc_v,ubc_w)

    ! deallocate bottom and top wkspace for pressure solver
    deallocate(bottom_wk_press,top_wk_press)

    ! finalize the MPI engine
    call decomp_info_finalize(dcpg)
    call decomp_info_finalize(dcp)
    call decomp_2d_finalize()

    call MPI_barrier(MPI_COMM_WORLD, ierror)

    return
  end subroutine clean_system

end module system_variables
