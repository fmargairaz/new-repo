!------------------------------------------------------------------------------
!!  module: decomp_2d_fft_2d
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
!!
!! - description: 
!!  - this code uese the FFTW (version 3.x)
!!
!------------------------------------------------------------------------------

module decomp_2d_fft_2d
  
  use decomp_2d  ! 2D decomposition module
  
  implicit none
  
  include "fftw3.f"
  
  ! Make everything private unless declared public
  private        
  
  ! engine-specific global variables
  ! integer, save :: plan_type = FFTW_MEASURE
  integer, save :: plan_type = FFTW_ESTIMATE

  ! FFTW plans
  integer*8, save :: plan_fx,plan_bx,plan_fy,plan_by
  integer*8, save :: plan_fxb,plan_bxb,plan_fyb,plan_byb

  ! This file contains common code shared by all FFT engines
  integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
  integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1
      
  ! The libary can only be initialised once
  logical, save :: initialised = .false. 

  ! Global size of the FFT
  integer, save :: nx_fft,ny_fft,nx_fftb,ny_fftb,nz_fft

  ! 2D processor grid
  integer, save, dimension(2) :: dims

  ! Decomposition objects
  TYPE(DECOMP_INFO), save :: ph  ! physical space
  TYPE(DECOMP_INFO), save :: sp,spg  ! spectral space

  TYPE(DECOMP_INFO), save :: phb  ! physical space
  TYPE(DECOMP_INFO), save :: spb,spgb  ! spectral space

  ! Workspace to store the intermediate Y-pencil data
  complex(mytype), allocatable, dimension(:,:,:) :: wk_4x,wk_4xb

  public :: decomp_2d_fft_2d_init,decomp_2d_fft_2d_finalize, &
       execute_decomp_2d_fft_2d,execute_decomp_2d_fft_2d_big, &
       decomp_2d_fft_2d_get_size,decomp_2d_fft_2d_get_sp_info

  ! Declare generic interfaces to handle different inputs
  interface decomp_2d_fft_2d_init
     module procedure fft_2d_init_noarg
     module procedure fft_2d_init_general
  end interface
  
  interface execute_decomp_2d_fft_2d
     module procedure fft_2d_r2c
     module procedure fft_2d_c2r
  end interface

  interface execute_decomp_2d_fft_2d_big
     module procedure fft_2d_r2c_big
     module procedure fft_2d_c2r_big
  end interface
  
contains !=======================================================================   

  !------------------------------------------------------------------------------
  !! subroutine: initialise fft_2d module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------

  subroutine fft_2d_init_noarg
    
    implicit none
    
    call fft_2d_init_general(nx_global, ny_global, nz_global)
        
    return
  end subroutine fft_2d_init_noarg

  !------------------------------------------------------------------------------
  !! subroutine: initialise fft_2d module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine fft_2d_init_general(nx, ny, nz)

    implicit none

    integer, intent(IN) :: nx,ny,nz

    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer :: status, errorcode, ierror

    if (initialised) then
       errorcode = 4
       call decomp_2d_abort(errorcode,&
            'FFT library should only be initialised once')
    end if
    
    nx_fft = nx
    ny_fft = ny

    nx_fftb = 3*nx/2
    ny_fftb = 3*ny/2

    nz_fft = nz

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, ierror)

    ! for c2r/r2c interface:
    ! if in physical space, a real array is of size: nx*ny*nz
    ! in spectral space, the complex array is of size: (nx/2+1)*ny*nz 
    call decomp_info_init(nx,ny,nz,ph)
    call decomp_info_init(nx/2+1,ny,nz,sp)
    call decomp_info_init(nx/2+1,ny,nz+2*dims(2),spg)

    call decomp_info_init(3*nx/2,3*ny/2,nz,phb)
    call decomp_info_init(3*nx/4+1,3*ny/2,nz,spb)
    call decomp_info_init(3*nx/4+1,3*ny/2,nz+2*dims(2),spgb)
    
    allocate(wk_4x(sp%xsz(1),sp%xsz(2),0:sp%xsz(3)+1), STAT=status)    
    allocate(wk_4xb(spb%xsz(1),spb%xsz(2),0:spb%xsz(3)+1), STAT=status) 

    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode,&
            'Out of memory when initialising FFT')
    end if

    wk_4x = (0._mytype,0._mytype)
    wk_4xb = (0._mytype,0._mytype)

    call init_fft_engine
    
    initialised = .true.
    
    return
  end subroutine fft_2d_init_general

  !------------------------------------------------------------------------------
  !! subroutine: finalise 2d FFT module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine decomp_2d_fft_2d_finalize
    
    implicit none

    call decomp_info_finalize(ph)
    call decomp_info_finalize(sp)
    call decomp_info_finalize(spg)
    deallocate(wk_4x)

    call decomp_info_finalize(phb)
    call decomp_info_finalize(spb)
    call decomp_info_finalize(spgb)
    deallocate(wk_4xb)

    call finalize_fft_engine

    initialised = .false.

    return
  end subroutine decomp_2d_fft_2d_finalize

  !------------------------------------------------------------------------------
  !! subroutine: gives size of 2d FFT variables
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - return the size, starting/ending index of the distributed array 
  !!    whose global size is (nx/2+1)*ny*nz, for defining data structures
  !!    in r2c and c2r interfaces
  !!
  !------------------------------------------------------------------------------
  subroutine decomp_2d_fft_2d_get_size(istart, iend, isize)
    
    implicit none
    integer, dimension(3), intent(OUT) :: istart, iend, isize
    
    istart = sp%zst
    iend   = sp%zen
    isize  = sp%zsz
    
    return
  end subroutine decomp_2d_fft_2d_get_size

  !------------------------------------------------------------------------------
  !! subroutine: gives 2Decomp info for spectral space
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - return the sp
  !!
  !------------------------------------------------------------------------------
  subroutine decomp_2d_fft_2d_get_sp_info(dcp_out)
    implicit none

    TYPE(DECOMP_INFO),intent(OUT) :: dcp_out
    dcp_out = sp

    return
  end subroutine decomp_2d_fft_2d_get_sp_info

  !------------------------------------------------------------------------------
  !! subroutine: build FFT plan in y
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - return a FFTW3 plan for 1D c2c FFTs in Y direction
  !!
  !------------------------------------------------------------------------------
  subroutine c2c_1d_y_plan(plan1, decomp, isign)
    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), allocatable, dimension(:) :: a1

    allocate(a1(decomp%ysz(2)))

    call dfftw_plan_dft_1d(plan1,decomp%ysz(2),a1,a1,isign,plan_type);

    deallocate(a1)

    return
  end subroutine c2c_1d_y_plan

  !------------------------------------------------------------------------------
  !! subroutine: build FFT plan in x
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - return a FFTW3 plan for 1D r2c FFTs in X direction
  !!
  !------------------------------------------------------------------------------
  subroutine r2c_1d_x_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), allocatable, dimension(:) :: a1
    complex(mytype), allocatable, dimension(:) :: a2

    allocate(a1(decomp_ph%xsz(1)))
    allocate(a2(decomp_sp%xsz(1)))

    call dfftw_plan_dft_r2c_1d(plan1,decomp_ph%xsz(1),a1,a2,plan_type)
    
    deallocate(a1,a2)    

    return
  end subroutine r2c_1d_x_plan

  !------------------------------------------------------------------------------
  !! subroutine: build FFT plan in x
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - return a FFTW3 plan for 1D c2r FFTs in X direction
  !!
  !------------------------------------------------------------------------------
  subroutine c2r_1d_x_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), allocatable, dimension(:) :: a1
    real(mytype), allocatable, dimension(:) :: a2

    allocate(a1(decomp_sp%xsz(1)))
    allocate(a2(decomp_ph%xsz(1)))

    call dfftw_plan_dft_c2r_1d(plan1,decomp_ph%xsz(1),a1,a2,plan_type)

    deallocate(a1,a2)

    return
  end subroutine c2r_1d_x_plan

  !------------------------------------------------------------------------------
  !! subroutine: initialise FFT engine
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - this routine performs one-time initialisations for the FFT engine
  !!
  !------------------------------------------------------------------------------
  subroutine init_fft_engine

    implicit none

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the FFTW (version 3.x) engine *****'
       write(*,*) ' '
    end if

    ! For R2C/C2R tranforms
    call r2c_1d_x_plan(plan_fx, ph, sp)
    call c2c_1d_y_plan(plan_fy, sp, FFTW_FORWARD )
    call c2c_1d_y_plan(plan_by, sp, FFTW_BACKWARD)
    call c2r_1d_x_plan(plan_bx, sp, ph)

    ! For R2C/C2R tranforms
    call r2c_1d_x_plan(plan_fxb, phb, spb)
    call c2c_1d_y_plan(plan_fyb, spb, FFTW_FORWARD )
    call c2c_1d_y_plan(plan_byb, spb, FFTW_BACKWARD)
    call c2r_1d_x_plan(plan_bxb, spb, phb)

    return
  end subroutine init_fft_engine
  
  !------------------------------------------------------------------------------
  !! subroutine: finalise FFT engine
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:t
  !!  - this routine performs one-time finalisations for the FFT engine
  !!
  !------------------------------------------------------------------------------
  subroutine finalize_fft_engine

    implicit none

    call dfftw_destroy_plan(plan_fx)
    call dfftw_destroy_plan(plan_bx)
    call dfftw_destroy_plan(plan_fy)
    call dfftw_destroy_plan(plan_by)

    call dfftw_destroy_plan(plan_fxb)
    call dfftw_destroy_plan(plan_bxb)
    call dfftw_destroy_plan(plan_fyb)
    call dfftw_destroy_plan(plan_byb)

    return
  end subroutine finalize_fft_engine

  !------------------------------------------------------------------------------
  !! subroutine: compute multiple 1d fft in y
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - c2c transform, multiple 1D FFTs in y direction
  !!
  !------------------------------------------------------------------------------
  subroutine c2c_1d_y(inout,plan1,decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer*8, intent(IN) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    complex(mytype), dimension(decomp%ysz(2)) :: tmp
    integer :: i,k

    !$OMP PARALLEL DO PRIVATE(TMP)
    do k=1,decomp%ysz(3)
       do i=1,decomp%ysz(1)
          tmp=inout(i,:,k)
          call dfftw_execute_dft(plan1, tmp, tmp)
          inout(i,:,k)=tmp
       enddo
    enddo
    !$OMP END PARALLEL DO

    return
  end subroutine c2c_1d_y

  !------------------------------------------------------------------------------
  !! subroutine: compute multiple 1d fft in x
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - r2c transform, multiple 1D FFTs in x direction
  !!
  !------------------------------------------------------------------------------ 
  subroutine r2c_1d_x(input,output,plan1,decomp)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output
    integer*8, intent(IN) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: k,j

    !$OMP PARALLEL DO
    do k=1,decomp%xsz(3)
       do j=1,decomp%xsz(2) 
          call dfftw_execute_dft_r2c(plan1,input(:,j,k),output(:,j,k))
       enddo
    enddo
    !$OMP END PARALLEL DO
 
    return
  end subroutine r2c_1d_x
  
  !------------------------------------------------------------------------------
  !! subroutine: compute multiple 1d fft in x
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - c2r transform, multiple 1D FFTs in x direction
  !!
  !------------------------------------------------------------------------------ 
  subroutine c2r_1d_x(input,output,plan1,decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output
    integer*8, intent(IN) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: k,j

    !$OMP PARALLEL DO
    do k=1,decomp%xsz(3)
       do j=1,decomp%xsz(2)
          call dfftw_execute_dft_c2r(plan1,input(:,j,k),output(:,j,k))
       enddo
    enddo
    !$OMP END PARALLEL DO
 
    return

  end subroutine c2r_1d_x
  
  !------------------------------------------------------------------------------
  !! subroutine: compute multiple 2d FFT
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - forward 2d FFT based on X-prencil
  !!  - input: X-prencil/physical space
  !!  - output: Y-prencil/specrtal space
  !!
  !------------------------------------------------------------------------------ 
  subroutine fft_2d_r2c(in_r,out_c,GFLAG)
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c
    logical,intent(IN) :: GFLAG

    if(GFLAG)then
       ! ===== 1D FFTs in X =====
       call r2c_1d_x(in_r,wk_4x,plan_fx,spg)
       
       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_x_to_y(wk_4x,out_c,spg)
          call c2c_1d_y(out_c,plan_fy,spg)
       else
          out_c = wk_4x
          call c2c_1d_y(out_c,plan_fy,spg)
       end if
    else   
        ! ===== 1D FFTs in X =====
       call r2c_1d_x(in_r,wk_4x(:,:,1:sp%xsz(3)),plan_fx,sp)
       
       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_x_to_y(wk_4x(:,:,1:sp%xsz(3)),out_c,sp)
          call c2c_1d_y(out_c,plan_fy,sp)
       else
          out_c = wk_4x
          call c2c_1d_y(out_c,plan_fy,sp)
       end if
    endif

    return
  end subroutine fft_2d_r2c

  !------------------------------------------------------------------------------
  !! subroutine: compute multiple 2d FFT
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - backward 2d FFT based on Y-prencil
  !!  - input: Y-prencil/spectral space
  !!  - output: X-prencil/physical space
  !!
  !------------------------------------------------------------------------------ 
  subroutine fft_2d_c2r(in_c,out_r,GFLAG)
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r
    logical,intent(IN) :: GFLAG

    if(GFLAG)then
       ! ===== 1D FFTs in Y =====
       call c2c_1d_y(in_c,plan_by,spg)
       
       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(in_c,wk_4x,spg)
          call c2r_1d_x(wk_4x,out_r,plan_bx,spg)
       else
          call c2r_1d_x(in_c,out_r,plan_bx,spg)
       end if
    else
       ! ===== 1D FFTs in Y =====
       call c2c_1d_y(in_c,plan_by,sp)
       
       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(in_c,wk_4x(:,:,1:sp%xsz(3)),sp)
          call c2r_1d_x(wk_4x(:,:,1:sp%xsz(3)),out_r,plan_bx,sp)
       else
          call c2r_1d_x(in_c,out_r,plan_bx,sp)
       end if
    endif
    return
  end subroutine fft_2d_c2r
  
  !------------------------------------------------------------------------------
  !! subroutine: compute multiple 2d FFT for big
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - forward 2d FFT based on X-prencil
  !!  - input: X-prencil/physical space
  !!  - output: Y-prencil/specrtal space
  !!
  !------------------------------------------------------------------------------ 
  subroutine fft_2d_r2c_big(in_r,out_c,GFLAG)
    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c
    logical,intent(IN) :: GFLAG

    if(GFLAG)then
       ! ===== 1D FFTs in X =====
       call r2c_1d_x(in_r,wk_4xb,plan_fxb,spgb)
       
       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_x_to_y(wk_4xb,out_c,spgb)
          call c2c_1d_y(out_c,plan_fyb,spgb)
       else
          out_c = wk_4xb
          call c2c_1d_y(out_c,plan_fyb,spgb)
       end if
    else   
        ! ===== 1D FFTs in X =====
       call r2c_1d_x(in_r,wk_4xb(:,:,1:sp%xsz(3)),plan_fxb,spb)
       
       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_x_to_y(wk_4xb(:,:,1:sp%xsz(3)),out_c,spb)
          call c2c_1d_y(out_c,plan_fyb,spb)
       else
          out_c = wk_4xb
          call c2c_1d_y(out_c,plan_fyb,spb)
       end if
    endif

    return
  end subroutine fft_2d_r2c_big

  !------------------------------------------------------------------------------
  !! subroutine: compute multiple 2d FFT for big 
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!  - backward 2d FFT based on Y-prencil
  !!  - input: Y-prencil/spectral space
  !!  - output: X-prencil/physical space
  !!
  !------------------------------------------------------------------------------ 
  subroutine fft_2d_c2r_big(in_c,out_r,GFLAG)
    implicit none
    
    complex(mytype), dimension(:,:,:), intent(INOUT) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r
    logical,intent(IN) :: GFLAG

    if(GFLAG)then
       ! ===== 1D FFTs in Y =====
       call c2c_1d_y(in_c,plan_byb,spgb)
       
       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(in_c,wk_4xb,spgb)
          call c2r_1d_x(wk_4x,out_r,plan_bxb,spgb)
       else
          call c2r_1d_x(in_c,out_r,plan_bxb,spgb)
       end if
    else
       ! ===== 1D FFTs in Y =====
       call c2c_1d_y(in_c,plan_byb,spb)
       
       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(in_c,wk_4xb(:,:,1:sp%xsz(3)),spb)
          call c2r_1d_x(wk_4xb(:,:,1:sp%xsz(3)),out_r,plan_bxb,spb)
       else
          call c2r_1d_x(in_c,out_r,plan_bxb,spb)
       end if
    endif
    return
  end subroutine fft_2d_c2r_big

end module decomp_2d_fft_2d
