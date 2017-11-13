!------------------------------------------------------------------------------
!!  module: treat the aliasing errors
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------

module dealiasing_module

  use parameters_IO
  use fft_engine  
  use decomp_2d
  use decomp_2d_custom
  use utility_tools

  implicit none

  private
  TYPE(DECOMP_INFO),save :: dcp,dcpg
  TYPE(DECOMP_INFO),save :: dcp_big
  TYPE(DECOMP_INFO),save :: dcp_hyb
  
  real(rprec),allocatable,dimension(:,:,:) :: wk_hyb_x,wk_hyb_y
  real(rprec),allocatable,dimension(:,:,:) :: wk1_big
  
  public :: dealiasing_module_init,dealiasing_module_finalize
  public :: dealias_forward,dealias_backward
contains !=======================================================================

!------------------------------------------------------------------------------
  !! subroutine: initialise scalars module
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine dealiasing_module_init(dcp_main,dcp_ghost)
    implicit none
    
        TYPE(DECOMP_INFO),intent(in) :: dcp_main,dcp_ghost
    integer :: status, errorcode

    ! import dcp for main MPI engine
    dcp = dcp_main
    dcpg = dcp_ghost

    call decomp_info_init(nx2,ny2,nz,dcp_big)
    call decomp_info_init(nx2,ny,nz+2*p_on_z,dcp_hyb)
    
    ! workspaces for (3/2nx,ny,nz) variables
    allocate(wk_hyb_x(dcp_hyb%xsz(1),dcp_hyb%xsz(2),0:dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising convective module workspace variables')
    endif

    allocate(wk_hyb_y(dcp_hyb%ysz(1),dcp_hyb%ysz(2),0:dcp%ysz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising convective module workspace variables')
    endif
    
    ! workspace for (3/2nx,3/2ny,nz) variables on 
    allocate(wk1_big(dcp_big%ysz(1),dcp_big%ysz(2),0:dcp_big%ysz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising convective module workspace variables')
    endif

    return
  end subroutine dealiasing_module_init
  
  !------------------------------------------------------------------------------
  !! subroutine: finalise scalars module
  !------------------------------------------------------------------------------=
  subroutine dealiasing_module_finalize()
    implicit none

    deallocate(wk_hyb_x,wk_hyb_y)
    deallocate(wk1_big)

    call decomp_info_finalize(dcp_big)
    call decomp_info_finalize(dcp_hyb)
    
    return
  end subroutine dealiasing_module_finalize

  !------------------------------------------------------------------------------
  !! subroutine: transpose data to big Y-pencil
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz  fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - both cc and cc_big do need ghost cells support
  !!
  !------------------------------------------------------------------------------
  subroutine dealias_forward(cc_big,cc) 
    implicit none
 
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(in) :: cc
    real(rprec),dimension(dcp_big%ysz(1),dcp_big%ysz(2),0:dcp_big%ysz(3)+1),intent(inout) :: cc_big

    complex(rprec),allocatable,dimension(:) :: scr,scr_big
    real(rprec) :: const
    integer :: i,j,k
    
    const=1.0_rprec/(nx2*ny2)

    !$omp parallel private(scr,scr_big) 
    allocate(scr(nx/2+1))
    allocate(scr_big(nx2/2+1))

    !$omp do
    do k=0,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          call dfftw_execute_dft_r2c(planx,cc(:,j,k),scr)
          scr_big=(0.0_rprec,0.0_rprec)
          scr_big(1:lhx-1)=scr(1:lhx-1)
          call dfftw_execute_dft_c2r(planx2_big,scr_big,wk_hyb_x(:,j,k))
       enddo
    enddo
    !$omp end do

    deallocate(scr,scr_big)

    !$omp barrier
    !$omp single
    call transpose_x_to_y(wk_hyb_x,wk_hyb_y,dcp_hyb)
    !$omp end single

    allocate(scr(ny/2+1))
    allocate(scr_big(ny2/2+1))
    
    !$omp do
    do k=0,dcp_big%ysz(3)+1
       do i=1,dcp_big%ysz(1)
          call dfftw_execute_dft_r2c(plany,wk_hyb_y(i,:,k),scr)
          scr_big=(0.0_rprec,0.0_rprec)
          scr_big(1:lhy-1)=scr(1:lhy-1)
          call dfftw_execute_dft_c2r(plany2_big,scr_big,cc_big(i,:,k))
       enddo
    enddo
    !$omp end do

    deallocate(scr,scr_big)
    !$omp end parallel

    cc_big=const*cc_big
    
    return
  end subroutine dealias_forward

  !------------------------------------------------------------------------------
  !! subroutine: dealiase the data
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - cc_big does not need ghost cells support
  !!  - cc does need ghost cells support
  !!
  !------------------------------------------------------------------------------
  subroutine dealias_backward(cc,cc_big)
    implicit none
    
    real(rprec),dimension(dcp_big%ysz(1),dcp_big%ysz(2),0:dcp_big%ysz(3)+1),intent(in) :: cc_big 
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(out) :: cc
    
    complex(rprec),allocatable,dimension(:) :: scr,scr_big
    real(rprec) :: const
    integer :: i,j,k

    const=1.0_rprec/(nx2*ny2)
    
    !$omp parallel private(scr,scr_big)
    allocate(scr(ny/2+1))
    allocate(scr_big(ny2/2+1))
    
    !$omp do
    do k=1,dcp%xsz(3)
       do i=1,dcp_big%ysz(1)
          call dfftw_execute_dft_r2c(plany_big,cc_big(i,:,k),scr_big)
          scr(1:lhy-1)=scr_big(1:lhy-1)
          scr(lhy)=(0.0_rprec,0.0_rprec)
          call dfftw_execute_dft_c2r(plany2,scr,wk_hyb_y(i,:,k))
       enddo
    enddo
    !$omp end do

    deallocate(scr,scr_big)

    !$omp barrier
    !$omp single
    call transpose_y_to_x(wk_hyb_y,wk_hyb_x,dcp_hyb)
    !$omp end single
    
    allocate(scr(nx/2+1))
    allocate(scr_big(nx2/2+1))
    
    !$omp do
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          call dfftw_execute_dft_r2c(planx_big,wk_hyb_x(:,j,k),scr_big)
          scr(1:lhx-1)=scr_big(1:lhx-1)
          scr(lhx)=(0.0_rprec,0.0_rprec)
          call dfftw_execute_dft_c2r(planx2,scr,cc(:,j,k))
       enddo
    enddo
    !$omp end do

    deallocate(scr,scr_big)
    !$omp end parallel

    cc=const*cc
    
    return
  end subroutine dealias_backward

end module dealiasing_module
