!------------------------------------------------------------------------------
!!  module: compute derivatives (gradients and divergences)
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------

module compute_derivatives

  use parameters_IO
  use fft_engine
  use decomp_2d
  use decomp_2d_custom
  use decomp_2d_fft_2d

  implicit none
  include 'mpif.h'

  !! - workspace for y-transposed variables
  real(rprec), allocatable, dimension(:,:,:) :: wk1_4y,wk2_4y
  
  TYPE(DECOMP_INFO),save :: dcp,dcpg,sp

  private
  
  public :: deriv_init, deriv_finalize
  public :: grad_velocity_filt,ddx,ddy,ddz_uv,ddz_w
  public :: div_stress
  
contains !=======================================================================
  
  !------------------------------------------------------------------------------
  !! subroutine: initialise derivatives module
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine deriv_init(dcp_main,dcp_ghost)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main,dcp_ghost
    
    integer :: status, errorcode

    dcp=dcp_main
    dcpg=dcp_ghost
    call decomp_2d_fft_2d_get_sp_info(sp)
    
    allocate(wk1_4y(dcp%ysz(1),dcp%ysz(2),0:dcp%ysz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising derivatives module workspace variables')
    end if

    allocate(wk2_4y(dcp%ysz(1),dcp%ysz(2),0:dcp%ysz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising derivatives module workspace variables')
    end if

    wk1_4y=0._rprec
    wk2_4y=0._rprec

  end subroutine deriv_init

  !------------------------------------------------------------------------------
  !! subroutine: finalise derivatives module
  !------------------------------------------------------------------------------
  subroutine deriv_finalize()
    
    implicit none
    
    deallocate(wk1_4y, wk2_4y)
  
  end subroutine deriv_finalize

  !------------------------------------------------------------------------------
  !! subroutine: compute flistered gradients in x&y
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine grad_velocity_filt(f,dfdx,dfdy)
    implicit none
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(INOUT) :: f,dfdx,dfdy
    complex(rprec),dimension(sp%ysz(1),sp%ysz(2),0:sp%ysz(3)+1) :: f_hat,ikxf,ikyf
    integer :: i,j,k
    real(rprec) :: const
    
    const=1._rprec/(nx*ny)
    call execute_decomp_2d_fft_2d(const*f,f_hat,.true.)

    ! kills the oddballs
    if((lhx.ge.sp%yst(1)).and.(lhx.le.sp%yen(1)))then
       !$omp parallel do
       do k=0,sp%ysz(3)+1
          f_hat(lhx-sp%yst(1)+1,:,k)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif
    if((lhy.ge.sp%yst(2)).and.(lhy.le.sp%yen(2)))then
       !$omp parallel do
       do k=0,sp%ysz(3)+1
          f_hat(:,lhy-sp%yst(2)+1,k)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif
    
    !$omp parallel do
    do k=0,dcp%ysz(3)+1
       do j=1,sp%ysz(2)
          do i=1,sp%ysz(1)
             ikxf(i,j,k)=eye*kkx(sp%yst(1)-1+i,sp%yst(2)-1+j)*f_hat(i,j,k)
             ikyf(i,j,k)=eye*kky(sp%yst(1)-1+i,sp%yst(2)-1+j)*f_hat(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    call execute_decomp_2d_fft_2d(ikxf,dfdx,.true.)
    call execute_decomp_2d_fft_2d(ikyf,dfdy,.true.)
    call execute_decomp_2d_fft_2d(f_hat,f,.true.)

    return
  end subroutine grad_velocity_filt

  !------------------------------------------------------------------------------
  !! subroutine: compute flistered gradients in x&y
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine grad_velocity_filt1d(f,dfdx,dfdy)
    implicit none
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(INOUT) :: f,dfdx,dfdy
    integer :: i,j,k,k_max
    
    !$omp parallel
    !$omp do
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          call compute_filt_dfdx(f(:,j,k),dfdx(:,j,k))
       enddo
    enddo
    !$omp end do

    !$omp barrier
    !$omp single
    call transpose_x_to_y(f,wk1_4y,dcpg)
    !$omp end single
    
    !compute ft in y to compute the derivative
    !$omp do 
    do k=1,dcp%ysz(3)
       do i=1,dcp%ysz(1)
          call compute_filt_dfdy(wk1_4y(i,:,k),wk2_4y(i,:,k))
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    call transpose_y_to_x(wk1_4y,f,dcpg)
    call transpose_y_to_x(wk2_4y,dfdy,dcpg)

    return
  end subroutine  grad_velocity_filt1d

  !------------------------------------------------------------------------------
  !! subroutine: compute filtered dfdx
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_filt_dfdx(f,dfdx)

    implicit none
    
    real(rprec),dimension(nx),intent(inout) :: f,dfdx
    complex(rprec),dimension(lhx) :: scr
    real(rprec) :: const

    !normalization constant
    const=1._rprec/nx
    
    !fft to go in wavenumber space (and compute there the derivatives to have high accuracy)
    call dfftw_execute_dft_r2c(planx,f*const,scr)
    
    !kills oddballs
    scr(lhx)=0._rprec
    
    !inverse transform to get filtred velocity (values back in physical space
    call dfftw_execute_dft_c2r(planx2,scr,f)
    
    !compute coefficients for pseudospectral derivative calculation
    scr=eye*kx*scr
    
    !inverse transform to get pseudospectral derivative (values back in physical space)
    call dfftw_execute_dft_c2r(planx2,scr,dfdx)
    
    return
  end subroutine compute_filt_dfdx
  
  
  !------------------------------------------------------------------------------
  !! subroutine: compute filtered dfdy
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_filt_dfdy(f,dfdy)
    
    implicit none

    real(rprec),dimension(ny),intent(inout) :: f,dfdy
    complex(rprec),dimension(lhy) :: scr
    real(rprec) :: const
    
    !normalization constant
    const=1._rprec/ny
    
    !fft to go in wavenumber space (and compute there the derivatives to have high accuracy)
    call dfftw_execute_dft_r2c(plany,f*const,scr)
    
    !kills oddballs
    scr(lhy)=0._rprec

    !inverse transform to get filtred velocity (values back in physical space)
    call dfftw_execute_dft_c2r(plany2,scr,f)
    
    !compute coefficients for pseudospectral derivative calculation
    scr=eye*ky*scr
    
    !inverse transform to get pseudospectral derivative (values back in physical space)
    call dfftw_execute_dft_c2r(plany2,scr,dfdy)
        
  end subroutine compute_filt_dfdy

  !------------------------------------------------------------------------------
  !! subroutine: compute derivative in x
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - compute the derivatives on a line along x
  !!  - this subroutine performs differentiation along the given direction adopting 
  !!  an out of place transform.
  !!  - remember that the input of c2r is anyway destroyed (scr)
  !!
  !------------------------------------------------------------------------------
  subroutine ddx(dfdx,f)

    implicit none
    
    real(rprec),dimension(:,:,:),intent(in) :: f
    real(rprec),dimension(:,:,:),intent(out) :: dfdx
    complex(rprec),dimension(lhx) :: scr
    real(rprec) :: const
    integer :: j,k,s3

    s3=size(f,3)
    
    ! normalization constant
    const=1._rprec/nx
    
    do k=1,s3
       do j=1,dcp%xsz(2)
          
          ! fft to go in wavenumber space (and compute there the derivatives to have high accuracy)
          call dfftw_execute_dft_r2c(planx,const*f(:,j,k),scr)
          
          ! kills oddballs
          scr(lhx)=0._rprec
          
          ! compute coefficients for pseudospectral derivative calculation
          scr=eye*kx*scr
          
          ! inverse transform to get pseudospectral derivative (values back in physical space)
          call dfftw_execute_dft_c2r(planx2,scr,dfdx(:,j,k))
       
       enddo
    enddo

    return
  end subroutine ddx
  
  
  !------------------------------------------------------------------------------
  !! subroutine: compute derivative in y
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - compute the derivatives on a line along x
  !!  - this subroutine performs differentiation along the given direction adopting 
  !!  an out of place transform.
  !!  - remember that the input of c2r is anyway destroyed (scr)
  !!
  !------------------------------------------------------------------------------
  subroutine ddy(dfdy,f)
    
    implicit none
    
    real(rprec),dimension(:,:,:),intent(in) :: f
    real(rprec),dimension(:,:,:),intent(out) :: dfdy
    real(rprec),dimension(dcp%ysz(1),dcp%ysz(2),dcp%ysz(3)) :: wk
    complex(rprec),dimension(lhy) :: scr
    real(rprec) :: const
    integer :: i,k,s3

    s3=size(f,3)
    
    ! normalization constant
    const=1._rprec/ny

    ! Y-pencil start ------------------------------------------------------------
    call transpose_x_to_y(const*f,wk,dcp)

    do k=1,s3
       do i=1,dcp%ysz(1)
          ! fft to go in wavenumber space (and compute there the derivatives to have high accuracy)
          call dfftw_execute_dft_r2c(plany,wk(i,:,k),scr)
          
          ! kills oddballs
          scr(lhy)=0._rprec
          
          ! compute coefficients for pseudospectral derivative calculation
          scr=eye*ky*scr
          
          ! inverse transform to get pseudospectral derivative (values back in physical space)
          call dfftw_execute_dft_c2r(plany2,scr,wk1_4y(i,:,k))
       enddo
    enddo

    call transpose_y_to_x(wk,dfdy,dcp)
    ! Y-pencil end --------------------------------------------------------------

    return
  end subroutine ddy

  
  !------------------------------------------------------------------------------
  !! subroutine: compute derivative in z for uvp nodes
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - input on uvp-nodes output on w-nodes
  !!
  !------------------------------------------------------------------------------
  subroutine ddz_uv(dfdz,f)

    implicit none
             
    !f and dfdx are real arrays as they don't need to go through the fft
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(in) :: f
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(out) :: dfdz
    
    integer :: i,j,k
    real(rprec) :: const

    ! vertical distance needed to calculate derivatives
    const=1._rprec/dz

    !$omp parallel do
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             dfdz(i,j,k)=const*(f(i,j,k)-f(i,j,k-1))
          enddo
       enddo
    enddo
    !$omp end parallel do

    return 
  end subroutine ddz_uv

  
  !------------------------------------------------------------------------------
  !! subroutine: compute derivative in z for w nodes
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - input on w-nodes output on uvp-nodes
  !!
  !------------------------------------------------------------------------------
  subroutine ddz_w(dfdz,f)
    
    implicit none
        
    ! f and dfdx are real arrays as they don't need to go through the fft
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(in) :: f
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(out) :: dfdz

    real(rprec) :: const
    integer :: i,j,k

    const=1._rprec/dz

    !$omp parallel do
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             dfdz(i,j,k)=const*(f(i,j,k+1)-f(i,j,k))
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    if(dcp%xen(3)==nz)then
       dfdz(:,:,dcp%xsz(3)+1)=0._rprec
    endif

    return
  end subroutine ddz_w

  !--------------------------------------------------------------------------------
  !! subroutine: compute divergence of stress
  !--------------------------------------------------------------------------------
  !!
  !! - marco giometto ( mgiometto@gmail.com )
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - add the divergence of stress to rhs
  !!  - provides divt, jz=1:nz-1
  !!  - stag=.true. for uvp nodes
  !!  - stag=.false. for w nodes
  !!
  !! - notes: 
  !!  - this way could parallelize this 3 derivative calculation at the
  !!    price of a little more memory, alternative use a scr for dtxdx
  !!    etc and update directly the rhs
  !!  - we compute the stress divergence from 1 to nz-1 since for uvp
  !!    the nz level is outside and since we prescribe a bc for w(nz)
  !!
  !--------------------------------------------------------------------------------
  subroutine div_stress(rhs,tx,ty,tz,stag,bottom_wk_press,top_wk_press)

    implicit none
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(inout) :: rhs
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(in) :: tx,ty,tz
    logical,intent(in) :: stag
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(out),optional :: bottom_wk_press,top_wk_press

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1) :: dtxdx,dtydy,dtzdz
    complex(rprec),dimension(sp%ysz(1),sp%ysz(2),0:sp%ysz(3)+1) :: f_hat,ikf

    integer :: i,j,k,k_min
    real(rprec) :: const
    integer ierror

    const=1._rprec/(nx*ny)

    ! derivative in x -----------------------------------------------------
    !call execute_decomp_2d_fft_2d(const*tx(:,:,1:dcp%xsz(3)),f_hat,.false.)
    call execute_decomp_2d_fft_2d(const*tx,f_hat,.true.)
    ! kills the oddballs
    if((lhx.ge.sp%yst(1)).and.(lhx.le.sp%yen(1)))then
       !$omp parallel do
       do k=0,sp%ysz(3)+1
          f_hat(lhx-sp%yst(1)+1,:,k)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif
    if((lhy.ge.sp%yst(2)).and.(lhy.le.sp%yen(2)))then
       !$omp parallel do
       do k=0,sp%ysz(3)+1
          f_hat(:,lhy-sp%yst(2)+1,:)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif
    
    !$omp parallel do
    do k=0,dcp%ysz(3)+1
       do j=1,sp%ysz(2)
          do i=1,sp%ysz(1)
             ikf(i,j,k)=eye*kkx(sp%yst(1)-1+i,sp%yst(2)-1+j)*f_hat(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    call execute_decomp_2d_fft_2d(ikf,dtxdx,.true.)

    ! derivative in y -----------------------------------------------------
    !call execute_decomp_2d_fft_2d(const*ty(:,:,1:dcp%xsz(3)),f_hat,.false.)
    call execute_decomp_2d_fft_2d(const*ty,f_hat,.true.)

    ! kills the oddballs
    if((lhx.ge.sp%yst(1)).and.(lhx.le.sp%yen(1)))then
       !$omp parallel do
       do k=0,sp%ysz(3)+1
          f_hat(lhx-sp%yst(1)+1,:,k)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif
    if((lhy.ge.sp%yst(2)).and.(lhy.le.sp%yen(2)))then
       !$omp parallel do
       do k=0,sp%ysz(3)+1
          f_hat(:,lhy-sp%yst(2)+1,:)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif

    !$omp parallel do
    do k=0,dcp%ysz(3)+1
       do j=1,sp%ysz(2)
          do i=1,sp%ysz(1)
             ikf(i,j,k)=eye*kky(sp%yst(1)-1+i,sp%yst(2)-1+j)*f_hat(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    call execute_decomp_2d_fft_2d(ikf,dtydy,.true.)
    
    ! derivative in z ----------------------------------------------------- 
    if(stag)then
       call ddz_w(dtzdz,tz)
    else
       call ddz_uv(dtzdz,tz)
    endif

    ! compute the divergence ----------------------------------------------
    !$omp parallel do
    do k=1,dcp%xsz(3)
       rhs(:,:,k) = rhs(:,:,k)+dtxdx(:,:,k)+dtydy(:,:,k)+dtzdz(:,:,k)
    enddo
    !$omp end parallel do
    

    ! save dtztz for pressure solver
    if(present(bottom_wk_press))then
       if((.not.stag).and.(dcp%xst(3)==1))then
          bottom_wk_press(:,:) = dtxdx(:,:,1)+dtydy(:,:,1)+dtzdz(:,:,1)
       else
          bottom_wk_press(:,:) = BOGUS
       endif
    endif
    
    if(present(top_wk_press))then  
       if((.not.stag).and.(dcp%xen(3)==nz))then
          top_wk_press(:,:) = dtxdx(:,:,dcp%xsz(3)+1)+dtydy(:,:,dcp%xsz(3)+1)
       else
          top_wk_press(:,:) = BOGUS
       endif
    endif

    return
  end subroutine div_stress
  
end module compute_derivatives

