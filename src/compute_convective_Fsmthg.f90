!------------------------------------------------------------------------------
!! module: compute convective terms
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!  - use the 2/3 dealiasing rule
!!  
!!
!------------------------------------------------------------------------------

module compute_convective_Fsmthg
  
  use parameters_IO
  use fft_engine
  use decomp_2d
  use decomp_2d_custom

  !! - global variables
  use system_variables, only : u,v,w,dudy,dudz,dvdx,dvdz,dwdx,dwdy

  implicit none
  include 'mpif.h'
  
  private
  TYPE(DECOMP_INFO),save :: dcp,dcpg
  TYPE(DECOMP_INFO),save :: dcp_big
  TYPE(DECOMP_INFO),save :: dcp_hyb
  
  real(rprec),allocatable,dimension(:,:,:) :: wk_x,wk_y
  real(rprec),allocatable,dimension(:,:,:) :: wk1
  
  real(rprec),allocatable,dimension(:) :: Fsmthg_x,Fsmthg_y
  integer :: ierr

  public :: conv_init,conv_finalize
  public :: compute_dealias_convec

contains !=======================================================================
  
  !------------------------------------------------------------------------------
  !! subroutine: initiate convective module
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - need dcp_info for 2/3 rule (to allocate 2Decomp workspace)
  !!
  !------------------------------------------------------------------------------
  subroutine conv_init(dcp_main,dcp_ghost)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main,dcp_ghost
    integer :: i,j
    integer :: status, errorcode

    ! import dcp for main MPI engine
    dcp = dcp_main
    dcpg = dcp_ghost
    call decomp_info_init(nx2, ny2, nz, dcp_big)
    call decomp_info_init(nx2, ny, nz+2*p_on_z, dcp_hyb)
    
    ! workspaces for (nx,ny,nz) variables
    allocate(wk_x(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising convective module workspace variables')
    endif

    allocate(wk_y(dcp%ysz(1),dcp%ysz(2),0:dcp%ysz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising convective module workspace variables')
    endif
    
    ! workspace for (nx,ny,nz) variables on 
    allocate(wk1(dcp%ysz(1),dcp%ysz(2),0:dcp%ysz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising convective module workspace variables')
    endif

    ! Fourier smoothing functons for 
    allocate(Fsmthg_x(lhx), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising convective module variables')
    endif
    do i=1,(lhx)
       Fsmthg_x(i) = dexp(-36.0_rprec*(real(i-1,rprec)/real(lhx-1,rprec))**36.0_rprec)
    enddo
    Fsmthg_x(lhx) = 0._rprec
    
    allocate(Fsmthg_y(lhy), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising convective module variables')
    endif
    do j=1,lhy
       Fsmthg_y(j) = dexp(-36.0_rprec*(real(j-1,rprec)/real(lhy-1,rprec))**36.0_rprec)
    enddo
    Fsmthg_y(lhy) = 0._rprec

    return
  end subroutine conv_init

  !------------------------------------------------------------------------------
  !! subroutine: finilise convective module
  !------------------------------------------------------------------------------
  subroutine conv_finalize() 
    implicit none

    deallocate(wk_x,wk_y)
    deallocate(wk1)
    
    deallocate(Fsmthg_x,Fsmthg_y)

    call decomp_info_finalize(dcp_big)
    call decomp_info_finalize(dcp_hyb)

    return
  end subroutine conv_finalize

  !------------------------------------------------------------------------------
  !! subroutine: compute aliased convective terms
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - this routine uses ghost cell to impose BC on u and v at the bottom
  !!
  !------------------------------------------------------------------------------
  subroutine compute_dealias_convec(cx,cy,cz) 
    implicit none
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(inout) :: cx,cy,cz
    real(rprec),dimension(dcp%ysz(1),dcp%ysz(2),0:dcp%ysz(3)+1) :: u1,u2,u3
    real(rprec),dimension(dcp%ysz(1),dcp%ysz(2),0:dcp%ysz(3)+1) :: vort1,vort2,vort3

    real(rprec) :: const
    integer :: i,j,k,k_min,k_max
    integer :: ierror

    if(dcp%xst(3)==1)then
       k_min = 0
    else
       k_min = 1
    endif

    ! compute truncated velocity
    const=1.0_rprec/(nx*ny)

    !$omp parallel do
    do k=0,dcp%xsz(3)+1
       cx(:,:,k) = -const*u(:,:,k)
       cy(:,:,k) = -const*v(:,:,k)
       cz(:,:,k) = -const*w(:,:,k)
    enddo
    !$omp end parallel do
    
    call transpose_to_truncate(u1,cx)
    call transpose_to_truncate(u2,cy)
    call transpose_to_truncate(u3,cz)
    
    ! compute truncated vorticity
    !$omp parallel do
    do k=0,dcp%xsz(3)+1
       cx(:,:,k) = const*(dwdy(:,:,k)-dvdz(:,:,k))
       cy(:,:,k) = const*(dudz(:,:,k)-dwdx(:,:,k))
       cz(:,:,k) = const*(dvdx(:,:,k)-dudy(:,:,k))
    enddo
    !$omp end parallel do    

#ifndef DNS_MOD
    ! special treatment k=1 cause vort1,vort2 is on uvp (vort3 always uvp)
    if((dcp%xst(3)==1).and.(lbc%special.eq.'wall_law'))then
       cx(:,:,0) = 0.0_rprec
       cy(:,:,0) = 0.0_rprec

       cx(:,:,1) = const*(0.5_rprec*(dwdy(:,:,1)+dwdy(:,:,2))-dvdz(:,:,1))      
       cy(:,:,1) = const*(dudz(:,:,1)-0.5_rprec*(dwdx(:,:,1)+dwdx(:,:,2)))
    endif
#endif

    call transpose_to_truncate(vort1,cx)
    call transpose_to_truncate(vort2,cy)
    call transpose_to_truncate(vort3,cz)

    ! update ghost cell
    call update_ghost(u1,dcp%ysz(1),dcp%ysz(2),dcp%ysz(3),dcp)
    call update_ghost(u2,dcp%ysz(1),dcp%ysz(2),dcp%ysz(3),dcp)
    call update_ghost(u3,dcp%ysz(1),dcp%ysz(2),dcp%ysz(3),dcp)
    call update_ghost(vort1,dcp%ysz(1),dcp%ysz(2),dcp%ysz(3),dcp)
    call update_ghost(vort2,dcp%ysz(1),dcp%ysz(2),dcp%ysz(3),dcp)
    call update_ghost(vort3,dcp%ysz(1),dcp%ysz(2),dcp%ysz(3),dcp)

    ! compute convective terms in each direction
    const=1.0_rprec/(nx*ny)

    ! compute x-term ====================================================
    ! uvp-nodes so 1,Nz-1
    if(dcp%xst(3)==1)then
       wk1(:,:,1) = const*(-u2(:,:,1)*vort3(:,:,1) + &
            0.5_rprec*u3(:,:,2)*vort2(:,:,2))
       k_min = 2
    else
       k_min = 1
    endif
    !$omp parallel do
    do k=k_min,dcp%ysz(3)
       do j=1,dcp%ysz(2)
          do i=1,dcp%ysz(1)
             wk1(i,j,k) = const*(-u2(i,j,k)*vort3(i,j,k) + &
                  0.5_rprec*(u3(i,j,k+1)*vort2(i,j,k+1) + &
                  u3(i,j,k)*vort2(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel do

    call dealias(cx,wk1)

    ! compute y-term ====================================================
    ! uvp-nodes so 1,Nz-1
     if(dcp%xst(3)==1)then
       wk1(:,:,1) = const*(u1(:,:,1)*vort3(:,:,1) + &
            0.5_rprec*u3(:,:,2)*(-vort1(:,:,2)))
       k_min = 2
    else
       k_min = 1
    endif
    !$omp parallel do
    do k=k_min,dcp%ysz(3)
       do j=1,dcp%ysz(2)
          do i=1,dcp%ysz(1)
             wk1(i,j,k) = const*(u1(i,j,k)*vort3(i,j,k) + &
                  0.5_rprec*(u3(i,j,k)*(-vort1(i,j,k)) + &
                  u3(i,j,k+1)*(-vort1(i,j,k+1))))
          enddo
       enddo
    enddo
    !$omp end parallel do

    call dealias(cy,wk1)

    ! compute z-term ====================================================
    ! w-nodes so 1,Nz
    if(dcp%xst(3)==1)then
       wk1(:,:,1) = 0.0_rprec
       k_min = 2
    else
       k_min = 1
    endif
    !$omp parallel do
    do k=k_min,dcp%ysz(3)
       do j=1,dcp%ysz(2)
          do i=1,dcp%ysz(1)
             wk1(i,j,k) = const*0.5_rprec*(&
                  (u1(i,j,k)+u1(i,j,k-1))*(-vort2(i,j,k)) + &
                  (u2(i,j,k)+u2(i,j,k-1))*(vort1(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel do

    call dealias(cz,wk1)

    return
  end subroutine compute_dealias_convec
  
  !------------------------------------------------------------------------------
  !! subroutine: transpose data to big Y-pencil
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!  - both cc and cc_big do need ghost cells support
  !!
  !------------------------------------------------------------------------------
  subroutine transpose_to_truncate(ccout,ccin)
    implicit none
 
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(in) :: ccin
    real(rprec),dimension(dcp%ysz(1),dcp%ysz(2),0:dcp%ysz(3)+1),intent(inout) :: ccout

    complex(rprec),allocatable,dimension(:) :: scr

    integer :: i,j,k
    integer :: tt
    
    !$omp parallel private(scr,tt) 
    allocate(scr(nx/2+1))
    tt=ceiling(real(nx)/3._rprec)

    !$omp do
    do k=0,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          call dfftw_execute_dft_r2c(planx,ccin(:,j,k),scr)
          do i=1,lhx
             scr(i)=scr(i)*Fsmthg_x(i)
          enddo
          call dfftw_execute_dft_c2r(planx2,scr,wk_x(:,j,k))
       enddo
    enddo
    !$omp end do

    deallocate(scr)

    !$omp barrier
    !$omp single
    call transpose_x_to_y(wk_x,wk_y,dcpg)
    !$omp end single

    allocate(scr(ny/2+1))
    tt=ceiling(real(ny)/3._rprec)
    
    !$omp do
    do k=0,dcp%ysz(3)+1
       do i=1,dcp%ysz(1)
          call dfftw_execute_dft_r2c(plany,wk_y(i,:,k),scr)
          do j=1,lhy
             scr(j)=scr(j)*Fsmthg_y(j)
          enddo
          call dfftw_execute_dft_c2r(plany2,scr,ccout(i,:,k))
       enddo
    enddo
    !$omp end do

    deallocate(scr)
    !$omp end parallel

    return
  end subroutine transpose_to_truncate

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
  subroutine dealias(ccout,ccin)
    implicit none
    
    real(rprec),dimension(dcp%ysz(1),dcp%ysz(2),0:dcp%ysz(3)+1),intent(in) :: ccin 
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(out) :: ccout

    complex(rprec),allocatable,dimension(:) :: scr

    integer :: i,j,k
    integer :: tt

    !$omp parallel private(scr)
    allocate(scr(ny/2+1))
    tt=ceiling(real(ny)/3._rprec)
 
    !$omp do
    do k=1,dcp%ysz(3)
       do i=1,dcp%ysz(1)
          call dfftw_execute_dft_r2c(plany,ccin(i,:,k),scr)
          do j=1,lhy
             scr(j)=scr(j)*Fsmthg_y(j)
          enddo
          call dfftw_execute_dft_c2r(plany2,scr,wk_y(i,:,k))
       enddo
    enddo
    !$omp end do

    deallocate(scr)

    !$omp barrier
    !$omp single
    call transpose_y_to_x(wk_y,wk_x,dcpg)
    !$omp end single
    
    allocate(scr(nx/2+1))
    tt=ceiling(real(nx)/3._rprec)

    !$omp do
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          call dfftw_execute_dft_r2c(planx,wk_x(:,j,k),scr)
          do i=1,lhx
             scr(i)=scr(i)*Fsmthg_x(i)
          enddo
          call dfftw_execute_dft_c2r(planx2,scr,ccout(:,j,k))
       enddo
    enddo
    !$omp end do

    deallocate(scr)
    !$omp end parallel
    
    return
  end subroutine dealias

end module compute_convective_Fsmthg
