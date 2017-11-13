!------------------------------------------------------------------------------
!! module: compute convective terms
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!  - use the 3/2 dealiasing rule
!!  - pencil decomposition: variables called _big are in the Y-pencil
!!  - vatiable size: _hyb are (3/2nx,ny,nz), _big are (3/2nx,3/2ny,nz)
!!
!------------------------------------------------------------------------------

module compute_convective_padd

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
  
  real(rprec),allocatable,dimension(:,:,:) :: wk_hyb_x,wk_hyb_y
  real(rprec),allocatable,dimension(:,:,:) :: wk1_big
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
  !!  - need dcp_info for 3/2 rule (to allocate 2Decomp workspace)
  !!
  !------------------------------------------------------------------------------
  subroutine conv_init(dcp_main,dcp_ghost)
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
  end subroutine conv_init

  !------------------------------------------------------------------------------
  !! subroutine: finilise convective module
  !------------------------------------------------------------------------------
  subroutine conv_finalize() 
    implicit none

    deallocate(wk_hyb_x,wk_hyb_y)
    deallocate(wk1_big)

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
    real(rprec),dimension(dcp_big%ysz(1),dcp_big%ysz(2),0:dcp_big%ysz(3)+1) :: u_big,v_big,w_big
    real(rprec),dimension(dcp_big%ysz(1),dcp_big%ysz(2),0:dcp_big%ysz(3)+1) :: vort1_big,vort2_big,vort3_big

    real(rprec) :: const
    integer :: i,j,k,k_min
    integer :: ierror

    if(dcp%xst(3)==1)then
       k_min = 0
    else
       k_min = 1
    endif

    ! compute padded velocity
    const=1.0_rprec/(nx*ny)

    !$omp parallel do
    do k=0,dcp%xsz(3)+1
       cx(:,:,k) = -const*u(:,:,k)
       cy(:,:,k) = -const*v(:,:,k)
       cz(:,:,k) = -const*w(:,:,k)
    enddo
    !$omp end parallel do
    
    call transpose_to_big(u_big,cx)
    call transpose_to_big(v_big,cy)
    call transpose_to_big(w_big,cz)
    
    ! compute padded vorticity
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

    call transpose_to_big(vort1_big,cx)
    call transpose_to_big(vort2_big,cy)
    call transpose_to_big(vort3_big,cz)

    ! update ghost cell
    call update_ghost(u_big,dcp_big%ysz(1),dcp_big%ysz(2),dcp_big%ysz(3),dcp_big)
    call update_ghost(v_big,dcp_big%ysz(1),dcp_big%ysz(2),dcp_big%ysz(3),dcp_big)
    call update_ghost(w_big,dcp_big%ysz(1),dcp_big%ysz(2),dcp_big%ysz(3),dcp_big)
    call update_ghost(vort1_big,dcp_big%ysz(1),dcp_big%ysz(2),dcp_big%ysz(3),dcp_big)
    call update_ghost(vort2_big,dcp_big%ysz(1),dcp_big%ysz(2),dcp_big%ysz(3),dcp_big)
    call update_ghost(vort3_big,dcp_big%ysz(1),dcp_big%ysz(2),dcp_big%ysz(3),dcp_big)

    ! compute convective terms in each direction
    const=1.0_rprec/(nx2*ny2)

    ! compute x-term ====================================================
    ! uvp-nodes so 1,Nz-1
    if(dcp%xst(3)==1)then
       wk1_big(:,:,1) = const*(-v_big(:,:,1)*vort3_big(:,:,1) + &
            0.5_rprec*w_big(:,:,2)*vort2_big(:,:,2))
       k_min = 2
    else
       k_min = 1
    endif

    !$omp parallel do
    do k=k_min,dcp%xsz(3)
       do j=1,dcp_big%ysz(2)
          do i=1,dcp_big%ysz(1)
             wk1_big(i,j,k) = const*(-v_big(i,j,k)*vort3_big(i,j,k) + &
                  0.5_rprec*(w_big(i,j,k+1)*vort2_big(i,j,k+1) + &
                  w_big(i,j,k)*vort2_big(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel do

    call dealias(cx,wk1_big)

    ! compute y-term ====================================================
    ! uvp-nodes so 1,Nz-1
    if(dcp%xst(3)==1)then
       wk1_big(:,:,1) = const*(u_big(:,:,1)*vort3_big(:,:,1) + &
            0.5_rprec*w_big(:,:,2)*(-vort1_big(:,:,2)))
       k_min = 2
    else
       k_min = 1
    endif
    !$omp parallel do
    do k=k_min,dcp%xsz(3)
       do j=1,dcp_big%ysz(2)
          do i=1,dcp_big%ysz(1)
             wk1_big(i,j,k) = const*(u_big(i,j,k)*vort3_big(i,j,k) + &
                  0.5_rprec*(w_big(i,j,k)*(-vort1_big(i,j,k)) + &
                  w_big(i,j,k+1)*(-vort1_big(i,j,k+1))))
          enddo
       enddo
    enddo
    !$omp end parallel do

    call dealias(cy,wk1_big)

    ! compute z-term ====================================================
    ! w-nodes so 1,Nz
    if(dcp%xst(3)==1)then
       wk1_big(:,:,1) = 0.0_rprec
       k_min = 2
    else
       k_min = 1
    endif
    !$omp parallel do
    do k=k_min,dcp%xsz(3)
       do j=1,dcp_big%ysz(2)
          do i=1,dcp_big%ysz(1)
             wk1_big(i,j,k) = const*0.5_rprec*(&
                  (u_big(i,j,k)+u_big(i,j,k-1))*(-vort2_big(i,j,k)) + &
                  (v_big(i,j,k)+v_big(i,j,k-1))*(vort1_big(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel do

    call dealias(cz,wk1_big)

    !--bottom level is not valid
    cx(:,:,0) = BOGUS
    cy(:,:,0) = BOGUS
    cz(:,:,0) = BOGUS
    
    !--top level is not valid
    cx(:,:,dcp%xsz(3)+1) = BOGUS
    cy(:,:,dcp%xsz(3)+1) = BOGUS
    cz(:,:,dcp%xsz(3)+1) = BOGUS 
    
    return
  end subroutine compute_dealias_convec
  
  !------------------------------------------------------------------------------
  !! subroutine: transpose data to big Y-pencil
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - both cc and cc_big do need ghost cells support
  !!
  !------------------------------------------------------------------------------
  subroutine transpose_to_big(cc_big,cc) 
    implicit none
 
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(in) :: cc
    real(rprec),dimension(dcp_big%ysz(1),dcp_big%ysz(2),0:dcp_big%ysz(3)+1),intent(inout) :: cc_big

    complex(rprec),allocatable,dimension(:) :: scr,scr_big

    integer :: i,j,k
    
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

    return
  end subroutine transpose_to_big

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
  subroutine dealias(cc,cc_big)
    implicit none
    
    real(rprec),dimension(dcp_big%ysz(1),dcp_big%ysz(2),0:dcp_big%ysz(3)+1),intent(in) :: cc_big 
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(out) :: cc
    

    complex(rprec),allocatable,dimension(:) :: scr,scr_big

    integer :: i,j,k
    
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
    
    return
  end subroutine dealias

end module compute_convective_padd
