!------------------------------------------------------------------------------
!! module: compute convective terms
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
!!
!! - description:
!!  - use the 2/3 dealiasing rule
!!  
!!
!------------------------------------------------------------------------------

module compute_convective_nodlsg
  
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
  integer :: ierr

  public :: conv_init,conv_finalize
  public :: compute_dealias_convec

contains !=======================================================================
  
  !------------------------------------------------------------------------------
  !! subroutine: initiate convective module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!  - need dcp_info for 2/3 rule (to allocate 2Decomp workspace)
  !!
  !------------------------------------------------------------------------------
  subroutine conv_init(dcp_main,dcp_ghost)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main,dcp_ghost
    integer :: status, errorcode

    ! import dcp for main MPI engine
    dcp = dcp_main
    dcpg = dcp_ghost
   
    return
  end subroutine conv_init

  !------------------------------------------------------------------------------
  !! subroutine: finilise convective module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine conv_finalize() 
    implicit none
    
    return
  end subroutine conv_finalize

  !------------------------------------------------------------------------------
  !! subroutine: compute aliased convective terms
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
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
       u1(:,:,k) = -const*u(:,:,k)
       u2(:,:,k) = -const*v(:,:,k)
       u3(:,:,k) = -const*w(:,:,k)
    enddo
    !$omp end parallel do
    
    ! compute truncated vorticity
    !$omp parallel do
    do k=0,dcp%xsz(3)+1
       vort1(:,:,k) = const*(dwdy(:,:,k)-dvdz(:,:,k))
       vort2(:,:,k) = const*(dudz(:,:,k)-dwdx(:,:,k))
       vort3(:,:,k) = const*(dvdx(:,:,k)-dudy(:,:,k))
    enddo
    !$omp end parallel do    

#ifndef DNS_MOD
    ! special treatment k=1 cause vort1,vort2 is on uvp (vort3 always uvp)
    if((dcp%xst(3)==1).and.(lbc%special.eq.'wall_law'))then
       vort1(:,:,0) = 0.0_rprec
       vort2(:,:,0) = 0.0_rprec

       vort1(:,:,1) = const*(0.5_rprec*(dwdy(:,:,1)+dwdy(:,:,2))-dvdz(:,:,1))      
       vort2(:,:,1) = const*(dudz(:,:,1)-0.5_rprec*(dwdx(:,:,1)+dwdx(:,:,2)))
    endif
#endif

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
       cx(:,:,1) = const*(-u2(:,:,1)*vort3(:,:,1) + &
            0.5_rprec*u3(:,:,2)*vort2(:,:,2))
       k_min = 2
    else
       k_min = 1
    endif
    !$omp parallel do
    do k=k_min,dcp%ysz(3)
       do j=1,dcp%ysz(2)
          do i=1,dcp%ysz(1)
             cx(i,j,k) = const*(-u2(i,j,k)*vort3(i,j,k) + &
                  0.5_rprec*(u3(i,j,k+1)*vort2(i,j,k+1) + &
                  u3(i,j,k)*vort2(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! compute y-term ====================================================
    ! uvp-nodes so 1,Nz-1
     if(dcp%xst(3)==1)then
       cy(:,:,1) = const*(u1(:,:,1)*vort3(:,:,1) + &
            0.5_rprec*u3(:,:,2)*(-vort1(:,:,2)))
       k_min = 2
    else
       k_min = 1
    endif
    !$omp parallel do
    do k=k_min,dcp%ysz(3)
       do j=1,dcp%ysz(2)
          do i=1,dcp%ysz(1)
             cy(i,j,k) = const*(u1(i,j,k)*vort3(i,j,k) + &
                  0.5_rprec*(u3(i,j,k)*(-vort1(i,j,k)) + &
                  u3(i,j,k+1)*(-vort1(i,j,k+1))))
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! compute z-term ====================================================
    ! w-nodes so 1,Nz
    if(dcp%xst(3)==1)then
       cz(:,:,1) = 0.0_rprec
       k_min = 2
    else
       k_min = 1
    endif
    !$omp parallel do
    do k=k_min,dcp%ysz(3)
       do j=1,dcp%ysz(2)
          do i=1,dcp%ysz(1)
             cz(i,j,k) = const*0.5_rprec*(&
                  (u1(i,j,k)+u1(i,j,k-1))*(-vort2(i,j,k)) + &
                  (u2(i,j,k)+u2(i,j,k-1))*(vort1(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine compute_dealias_convec

end module compute_convective_nodlsg
