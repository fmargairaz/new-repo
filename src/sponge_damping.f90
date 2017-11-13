!------------------------------------------------------------------------------
!!  module: sponge damping
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------

module sponge_damping

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use utility_tools

  !! - global variables
  use system_variables, only : u,v,w,rhsx,rhsy,rhsz

  implicit none
  include 'mpif.h'

  private
  TYPE(DECOMP_INFO),save :: dcp
  real(rprec),allocatable,dimension(:),save :: sponge

  public :: sponge_damping_init,sponge_damping_finalize
  public :: add_sponge_damping

contains !=======================================================================
  
  !------------------------------------------------------------------------------
  !! subroutine: initialise sponge damping module 
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine sponge_damping_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    real(rprec) :: sponge_layer_fraction=0._rprec!To add in parameter
    real(rprec) :: z,z_d,cfrdmp
    integer :: k
    integer :: status, errorcode

    dcp = dcp_main

    allocate(sponge(dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising splonge variable')
    end if

    sponge=0._rprec

    ! For the damping_method = 2 the sponge is computed at the nodes of w (not in u,v,p nodes) 
    ! Gerard Cortina 13/05/15
    z_d = sponge_layer_fraction*lz
    cfrdmp = 3.9_rprec
    
    do k=1,dcp%xsz(3)+1
       z=(dcp%xst(3)+k-2._rprec)*dz
       !z=(dcp%xst(3)-1+k-0.5_rprec)*dz*H0
       
       if((z.ge.z_d).and.(z.le.lz))then
          sponge(k)=cfrdmp*0.5_rprec*(1._rprec-cos(pi*(z-z_d)/(lz-z_d)))
       else
          sponge(k)=0._rprec
       end if
    end do
    
    if(dcp%xen(3)==nz)then
       sponge(dcp%xsz(3)+1)=sponge(dcp%xsz(3))
    end if

    sponge=0._rprec

!!$
!!$
!!$
!!$    elseif (damping_method==1) then
!!$       ! sets relaxation term to vertical momentum equation in top quarter
!!$       ! of domain, relaxation timescale on top is 50s with a factor of 5 if we
!!$       ! had Nz=40 for the layers 40...31, Nieuwstadt et al. 1991, turbulent shear 
!!$       ! flows
!!$
!!$       $if ($MPI)
!!$       !--the non-MPI recursive form in inconvenient, so instead replace with
!!$       !  analytically evaluated version (should be this way anyway...)
!!$       sponge=0._rprec
!!$       !  factor=9._rprec/(nz_global - 3*nz_global/4 + 1)
!!$       factor=9._rprec/(nz_tot - 3*nz_tot/4 + 1)
!!$       sponge_top = z_i / (50._rprec * u_star)
!!$       do k = 1, nz
!!$          k_global = k + coord * nz
!!$          if (k_global > 3*nz_global/4 + 1) then
!!$             sponge(k) = sponge_top * 5._rprec**((k_global-nz_global) * factor)
!!$          end if
!!$       end do
!!$       $else
!!$       sponge=0._rprec
!!$       factor=9._rprec/(nz-3*nz/4+1)
!!$       sponge(nz)=z_i/(50._rprec*u_star)
!!$       do k=nz-1,3*nz/4+1,-1
!!$          sponge(k)=sponge(k+1)/5._rprec**factor
!!$       end do
!!$       $endif
!!$
!!$    end if
    
    return
  end subroutine sponge_damping_init
  
  !------------------------------------------------------------------------------
  !! subroutine: finalise sponge damping module
  !------------------------------------------------------------------------------
  subroutine sponge_damping_finalize()
    
    implicit none
    
    deallocate(sponge)
  
  end subroutine sponge_damping_finalize

  !------------------------------------------------------------------------------
  !! subroutine: compute planar average 
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - add damping terms to the momentum RHS
  !!
  !------------------------------------------------------------------------------
  subroutine add_sponge_damping()

    real(rprec),dimension(dcp%xsz(3)) :: u_avg,v_avg,w_avg
    integer :: i,j,k

    call compute_planar_average(u_avg,u(:,:,1:dcp%xsz(3)),dcp)
    call compute_planar_average(v_avg,v(:,:,1:dcp%xsz(3)),dcp)
    call compute_planar_average(w_avg,w(:,:,1:dcp%xsz(3)),dcp)

    !$omp parallel do
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             RHSx(i,j,k)=RHSx(i,j,k)-0.5*(sponge(k)+sponge(k+1))*(u(i,j,k)-u_avg(k))
             RHSy(i,j,k)=RHSy(i,j,k)-0.5*(sponge(k)+sponge(k+1))*(v(i,j,k)-v_avg(k))
             RHSz(i,j,k)=RHSz(i,j,k)-sponge(k)*(w(i,j,k)-w_avg(k))
          enddo
       enddo
    enddo
    !$omp end parallel do
     
     return
   end subroutine add_sponge_damping
   
 end module sponge_damping
