!------------------------------------------------------------------------------
!!  module: compute scalar
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!! - travis morrison (tjmorrison635@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------

module compute_scalars

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use utility_tools
  use dealiasing_module

  !! - global variables
  use system_variables, only : u,v,w
  use compute_derivatives
  use compute_sgs, only : nu_t   

  implicit none

  private
  TYPE(DECOMP_INFO),save :: dcp
  
  !real(rprec),dimension(:,:,:),allocatable,public :: cs_opt2,nu_t
  !real(rprec),parameter::n_dmp=2._rprec

  public :: scalars_module_init,scalars_module_finalize

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: initialise scalars module
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine scalars_module_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main

    !allocate(cs_opt2(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    !if (status /= 0) then
    !   errorcode = 3
    !   call decomp_2d_abort(errorcode, &
    !        'Out of memory when initialising sgs module variables')
    !end if

    return
  end subroutine scalars_module_init
  
  !------------------------------------------------------------------------------
  !! subroutine: finalise scalars module
  !------------------------------------------------------------------------------=
  subroutine scalars_module_finalize()
    
    implicit none
    
    !deallocate(cs_opt2, nu_t)
    
  end subroutine scalars_module_finalize

  !------------------------------------------------------------------------------
  !! subroutine: compute scalars feedback - Boussinesq approximation
  !------------------------------------------------------------------------------
  !! - description:
  !!  - beta is stored on W nodes, but Theta is on UVP nodes
  !!
  !------------------------------------------------------------------------------
  subroutine compute_scalars_main()
    implicit none

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1) ::  scalar
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)) ::  rhs1,rhs0

    call compute_rhs_scalars()
    call step_scalars(scalar,rhs1,rhs0)
    call compute_beta_scalars()

    return
  end subroutine compute_scalars_main
  
  !------------------------------------------------------------------------------
  !! subroutine: compute rhs for scalar
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_rhs_scalars()
    implicit none

    return
  end subroutine compute_rhs_scalars

  !------------------------------------------------------------------------------
  !! subroutine: compute time intergation of scalar
  !------------------------------------------------------------------------------
  !! - description:
  !!  - using second order Adams-Bashforth scheme
  !!
  !------------------------------------------------------------------------------
  subroutine step_scalars(scalar,rhs1,rhs0)
    implicit none
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(INOUT) ::  scalar
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)),intent(IN) ::  rhs1,rhs0

    integer :: i,j,k

    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             scalar(i,j,k)=scalar(i,j,k)+dt*(1.5_rprec*rhs1(i,j,k)-0.5_rprec*rhs0(i,j,k))
          enddo
       enddo
    enddo

    ! update the ghost cell
    call update_ghost(scalar,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)

    !below and above the domaine
    if(dcp%xst(3)==1)then
       scalar(:,:,0)=0._rprec
    endif
    if(dcp%xen(3)==nz)then
       scalar(:,:,dcp%xsz(3)+1)=scalar(:,:,dcp%xsz(3))
    endif
    
    return
  end subroutine step_scalars

  !------------------------------------------------------------------------------
  !! subroutine: compute beta for scalar
  !------------------------------------------------------------------------------
  !! - description:
  !!  - beta is stored on W nodes, but Theta is on UVP nodes
  !------------------------------------------------------------------------------
  subroutine compute_beta_scalars()
    implicit none

    return
  end subroutine compute_beta_scalars

end module compute_scalars
