!------------------------------------------------------------------------------
!!  module: compute wall law
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
!!
!! - description:
!!
!------------------------------------------------------------------------------

module compute_wall_law

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use testfilter
  use utility_tools

  !! - global variables
  use system_variables, only : u,v,w,dudz,dvdz,txz,tyz

  implicit none
  include 'mpif.h'

  private
  TYPE(DECOMP_INFO),save :: dcp
  
  real(rprec),dimension(:,:),allocatable,public :: lbc_z0,lbc_sc_z0,u_star_avg

  public :: wall_init,wall_finalize
  public :: wallstress,obukhov

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: initialise wall module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 02/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine wall_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main

    allocate(lbc_z0(dcp%xsz(1),dcp%xsz(2)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising wall module variables')
    end if
    
    allocate(lbc_sc_z0(dcp%xsz(1),dcp%xsz(2)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising wall module variables')
    end if
    
    allocate(u_star_avg(dcp%xsz(1),dcp%xsz(2)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising wall module variables')
    end if
    
    lbc_z0=0._rprec
    lbc_sc_z0=0._rprec
    u_star_avg=0._rprec
    
    return
  end subroutine wall_init

  !------------------------------------------------------------------------------
  !! subroutine: finalise wall module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 02/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine wall_finalize()
    implicit none
    if(dcp%xst(3)==1)then
       deallocate(lbc_z0,lbc_sc_z0,u_star_avg)
    endif
    return
  end subroutine wall_finalize

  !------------------------------------------------------------------------------
  !! subroutine: compute law of the wall
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margaira@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine wallstress(bc_u,bc_v,wl_type)
    implicit none

    character(*),intent(in) :: wl_type
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(INOUT) :: bc_u,bc_v

    real(rprec) :: u1_s,v1_s,u1_sum,v1_sum,const
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)) :: u1,v1,u_avg
    integer:: i,j,k
    integer :: errorcode,ierror

    !compute u1,v1,u_avg
    if(wl_type .eq. "wall_homog")then
       ! average way ======================================
       call get_layer_average(u1,u,1,dcp)
       call get_layer_average(v1,v,1,dcp)
       
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             u_avg(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
          enddo
       enddo
       
    elseif(wl_type .eq. "wall_inst")then
       ! instantaneous way ======================================
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             u1(i,j)=u(i,j,1)
             v1(i,j)=v(i,j,1)
             u_avg(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
          enddo
       enddo
    elseif(wl_type .eq. "wall_filter")then
       call get_layer_filter(u1,u,1,G_test,dcp)
       call get_layer_filter(v1,v,1,G_test,dcp)
       ! compute u_avg
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)                
             u_avg(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
          end do
       end do
    else
       errorcode = 1
       call decomp_2d_abort(errorcode,'Invalid wall law type')
    endif

    ! compute ustar_avg
    !call obukhov(u_star_avg,u1,v1,u_avg)

    ! add BC to nodes @dz/2
    do j=1,dcp%xsz(2)
       do i=1,dcp%xsz(1)
          ! check for nan values in case of 0 velocity
          if(u_avg(i,j).eq.0.0_rprec)then
             dudz(i,j,1)=0.0_rprec
             dvdz(i,j,1)=0.0_rprec
             
             !using ghost cell k=0 to impose Newmann bc
             bc_u(i,j)=0.0_rprec
             bc_v(i,j)=0.0_rprec

             txz(i,j,1)=0.0_rprec
             tyz(i,j,1)=0.0_rprec
          else
             !using ghost cell k=0 to impose Newmann bc
             bc_u(i,j)=u_star_avg(i,j)/(0.5_rprec*dz*KvonK)*u1(i,j)/u_avg(i,j)
             bc_v(i,j)=u_star_avg(i,j)/(0.5_rprec*dz*KvonK)*v1(i,j)/u_avg(i,j)
             
             dudz(i,j,1)=bc_u(i,j)
             dvdz(i,j,1)=bc_v(i,j)
             
             txz(i,j,1)=(u_star_avg(i,j)**2)*u1(i,j)/u_avg(i,j)
             tyz(i,j,1)=(u_star_avg(i,j)**2)*v1(i,j)/u_avg(i,j)
          endif
       enddo
    enddo

    return
  end subroutine wallstress


  !!------------------------------------------------------------------------------
  !!  subroutine: obukhov
  !!------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!   - add phyisics afterwards 
  !!
  !!------------------------------------------------------------------------------
  subroutine obukhov(ustar_avg,wl_type)
    implicit none
    
    character(*),intent(in) :: wl_type
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(out) :: ustar_avg

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)) :: u1,v1,u_avg
    real(rprec) :: u1_s,v1_s,u1_sum,v1_sum,const
    integer:: i,j,k
    integer :: errorcode,ierror

    !compute u1,v1,u_avg
    if(wl_type .eq. "wall_homog")then
       ! average way ======================================
       call get_layer_average(u1,u,1,dcp)
       call get_layer_average(v1,v,1,dcp)
       
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             u_avg(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
          enddo
       enddo

    elseif(wl_type .eq. "wall_inst")then
       ! instantaneous way ======================================
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             u1(i,j)=u(i,j,1)
             v1(i,j)=v(i,j,1)
             u_avg(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
          enddo
       enddo
    elseif(wl_type .eq. "wall_filter")then
       call get_layer_filter(u1,u,1,G_test,dcp)
       call get_layer_filter(v1,v,1,G_test,dcp)
       
       ! compute u_avg
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)                
             u_avg(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
          end do
       end do
    else
       errorcode = 1
       call decomp_2d_abort(errorcode,'Invalid wall law type')
    endif

    ! neutral case ======================================
    !compute ustar_avg
    do j=1,dcp%xsz(2)
       do i=1,dcp%xsz(1)
          ustar_avg(i,j) = KvonK*u_avg(i,j)/dlog((0.5_rprec*dz)/lbc_z0(i,j)) 
       enddo
    enddo

    return
  end subroutine obukhov

end module compute_wall_law
