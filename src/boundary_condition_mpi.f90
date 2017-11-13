!------------------------------------------------------------------------------
!!  module: impose bounrady conditions
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
!!
!! - description:
!!
!------------------------------------------------------------------------------

module boundary_conditions

  use decomp_2d
  use parameters_IO

  implicit none
  include 'mpif.h'

  private

  public :: set_bc, set_bc_surface, set_lbc_z0, set_lbc_sc_z0
  public :: compute_bc
  
contains !=======================================================================
  
  !--------------------------------------------------------------------------------
  !! subroutine: set_lbc()
  !--------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 05/11/2013
  !!
  !--------------------------------------------------------------------------------
  subroutine set_bc(pos,var,bc_type,bc,stag,dcp)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp

    real(rprec),dimension(:,:),intent(in) :: bc
    logical,intent(in) :: stag
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(INOUT) :: var
    character(*),intent(in) :: pos,bc_type
    integer :: jx,jy
    
    if(pos.eq.'bot')then
       if(stag)then
          if(bc_type .eq. "drc") var(:,:,0)=2.0_rprec*bc(:,:)-var(:,:,1)
          if(bc_type .eq. "neu") var(:,:,0)=var(:,:,1)-bc(:,:)*dz
       else
          if(bc_type .eq. "drc") var(:,:,1)=bc(:,:)
          if(bc_type .eq. "neu") write(*,*) 'do not need'
       endif
 
    elseif(pos.eq.'top')then
       if(stag)then
          if(bc_type .eq. "drc") var(:,:,dcp%xsz(3)+1)=2.0_rprec*bc(:,:)-var(:,:,dcp%xsz(3))
          if(bc_type .eq. "neu") var(:,:,dcp%xsz(3)+1)=var(:,:,dcp%xsz(3))+bc(:,:)*dz
       else
          if(bc_type .eq. "drc") var(:,:,dcp%xsz(3)+1)=bc(:,:)
          if(bc_type .eq. "neu") write(*,*) 'do not need'
       endif
    endif
    
  end subroutine set_bc


  !--------------------------------------------------------------------------------
  !! subroutine: compute_lbc()
  !--------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 05/11/2013
  !!
  !--------------------------------------------------------------------------------
  subroutine compute_bc(pos,bc_u,bc_v,bc_w,dcp)
    implicit none
    TYPE(DECOMP_INFO),intent(in) :: dcp

    real(rprec),dimension(:,:),intent(INOUT) :: bc_u,bc_v,bc_w
    character(*),intent(in) :: pos

    if(pos.eq.'bot')then
       bc_u(:,:)=lbc%u
       bc_v(:,:)=lbc%v
       bc_w(:,:)=lbc%v
    elseif(pos.eq.'top')then
       bc_u(:,:)=ubc%u
       bc_v(:,:)=ubc%v
       bc_w(:,:)=ubc%w
    endif

    return
  end subroutine compute_bc


  !------------------------------------------------------------------------------
  !! subroutine: 
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine set_bc_surface(lbc_z0,dcp)
    implicit none
    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(OUT) :: lbc_z0

    real(rprec),dimension(nx,ny) :: z0_in
    integer :: i,j,reclen,ierror
    
    if(surface_files)then
       if(nrank==0)then
          inquire(iolength=reclen) z0_in
          open(98,file='./input/surface_z0.dat',form='unformatted',access='direct',recl=reclen)
          read(98,rec=1) z0_in
          close(98)
       endif
       
       call MPI_Bcast(z0_in,nx*ny,real_type,0,MPI_COMM_WORLD, ierror);
   
       ! Non-dimensionalize
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             lbc_z0(i,j)=z0_in(dcp%xst(1)-1+i,dcp%xst(2)-1+j)/zi
          enddo
       enddo
    else
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             lbc_z0(i,j)=lbc_dz0/zi
          enddo
       enddo
    endif

    return
  end subroutine set_bc_surface
  
  
  !!--------------------------------------------------------------------------------
  !!  subroutine: set_lbc_z0()
  !!--------------------------------------------------------------------------------
  !!
  !!  last checked/modified: marco giometto ( mgiometto@gmail.com ) on 16/10/2012
  !!
  !!  its done like this so that can be nicely connected with remote or diff values
  !!
  !!--------------------------------------------------------------------------------
  subroutine set_lbc_z0(lbc_z0)
    implicit none
    
    real(rprec),dimension(:,:),intent(inout) :: lbc_z0 
    
    lbc_z0(:,:)=lbc_dz0
    
    return
  end subroutine set_lbc_z0
  
  !!--------------------------------------------------------------------------------
  !!  subroutine set_lbc_sc_z0()
  !!--------------------------------------------------------------------------------
  !!
  !!  last checked/modified: marco giometto ( mgiometto@gmail.com ) on 16/10/2012
  !!
  !!  its done like this so that can be nicely connected with remote or diff values
  !!
  !!--------------------------------------------------------------------------------
  subroutine set_lbc_sc_z0(lbc_sc_z0)
    implicit none
    
    real(rprec),dimension(:,:),intent(inout) :: lbc_sc_z0
    
    lbc_sc_z0=lbc_dz0
    
    return
  end subroutine set_lbc_sc_z0
  
end module boundary_conditions
