module bc_mod

  use decomp_2d
  use parameters_IO

  implicit none

  private

  public :: set_bc
  
  !---------------------------------------------------------------------------------
contains
  
  
  
  !!--------------------------------------------------------------------------------
  !!  object: subroutine set_lbc()
  !!--------------------------------------------------------------------------------
  !!
  !!  last checked/modified: marco giometto ( mgiometto@gmail.com ) on 05/11/2013
  !!
  !!--------------------------------------------------------------------------------
  subroutine set_bc(ind,var,bc_type,bc,stag,decomp)
    
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: decomp

    real(rprec),intent(in) :: bc
    integer,intent(in) :: ind,stag
    real(rprec),dimension(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)),intent(inout) :: var
    character(*),intent(in) :: bc_type
    integer :: jx,jy
    
    if(ind.eq.0)then
                     
       if(stag.eq.1)then
          if(bc_type .eq. "drc") var(:,:,1)=2.0_rprec*bc-var(:,:,1)
          if(bc_type .eq. "neu") var(:,:,1)=var(:,:,1)-bc*dz
       elseif(stag.eq.0)then
          if(bc_type .eq. "drc") var(:,:,1)=bc
          if(bc_type .eq. "neu") print*,'do not need'
       endif
       
    elseif(ind.eq.1)then
       
       if(stag.eq.1)then
          if(bc_type .eq. "drc") var(:,:,decomp%xsz(3))=2.0_rprec*bc-var(:,:,decomp%xsz(3)-1)
          if(bc_type .eq. "neu") var(:,:,decomp%xsz(3))=var(:,:,decomp%xsz(3)-1)+bc*dz
       elseif(stag.eq.0)then
          if(bc_type .eq. "drc") var(:,:,decomp%xsz(3))=bc
          if(bc_type .eq. "neu") print*,'do not need'
       endif
    endif
    
  end subroutine set_bc
  
end module bc_mod
