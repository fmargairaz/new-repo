!------------------------------------------------------------------------------
!!  module: compute sub grid stress 
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
!!
!! - description:
!!  - this module contains all the routines used for computing sgs
!!  - cs_opt2 & nu_t are definded on w-nodes
!!
!------------------------------------------------------------------------------

module compute_dns_stress

  use decomp_2d
  use parameters_IO, only : rprec

  !! - global variables
  use system_variables, only : dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,&
       txx,tyy,tzz,txy,txz,tyz

  implicit none

  private

  public :: dns_stress
   
contains !=======================================================================
  
  !!--------------------------------------------------------------------------------
  !! subroutine: dns_stress
  !!--------------------------------------------------------------------------------
  !!
  !!  last checked/modified: marco giometto ( mgiometto@gmail.com ) on 26/03/2013
  !!
  !!--------------------------------------------------------------------------------
  subroutine dns_stress (dcp,Re,lbc_special)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp  
    real(rprec),intent(in) :: Re
    character(*),intent(in) :: lbc_special
    
    integer :: i,j,k,k_min
    
    ! txx,tyy,tzz,txy on uvp-nodes & txz,tyz on w-nodes
    if((dcp%xst(3)==1).and.(lbc_special.eq.'wall_law'))then
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1) 
             txx(i,j,1)=2.0_rprec*dudx(i,j,1)/Re
             tyy(i,j,1)=2.0_rprec*dvdy(i,j,1)/Re
             tzz(i,j,1)=2.0_rprec*dwdz(i,j,1)/Re
             txy(i,j,1)=(dudy(i,j,k)+dvdx(i,j,1))/Re
          enddo
       enddo
       k_min=2
    else
       k_min=1
    end if

    !$omp parallel do 
    do k=k_min,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1) 
             txx(i,j,k)=2.0_rprec*dudx(i,j,k)/Re
             tyy(i,j,k)=2.0_rprec*dvdy(i,j,k)/Re
             tzz(i,j,k)=2.0_rprec*dwdz(i,j,k)/Re
             txy(i,j,k)=(dudy(i,j,k)+dvdx(i,j,k))/Re
             txz(i,j,k)=(dudz(i,j,k)+dwdx(i,j,k))/Re
             tyz(i,j,k)=(dvdz(i,j,k)+dwdy(i,j,k))/Re
          enddo
       enddo
    enddo
   !$omp end parallel do

    return
  end subroutine dns_stress

end module compute_dns_stress
