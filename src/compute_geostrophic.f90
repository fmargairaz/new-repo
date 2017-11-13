!------------------------------------------------------------------------------
!!  module: compute geostrophic
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 29/10/2015
!!
!! - description:
!!
!------------------------------------------------------------------------------

module compute_geostrophic
  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  
  real(rprec),parameter :: t_period = 72.0_rprec*3600.0_rprec
  real(rprec),parameter :: up_limit = 6.0_rprec*3600.0_rprec 
  
  private
  public :: compute_geostrophic_correct
contains

  subroutine compute_geostrophic_correct(jt_total)
    implicit none
    
    integer,intent(IN) :: jt_total
    
    real(rprec) :: ug_new,vg_new
    real(rprec) :: curr_time

    curr_time=real((jt_total-geostr_start)*dt_dim,rprec)

    ug_new=sqrt(ug**2+vg**2)*cos(2.0_rprec*pi*curr_time/t_period)
    vg_new=sqrt(ug**2+vg**2)*sin(2.0_rprec*pi*curr_time/t_period)

    if(curr_time.ge.up_limit)then
       ug_new=sqrt(ug**2+vg**2)*cos((2.0_rprec*pi*up_limit)/t_period)
       vg_new=sqrt(ug**2+vg**2)*sin((2.0_rprec*pi*up_limit)/t_period)
    endif

    ug=ug_new
    vg=vg_new
  end subroutine compute_geostrophic_correct

end module compute_geostrophic
