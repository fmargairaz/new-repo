!-------------------------------------
!fftw (3.2.x versions)
!---------------------------------------------------------------------------
module fft_mod

  use parameters_IO
  
  implicit none

  include "fftw3.f"
  
  !private

  !public :: kx,ky,k2,eye,forw,back,
  !public :: init_fft
  
  integer*8 :: forw,back                        !as in the manual
  real(rprec), allocatable, dimension(:) :: kx
  real(rprec), allocatable, dimension(:) :: ky
  integer*8 :: planx,plany,planx2,plany2
  
contains
  
  !!-----------------------------------------------------------------------------------
  !!  object: subroutine init_fft ()
  !!-----------------------------------------------------------------------------------
  !!
  !!-----------------------------------------------------------------------------------
  subroutine fft_init()
    
    
    implicit none
    
    real(rprec),dimension(nx) :: inx
    real(rprec),dimension(ny) :: iny
    double complex,dimension(lhx) :: outx
    double complex,dimension(lhy) :: outy
    
    !formulate the fft plans (forw, n0, n1, double in, complex out, flags)
    call dfftw_plan_dft_r2c_1d(planx,nx,inx,outx,fftw_estimate)
    call dfftw_plan_dft_r2c_1d(plany,ny,iny,outy,fftw_estimate)
    
    call dfftw_plan_dft_c2r_1d(planx2,nx,outx,inx,fftw_estimate)
    call dfftw_plan_dft_c2r_1d(plany2,ny,outy,iny,fftw_estimate)

    call init_wavenumber()
    
  end subroutine fft_init
  
  
  
  !!-----------------------------------------------------------------------------------
  !!  object: subroutine init_wavenumber ()
  !!-----------------------------------------------------------------------------------
  !!
  !!-----------------------------------------------------------------------------------
  subroutine init_wavenumber ()
    
    implicit none    

    integer :: jx,jy

    allocate(kx(lhx))
    allocate(ky(lhy))

    !fills kx starting from jx=0 to jx=nx/2 (the freq nx/2+1 is filled with 0 later)
    do jx=1,lhx
       kx(jx)=real(jx-1,rprec)
    enddo
    
    do jy=1,lhy
       ky(jy)=real(jy-1,rprec)
    enddo
    
    !eliminates the nyquist value to calculate a correct first derivative later
    kx(lhx)=0._rprec
    ky(lhy)=0._rprec
    
    !renormalization for the aspect ratio change
    kx=2._rprec*pi/lx*kx
    ky=2._rprec*pi/ly*ky
    
  end subroutine init_wavenumber
  
  
end module fft_mod
