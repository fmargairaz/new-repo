!------------------------------------------------------------------------------
!!  module: FFT engine
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
!!
!! - description: 
!!  - this code uese the FFTW (version 3.x)
!!
!------------------------------------------------------------------------------
module fft_engine

  use parameters_IO
  
  implicit none

  include "fftw3.f"
  
  private
  
  ! wavelength in x and y
  real(rprec),allocatable,dimension(:),public :: kx,ky
  real(rprec),allocatable,dimension(:,:),public :: kkx,kky,kk2

  ! fftw plans for nx & ny variables 
  integer*8,public :: planx,plany,planx2,plany2
  
  ! fftw plans for 3/2nx & 3/2ny varilabes
  integer*8,public :: planx_big,plany_big,planx2_big,plany2_big
  
  ! fftw plans for 2d varilabes
  integer*8,public :: plan2d_f,plan2d_b
  
  public :: fft_init

contains !=======================================================================
  
  !-----------------------------------------------------------------------------------
  !! subroutine: init_fft ()
  !-----------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
  !!
  !! - description:
  !!
  !-----------------------------------------------------------------------------------
  subroutine fft_init()
    
    implicit none

    real(rprec),dimension(nx) :: inx
    real(rprec),dimension(ny) :: iny
    complex(rprec),dimension(lhx) :: outx
    complex(rprec),dimension(lhy) :: outy
    
    real(rprec),dimension(nx2) :: inx_big
    real(rprec),dimension(ny2) :: iny_big
    complex(rprec),dimension(nx2/2+1) :: outx_big
    complex(rprec),dimension(ny2/2+1) :: outy_big

    real(rprec),dimension(nx,ny) :: in2d
    complex(rprec),dimension(lhx,ny) :: out2d

    !formulate the fft plans (forw, n0, n1, double in, complex out, flags)
    call dfftw_plan_dft_r2c_1d(planx,nx,inx,outx,fftw_estimate)
    call dfftw_plan_dft_r2c_1d(plany,ny,iny,outy,fftw_estimate)
    
    call dfftw_plan_dft_c2r_1d(planx2,nx,outx,inx,fftw_estimate)
    call dfftw_plan_dft_c2r_1d(plany2,ny,outy,iny,fftw_estimate)

    call dfftw_plan_dft_r2c_1d(planx_big,nx2,inx_big,outx_big,fftw_estimate)
    call dfftw_plan_dft_r2c_1d(plany_big,ny2,iny_big,outy_big,fftw_estimate)
    
    call dfftw_plan_dft_c2r_1d(planx2_big,nx2,outx_big,inx_big,fftw_estimate)
    call dfftw_plan_dft_c2r_1d(plany2_big,ny2,outy_big,iny_big,fftw_estimate)
    
    call dfftw_plan_dft_r2c_2d(plan2d_f,nx,ny,in2d,out2d,fftw_estimate)
    call dfftw_plan_dft_c2r_2d(plan2d_b,nx,ny,out2d,in2d,fftw_estimate)
        
    call build_wavenumbers()
    
    return
  end subroutine fft_init
  
  !-----------------------------------------------------------------------------------
  !! subroutine: build wavenumbers arrays in sprectral space
  !-----------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/06/2014
  !!
  !! - description:
  !!  - kx&ky in sprectral space
  !!  - adds at nyquist freq values of kx, ky = 0 for the derivatives  
  !!  - renormalizes the wavenumbers to take into account the real domain dimention
  !!
  !-----------------------------------------------------------------------------------
  subroutine build_wavenumbers ()
    
    implicit none    

    integer :: status, errorcode
    integer i,j

    ! = 1D wavenumber ==========
    allocate(kx(lhx));allocate(ky(lhy))
    kx=0._rprec;ky=0._rprec;

    !fills kx starting from jx=0 to jx=nx/2 (the freq nx/2+1 is filled with 0 later)
    do i=1,lhx
       kx(i)=real(i-1,rprec)
    enddo
    do j=1,lhy
       ky(j)=real(j-1,rprec)
    enddo
    !eliminates the nyquist value to calculate a correct first derivative later
    kx(lhx)=0._rprec
    ky(lhy)=0._rprec
    !renormalization for the aspect ratio change
    kx=2._rprec*pi/lx*kx
    ky=2._rprec*pi/ly*ky

    ! = 2D wavenumber ==========
    allocate(kkx(lhx,ny));allocate(kky(lhx,ny));allocate(kk2(lhx,ny))
    kkx=0._rprec;kky=0._rprec;kk2=0._rprec

    !fills kx starting from jx=0 to jx=nx/2 (the freq nx/2+1 is filled with 0 later)
    do i=1,lhx
       kkx(i,:) = real(i-1,kind=rprec)
    end do
    do j=1,ny
       kky(:,j) = real(modulo(j - 1 + ny/2,ny) - ny/2, kind=rprec)
    end do

    ! eliminates the nyquist value to calculate a correct first derivative later
    kkx(lhx,:)=0._rprec
    kky(lhx,:)=0._rprec
    kkx(:,lhy)=0._rprec
    kky(:,lhy)=0._rprec
    ! renormalization for the aspect ratio change
    kkx=2._rprec*pi/lx*kkx
    kky=2._rprec*pi/ly*kky

    kk2 = kkx*kkx + kky*kky
    return
  end subroutine build_wavenumbers

end module fft_engine
