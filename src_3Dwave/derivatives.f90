module derivatives_mod

  use parameters_IO
  use fft_mod
  use decomp_2d

  implicit none
  
  private
  
  public :: ddz_uv_4x,ddz_w_4x, ddx, ddy
  
  !##############################################################################

contains
  
  
  !!------------------------------------------------------------------------------
  !!  object: subroutine ddz_uv (dfdz, f)
  !!------------------------------------------------------------------------------
  !!
  !!  last checked/modified: marco giometto ( mgiometto@gmail.com ), day 05/11/2013
  !!
  !!  description:
  !!
  !!------------------------------------------------------------------------------
  subroutine ddz_uv_4x(dfdz,f,decomp)

    implicit none
             
    TYPE(DECOMP_INFO),intent(in) :: decomp

    !f and dfdx are real arrays as they don't need to go through the fft
    real(rprec),dimension(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)),intent(in) :: f
    real(rprec),dimension(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)),intent(inout) :: dfdz
    
    real(rprec), allocatable, dimension(:,:,:) :: f_halo
  
    integer :: k
    real(rprec) :: const
    
    call update_halo(f,f_halo,level=1)

    !vertical distance needed to calculate derivatives
    const=1._rprec/dz
    
    !operate from 0-nz uvp --> 0-nz w
    !$OMP PARALLEL DO
    do k=1,decomp%xsz(3)
       dfdz(:,:,k)=const*(f_halo(1:decomp%xsz(1),1:decomp%xsz(2),k)-f_halo(1:decomp%xsz(1),1:decomp%xsz(2),k-1))
       !uvp i take adv of level 0
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine ddz_uv_4x

  
  !!------------------------------------------------------------------------------
  !!  object: subroutine ddz_w (dfdz, f)
  !!------------------------------------------------------------------------------
  !!
  !!  last checked/modified: marco giometto ( mgiometto@gmail.com ), day 05/11/2013
  !!
  !!  description:
  !!
  !!------------------------------------------------------------------------------
  subroutine ddz_w_4x (dfdz,f,decomp)
    
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: decomp
    
    !f and dfdx are real arrays as they don't need to go through the fft
    real(rprec), dimension(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)), intent(in) :: f
    real(rprec), dimension(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)), intent(out) :: dfdz

    real(rprec), allocatable, dimension(:,:,:) :: f_halo

    real(rprec) :: const
    integer :: k
    
    call update_halo(f,f_halo,level=1)
    
    const=1.0_rprec/dz

    !$OMP PARALLEL DO
    do k=1,decomp%xsz(3)   !1 not 0 cause anyway the 0 for uvp is below comp. domain
       dfdz(:,:,k)=const*(f_halo(1:decomp%xsz(1),1:decomp%xsz(2),k+1)-f_halo(1:decomp%xsz(1),1:decomp%xsz(2),k))
    enddo
    !$OMP END PARALLEL DO
    
  end subroutine ddz_w_4x
  
  
  !!------------------------------------------------------------------------------
  !!  object: subroutine ddx (dfdx, f)
  !!------------------------------------------------------------------------------
  !!
  !!  last modified: marco giometto ( mgiometto@gmail.com ), day 05/11/2013
  !!
  !!  description:
  !!
  !!  this subroutine performs differentiation along the given direction adopting 
  !!  an out of place transform.
  !!  remember that the input of c2r is anyway destroyed (scr)
  !!
  !!------------------------------------------------------------------------------
  subroutine ddx(dfdx,f)

    implicit none
    
    real(rprec),dimension(nx),intent(in) :: f
    real(rprec),dimension(nx),intent(out) :: dfdx
    complex(rprec),dimension(lhx) :: scr
    real(rprec) :: const
    complex, parameter :: ii =(0.0,1.0)

    !normalization constant
    const=1._rprec/nx
    
    !fft to go in wavenumber space (and compute there the derivatives to have high accuracy)
    call dfftw_execute_dft_r2c(planx,f*const,scr)
    
    !kills oddballs
    scr(lhx)=0._rprec
    
    !compute coefficients for pseudospectral derivative calculation
    scr=ii*kx*scr
    
    !inverse transform to get pseudospectral derivative (values back in physical space)
    call dfftw_execute_dft_c2r(planx2,scr,dfdx)
    
    
  end subroutine ddx
  
  
  !!------------------------------------------------------------------------------
  !!  object: subroutine ddy (dfdy, f)
  !!------------------------------------------------------------------------------
  !!
  !!  last modified: marco giometto ( mgiometto@gmail.com ), day 05/11/2013
  !!
  !!  description:
  !!
  !!  this subroutine performs differentiation along the given direction adopting 
  !!  an out of place transform.
  !!  remember that the input of c2r is anyway destroyed (scr)
  !!
  !!------------------------------------------------------------------------------
  subroutine ddy(dfdy,f)
    
    implicit none
    
    real(rprec),dimension(ny),intent(in) :: f
    real(rprec),dimension(ny),intent(out) :: dfdy
    complex(rprec),dimension(lhy) :: scr
    real(rprec) :: const
    complex, parameter :: ii =(0.0,1.0)
    
    !normalization constant
    const=1._rprec/ny
    
    !fft to go in wavenumber space (and compute there the derivatives to have high accuracy)
    call dfftw_execute_dft_r2c(plany,f*const,scr)
    
    !kills oddballs
    scr(lhy)=0._rprec
    
    !compute coefficients for pseudospectral derivative calculation
    scr=ii*ky*scr
    
    !inverse transform to get pseudospectral derivative (values back in physical space)
    call dfftw_execute_dft_c2r(plany2,scr,dfdy)
        
  end subroutine ddy
  
  
end module derivatives_mod

