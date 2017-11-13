!------------------------------------------------------------------------------
!!  module: test filter
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/05/2014
!!
!! - description: 
!!
!------------------------------------------------------------------------------
module testfilter

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use decomp_2d_fft_2d
  use fft_engine

  implicit none

  include "fftw3.f"

  private
   TYPE(DECOMP_INFO),save :: dcp,sp
  
  integer,parameter,public::filter_size=1
  real(rprec),allocatable,dimension(:,:),public :: G_test,G_test_test

  public :: test_filter_init,test_filter_finalize,test_filter_create
  public :: test_filter,test_filter_layer

contains !=======================================================================
  
  !------------------------------------------------------------------------------
  !! subroutine: initialisation of test filter kernel function
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 18/06/2014
  !!
  !! - description:
  !!  - note: this filters in-place, so input is ruined
  !!
  !------------------------------------------------------------------------------
  subroutine test_filter_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main
    call decomp_2d_fft_2d_get_sp_info(sp)
    
    allocate(G_test(lhx,ny), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising test fliter module variables')
    end if
    
    allocate(G_test_test(lhx,ny), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising test fliter module variables')
    end if

    return
  end subroutine test_filter_init

  !------------------------------------------------------------------------------
  !! subroutine: initialisation test filter module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 18/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine test_filter_finalize()
    implicit none

    deallocate(G_test)

    if(allocated(G_test_test))then
       deallocate(G_test_test)
    endif
    
    return
  end subroutine test_filter_finalize


  !------------------------------------------------------------------------------
  !! subroutine: finialise test filter module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 18/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine test_filter_create(alpha,G)
    implicit none
    real(rprec),intent(IN) :: alpha
    real(rprec),dimension(lhx,ny),intent(INOUT) :: G
    
    real(rprec) delta, kc2

    integer :: ifilter=1,model=1 ! to add to parameter

    G=1._rprec/(nx*ny)  ! normalization for the forward FFT
    delta = alpha*sqrt(dx*dy)  ! "2d-delta", not full 3d one
    
    if(ifilter==1) then ! spectral cutoff filter
       !if (model==6.OR.model==7) then
       !   print *, 'Use Gaussian or Top-hat filter for mixed models'
       !   stop
       !endif
       kc2=(pi/(delta))**2
       where (real(kk2) >= kc2) G = 0._rprec
    else if(ifilter==2) then ! Gaussian filter
       G=exp(-(2._rprec*delta)**2*kk2/(4._rprec*6._rprec))*G       
    else if(ifilter==3) then ! Top-hat (Box) filter
       G=(sin(kkx*delta/2._rprec)*sin(kky*delta/2._rprec)+1E-8)/&
            (kkx*delta/2._rprec*kky*delta/2._rprec+1E-8)*G
    endif
    
    ! since our k2 has zero at Nyquist, we have to do this by hand
    G(lhx,:) = 0._rprec
    G(:,ny/2+1) = 0._rprec

  end subroutine test_filter_create

  !------------------------------------------------------------------------------
  !! subroutine: test filter
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 18/06/2014
  !!
  !! - description:
  !!  - note: this filters in-place, so input is ruined
  !!  - update needed better management of the upper ghost cell
  !!
  !------------------------------------------------------------------------------
  subroutine test_filter(f,G)
    implicit none

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(INOUT) :: f
    real(rprec),dimension(lhx,ny),intent(IN) :: G
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1) :: tmp
    complex(rprec),dimension(sp%ysz(1),sp%ysz(2),0:sp%ysz(3)+1) :: f_hat
    integer :: i,j,k
    
    tmp=0._rprec
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1) 
             tmp(i,j,k)=f(i,j,k);
          enddo
       enddo
    enddo

    call execute_decomp_2d_fft_2d(tmp,f_hat,.true.) 

    do k=1,sp%ysz(3)+1
       do j=1,sp%ysz(2)
          do i=1,sp%ysz(1)
             f_hat(i,j,k)=f_hat(i,j,k)*G(sp%yst(1)-1+i,sp%yst(2)-1+j)
          enddo
       enddo
    enddo
    
    call execute_decomp_2d_fft_2d(f_hat,tmp,.true.)
    
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1) 
             f(i,j,k)=tmp(i,j,k)
          enddo
       enddo
    enddo

    return
  end subroutine test_filter

  !------------------------------------------------------------------------------
  !! subroutine: test filter layer
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 18/06/2014
  !!
  !! - description:
  !!  - the normalization is "in" G_test
  !!  - note: this filters in-place, so input is ruined
  !!
  !------------------------------------------------------------------------------
  subroutine test_filter_layer(f,G)
    implicit none

    real(rprec), dimension(nx,ny),intent(INOUT) :: f
    real(rprec), dimension(lhx,ny),intent(IN) :: G

    complex(rprec),dimension(lhx,ny) :: scr
    
    call dfftw_execute_dft_r2c(plan2d_f,f,scr)
    scr = G*scr
    call dfftw_execute_dft_c2r(plan2d_b,scr,f)
    
    return
  end subroutine test_filter_layer

end module testfilter
