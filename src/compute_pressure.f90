!------------------------------------------------------------------------------
!! module: compute pressure
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 01/07/2014
!!
!! - description:
!!
!------------------------------------------------------------------------------

module compute_pressure

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use decomp_2d_fft_2d
  use fft_engine
  use utility_tools

  !! - global variables
  use system_variables, only : u,v,w,p,dpdx,dpdy,dpdz,bottom_wk_press,top_wk_press

  implicit none
  include 'mpif.h'

  private
  TYPE(DECOMP_INFO),save :: dcp,sp
  
  real(rprec),allocatable,dimension(:,:,:) :: a,b,c

  public :: press_init,press_finalize
  public :: press_stag_array

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: initialise the pressure module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 23/06/2014
  !!
  !! - description:
  !!  - kkx and kky are in global indexes, size (nx/2+1,ny)
  !!  - a,b,c upper-,main-,lower-diag of Poisson solver matrix are in local 
  !!    indexes, size (nx/2+1,ny,nz+1)
  !!
  !------------------------------------------------------------------------------
  subroutine press_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main
    integer :: status, errorcode

    dcp = dcp_main
    call decomp_2d_fft_2d_get_sp_info(sp)

    allocate(a(sp%zsz(1),sp%zsz(2),sp%zsz(3)+2), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising pressure module workspace variables')
    end if

    allocate(b(sp%zsz(1),sp%zsz(2),sp%zsz(3)+2), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising pressure module workspace variables')
    end if

    allocate(c(sp%zsz(1),sp%zsz(2),sp%zsz(3)+2), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising pressure module workspace variables')
    end if
    
    a=0._rprec
    b=0._rprec
    c=0._rprec

    call build_matrix()
    
    return
  end subroutine press_init

  !------------------------------------------------------------------------------
  !! subroutine: finalise the pressure module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 23/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine press_finalize()
    implicit none

    deallocate(a,b,c)
    
    return
  end subroutine press_finalize

  !------------------------------------------------------------------------------
  !! subroutine: compute the pressure and its gradient
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/06/2014
  !!
  !! - description:
  !!  - input: tadv1 time adv coeff
  !!  - global variables: p,dpdx,dpdy on uvp-nodes k=1,...,Nz-1
  !!  - global variables: dpdz on w-nodes k=1,...,Nz-1
  !!
  !------------------------------------------------------------------------------
  subroutine press_stag_array(tadv1)

    real(rprec),intent(IN) :: tadv1
    
    integer ::  i,j,k,k_min,k_max
    real(rprec) :: const
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1) :: rH_x,rH_y,rH_z
    complex(rprec),dimension(sp%ysz(1),sp%ysz(2),0:sp%ysz(3)+1) :: H_x,H_y,H_z
    complex(rprec),dimension(sp%ysz(1),sp%ysz(2),0:sp%ysz(3)+1) :: wk_press
    complex(rprec),dimension(lhx,ny) :: bot_layer,top_layer
    integer :: icount,root,ierror

    wk_press=(0._rprec,0._rprec)

    ! on X-pencil ========================================

    const=1._rprec/(nx*ny)

    ! x&y on uvp-nodes and z on w-nodes
    !$omp parallel do
    do k=1,dcp%xsz(3)
       rH_x(:,:,k) = const/tadv1*(u(:,:,k)/dt)
       rH_y(:,:,k) = const/tadv1*(v(:,:,k)/dt)
       rH_z(:,:,k) = const/tadv1*(w(:,:,k)/dt)
    enddo
    !$omp end parallel do

    rH_x(:,:,0) = 0._rprec
    rH_y(:,:,0) = 0._rprec
    rH_z(:,:,0) = 0._rprec
    rH_x(:,:,dcp%xsz(3)+1) = 0._rprec
    rH_y(:,:,dcp%xsz(3)+1) = 0._rprec
    rH_z(:,:,dcp%xsz(3)+1) = 0._rprec

    ! bottom and top layer
    if(dcp%xst(3)==1)then
       bottom_wk_press=const*bottom_wk_press
       call get_layer_fft_x(bot_layer,bottom_wk_press,dcp)
    endif
    icount = lhx*ny
    root = 0
    call MPI_BCAST(bot_layer,icount,MPI_DOUBLE_COMPLEX,root,MPI_COMM_WORLD,ierror)
   
    if(dcp%xen(3)==nz)then 
       call get_layer_fft_x(top_layer,const*top_wk_press,dcp)
    endif
    icount = (nx/2+1)*ny
    root = dims(2)-1
    call MPI_BCAST(top_layer,icount,MPI_DOUBLE_COMPLEX,root,MPI_COMM_WORLD,ierror)
    
    ! fft transpose data on Y-pencil
    call execute_decomp_2d_fft_2d(rH_x,H_x,.true.)
    call execute_decomp_2d_fft_2d(rH_y,H_y,.true.)
    call execute_decomp_2d_fft_2d(rH_z,H_z,.true.)    
    
    ! on Y-pencil ========================================

    ! kill all the oddballs (on x and y)
    if((lhx.ge.sp%yst(1)).and.(lhx.le.sp%yen(1)))then
       !$omp parallel do
       do k=1,sp%ysz(3)
          H_x(lhx-sp%yst(1)+1,:,k)=(0._rprec,0._rprec)
          H_y(lhx-sp%yst(1)+1,:,k)=(0._rprec,0._rprec)
          H_z(lhx-sp%yst(1)+1,:,k)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif
    if((lhy.ge.sp%yst(2)).and.(lhy.le.sp%yen(2)))then
       !$omp parallel do
       do k=1,sp%ysz(3)
          H_x(:,lhy-sp%yst(2)+1,k)=(0._rprec,0._rprec)
          H_y(:,lhy-sp%yst(2)+1,k)=(0._rprec,0._rprec)
          H_z(:,lhy-sp%yst(2)+1,k)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif
    
    call update_ghost(H_z,sp%ysz(1),sp%ysz(2),sp%ysz(3),sp)

    ! bottom and top layer
    bot_layer(lhx,:)=(0._rprec,0._rprec)
    bot_layer(:,lhy)=(0._rprec,0._rprec)
    bot_layer=bot_layer*dz

    top_layer(lhx,:)=(0._rprec,0._rprec)
    top_layer(:,lhy)=(0._rprec,0._rprec)
    top_layer=top_layer*dz

    !$omp parallel do
    do k=1,dcp%ysz(3)
       do j=1,sp%ysz(2) 
          do i=1,sp%ysz(1)
             wk_press(i,j,k) = eye * (&
                  kkx(sp%yst(1)-1+i,sp%yst(2)-1+j)*H_x(i,j,k)+&
                  kky(sp%yst(1)-1+i,sp%yst(2)-1+j)*H_y(i,j,k))+&
                  (H_z(i,j,k+1)-H_z(i,j,k))/dz
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! solve the poisson equation
    call poisson_solver(wk_press,bot_layer,top_layer)
    call update_ghost(wk_press,sp%ysz(1),sp%ysz(2),sp%ysz(3),sp)

    ! kills the oddballs
    if((lhx.ge.sp%yst(1)).and.(lhx.le.sp%yen(1)))then
       !$omp parallel do
       do k=0,sp%ysz(3)+1
          wk_press(lhx-sp%yst(1)+1,:,k)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif
    if((lhy.ge.sp%yst(2)).and.(lhy.le.sp%yen(2)))then
       !$omp parallel do
       do k=0,sp%ysz(3)+1
          wk_press(:,lhy-sp%yst(2)+1,k)=(0._rprec,0._rprec)
       enddo
       !$omp end parallel do
    endif

    ! compute the pressure gradient
    call compute_press_grad(p,dpdx,dpdy,dpdz,wk_press)

    return
  end subroutine press_stag_array
  
  !------------------------------------------------------------------------------
  !! subroutine: Poisson solver
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 25/06/2014
  !!
  !! - description:
  !!  - input&output in the same variable
  !!  - input: RHS in Y-pencil
  !!  - output: p_hat in Y-pencil
  !!
  !------------------------------------------------------------------------------
  subroutine poisson_solver(wk_press,bot_layer,top_layer)
    implicit none
    complex(rprec),dimension(sp%ysz(1),sp%ysz(2),0:sp%ysz(3)+1),intent(INOUT) :: wk_press
    complex(rprec),dimension(lhx,ny),intent(IN) :: bot_layer,top_layer
    complex(rprec),dimension(sp%zsz(1),sp%zsz(2),sp%zsz(3)+2) :: p_hat_z,RHS_col_z
    complex(rprec),dimension(lhx,ny) :: wk_layer

    integer :: i,j,k

    call transpose_y_to_z(wk_press(:,:,1:sp%ysz(3)),RHS_col_z(:,:,2:sp%zsz(3)+1),sp)
    !RHS_col_z(:,:,2:sp%zsz(3)+1)=wk_press(:,:,1:sp%ysz(3))

    ! bottom layer need special treatment
    ! -- Z-pencil ---------------------------------
    do j=1,sp%zsz(2)
       do i=1,sp%zsz(1)
          RHS_col_z(i,j,1) = bot_layer(sp%zst(1)-1+i,sp%zst(2)-1+j)
       enddo
    enddo
    
    ! top layer need special treatment
    ! -- Z-pencil --------------------------------- 
    do j=1,sp%zsz(2)
       do i=1,sp%zsz(1)
          RHS_col_z(i,j,sp%zsz(3)+2) = top_layer(sp%zst(1)-1+i,sp%zst(2)-1+j)
       enddo
    enddo

    if((sp%zst(1)==1).and.(sp%zst(2)==1))then
       RHS_col_z(1,1,1) = (0._rprec,0._rprec)
    endif

    call tridiag_solver_parallel(p_hat_z,RHS_col_z)

    call transpose_z_to_y(p_hat_z(:,:,2:sp%zsz(3)+1),wk_press(:,:,1:sp%ysz(3)),sp)
    !wk_press(:,:,1:sp%ysz(3))=p_hat_z(:,:,2:sp%zsz(3)+1)    

    !need to move p_hat_z(:,:,1) -> wk_press(:,:,0)
    call gather_layer_spz(wk_layer,p_hat_z(:,:,1),sp)

    if(dcp%xst(3)==1)then
       do j=1,sp%ysz(2) 
          do i=1,sp%ysz(1)
             wk_press(i,j,0) = wk_layer(sp%yst(1)-1+i,sp%yst(2)-1+j)
          enddo
       enddo
    endif

    return
  end subroutine poisson_solver
  
  !------------------------------------------------------------------------------
  !! subroutine: compute pressure gradient
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/06/2014
  !!
  !! - description:
  !!  - input: p_hat in Y-pencil
  !!  - output: p,dpdx,dpdy,dpdz in X-pencil
  !!  - ghostcell are up-to-date for pressure
  !!  - no ghostcell for gradient
  !!
  !------------------------------------------------------------------------------
  subroutine compute_press_grad(p,dpdx,dpdy,dpdz,p_hat)

    complex(rprec),dimension(sp%ysz(1),sp%ysz(2),0:sp%ysz(3)+1),intent(INOUT) :: p_hat
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(OUT) :: p
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)),intent(OUT) :: dpdx,dpdy,dpdz
    
    complex(rprec),dimension(sp%ysz(1),sp%ysz(2),sp%ysz(3)) ::ikxp,ikyp
    
    real(rprec) :: const
    integer i,j,k
        
    ! on Y-pencil ========================================
        
    ! compute derivatives along x&y in spectral space
    !$omp parallel do private(i,j)
    do k=1,dcp%ysz(3)
       do j=1,sp%ysz(2)
          do i=1,sp%ysz(1)
             ikxp(i,j,k)=eye*kkx(sp%yst(1)-1+i,sp%yst(2)-1+j)*p_hat(i,j,k)
             ikyp(i,j,k)=eye*kky(sp%yst(1)-1+i,sp%yst(2)-1+j)*p_hat(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! ifft transpose data on X-pencil
    call execute_decomp_2d_fft_2d(ikxp,dpdx,.false.)
    call execute_decomp_2d_fft_2d(ikyp,dpdy,.false.)
    call execute_decomp_2d_fft_2d(p_hat,p,.true.)

    ! on X-pencil ========================================
    
    ! set layers above domain
    if(dcp%xen(3)==nz)then
       p(:,:,dcp%xsz(3)+1)=BOGUS
       !p(:,:,dcp%xsz(3)+1)=p(:,:,dcp%xsz(3))
       !p(:,:,dcp%xsz(3)+1)=BOGUS
    endif

    const=1._rprec/dz

    !$omp parallel do private(i,j)
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2) 
          do i=1,dcp%xsz(1) 
             dpdz(i,j,k) = const*(p(i,j,k)-p(i,j,k-1))
          enddo
       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine compute_press_grad

  !------------------------------------------------------------------------------
  !! subroutine: tridiagonal solver
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 24/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine tridiag_solver(x,r)
    implicit none
    
    complex(rprec),dimension(sp%zsz(1),sp%zsz(2),sp%zsz(3)+2),intent(IN) :: r
    complex(rprec),dimension(sp%zsz(1),sp%zsz(2),sp%zsz(3)+2),intent(OUT):: x
  
    integer :: i,j,k
    integer :: errorcode
    real(rprec),dimension(sp%zsz(1),sp%zsz(2)) :: bet
    real(rprec),dimension(sp%zsz(1),sp%zsz(2),nz+2) :: gam
    
    ! want to skip ny/2+1 and 1, 1
    do j = 1,sp%zsz(2)
       do i = 1,sp%zsz(1)
          bet(i,j) = b(i,j,1)
          x(i,j,1) = r(i,j,1)/bet(i,j)
       enddo
    enddo
    
    ! forward step
    do k = 2,sp%zsz(3)+2
       
       do j = 1,sp%zsz(2)
          if(sp%zst(2)-1+j==lhy) cycle
          do i=1,sp%zsz(1)
             if(sp%zst(1)-1+i==lhx) cycle
             gam(i,j,k)=c(i,j,k-1)/bet(i,j)
             bet(i,j)=b(i,j,k)-a(i,j,k)*gam(i,j,k)
             if (bet(i,j) == 0._rprec) then
                write (*, *) 'tridag_array failed at (i,j,k)=',i,j,k
                write(*,*) gam(i,j,k),bet(i,j)
                write(*,*) a(i,j,k),b(i,j,k),c(i,j,k)
                errorcode = 9
                call decomp_2d_abort(errorcode, &
                     'tridiag solver: error bet(i,j)=0, dividing by 0!') 
             endif
             x(i,j,k) = (r(i,j,k)-a(i,j,k)*x(i,j,k-1))/bet(i,j)
          enddo
       enddo
    enddo

    ! backward substitution
    do k=sp%zsz(3)+1,2,-1
       do j=1,sp%zsz(2)
          if(sp%zst(2)-1+j==lhy) cycle
          do i=1,sp%zsz(1)
             If(sp%zst(1)-1+i==lhx) cycle
             x(i,j,k) = x(i,j,k)-gam(i,j,k+1)*x(i,j,k+1)          
          enddo
       enddo
    enddo

    return
  end subroutine tridiag_solver

  !------------------------------------------------------------------------------
  !! subroutine: tridiagonal solver in parallel
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 02/07/2014
  !!
  !! - description:
  !!  - i,j are independent because of Fourier decomposition
  !!  - parallelised on j not on i because of memory alignement
  !------------------------------------------------------------------------------
  subroutine tridiag_solver_parallel(x,r)
    implicit none
    
    complex(rprec),dimension(sp%zsz(1),sp%zsz(2),sp%zsz(3)+2),intent(IN) :: r
    complex(rprec),dimension(sp%zsz(1),sp%zsz(2),sp%zsz(3)+2),intent(OUT):: x
  
    integer :: i,j,k
    integer :: errorcode
    real(rprec),dimension(sp%zsz(1)) :: bet
    real(rprec),dimension(sp%zsz(1),nz+2) :: gam
    
    ! want to skip ny/2+1 and 1, 1

    !$omp parallel private(bet,gam,i,k)
    !$omp do
    do j=1,sp%zsz(2)
       if(sp%zst(2)-1+j==lhy) cycle
       do i=1,sp%zsz(1)
          if(sp%zst(1)-1+i==lhx) cycle
          bet(i) = b(i,j,1)
          x(i,j,1) = r(i,j,1)/bet(i)
       enddo
    
       ! forward step
       do k=2,sp%zsz(3)+2
          do i=1,sp%zsz(1)
             if(sp%zst(1)-1+i==lhx) cycle
             gam(i,k) = c(i,j,k-1)/bet(i)
             bet(i) = b(i,j,k)-a(i,j,k)*gam(i,k)
             x(i,j,k) = (r(i,j,k)-a(i,j,k)*x(i,j,k-1))/bet(i)
          enddo
       enddo

       ! backward substitution
       do k=sp%zsz(3)+1,1,-1
          do i=1,sp%zsz(1)
             if(sp%zst(1)-1+i==lhx) cycle
             x(i,j,k) = x(i,j,k)-gam(i,k+1)*x(i,j,k+1)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine tridiag_solver_parallel
  
  !------------------------------------------------------------------------------
  !! subroutine: build wavenumbers arrays in sprectral space
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/06/2014
  !!
  !! - description:
  !!  - kx&ky in sprectral space
  !!  - adds at nyquist freq values of kx, ky = 0 for the derivatives  
  !!  - renormalizes the wavenumbers to take into account the real domain dimention
  !!
  !------------------------------------------------------------------------------
  subroutine build_wavenumbers()
    implicit none
    
    integer i,j
    
    do i=1,lhx
       kkx(i,:) = real(i-1,kind=rprec)
    end do
    
    do j=1,ny
       kky(:,j) = real(modulo(j-1+ny/2,ny)-ny/2,kind=rprec)
    end do
    
    ! eliminates the nyquist value to calculate a correct first derivative later
    kkx(lhx,:)=0._rprec
    kky(lhx,:)=0._rprec
    
    kkx(:,lhy)=0._rprec
    kky(:,lhy)=0._rprec
    
    ! renormalization for the aspect ratio change
    kkx=2._rprec*pi/lx*kkx
    kky=2._rprec*pi/ly*kky
    
    return
  end subroutine build_wavenumbers
  
  !------------------------------------------------------------------------------
  !! subroutine: build diagonales of Poisson solver matrix
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 20/06/2014
  !!
  !! - description:
  !!  - equation (-kx^2-ky^2+dzdz) P = RHS 
  !!
  !------------------------------------------------------------------------------
  subroutine build_matrix()
    implicit none
    
    integer i,j,k,k_min,k_max
    integer errorcode
    real(rprec) bet,gam

    ! first component of Fourier its the 0 frequency
    a(:,:,1) = BOGUS !0._rprec
    b(:,:,1) = -1._rprec
    c(:,:,1) = 1._rprec

    a(:,:,sp%zsz(3)+2) = -1._rprec
    b(:,:,sp%zsz(3)+2) = 1._rprec
    c(:,:,sp%zsz(3)+2) = BOGUS!0._rprec

    !$omp parallel private(gam,bet)
    !$omp do
    do k=2,sp%zsz(3)+1
       do j=1,sp%zsz(2)
          !if(sp%yst(2)-1+j==lhy) cycle
          do i=1,sp%zsz(1)
             !if(sp%zst(1)-1+i==lhx) cycle
             a(i,j,k) = 1._rprec/(dz**2)
             b(i,j,k) = -(kkx(sp%zst(1)-1+i,sp%zst(2)-1+j)**2 & 
                  + kky(sp%zst(1)-1+i,sp%zst(2)-1+j)**2 & 
                  + 2._rprec/(dz**2))
             c(i,j,k) = 1._rprec/(dz**2)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    if((sp%zst(1)==1).and.(sp%zst(2)==1))then
       a(1, 1, 1) = 0._rprec
       b(1, 1, 1) = 1._rprec
       c(1, 1, 1) = 0._rprec
    endif

    ! = check matrix inversability ================================
    ! !$omp do
    do j=1,sp%zsz(2)
       do i=1,sp%zsz(1)        
          if (b(i,j,1)==0._rprec) then ! = ABORT ===================
             errorcode = 9
             call decomp_2d_abort(errorcode, &
                  'matric solver: error b(i,j,1)=0 while building')
          end if
       enddo
    enddo
    ! !$omp end do

    ! !$omp do
    bet=2;
    do k=2,sp%zsz(3)+1
       do j=1,sp%zsz(2)
          !if(sp%yst(2)-1+j==lhy) cycle
          do i=1,sp%zsz(1)
             !if(sp%zst(1)-1+i==lhx) cycle
             gam=c(i,j,k-1)/bet
             bet=b(i,j,k)-a(i,j,k)*gam
             if (bet == 0._rprec) then ! = ABORT ===================
                write (*, *) 'tridag solver failed at (i,j,k)=',i,j,k
                errorcode = 9
                call decomp_2d_abort(errorcode, &
                     'tridiag solver: error beta=0, dividing by 0!') 
             endif
          enddo
       enddo
    enddo
    ! !$omp end do
    ! !$omp end parallel

    return
  end subroutine build_matrix

end module compute_pressure
