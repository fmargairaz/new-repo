!------------------------------------------------------------------------------
!!  module: compute Lagrangian Scale Dependent dynamic model
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------

module compute_LagScaleDep

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use utility_tools
  use testfilter

  !! - global variables
  use system_variables, only : u,v,w
  implicit none

  private
  TYPE(DECOMP_INFO),save :: dcp

  logical,parameter :: DEBUG=.false.
  real(rprec),parameter :: tf1=2._rprec,tf1_2=4._rprec
  real(rprec),parameter :: tf2=4._rprec,tf2_2=16._rprec
  real(rprec),parameter :: opftime=1.5_rprec

  logical :: F_coeff_init=.false.

  real(rprec),allocatable,dimension(:,:,:),public :: F_LM,F_MM,F_QN,F_NN
  real(rprec),allocatable,dimension(:,:,:),public :: Beta_mum

  real(rprec),dimension(:,:,:),allocatable :: u_lag,v_lag,w_lag

  real(rprec),dimension(:,:,:),allocatable :: S,S_bar,S_hat
  real(rprec),dimension(:,:,:),allocatable :: u_bar,v_bar,w_bar,u_hat,v_hat,w_hat
  
  real(rprec),parameter::n_dmp=2._rprec

  public :: LagScaleDep_init
  public :: compute_LagScaleDep_coeff

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: initialise sgs module
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine LagScaleDep_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main

    if(Lagran_init==.true.) F_coeff_init=.true.

    ! F_LM,F_MM,F_QN,F_NN,Beta_mum
    allocate(F_LM(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(F_MM(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(F_QN(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(F_NN(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(Beta_mum(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if

    ! u_lag,v_lag,w_lag
    allocate(u_lag(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(v_lag(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(w_lag(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    u_lag = 0._rprec;v_lag = 0._rprec;w_lag = 0._rprec

    ! S,S_bar,S_hat
    allocate(S(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(S_bar(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(S_hat(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if

    !u_bar,v_bar,w_bar,u_hat,v_hat,w_hat
    allocate(u_bar(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(v_bar(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(w_bar(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(u_hat(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(v_hat(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    allocate(w_hat(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising LagScDep module variables')
    end if
    
    return
  end subroutine LagScaleDep_init

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  subroutine LagScaleDep_finalize()
    implicit none

    deallocate(F_LM,F_MM,F_QN,F_NN,Beta_mum)
    deallocate(u_lag,v_lag,w_lag)
    deallocate(S,S_bar,S_hat,u_bar,v_bar,w_bar,u_hat,v_hat,w_hat)
  
    return
  end subroutine LagScaleDep_finalize
  
  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_LagScaleDep_coeff(Cs_opt2,s11,s12,s22,s33,s13,s23,delta,jt,jtt)
    implicit none
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(OUT) :: Cs_opt2
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(IN) :: s11,s12,s22,s33,s13,s23
    real(rprec), intent(IN) :: delta
    integer, intent(IN) :: jt,jtt

    integer,dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: iaddx,jaddy,kaddz
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: comp_x,comp_y,comp_z
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: LM,MM,QN,NN
    
    real(rprec) :: powcoeff
    integer :: i,j,k

    if(jtt.ge.dyn_init)then
       if(DEBUG) write(*,*) nrank,jtt,'add to u_lag'
       !$omp parallel do
       do k=0,dcp%xsz(3)+1
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1)
                u_lag(i,j,k)=u_lag(i,j,k)+u(i,j,k)
                v_lag(i,j,k)=v_lag(i,j,k)+v(i,j,k)
                w_lag(i,j,k)=w_lag(i,j,k)+w(i,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do       
    endif
    
    if(jtt.lt.dyn_init)then
       cs_opt2(:,:,:) = 0.03_rprec
    elseif((jtt.ge.dyn_init).and.(mod(jtt,cs_count)==0))then
       if(DEBUG) write(*,*) nrank,jtt,'compute_cs'
       !$omp parallel do
       do k=1,dcp%xsz(3)+1
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1) 
                S(i,j,k) = sqrt(2._rprec*(s11(i,j,k)**2 + s22(i,j,k)**2 + s33(i,j,k)**2 + &
                     2._rprec*(s12(i,j,k)**2 + s13(i,j,k)**2 + s23(i,j,k)**2)))
             enddo
          enddo
       enddo
       !$omp end parallel do

       ! 2Delta - LM and MM
       call compute_dynamic_mum(LM,MM,S_bar,u_bar,v_bar,w_bar,S,S11,S12,S22,S33,S13,S23,tf1_2,delta,G_test)
       ! 4Delta - QN and NN
       call compute_dynamic_mum(QN,NN,S_hat,u_hat,v_hat,w_hat,S,S11,S12,S22,S33,S13,S23,tf2_2,delta,G_test_test)

       if(F_coeff_init==.true.)then
          if(DEBUG) write(*,*) nrank,jtt,'compute F_coeff'
          ! compute displacement and weigh for Lagrangian interpolation.
          call compute_interpolag_weight(iaddx,jaddy,kaddz,comp_x,comp_y,comp_z) 
          ! compute F_LM,F_MM,F_QN,F_NN
          call compute_interpolag_coeff(F_LM,iaddx,jaddy,kaddz,comp_x,comp_y,comp_z)
          call compute_interpolag_coeff(F_MM,iaddx,jaddy,kaddz,comp_x,comp_y,comp_z)
          call compute_interpolag_coeff(F_QN,iaddx,jaddy,kaddz,comp_x,comp_y,comp_z)
          call compute_interpolag_coeff(F_NN,iaddx,jaddy,kaddz,comp_x,comp_y,comp_z)
       else
          if(DEBUG) write(*,*) nrank,jtt,'initialise F_coeff'
          !$omp parallel do
          do k=1,dcp%xsz(3)+1
             do j=1,dcp%xsz(2)
                do i=1,dcp%xsz(1)
                   F_MM(i,j,k)=MM(i,j,k);F_LM(i,j,k)=0.03_rprec*MM(i,j,k)
                   F_NN(i,j,k)=NN(i,j,k);F_QN(i,j,k)=0.03_rprec*NN(i,j,k)
                enddo
             enddo
          enddo
          !$omp end parallel do
          F_coeff_init=.true.
       end if
       
       powcoeff = -1._rprec/8._rprec !for momentum
       call update_interpolag_coeff(F_LM,F_MM,LM,MM,powcoeff,delta) ! 2Delta - LM and MM
       call update_interpolag_coeff(F_QN,F_NN,QN,NN,powcoeff,delta) ! 4Delta - QN and NN
       
       call compute_dynamic_coeff(Cs_opt2,Beta_mum,F_LM,F_MM,F_QN,F_NN) !compute Cs2

       !write(*,*) nrank,Cs_opt2(1,1,:)

       ! scalars
       if(Lagran_init==.false.) Lagran_init=.true.

       ! reset the veoliciy averages
       if(DEBUG) write(*,*) nrank,jtt,'u_lag zeroed'
       u_lag=0._rprec;v_lag=0._rprec;w_lag=0._rprec
    else
       !nothing to do
    endif
    
    return
  end subroutine compute_LagScaleDep_coeff
  
  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_interpolag_weight(iaddx,jaddy,kaddz,comp_x,comp_y,comp_z)
    implicit none
    
    integer,dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(OUT) :: iaddx,jaddy,kaddz
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(OUT) :: comp_x,comp_y,comp_z
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: xp,yp,zp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1) :: tmp

    integer :: i,j,k,k_min
    integer :: addx,addy,addz


    ! interpolate u_lag on w_nodes
    tmp = u_lag/real(cs_count,rprec)    
    if(dcp%xst(3)==1)then
       u_lag(:,:,1) = tmp(:,:,1)
       k_min = 2
    else
       k_min = 1
    endif
    do k=k_min,dcp%xsz(3)+1
       u_lag (:,:,k) = 0.5_rprec*(tmp(:,:,k)+tmp(:,:,k-1))
    enddo

    ! interpolate v_lag on w_nodes
    tmp = v_lag/real(cs_count,rprec)
    if(dcp%xst(3)==1)then
       v_lag(:,:,1) = tmp(:,:,1)
       k_min = 2
    else
       k_min = 1
    endif
    do k=k_min,dcp%xsz(3)+1
       v_lag(:,:,k) = 0.5_rprec*(tmp(:,:,k)+tmp(:,:,k-1))
    enddo

    w_lag = w_lag/real(cs_count,rprec)
    if(dcp%xst(3)==1)then
       w_lag(:,:,1) = 0.25_rprec*w_lag(:,:,2)
    end if

    ! Computes the 3-D inverse displacement arrays that describe the location 
    !  where the point was at the previous step
    !$omp parallel do
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             xp(i,j,k) = -u_lag(i,j,k)*dt*real(cs_count,rprec)/dx  
             yp(i,j,k) = -v_lag(i,j,k)*dt*real(cs_count,rprec)/dy   
             zp(i,j,k) = -w_lag(i,j,k)*dt*real(cs_count,rprec)/dz
          enddo
       enddo
    enddo
    !$omp end parallel do

    ! Because first plane is on u,v,p nodes this corrects for the fact that the first cell
    !  in the z direction has height dz/2 it doubles the zp fraction if this fraction relates 
    !  to the cell beneath it
    if(dcp%xst(3)==1)then
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1) 
             zp(i,j,2) = zp(i,j,2) + min(zp(i,j,2),0._rprec)
             zp(i,j,1) = zp(i,j,1) + max(zp(i,j,1),0._rprec)
          end do
       end do
       zp(:,:,2) = min(1._rprec,zp(:,:,2))
       zp(:,:,1) = min(1._rprec,zp(:,:,2))
       zp(:,:,2) = max(-1._rprec,zp(:,:,2))
       zp(:,:,1) = max(-1._rprec,zp(:,:,2))
    end if

    !if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
    !   if(mod(jt,p_count).eq.0) print*,'Lagrangian CFL condition= ',  &
    !        maxval(abs(xp(1:nx, :, 1:nz)))
    !end if
    
    ! The are the values to add to the indices i,j,k= +1 or -1 depending on what cube 
    !  should be used for interpolation.
    ! -> Computes the relative weights given to F_** in the cube depending on point location
    
    !$omp parallel do private(addx,addy,addz)
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             ! the are the values to add to the indices i,j,k= +1 or -1
             ! depending on what cube should be used for interpolation
             addx = int(sign(1._rprec,xp(i,j,k)))
             addy = int(sign(1._rprec,yp(i,j,k)))
             addz = int(sign(1._rprec,zp(i,j,k)))
             iaddx(i,j,k) = i + addx 
             jaddy(i,j,k) = j + addy
             kaddz(i,j,k) = k + addz
             comp_x(i,j,k) = abs(xp(i,j,k))
             comp_y(i,j,k) = abs(yp(i,j,k))
             comp_z(i,j,k) = abs(zp(i,j,k))
          enddo
       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine compute_interpolag_weight
    
  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!    The interpolation works in such a way that takes the values of the F_** function
  !!    around the actual point (i,j,k), and multiplies it with a weight function. This weight function
  !!    is the result of integrating over a streamline considering Taylor's Hypothesis. (Depending on how
  !!    far is your previous function value, the weight is bigger or smaller. The weight is based in spatial
  !!    distance).
  !!
  !------------------------------------------------------------------------------
  subroutine compute_interpolag_coeff(F,iaddx,jaddy,kaddz,comp_x,comp_y,comp_z)
    implicit none

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(INOUT) :: F
    integer,dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(IN) :: iaddx,jaddy,kaddz
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(IN) :: comp_x,comp_y,comp_z

    real(rprec),allocatable,dimension(:,:,:) :: FF

    integer :: i,j,k,ax,ay,az
    real(rprec) :: cx,cy,cz,fx,fy,fz
    
    call update_halo3(F,FF,2,dcp)

    if(dcp%xst(3)==1)then
       FF(:,:,0) = FF(:,:,1)
       FF(:,:,-1) = BOGUS
    endif
    if(dcp%xen(3)==nz)then
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             FF(i,j,dcp%xsz(3)+1) = F(i,j,dcp%xsz(3)+1)
             FF(i,j,dcp%xsz(3)+2) = F(i,j,dcp%xsz(3)+1)
          enddo
       enddo
    endif
    
    !$omp parallel do private(ax,ay,az,cx,cy,cz,fx,fy,fz)
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             
             ax=iaddx(i,j,k);ay=jaddy(i,j,k);az=kaddz(i,j,k);
             cx=comp_x(i,j,k);cy=comp_y(i,j,k);cz=comp_z(i,j,k);
             fx=1._rprec-cx;fy=1._rprec-cy;fz=1._rprec-cz;

             F(i,j,k)=fx*fy*(FF(i,j,k)*fz+FF(i,j,az)*cz) + fx*cy*(FF(i,ay,k)*fz+FF(i,ay,az)*cz)&
                   + cx*fy*(FF(ax,j,k)*fz+FF(ax,j,az)*cz) + cx*cy*(FF(ax,ay,k)*fz+FF(ax,ay,az)*cz)

          enddo
       enddo
    enddo
    !$omp end parallel do

    

    return
  end subroutine compute_interpolag_coeff

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_dynamic_mum(AB,BB,S_ft,u_ft,v_ft,w_ft,S,S11,S12,S22,S33,S13,S23,tf_2,delta,G)
    implicit none

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(OUT) :: AB,BB,S_ft,u_ft,v_ft,w_ft
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(IN) :: S,S11,S12,S22,S33,S13,S23
    real(rprec),intent(IN) :: tf_2,delta
    real(rprec),dimension(:,:),intent(IN) :: G

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: A11,A12,A13,A22,A23,A33
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: B11,B12,B13,B22,B23,B33
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: S11_ft,S12_ft,S13_ft,S22_ft,S23_ft,S33_ft
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: SS11_ft,SS12_ft,SS13_ft,SS22_ft,SS23_ft,SS33_ft
    
    real(rprec) :: const
    integer :: i,j,k,k_min

    const = 2._rprec*delta**2

    if(dcp%xst(3)==1)then
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             u_ft(i,j,1) = u(i,j,1)! first uv-node
             v_ft(i,j,1) = v(i,j,1)! first uv-node
             w_ft(i,j,1) = 0.25_rprec*w(i,j,2)! first uv-node  
          enddo
       enddo
       k_min=2
    else
       k_min=1
    endif
    
    !$omp parallel
    !$omp do
    do k=k_min,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             u_ft(i,j,k) = 0.5_rprec*(u(i,j,k) + u(i,j,k-1))
             v_ft(i,j,k) = 0.5_rprec*(v(i,j,k) + v(i,j,k-1))
             w_ft(i,j,k) = w(i,j,k)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    A11=u_ft*u_ft
    A12=u_ft*v_ft
    A13=u_ft*w_ft
    A23=v_ft*w_ft
    A22=v_ft*v_ft
    A33=w_ft*w_ft

    call test_filter(u_ft,G)
    call test_filter(v_ft,G)
    call test_filter(w_ft,G)
    
    call test_filter(A11,G)
    A11 = A11 - u_ft*u_ft
    call test_filter(A12,G)
    A12 = A12 - u_ft*v_ft
    call test_filter(A13,G)
    A13 = A13 - u_ft*w_ft
    call test_filter(A22,G)
    A22 = A22 - v_ft*v_ft
    call test_filter(A23,G)
    A23 = A23 - v_ft*w_ft
    call test_filter(A33,G)
    A33 = A33 - w_ft*w_ft

    ! S_ij already on w-nodes
    S11_ft = S11  
    S12_ft = S12  
    S13_ft = S13  
    S22_ft = S22  
    S23_ft = S23  
    S33_ft = S33  

    call test_filter(S11_ft,G)
    call test_filter(S12_ft,G)
    call test_filter(S13_ft,G)
    call test_filter(S22_ft,G)
    call test_filter(S23_ft,G)
    call test_filter(S33_ft,G)

    S_ft = sqrt(2._rprec*(S11_ft**2 + S22_ft**2 + S33_ft**2 +&
          2._rprec*(S12_ft**2 + S13_ft**2 + S23_ft**2)))

    SS11_ft = S*S11
    SS12_ft = S*S12
    SS13_ft = S*S13
    SS22_ft = S*S22
    SS23_ft = S*S23
    SS33_ft = S*S33

    call test_filter(SS11_ft,G)
    call test_filter(SS12_ft,G)
    call test_filter(SS13_ft,G)
    call test_filter(SS22_ft,G)
    call test_filter(SS23_ft,G)
    call test_filter(SS33_ft,G)     

    B11 = const*(SS11_ft - tf_2*S_ft*S11_ft)
    B12 = const*(SS12_ft - tf_2*S_ft*S12_ft)
    B13 = const*(SS13_ft - tf_2*S_ft*S13_ft)
    B22 = const*(SS22_ft - tf_2*S_ft*S22_ft)
    B23 = const*(SS23_ft - tf_2*S_ft*S23_ft)
    B33 = const*(SS33_ft - tf_2*S_ft*S33_ft)

    AB = A11*B11+A22*B22+A33*B33+2._rprec*(A12*B12+A13*B13+A23*B23)
    BB = B11**2+B22**2+B33**2+2._rprec*(B12**2+B13**2+B23**2)

  end subroutine compute_dynamic_mum

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine update_interpolag_coeff(F_AB,F_BB,AB,BB,powcoeff,delta)
    implicit none

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(INOUT) :: F_AB,F_BB
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(IN) :: AB,BB
    real(rprec),intent(IN) :: powcoeff,delta

    real(rprec) :: opftdelta,lagran_dt,tmp,dumfac,epsi
    integer :: i,j,k

    opftdelta = opftime*delta
    lagran_dt=dt*real(cs_count,kind=rprec)

    !$omp parallel do private(tmp,dumfac,epsi)
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)   
             tmp=max(real(F_AB(i,j,k)*F_BB(i,j,k)),real(1E-24))
             tmp=opftdelta*tmp**powcoeff
             tmp=max(real(1E-24),real(tmp))
             dumfac=lagran_dt/tmp
             epsi=dumfac/(1._rprec+dumfac)

             F_AB(i,j,k)=(epsi*AB(i,j,k) + (1._rprec-epsi)*F_AB(i,j,k))
             F_BB(i,j,k)=(epsi*BB(i,j,k) + (1._rprec-epsi)*F_BB(i,j,k))
    
             F_AB(i,j,k)= max(real(1E-24),real(F_AB(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine update_interpolag_coeff

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_dynamic_coeff(Cs_,Beta_,F_AB,F_BB,F_CD,F_DD)
    implicit none
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(OUT) :: Beta_
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(OUT) :: Cs_
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(IN) :: F_AB,F_BB,F_CD,F_DD
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: Cs_opt2_2d,Cs_opt2_4d
    real(rprec),dimension(0:dcp%xsz(3)+1) :: cs_pavg
    real(rprec) :: Betaclip
    integer :: i,j,k

    !$omp parallel do
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             Cs_opt2_2d(i,j,k)=F_AB(i,j,k)/F_BB(i,j,k)
             Cs_opt2_2d(i,j,k)=max(real(1E-24),real(Cs_opt2_2d(i,j,k)))

             Cs_opt2_4d(i,j,k)=F_CD(i,j,k)/F_DD(i,j,k)
             Cs_opt2_4d(i,j,k)=max(real(1E-24),real(Cs_opt2_4d(i,j,k)))
             
             Beta_(i,j,k)=(Cs_opt2_4d(i,j,k)/Cs_opt2_2d(i,j,k))**(log(tf1)/(log(tf2)-log(tf1)))
             Beta_(i,j,k)=max(real(Beta_(i,j,k)),real(1._rprec/(tf1*tf2)))
          enddo
       enddo
    enddo
    !$omp end parallel do

    
    if(dcp%xen(3)==nz)then
       Beta_(:,:,dcp%xsz(3)+1)=1._rprec 
    end if
    
    !$omp parallel do
    do k=1,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             Cs_(i,j,k)=Cs_opt2_2d(i,j,k)/Beta_(i,j,k)
             Cs_(i,j,k)=max(real(1E-24),real(Cs_(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel do
    if(DEBUG)then
       call compute_planar_average_wghost(cs_pavg,CS_,dcp)
       if(dcp%xst(2)==1)then
          do k=1,dcp%xsz(3)+1
             write(*,*) nrank,dcp%xst(3)+k-1,sqrt(cs_pavg(k))
          enddo
       endif
    endif
    !Beta_avg(jz)=sum(Cs_opt2_2d(1:nx,1:ny,jz))/sum(Cs_opt2(1:nx,1:ny,jz))
    !Betaclip_avg(jz)=sum(Cs_opt2_4d(1:nx,1:ny,jz))/sum(Cs_opt2_2d(1:nx,1:ny,jz))
    
    return
  end subroutine compute_dynamic_coeff

end module compute_LagScaleDep
