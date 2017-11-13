!------------------------------------------------------------------------------
!!  module: compute dynamic Smagorinski model
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------

module compute_DynSmag

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use utility_tools
  use testfilter

  !! - global variables
  use system_variables, only : u,v,w
  implicit none

  private
  TYPE(DECOMP_INFO),save :: dcp,dcpg

  public :: DynSmag_init
  public :: compute_DynSmag_coeff

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
  subroutine DynSmag_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main
            
    return
  end subroutine DynSmag_init
  
  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_DynSmag_coeff(Cs_opt2,s11,s12,s22,s33,s13,s23,delta,jt,jtt)
    implicit none
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(OUT) :: Cs_opt2
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(IN) :: s11,s12,s22,s33,s13,s23
    real(rprec), intent(IN) :: delta
    integer, intent(IN) :: jt,jtt

    real(rprec),dimension(dcp%xsz(3)+1) :: CS_Z
    integer :: k

    if(jtt.lt.dyn_init)then
       cs_opt2(:,:,:) = 0.03_rprec
    elseif((jtt.ge.dyn_init).and.(mod(jt,cs_count)==0))then
       call compute_dynamic_mum(CS_Z,S11,S12,S22,S33,S13,S23,delta)

       do k=1,dcp%xsz(3)+1
          Cs_opt2(:,:,k)=CS_Z(k)
       enddo
    endif
    
    return
  end subroutine compute_DynSmag_coeff

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_dynamic_mum(CS_Z,S11,S12,S22,S33,S13,S23,delta)
    implicit none

    real(rprec),dimension(dcp%xsz(3)+1),intent(OUT) :: CS_Z
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(IN) :: S11,S12,S22,S33,S13,S23
    real(rprec),intent(IN) :: delta

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: S,S_ft,u_ft,v_ft,w_ft
    
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: L11,L12,L13,L22,L23,L33
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: M11,M12,M13,M22,M23,M33
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: S11_ft,S12_ft,S13_ft,S22_ft,S23_ft,S33_ft
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: SS11_ft,SS12_ft,SS13_ft,SS22_ft,SS23_ft,SS33_ft

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1) :: LM,MM
    real(rprec),dimension(0:dcp%xsz(3)+1) :: num,den
    
    real(rprec) :: const
    integer :: i,j,k,k_min
    
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

    L11=u_ft*u_ft
    L12=u_ft*v_ft
    L13=u_ft*w_ft
    L22=v_ft*v_ft
    L23=v_ft*w_ft
    L33=w_ft*w_ft

    !calculating the second part for the Germano identity
    call test_filter(u_ft,G_test)
    call test_filter(v_ft,G_test)
    call test_filter(w_ft,G_test)

    !calculating the L_ij term for the Germano identity
    call test_filter(L11,G_test)
    L11 = L11 - u_ft*u_ft
    call test_filter(L12,G_test)
    L12 = L12 - u_ft*v_ft
    call test_filter(L13,G_test)
    L13 = L13 - u_ft*w_ft
    call test_filter(L22,G_test)
    L22 = L22 - v_ft*v_ft
    call test_filter(L23,G_test)
    L23 = L23 - v_ft*w_ft
    call test_filter(L33,G_test)
    L33 = L33 - w_ft*w_ft

    !calculating the |S| doubly filtered (w-nodes)!
    S = sqrt(2._rprec*(s11**2 + s22**2 + s33**2 + &
         2._rprec*(s12**2 + s13**2 + s23**2)))
       
    ! S_ij already on w-nodes
    S11_ft = S11  
    S12_ft = S12  
    S13_ft = S13  
    S22_ft = S22  
    S23_ft = S23  
    S33_ft = S33  

    !initializing the S_ij doubly filtered (w-nodes)
    call test_filter(S11_ft,G_test)
    call test_filter(S12_ft,G_test)
    call test_filter(S13_ft,G_test)
    call test_filter(S22_ft,G_test)
    call test_filter(S23_ft,G_test)
    call test_filter(S33_ft,G_test)

    !calculating the |S| doubly filtered (w-nodes)
    S_ft = sqrt(2._rprec*(S11_ft**2 + S22_ft**2 + S33_ft**2 +&
          2._rprec*(S12_ft**2 + S13_ft**2 + S23_ft**2)))

    !calculating the product of |S| and S_ij not test-filtered
    SS11_ft = S*S11
    SS12_ft = S*S12
    SS13_ft = S*S13
    SS22_ft = S*S22
    SS23_ft = S*S23
    SS33_ft = S*S33

    !test-filtering the product of |S| and S_ij (first part of M_ij)
    call test_filter(SS11_ft,G_test)
    call test_filter(SS12_ft,G_test)
    call test_filter(SS13_ft,G_test)
    call test_filter(SS22_ft,G_test)
    call test_filter(SS23_ft,G_test)
    call test_filter(SS33_ft,G_test)
    
    !create the M_ij with Smagorinsky coefficient taken away (scale similarity)
    const = 2._rprec*delta**2
    M11 = const*(SS11_ft - 4._rprec*S_ft*S11_ft)
    M12 = const*(SS12_ft - 4._rprec*S_ft*S12_ft)
    M13 = const*(SS13_ft - 4._rprec*S_ft*S13_ft)
    M22 = const*(SS22_ft - 4._rprec*S_ft*S22_ft)
    M23 = const*(SS23_ft - 4._rprec*S_ft*S23_ft)
    M33 = const*(SS33_ft - 4._rprec*S_ft*S33_ft)
    
    !$omp parallel
    !$omp do
    do k=1,dcp%xsz(3)+1
       LM(:,:,k) = L11(:,:,k)*M11(:,:,k)+L22(:,:,k)*M22(:,:,k)+L33(:,:,k)*M33(:,:,k)&
            +2._rprec*(L12(:,:,k)*M12(:,:,k)+L13(:,:,k)*M13(:,:,k)+L23(:,:,k)*M23(:,:,k))
       MM(:,:,k) = M11(:,:,k)**2+M22(:,:,k)**2+M33(:,:,k)**2+&
            2._rprec*(M12(:,:,k)**2+M13(:,:,k)**2+M23(:,:,k)**2)
    enddo
    !$omp end do
    !$omp end parallel
    
    call compute_planar_average_wghost(num,LM,dcp)
    call compute_planar_average_wghost(den,MM,dcp)

    !calculating the Smagorinsky coefficient by plane averaging
    
    !$omp parallel
    !$omp do
    do k=1,dcp%xsz(3)+1
       if(den(k) .ne. 0.0_rprec)then
          CS_Z(k)=num(k)/den(k)
       else
          CS_Z(k)=0.0_rprec
       endif
       CS_Z(k) = max(0._rprec, CS_Z(k))  
    enddo
    !$omp end do
    !$omp end parallel
    
  end subroutine compute_dynamic_mum

end module compute_DynSmag
