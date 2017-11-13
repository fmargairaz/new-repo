!------------------------------------------------------------------------------
!!  module: compute sub grid stress 
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!  - this module contains all the routines used for computing sgs
!!  - cs_opt2 & nu_t are definded on w-nodes
!!
!------------------------------------------------------------------------------

module compute_sgs

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom

  !! - global variables
  use system_variables, only : dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,&
       txx,tyy,tzz,txy,txz,tyz

  use compute_DynSmag
  use compute_LagScaleDep

  implicit none

  private
  TYPE(DECOMP_INFO),save :: dcp

  real(rprec),dimension(:,:,:),allocatable,public :: cs_opt2,nu_t
  real(rprec),parameter::n_dmp=2._rprec

  public :: sgs_init,sgs_finalize
  public :: sgs_stress

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
  subroutine sgs_init(dcp_main)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main

    allocate(cs_opt2(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising sgs module variables')
    end if

    allocate(nu_t(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising sgs module variables')
    end if

    cs_opt2=0._rprec
    nu_t=0._rprec

    if(sgs_model .eq. "dyn_Smag")then
       call DynSmag_init(dcp)
    elseif(sgs_model .eq. "dyn_LagScaleDep")then
       call LagScaleDep_init(dcp)
    end if

    return
  end subroutine sgs_init

  !------------------------------------------------------------------------------
  !! subroutine: finalise sgs module
  !------------------------------------------------------------------------------=
  subroutine sgs_finalize()

    implicit none

    deallocate(cs_opt2, nu_t)

  end subroutine sgs_finalize

  !------------------------------------------------------------------------------
  !! subroutine: compute sgs stress
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - txx,txy,tyy,tzz on uvp-nodes except at bottom
  !!  - txz,tyz on w-nodes
  !!
  !------------------------------------------------------------------------------
  subroutine sgs_stress(jt,jtt)
    implicit none

    integer, intent(IN):: jt,jtt

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1) :: s11,s12,s22,s33,s13,s23

    real(rprec),dimension(dcp%xsz(3)+1) :: l

    real(rprec) :: delta,stmp,filter_size,nu_m,z
    integer :: i,j,k,k_min

    real(rprec) :: const
    real(rprec),save :: lagr_dt=0.0_rprec
    logical :: check=.true.

    ! molecular viscosity
    nu_m=1._rprec/Re

    if(sgs_model .eq. "dns")then
       const=2._rprec/Re
       !$omp parallel do
       do k=1,dcp%xsz(3)
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1)
                txx(i,j,k)=const*dudx(i,j,k)
                txy(i,j,k)=const*0.5_rprec*(dudy(i,j,k)+dvdx(i,j,k))
                tyy(i,j,k)=const*dvdy(i,j,k)
                tzz(i,j,k)=const*dwdz(i,j,k)

                txz(i,j,k)=const*0.5_rprec*(dudz(i,j,k)+dwdx(i,j,k))
                tyz(i,j,k)=const*0.5_rprec*(dvdz(i,j,k)+dwdy(i,j,k))
             enddo
          enddo
       enddo
       !$omp end parallel do       
    else
       ! filter size - nondimensional
       filter_size=1._rprec
       delta=filter_size*(dx*dy*dz)**(1._rprec/3._rprec)  
       l(:)=delta

       call compute_sij(dcp,s11,s12,s13,s23,s22,s33)

       if(sgs_model .eq. "cst_Smag")then
          ! wall damping
          if(dcp%xst(3)==1)then
             z=0.5_rprec*dz
             l(1)=(Cs0**(n_dmp)*(KvonK*z)**(-n_dmp)+(delta)**(-n_dmp))**(-1._rprec/n_dmp)
             k_min = 2
          else
             k_min = 1
          endif

          do k=k_min,dcp%xsz(3)+1
             z=(dcp%xst(3)+k-2._rprec)*dz
             ! z's nondimensional, l here is on w-nodes, except at bottom
             l(k)=(Cs0**(n_dmp)*(KvonK*z)**(-n_dmp)+(delta)**(-n_dmp))**(-1._rprec/n_dmp)
          end do
          !define mixing lenght
          ! TODO if not using CstSmag & wallfunction l(:)=delta

          ! Smagorinsky coefficient
          cs_opt2(:,:,:) = Cs0**2
       elseif(sgs_model .eq. "dyn_Smag")then
          call compute_DynSmag_coeff(Cs_opt2,s11,s12,s22,s33,s13,s23,delta,jt,jtt)
       elseif(sgs_model .eq. "dyn_LagScaleDep")then
          call compute_LagScaleDep_coeff(Cs_opt2,s11,s12,s22,s33,s13,s23,delta,jt,jtt)
       else
          call decomp_2d_abort(1,'Invalid SGS model')
       endif
       ! eddy viscosity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !$omp parallel do private(stmp)
       do k=1,dcp%xsz(3)+1
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1) 
                stmp = sqrt(2._rprec*(s11(i,j,k)**2 + s22(i,j,k)**2 + s33(i,j,k)**2 + &
                     2._rprec*(s12(i,j,k)**2 + s13(i,j,k)**2 + s23(i,j,k)**2)))
                !smagorinsky hypothesis: evaluate nu_t = c_s^2 l^2 |s|
                nu_t(i,j,k)=stmp*cs_opt2(i,j,k)*(l(k)**2)
             enddo
          enddo
       enddo
       !$omp end parallel do

       ! updata ghost cells for eddy viscosity
       !call update_ghost(nu_t,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)

       ! sgs stress !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! special treatment for bottom layer - compute stresses on uvp-nodes 
       if(dcp%xst(3)==1)then
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1)
                txx(i,j,1)=2._rprec*nu_t(i,j,1)*dudx(i,j,1)
                txy(i,j,1)=2._rprec*nu_t(i,j,1)*0.5_rprec*(dudy(i,j,1)+dvdx(i,j,1))
                tyy(i,j,1)=2._rprec*nu_t(i,j,1)*dvdy(i,j,1)
                tzz(i,j,1)=2._rprec*nu_t(i,j,1)*dwdz(i,j,1)
             enddo
          enddo
          k_min=2
       else
          k_min=1
       endif

       !$omp parallel private(const)

       ! here: t_ij & du_i/dx_j on uvp-nodes, nu_t on w-nodes
       !$omp do
       do k=k_min,dcp%xsz(3)
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1)
                const=0.5_rprec*(nu_t(i,j,k)+nu_t(i,j,k+1))
                txx(i,j,k)=2.0_rprec*const*dudx(i,j,k)     
                txy(i,j,k)=2.0_rprec*const*0.5_rprec*(dudy(i,j,k)+dvdx(i,j,k)) 
                tyy(i,j,k)=2.0_rprec*const*dvdy(i,j,k)
                tzz(i,j,k)=2.0_rprec*const*dwdz(i,j,k)
             enddo
          enddo
       enddo
       !$omp end do

       ! here: t_ij & s_ij on w_nodes
       !$omp do
       do k=k_min,dcp%xsz(3)+1
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1)
                txz(i,j,k)=2.0_rprec*nu_t(i,j,k)*s13(i,j,k)
                tyz(i,j,k)=2.0_rprec*nu_t(i,j,k)*s23(i,j,k)
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel

       ! at wall: assume dz(tzz)=0.0 (imposed using ghost cell below domain)
       if(dcp%xst(3)==1)then
          tzz(:,:,0)=tzz(:,:,1)
       endif

    endif

    return
  end subroutine sgs_stress

  !------------------------------------------------------------------------------
  !! subroutine: compute sij
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - sij do not need ghost cells
  !!  - computes sij on uvp node at jz=1, on w nodes for the others (jz>1)
  !!
  !------------------------------------------------------------------------------
  subroutine compute_sij(dcp,s11,s12,s13,s23,s22,s33)

    implicit none
    type(decomp_info),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1),intent(out) ::  s11,s12,s13,s23,s22,s33

    real(rprec) :: ux,uy,uz,vx,vy,vz,wx,wy,wz

    integer :: i,j,k,k_min

    ! first level has everything on uvp nodes
    if(dcp%xst(3)==1)then
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1) 

             !uvp-node @ dz/2
             ux=dudx(i,j,1)
             uy=dudy(i,j,1)
             if(lbc%special.eq.'wall_law')then
                uz=dudz(i,j,1)
                vz=dvdz(i,j,1)
             else
                uz=0.5_rprec*(dudz(i,j,1)+dudz(i,j,2))
                vz=0.5_rprec*(dvdz(i,j,1)+dvdz(i,j,2))
             endif
             vx=dvdx(i,j,1)                            
             vy=dvdy(i,j,1)       

             !uvp-node so i need to interp
             wx=0.5_rprec*(dwdx(i,j,1)+dwdx(i,j,2))   
             wy=0.5_rprec*(dwdy(i,j,1)+dwdy(i,j,2))
             wz=dwdz(i,j,1)

             !uvp-node @ dz/2
             s11(i,j,1)=ux                             
             s12(i,j,1)=0.5_rprec*(uy+vx)              
             s13(i,j,1)=0.5_rprec*(uz+wx)              
             s22(i,j,1)=vy                             
             s23(i,j,1)=0.5_rprec*(vz+wy)              
             s33(i,j,1)=wz                             

          enddo
       enddo
       k_min=2
    else
       k_min=1
    endif

    ! standard nodes (put everything on w, where Cs and nu_t are defined)
    !$omp parallel do private(ux,uy,uz,vx,vy,vz,wx,wy,wz)
    do k=k_min,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             ux=0.5_rprec*(dudx(i,j,k) + dudx(i,j,k-1))     
             uy=0.5_rprec*(dudy(i,j,k) + dudy(i,j,k-1))     
             uz=dudz(i,j,k)                                 
             vx=0.5_rprec*(dvdx(i,j,k) + dvdx(i,j,k-1))     
             vy=0.5_rprec*(dvdy(i,j,k) + dvdy(i,j,k-1))     
             vz=dvdz(i,j,k)                                 
             wx=dwdx(i,j,k)                                 
             wy=dwdy(i,j,k)                                 
             wz=0.5_rprec*(dwdz(i,j,k) + dwdz(i,j,k-1))     

             s11(i,j,k)=ux                                  
             s12(i,j,k)=0.5_rprec*(uy+vx)                   
             s13(i,j,k)=0.5_rprec*(uz+wx)                   
             s22(i,j,k)=vy                                  
             s23(i,j,k)=0.5_rprec*(vz+wy)                   
             s33(i,j,k)=wz                                  
          enddo
       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine compute_sij

end module compute_sgs
