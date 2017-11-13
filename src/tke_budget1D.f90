!------------------------------------------------------------------------------
!! module: 
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 21/05/2015
!!
!! - description:
!!
!------------------------------------------------------------------------------

module tke_budget1D
  use ifport
  
  use parameters_IO
  use decomp_2d
  use decomp_2d_custom

  !! - global variables
  use system_variables
  
  use compute_derivatives
  use compute_sgs
  use utility_tools
  
  implicit none
  include 'mpif.h'

  private
  integer,save :: ii
  character(256) :: oPath 
  logical, parameter :: s_flag=.false. ! to change in the future

  real(rprec),dimension(:),allocatable :: au,av,aw,ap,auu,avv,aww,auv,auw,avw,&
       adudx,adudy,adudz,advdx,advdy,advdz,adwdx,adwdy,adwdz,&
       adudx2,adudy2,adudz2,advdx2,advdy2,advdz2,adwdx2,adwdy2,adwdz2,  &
       atxx,atyy,atzz,atxy,atxz,atyz,&
       acs,abeta,&
       auuw,avvw,awww,&
       apw,apu,&
       adpdx,adpdy,adpdz,apdudx,apdvdy,apdwdz,&
       autxz,avtyz,awtzz,&
       adxx,adyy,adzz,adxy,adxz,adyz,&
       apdudz,apdwdx,auww,&
       adudxdwdx,adudydwdy,adudzdwdz

  real(rprec),dimension(:),allocatable :: as,ass,aus,avs,aws,aps,atsx,atsy,atsz,&
       adsdx,adsdy,adsdz,ads,abetas,apdsdz,asww,adsdxdwdx,adsdydwdy,adsdzdwdz
 
  public :: tke_budget_init, tke_budget_finalize
  public :: compute_tke_budget


contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: initialise tke budget module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 17/05/2015
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine tke_budget_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main
    integer :: actLaunch=1
    integer :: result,status, errorcode
    

    dcp = dcp_main
    ii=(actLaunch-1)*(nt/tke_period)+1
    oPath=trim(out_path)//'ta1_profiles/' 

    if(nrank==0)then
       !create output folder
       write(*,*) 'Creating folder for tke budget output...'
       result=system('mkdir '//trim(oPath))
       if(result.eq.0)then
          write(*,*) 'Created output folder!'
       else
          write(*,*) 'Output already exists!'
          write(*,*) 'Result=',Result
      endif
   endif

    allocate(au(dcp%xsz(3)));au=0._rprec;
    allocate(av(dcp%xsz(3)));av=0._rprec;
    allocate(aw(dcp%xsz(3)));aw=0._rprec;
    allocate(ap(dcp%xsz(3)));ap=0._rprec;

    allocate(auu(dcp%xsz(3)));auu=0._rprec;
    allocate(avv(dcp%xsz(3)));avv=0._rprec;
    allocate(aww(dcp%xsz(3)));aww=0._rprec;
    allocate(auv(dcp%xsz(3)));auv=0._rprec;
    allocate(auw(dcp%xsz(3)));auw=0._rprec;
    allocate(avw(dcp%xsz(3)));avw=0._rprec;

    allocate(adudx(dcp%xsz(3)));adudx=0._rprec;
    allocate(adudy(dcp%xsz(3)));adudy=0._rprec;
    allocate(adudz(dcp%xsz(3)));adudz=0._rprec;
    allocate(advdx(dcp%xsz(3)));advdx=0._rprec;
    allocate(advdy(dcp%xsz(3)));advdy=0._rprec;
    allocate(advdz(dcp%xsz(3)));advdz=0._rprec;
    allocate(adwdx(dcp%xsz(3)));adwdx=0._rprec;
    allocate(adwdy(dcp%xsz(3)));adwdy=0._rprec;
    allocate(adwdz(dcp%xsz(3)));adwdz=0._rprec;

    allocate(adudx2(dcp%xsz(3)));adudx2=0._rprec;
    allocate(adudy2(dcp%xsz(3)));adudy2=0._rprec;
    allocate(adudz2(dcp%xsz(3)));adudz2=0._rprec;
    allocate(advdx2(dcp%xsz(3)));advdx2=0._rprec;
    allocate(advdy2(dcp%xsz(3)));advdy2=0._rprec;
    allocate(advdz2(dcp%xsz(3)));advdz2=0._rprec;
    allocate(adwdx2(dcp%xsz(3)));adwdx2=0._rprec;
    allocate(adwdy2(dcp%xsz(3)));adwdy2=0._rprec;
    allocate(adwdz2(dcp%xsz(3)));adwdz2=0._rprec;

    allocate(atxx(dcp%xsz(3)));atxx=0._rprec;
    allocate(atyy(dcp%xsz(3)));atyy=0._rprec;
    allocate(atzz(dcp%xsz(3)));atzz=0._rprec;
    allocate(atxy(dcp%xsz(3)));atxy=0._rprec;
    allocate(atxz(dcp%xsz(3)));atxz=0._rprec;
    allocate(atyz(dcp%xsz(3)));atyz=0._rprec;
    
    allocate(acs(dcp%xsz(3)));acs=0._rprec;
    allocate(abeta(dcp%xsz(3)));abeta=0._rprec; 
    
    allocate(auuw(dcp%xsz(3)));auuw=0._rprec;
    allocate(avvw(dcp%xsz(3)));avvw=0._rprec; 
    allocate(awww(dcp%xsz(3)));awww=0._rprec;
    
    allocate(apw(dcp%xsz(3)));apw=0._rprec;
    allocate(apu(dcp%xsz(3)));apu=0._rprec;
    
    allocate(adpdx(dcp%xsz(3)));adpdx=0._rprec;
    allocate(adpdy(dcp%xsz(3)));adpdy=0._rprec;
    allocate(adpdz(dcp%xsz(3)));adpdz=0._rprec;
    allocate(apdudx(dcp%xsz(3)));apdudx=0._rprec;
    allocate(apdvdy(dcp%xsz(3)));apdvdy=0._rprec;
    allocate(apdwdz(dcp%xsz(3)));apdwdz=0._rprec;

    allocate(autxz(dcp%xsz(3)));autxz=0._rprec;
    allocate(avtyz(dcp%xsz(3)));avtyz=0._rprec;
    allocate(awtzz(dcp%xsz(3)));awtzz=0._rprec;
    
    allocate(adxx(dcp%xsz(3)));adxx=0._rprec;
    allocate(adyy(dcp%xsz(3)));adyy=0._rprec;
    allocate(adzz(dcp%xsz(3)));adzz=0._rprec;
    allocate(adxy(dcp%xsz(3)));adxy=0._rprec;
    allocate(adxz(dcp%xsz(3)));adxz=0._rprec;
    allocate(adyz(dcp%xsz(3)));adyz=0._rprec;
    
    allocate(apdudz(dcp%xsz(3)));apdudz=0._rprec;
    allocate(apdwdx(dcp%xsz(3)));apdwdx=0._rprec;
    allocate(auww(dcp%xsz(3)));auww=0._rprec;

    allocate(adudxdwdx(dcp%xsz(3)));adudxdwdx=0._rprec;
    allocate(adudydwdy(dcp%xsz(3)));adudydwdy=0._rprec;
    allocate(adudzdwdz(dcp%xsz(3)));adudzdwdz=0._rprec;

    return
  end subroutine tke_budget_init
  
  !------------------------------------------------------------------------------
  !! subroutine: finalise tke budget module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 18/05/2015
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine tke_budget_finalize()
    
    implicit none
    
    deallocate(au,av,aw,ap,auu,avv,aww,auv,auw,avw,&
       adudx,adudy,adudz,advdx,advdy,advdz,adwdx,adwdy,adwdz,&
       adudx2,adudy2,adudz2,advdx2,advdy2,advdz2,adwdx2,adwdy2,adwdz2,  &
       atxx,atyy,atzz,atxy,atxz,atyz,acs,abeta,auuw,avvw,awww,&
       apw,apu,adpdx,adpdy,adpdz,apdudx,apdvdy,apdwdz,autxz,avtyz,awtzz,&
       adxx,adyy,adzz,adxy,adxz,adyz,apdudz,apdwdx,auww,&
       adudxdwdx,adudydwdy,adudzdwdz)
       
  
  end subroutine tke_budget_finalize

  !!------------------------------------------------------------------------------
  !!  OBJECT: subroutine ComputeBudget()
  !!------------------------------------------------------------------------------
  !!
  !!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 26/07/2014
  !!
  !!  DESCRIPTION:
  !!
  !!  This subroutine extracts the needed statistics for the KAT simulation
  !!  all data are interpolated in w nodes
  !!
  !!
  !!  Production is OK
  !!  Buoyancy production is also OK
  !!  Pressure transport is OK
  !!  Transport by molecular diffusion and turbulence is OK
  !!  Dissipation is OK
  !!  
  !!------------------------------------------------------------------------------
  subroutine compute_tke_budget(jt)
    implicit none

    integer, intent(in) :: jt
    integer :: i,j,k
    !scratch useful for interpolation
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1) :: scr1,scr2,scr3,scr4  

    ! mean vars on uvp nodes
    call avg_var(u,au,'u',jt,ii)
    call avg_var(v,av,'v',jt,ii)
    call int_var(w,scr1); call avg_var(scr1,aw,'w',jt,ii)
    !if(s_flag) call avg_var(sc1,as,'s',jt,ii)
    ! mean derivatives
    call avg_var(dudx,adudx,'dudx',jt,ii)
    call avg_var(dudy,adudy,'dudy',jt,ii)
    call int_var(dudz,scr2); call avg_var(scr2,adudz,'dudz',jt,ii)
    call avg_var(dvdx,advdx,'dvdx',jt,ii)
    call avg_var(dvdy,advdy,'dvdy',jt,ii)
    call int_var(dvdz,scr2);
    call avg_var(scr2,advdz,'dvdz',jt,ii)
    call int_var(dwdx,scr2); call avg_var(scr2,adwdx,'dwdx',jt,ii)
    call int_var(dwdy,scr3); call avg_var(scr3,adwdy,'dwdy',jt,ii)
    call avg_var(dwdz,adwdz,'dwdz',jt,ii)
    if(s_flag)then
       !call avg_var(dsc1dx,adsdx,'dsdx',jt,ii)
       !call avg_var(dsc1dy,adsdy,'dsdy',jt,ii)
       !call int_var(dsc1dz,scr2); call avg_var(scr2,adsdz,'dsdz',jt,ii)
    endif

    ! res flux
    call avg_var(u*u,auu,'uu',jt,ii)
    call avg_var(v*v,avv,'vv',jt,ii)
    call avg_var(scr1*scr1,aww,'ww',jt,ii)
    call avg_var(u*v,auv,'uv',jt,ii)
    call avg_var(scr1*u,auw,'uw',jt,ii)
    call avg_var(scr1*v,avw,'vw',jt,ii)
    if(s_flag)then
       !call avg_var(sc1*sc1,ass,'ss',jt,ii)    
       !call avg_var(u*sc1,aus,'su',jt,ii)    
       !call avg_var(v*sc1,avs,'sv',jt,ii)
       !call avg_var(scr1*sc1,aws,'sw',jt,ii)
    endif
    
    !sgs flux or dns stress
    call avg_var(txx,atxx,'txx',jt,ii) 
    call avg_var(tyy,atyy,'tyy',jt,ii)
    call avg_var(tzz,atzz,'tzz',jt,ii)
    call avg_var(txy,atxy,'txy',jt,ii)
    call int_var(txz,scr1); call avg_var(scr1,atxz,'txz',jt,ii)
    call int_var(tyz,scr1); call avg_var(scr1,atyz,'tyz',jt,ii)
    if(s_flag)then
       !call avg_var(txsc1,atsx,'tsx',jt,ii)
       !call avg_var(tysc1,atsy,'tsy',jt,ii)
       !call int_var(tzsc1,scr1); call avg_var(scr1,atsz,'tsz',jt,ii)
    endif

    if(SGS)then
       scr1=cs_opt2
       call avg_var(scr1,acs,'cs',jt,ii)                     !les stuff not interpolated
       !call avg_var(BETA,abeta,'beta',ii)
       if(s_flag)then
          !scr1(:,:,1:nz)=ds_opt2_sc1(:,:,:)
          !scr1(:,:,0)=scr1(:,:,1)   !this way I take care of the fact that Cs is defined in UVP in jz=1
          !call avg_var(scr1,ads,'ds',jt,ii)                     !les stuff 
          !call avg_var(beta_sc1,abetas,'betas',jt,ii)
       endif
    endif

    !turb transport
    call int_var(w,scr1); call avg_var(u**2*scr1,auuw,'uuw',jt,ii)
    call avg_var(v**2*scr1,avvw,'vvw',jt,ii)
    call avg_var(scr1**3,awww,'www',jt,ii)
    call avg_var(u*scr1**2,auww,'uww',jt,ii)
    !if(s_flag) call avg_var(sc1*scr1**2,asww,'sww',jt,ii)
    
    !sgs transp
    call avg_var(scr1*tzz,awtzz,'wtzz',jt,ii)
    call int_var(txz,scr3); 
    if((dcp%xst(3)==1).and.(lbc%special.eq.'wall_law'))then
       scr3(:,:,1)=txz(:,:,1)
    endif
    call avg_var(scr3*u,autxz,'utxz',jt,ii)
    call int_var(tyz,scr2); 
    if((dcp%xst(3)==1).and.(lbc%special.eq.'wall_law'))then 
       scr2(:,:,1)=tyz(:,:,1)
    endif
    call avg_var(scr2*v,avtyz,'vtyz',jt,ii)

    !dissipation (sgs)
    call avg_var(dudx*txx,adxx,'dxx',jt,ii)
    call avg_var(dvdy*tyy,adyy,'dyy',jt,ii)
    call avg_var(dwdz*tzz,adzz,'dzz',jt,ii)
    call avg_var(0.5_rprec*(dudy+dvdx)*txy,adxy,'dxy',jt,ii)
    call int_var(dudz+dwdx,scr4); 
    if((dcp%xst(3)==1).and.(lbc%special.eq.'wall_law'))then
       scr4(:,:,1)=dudz(:,:,1)+0.5_rprec*(dwdx(:,:,1)+dwdx(:,:,2))
    endif
    call avg_var(0.5_rprec*scr4*scr3,adxz,'dxz',jt,ii)
    call int_var(dvdz+dwdy,scr1); 
    if((dcp%xst(3)==1).and.(lbc%special.eq.'wall_law'))then
       scr1(:,:,1)=dvdz(:,:,1)+0.5_rprec*(dwdy(:,:,1)+dwdy(:,:,2))
    endif
    call avg_var(0.5_rprec*scr2*scr1,adyz,'dyz',jt,ii)

    !dissipation
    call int_var(dudz,scr1); 
    call int_var(dvdz,scr2); 
    call int_var(dwdx,scr3); 
    call int_var(dwdy,scr4); 
    call avg_var(dudx*dudx,adudx2,'dudx2',jt,ii)
    call avg_var(dudy*dudy,adudy2,'dudy2',jt,ii)
    call avg_var(scr1*scr1,adudz2,'dudz2',jt,ii)
    call avg_var(dvdx*dvdx,advdx2,'dvdx2',jt,ii)
    call avg_var(dvdy*dvdy,advdy2,'dvdy2',jt,ii)
    call avg_var(scr2*scr2,advdz2,'dvdz2',jt,ii)
    call avg_var(scr3*scr3,adwdx2,'dwdx2',jt,ii)
    call avg_var(scr4*scr4,adwdy2,'dwdy2',jt,ii)
    call avg_var(dwdz*dwdz,adwdz2,'dwdz2',jt,ii)

    !dissipation (UW)
    call int_var(dwdx,scr1)
    call int_var(dwdy,scr2)
    call int_var(dudz,scr3)
    call avg_var(dudx*scr1,adudxdwdx,'dudxdwdx',jt,ii)
    call avg_var(dudy*scr2,adudydwdy,'dudydwdy',jt,ii)
    call avg_var(scr3*dwdz,adudzdwdz,'dudzdwdz',jt,ii)

    if(s_flag)then
       !dissipation (SW)
       !call int_var(dwdx,scr1)
       !call int_var(dwdy,scr2)
       !call int_var(dsc1dz,scr3)
       !call avg_var(dsc1dx*scr1,adsdxdwdx,'dsdxdwdx',jt,ii)
       !call avg_var(dsc1dy*scr2,adsdydwdy,'dsdydwdy',jt,ii)
       !call avg_var(scr3*dwdz,adsdzdwdz,'dsdzdwdz',jt,ii)
    endif

    !pressure stuff
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             scr1(i,j,k)=p(i,j,k)-0.5_rprec*(u(i,j,k)**2+v(i,j,k)**2+(0.5_rprec*(w(i,j,k)+w(i,j,k+1)))**2) 
          enddo
       enddo
    enddo
    call update_ghost(scr1,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)

    !derivatives
    call grad_velocity_filt(scr1,scr2,scr3)
    call avg_var(scr1,ap,'p',jt,ii)    !cleaned pressure
    call avg_var(scr2,adpdx,'dpdx',jt,ii)    !cleaned pressure
    call avg_var(scr3,adpdy,'dpdy',jt,ii)    !cleaned pressure
    call ddz_uv(scr2,scr1)
    call int_var(scr2,scr3); call avg_var(scr3,adpdz,'dpdz',jt,ii)    !cleaned pressure
    call int_var(w,scr2)
    call avg_var(scr1*scr2,apw,'pw',jt,ii)        !pressure transp tke
    call avg_var(scr1*u,apu,'pu',jt,ii)
    !call avg_var(scr1*v,apv,'pv',jt,ii)
    !if(s_flag) call avg_var(scr1*sc1,aps,'ps',jt,ii)
    
    call avg_var(scr1*dudx,apdudx,'pdudx',jt,ii)
    call avg_var(scr1*dvdy,apdvdy,'pdvdy',jt,ii)
    call avg_var(scr1*dwdz,apdwdz,'pdwdz',jt,ii)
    call int_var(dudz,scr2)
    call int_var(dwdx,scr3)
    call avg_var(scr1*scr3,apdwdx,'pdwdx',jt,ii)  !pressure transp tke
    call avg_var(scr1*scr2,apdudz,'pdudz',jt,ii)  !pressure transp tke
    if(s_flag)then
       !call int_var(dsc1dz,scr4)
       !call avg_var(scr1*scr4,apdsdz,'pdsdz',jt,ii)  !pressure transp tke
    endif
    !update record
    if(mod(jt,tke_period).eq.0)then
       ii=ii+1
    endif

  end subroutine compute_tke_budget


  !!------------------------------------------------------------------------------
  !!  OBJECT: subroutine avg_var()
  !!------------------------------------------------------------------------------
  subroutine avg_var(var,avar,str_v,jt,ii)
    implicit none

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(IN) :: var
    real(rprec),dimension(dcp%xsz(3)),intent(INOUT) :: avar
    character(*),intent(in) :: str_v
    integer,intent(in) :: jt,ii
    character(128) :: fname
    real(rprec),dimension(dcp%xsz(3)) :: plavg
    integer :: j,k
    real(rprec) :: n
    real(rprec),dimension(1:nz) :: avTot
    integer,dimension(dims(2)) :: rcount,rdispl
    integer :: ierror
    
    call compute_planar_average(plavg,var(:,:,1:dcp%xsz(3)),dcp)

    do k=1,dcp%xsz(3)
       avar(k)=avar(k)+plavg(k)
    enddo
    
    !save if index is reached
    if((mod(jt,tke_period).eq.0).and.(dcp%xst(2)==1))then
       avTot=1234567.0_rprec
       do k=1,dims(2)
          rcount(k) = dcp_sp_sizes(1,3,k); 
       enddo
       rdispl(1)=0
       do k=2,dims(2)
          rdispl(k) = rdispl(k-1)+rcount(k-1); 
       enddo
       call mpi_gatherv(avar(1),dcp%xsz(3),real_type,&
            avTot(1),rcount,rdispl,real_type,0,DECOMP_2D_VRTCL_X,ierror)
       if(nrank==0)then
          n=real(tke_period,rprec)
          write(fname,'(a,i0)') trim(oPath)//trim(str_v)//'.out'
          open(11,file=trim(fname),access='direct',recl=2*nz)
          write(11,rec=ii) avTot/n
          close(11)
       endif
       !reinitialize averages
       avar=0.0_rprec      
    endif

  end subroutine avg_var

  !!------------------------------------------------------------------------------
  !!  OBJECT: subroutine int_var()
  !!------------------------------------------------------------------------------
  !!
  !!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 20/01/2014
  !!
  !!  DESCRIPTION:
  !!
  !!  This subroutine interpolates from w to uvp
  !!
  !!------------------------------------------------------------------------------
  subroutine int_var(var,ivar)
    implicit none

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(in) :: var
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(out) :: ivar
    integer :: i,j,k

    !$omp parallel do
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2) 
          do i=1,dcp%xsz(1) 
             ivar(i,j,k)=0.5_rprec*(var(i,j,k)+var(i,j,k+1))
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    return
  end subroutine int_var


end module tke_budget1D
