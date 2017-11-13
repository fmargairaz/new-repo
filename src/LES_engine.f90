!------------------------------------------------------------------------------
!!  module: LES_engine
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------

module LES_engine

  use parameters_IO
  use system_variables

  use decomp_2d
  use decomp_2d_custom
  use decomp_2d_fft_2d
  use fft_engine
  use utility_tools

  use initial_conditions
  use boundary_conditions

  use compute_wall_law
  use compute_geostrophic
  use sponge_damping
  use testfilter

  use compute_derivatives
!choose dealiasing module
#ifdef DEA_TRUNC
  use compute_convective_trunc
#endif
#ifdef DEA_NODLSG
  use compute_convective_nodlsg
#endif
#ifdef DEA_FSMTHG
  use compute_convective_Fsmthg
#endif
#ifdef DEA_PADD
  use compute_convective_padd
#endif  

  use compute_pressure
  
#ifdef DNS_MOD
  use compute_dns_stress
#else
  use compute_sgs
#endif
  
  use running_diagnositcs
  use mod_IO_checkpoint  

  use mod_IO_main 
  use mod_IO_statistics
  use tke_budget1D

  implicit none
  include 'mpif.h'

  private
  public :: LES_main_engine

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: LES main engine
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine LES_main_engine()
    implicit none
    
    ! time loop parameters
    integer, parameter :: NaN_check=20
    logical, parameter :: profiling=.true.
    
    ! MPI variables
    integer :: me,nproc,ierror
    integer, dimension(2) :: dims

    ! indexes for spatial and time loops
    integer :: i,j,k
    integer :: k_min
    integer :: jt,jt_total

    ! physical time
    real(rprec) :: tt
    real(rprec) :: z

    ! other variables
    real(rprec) :: tadv1,tadv2

    ! profiling and timing variables
    real(rprec) :: time_grad=0._rprec
    real(rprec) :: time_stress=0._rprec
    real(rprec) :: time_rhs=0._rprec
    real(rprec) :: time_press=0._rprec
    real(rprec) :: time_loop=0._rprec
    real(rprec) :: time_out=0._rprec
    real(rprec) :: time_tmp,time_sub

    call MPI_comm_rank(MPI_COMM_WORLD, me, ierror)
    call MPI_comm_size(MPI_COMM_WORLD, nproc, ierror)  

    if(me==0)then
       write(*,*) '## Starting main engine'
    endif

    call read_parameters(me,nproc,ierror)
    ! variables to modifiy for restart
    jt=0;jt_total=0;tt=0._rprec
    if(restart_flag)then
       call read_parameters_restart(jt_total,me)
       tt=jt_total*dt
    endif

    ! initiate modules =================================================
    call create_system(me,nproc,dims,ierror)
    if(me==0)then
       write(*,*) '## system allocated - ready to initialze modules'
    endif

    call fft_init()
    call decomp_2d_fft_2d_init()
    call decomp_2d_partition_info(dcp)

    call deriv_init(dcp,dcpg)
    call conv_init(dcp,dcpg)
    call press_init(dcp)
    call wall_init(dcp)
    call sponge_damping_init(dcp)
    call test_filter_init(dcp)
#ifndef DNS_MOD
    call sgs_init(dcp)
#endif

    ! IO modules
    call running_diagnostics_init(dcp)
    call checkpoint_io_init(dcp)
    
    if(output) call io_init(dcp)
    if(stat_flag) call io_stat_init(dcp)
    if(tke_flag) call tke_budget_init(dcp)

    call test_filter_create(2._rprec,G_test)
    call test_filter_create(4._rprec,G_test_test)

    if(me==0)then
       write(*,*) '## Modules initialised'
    endif

    ! initial conditions =================================================
    call set_bc_surface(lbc_z0,dcp) !surface roughness for wall law
    if(.not.restart_flag)then
       call set_initial_conditions(ic_type,Re,u_star,pfx,pfy,dcp)
    else
       call checkpoint_in(jt_total)
    endif

    !using ghost cell k=0 to impose Newmann bc
    if(lbc%special.eq.'wall_law')then
       lbc%typ_u='neu'
       lbc%typ_v='neu'
    endif

    if(output) call output_loop(0,jt_total,nrank)    
    call compute_running_diagnistics(jt,jt_total,tt,.false.)
       
    if(DEBUG_TIMELOOP)then 
       write(*,'(A,I3,A)') 'DEBUG - system initialisation on node #',nrank,' <DONE>'
    endif

    ! not needed
    call MPI_barrier(MPI_COMM_WORLD, ierror)
    
    if(me==0)then
       write(*,*) '## Flow initialized - read for simulation'
       write(*,*) '=============================================================================='       
    endif

    ! START OF MAIN TIME LOOP #################################################
    do jt=1,nt

       time_tmp = mpi_wtime()

       ! advance time and step
       tt = tt+dt
       jt_total = jt_total+1

       !===================================================================
       ! geostrophic wind correction
       if((geostr_rot).and.(jt_total.ge.geostr_start))then
          call compute_geostrophic_correct(jt_total)
       endif

       !===================================================================
       ! compute the friction velocity at k=1
       if((dcp%xst(3)==1).and.(lbc%special.eq.'wall_law'))then
          call obukhov(u_star_avg,wall_model)
          call wallstress(lbc_u,lbc_v,wall_model)
       endif
       
       !===================================================================
       ! derivatives
       if(profiling) time_sub=mpi_wtime()
       call compute_vel_grad()
       if(profiling) time_grad=time_grad+mpi_wtime()-time_sub

       !===================================================================
       ! overwrite vel.deriv. + sc.deriv. + stresses at k=1
       ! surface model (it's always +1 actually)
       if((dcp%xst(3)==1).and.(lbc%special.eq.'wall_law'))then
          !call wallstress(lbc_u,lbc_v,wall_model)
       endif
       
       !===================================================================
       ! stress calculations
       if(profiling) time_sub=mpi_wtime()
#ifdef DNS_MOD
       call dns_stress(dcp,Re,lbc%special)
#else
       call sgs_stress(jt,jt_total)
#endif
       if(profiling) time_stress=time_stress+mpi_wtime()-time_sub

       ! update ghost cell on txz,tyz,tzz
       call update_ghost(txz,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_ghost(tyz,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_ghost(tzz,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       
       if(DEBUG_TIMELOOP)then 
          write(*,'(A,I3,A)') 'DEBUG - compuding stress tensor on node #',nrank,' <DONE>'
       endif
       
       !===================================================================
       ! build RHS for advection step

       if(profiling) time_sub=mpi_wtime()
       call compute_RHS_adv()
       if(profiling) time_rhs=time_rhs+mpi_wtime()-time_sub
       
       !===================================================================
       ! time advancement 
       
       ! choose the time scheme (Adam-Bashforth or Euler)
       if (jt_total == 1) then
          tadv1 = 1.0_rprec
          tadv2 = 0.0_rprec
       else
          tadv1 = 1.5_rprec
          tadv2 = 1.0_rprec-tadv1
       endif

       ! advection step
       !$omp parallel do 
       do k=1,dcp%xsz(3)
          u(:,:,k) = u(:,:,k)+dt*(tadv1*rhsx(:,:,k)+tadv2*rhsx_f(:,:,k)) 
          v(:,:,k) = v(:,:,k)+dt*(tadv1*rhsy(:,:,k)+tadv2*rhsy_f(:,:,k))
          w(:,:,k) = w(:,:,k)+dt*(tadv1*rhsz(:,:,k)+tadv2*rhsz_f(:,:,k))
       enddo
       !$omp end parallel do

       ! update ghost cells after advection step
       call update_ghost(u,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_ghost(v,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_ghost(w,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       
       ! impose zero "intermediate (in a fractional step sense)" velocity at 
       ! top and bottom layer (and below and above domain) 
       if(dcp%xst(3)==1)then
          w(:,:,0)=0._rprec
          !w(:,:,1)=0._rprec
       endif
       if(dcp%xen(3)==nz)then
          w(:,:,dcp%xsz(3)+1)=0._rprec   
       endif
       
       if(DEBUG_TIMELOOP)then 
          write(*,'(A,I3,A)') 'DEBUG - advection step on node #',nrank,' <DONE>'
       endif

       !===================================================================
       ! poisson eq.
       if(profiling) time_sub=mpi_wtime()
       call press_stag_array(tadv1)
       if(profiling) time_press=time_press+mpi_wtime()-time_sub

      if(DEBUG_TIMELOOP)then 
          write(*,'(A,I3,A)') 'DEBUG - compuding pressure on node #',nrank,' <DONE>'
       endif 

       !===================================================================
       ! projection step
       !$omp parallel 
       !$omp do 
       do k=1,dcp%xsz(3)
          u(:,:,k) = u(:,:,k) - dt*tadv1*dpdx(:,:,k)
          v(:,:,k) = v(:,:,k) - dt*tadv1*dpdy(:,:,k)
          w(:,:,k) = w(:,:,k) - dt*tadv1*dpdz(:,:,k)
       enddo
       !$omp end do

       ! update rhs with pressure gradient for rhs_f
       !$omp do
       do k=1,dcp%xsz(3)
          rhsx(:,:,k) = rhsx(:,:,k)-dpdx(:,:,k)
          rhsy(:,:,k) = rhsy(:,:,k)-dpdy(:,:,k)
          rhsz(:,:,k) = rhsz(:,:,k)-dpdz(:,:,k)
       enddo
       !$omp end do 
       !$omp end parallel 

       ! update ghost cells after projection step
       call update_ghost(u,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_ghost(v,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
       call update_ghost(w,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)

       ! top and bottom layer (and below and above domain) 
       if(dcp%xst(3)==1)then
          w(:,:,1)=0._rprec
       endif
       if(dcp%xen(3)==nz)then
          w(:,:,dcp%xsz(3)+1)=0._rprec   
       endif

       if(DEBUG_TIMELOOP)then 
          write(*,'(A,I3,A)') 'DEBUG - projection step on node #',nrank,' <DONE>'
       endif

       !===================================================================
       ! re-apply bc for output purposes
       if(dcp%xst(3)==1)then
          !call compute_bc('bot',lbc_u,lbc_v,lbc_w,dcp)
          
          call set_bc('bot',u,lbc%typ_u,lbc_u,.true.,dcp)
          call set_bc('bot',v,lbc%typ_v,lbc_v,.true.,dcp)
          call set_bc('bot',w,lbc%typ_w,lbc_w,.false.,dcp)
       endif
       if(dcp%xen(3)==nz)then
          !call compute_bc('top',ubc_u,ubc_v,ubc_w,dcp)
          
          call set_bc('top',u,ubc%typ_u,ubc_u,.true.,dcp)
          call set_bc('top',v,ubc%typ_v,ubc_v,.true.,dcp)
          call set_bc('top',w,ubc%typ_w,ubc_w,.false.,dcp)
       endif
       if(profiling) time_loop=time_loop+mpi_wtime()-time_tmp

       !===================================================================
       ! output & post time step chesks

       if(profiling) time_tmp=mpi_wtime()
       
       call compute_running_diagnistics(jt,jt_total,tt,.true.)
       
       ! check if velocity fields containes NaN's
       if(mod(jt,NaN_check).eq.0)then
          call check_NaN(jt)
       endif

       ! output stats
       if((stat_flag).and.(jt_total.gt.start_stats))then 
          call  io_stat_averages(jt,jt_total)
       endif
       
       ! spectra 
       if(spectra_flag)then
          call io_stat_spectra(jt,jt_total)
       endif

       ! tke budget
       if((tke_flag).and.(jt_total.gt.start_tke))then       
          call compute_tke_budget(jt)
       endif

       ! output
       if(output .and. mod(jt,step_out).eq.0)then
          call output_loop(jt,jt_total,nrank) 
       endif

       if(mod(jt,scheckpoint).eq.0)then
          call scheckpoint_out()
       endif
       if(mod(jt_total,checkpoint).eq.0)then
          call checkpoint_out(jt_total)
       endif

       if(profiling) time_out=time_out+mpi_wtime()-time_tmp

       if(DEBUG_TIMELOOP)then 
          write(*,'(A,I7,A,I3,A)') 'DEBUG - time step ',jt,' on node #',nrank,' <DONE>'
       endif

    enddo
    ! END OF MAIN TIME LOOP #################################################
    
    ! write checkpoint at teh end of the run (with check to avoid double output) 
    if(mod(jt_total,checkpoint).ne.0)then
       call checkpoint_out(jt_total)
    endif
    call write_parameters_restart(jt_total,nrank)

    if(output) call io_finalize()
    if(stat_flag) call io_stat_finalize()
    if(tke_flag) call tke_budget_finalize()

    call wall_finalize()
    call sponge_damping_finalize()
    call deriv_finalize()
    call conv_finalize()
    call press_finalize()

#ifndef DNS_MOD
    call sgs_finalize()
#endif
    call clean_system()
    
    if(DEBUG_TIMELOOP)then 
       write(*,'(A,I3,A)') 'DEBUG - system cleaning on node #',nrank,' <DONE>'
    endif

    if((nrank==0).and.(profiling==.true.))then
       write(*,'(A,F10.4,A)') ' wall time for time loop: ',time_loop,'s' 
       write(*,'(A,F10.4,A)') ' wall time for output: ',time_out,'s'
       
       open(99,file='./profiling.txt')
       write(99,'(I15.0)') p_on_y,p_on_z,nx,ny,nz,nt
       write(99,'(F15.7)') time_grad,time_stress,time_rhs,time_press,time_loop,time_out
    endif

    return
  end subroutine LES_main_engine

  !------------------------------------------------------------------------------
  !! subroutine: compute velocity gradients
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_vel_grad()
    implicit none

    !compute gradients in x&y
    call grad_velocity_filt(u,dudx,dudy)
    call grad_velocity_filt(v,dvdx,dvdy)
    call grad_velocity_filt(w,dwdx,dwdy)

    ! update ghost cells with filtered velocity
    !call update_ghost(u,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    !call update_ghost(v,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    !call update_ghost(w,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)

    !===================================================================
    ! boundary conditions
    if(dcp%xst(3)==1)then
       call set_bc('bot',u,lbc%typ_u,lbc_u,.true.,dcp)
       call set_bc('bot',v,lbc%typ_v,lbc_v,.true.,dcp)
       call set_bc('bot',w,lbc%typ_w,lbc_w,.false.,dcp)
    endif
    if(dcp%xen(3)==nz)then
       call set_bc('top',u,ubc%typ_u,ubc_u,.true.,dcp)
       call set_bc('top',v,ubc%typ_v,ubc_v,.true.,dcp)
       call set_bc('top',w,ubc%typ_w,ubc_w,.false.,dcp)
    endif
    
    !compute gradients in z
    call ddz_uv(dudz,u)
    call ddz_uv(dvdz,v)
    call ddz_w(dwdz,w)

    ! update ghost cells
    !call update_ghost(dudx,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    !call update_ghost(dudy,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(dudz,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)

    !call update_ghost(dvdx,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    !call update_ghost(dvdy,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(dvdz,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    
    !call update_ghost(dwdx,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    !call update_ghost(dwdy,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(dwdz,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)

    if(DEBUG_TIMELOOP)then 
       write(*,'(A,I3,A)') 'DEBUG - compuding velocity gradients on node #',nrank,' <DONE>'
    endif

    return
  end subroutine compute_vel_grad  

  !------------------------------------------------------------------------------
  !! subroutine: build the Right Hand Side for advection step
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_RHS_adv()
    implicit none

    real(rprec) :: ug_time_fact
    integer :: k,ierror
   
    !update old rhs
    !$omp parallel do
    do k=1,dcp%xsz(3)
       rhsx_f(:,:,k)=rhsx(:,:,k)
       rhsy_f(:,:,k)=rhsy(:,:,k)
       rhsz_f(:,:,k)=rhsz(:,:,k)
    enddo
    !$omp end parallel do

    !$omp parallel do
    do k=1,dcp%xsz(3)
       rhsx(:,:,k)=0._rprec
       rhsy(:,:,k)=0._rprec
       rhsz(:,:,k)=0._rprec
    enddo
    !$omp end parallel do
    
    ! add convective terms 
    call compute_dealias_convec(rhsx,rhsy,rhsz)
   
   if(DEBUG_TIMELOOP)then 
       write(*,'(A,I3,A)') 'DEBUG - compuding convective on node #',nrank,' <DONE>'
    endif
    

    ! add divergence of stress
    call div_stress(rhsx,txx,txy,txz,.true.)
    call div_stress(rhsy,txy,tyy,tyz,.true.)
    call div_stress(rhsz,txz,tyz,tzz,.false.,bottom_wk_press,top_wk_press)

    ! pressure gradient forcing
    !$omp parallel do 
    do k=1,dcp%xsz(3)
       rhsx(:,:,k)=rhsx(:,:,k)+pfx
       rhsy(:,:,k)=rhsy(:,:,k)+pfy
    enddo
    !$omp end parallel do

    ! add forcing
    if (coriolis_flag) then
       ug_time_fact = 1.d0
       !$omp parallel do
       do k=1,dcp%xsz(3)
          rhsx(:,:,k)=rhsx(:,:,k) + coriol*v(:,:,k) - ug_time_fact*coriol*vg;
          rhsy(:,:,k)=rhsy(:,:,k) - coriol*u(:,:,k) + ug_time_fact*coriol*ug;
       enddo
       !$omp end parallel do	
    end if
    
    ! add sponge damping to the RHS
    !call add_sponge_damping()

    if(DEBUG_TIMELOOP)then 
       write(*,'(A,I3,A)') 'DEBUG - compuding RHS adv. step on node #',nrank,' <DONE>'
    endif
  
    return
  end subroutine compute_RHS_adv

  !------------------------------------------------------------------------------
  !! subroutine: 
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine check_NaN(jt)
    implicit none

    integer, intent(IN) :: jt
    integer :: i,j,k
    logical :: got_NaN=.false.
    
    if(DEBUG_TIMELOOP) then
       write(*,'(A,I5,A,I3,A)') 'DEBUG - NaN check at jt = ',jt,' on node #',nrank,' <START>'
    endif
    
    !$omp parallel do private(got_NaN)
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             if(u(i,j,k).ne.u(i,j,k)) got_NaN = .true.
             if(v(i,j,k).ne.v(i,j,k)) got_NaN = .true.
             if(w(i,j,k).ne.w(i,j,k)) got_NaN = .true.
             
             if(got_NaN)then
                write(*,'(A,I4,I4,I4,A,I6,A,I5)') 'ABORT - got a NaN at (',i,j,k,') at jt=',jt,' - on process',nrank
                call decomp_2d_abort(9,'I got a NaN')
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    if(DEBUG_TIMELOOP) then
       write(*,'(A,I5,A,I3,A)') 'DEBUG - NaN check at jt = ',jt,' on node #',nrank,' <END>'
    endif
  
    return
  end subroutine check_NaN
  
end module LES_engine
