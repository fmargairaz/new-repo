!------------------------------------------------------------------------------
!!  module: running diagnositcs
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 31/01/2015
!!
!! - description:
!!
!------------------------------------------------------------------------------

module running_diagnositcs

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom

  !! - global variables
  use system_variables, only : u,v,w,dudx,dvdy,dwdz
  use compute_wall_law, only : u_star_avg

  implicit none
  include 'mpif.h'

  private
  TYPE(DECOMP_INFO),save :: dcp

  integer, parameter :: buffer_size=100
  integer, parameter :: running_diag_promt=50
  integer,save :: buffer_pos

  integer,dimension(buffer_size),save :: buffer_jt
  real(rprec),dimension(buffer_size),save :: buffer_tt
  real(rprec),dimension(buffer_size),save :: buffer_cfl,buffer_mke,buffer_u_star
  real(rprec),dimension(buffer_size),save :: buffer_divu

  public :: running_diagnostics_init,compute_running_diagnistics

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: initialise running diagnostics  module
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 02/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine running_diagnostics_init(dcp_main)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main

    buffer_jt=0 
    buffer_tt=0._rprec
    buffer_cfl=0._rprec
    buffer_mke=0._rprec
    buffer_u_star=0._rprec
    buffer_divu=0._rprec

    buffer_pos=0

    return
  end subroutine running_diagnostics_init


  !------------------------------------------------------------------------------
  !! subroutine: compute running diagnostic 
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 31/01/2015
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine compute_running_diagnistics(jt,jt_total,tt,flag_add2buffer)

    integer, intent(IN) :: jt,jt_total
    real(rprec),intent(IN) :: tt
    logical,intent(IN) :: flag_add2buffer
    real(rprec) :: cfl,mke,u_star_diag,rmsdivu
    integer :: n

    call compute_CFL(cfl)
    call compute_MKE(mke)
    call compute_u_star_diag(u_star_diag)
    call compute_rms_divu(rmsdivu)
    
    if(jt==0 .and. nrank==0)then 
       write(*,*) '## running diagnostics initilized - flow summary'
       write(*,'(A,I10,A,E11.4,A,E11.4,A,E11.4)') &
            '    jt_t=',jt_total,', t=',tt, ', CFL=',cfl,', MKE=',mke              
       
    elseif(mod(jt,running_diag_promt)==0 .and. nrank==0)then    
       write(*,'(A,I10,A,I10,A,E11.4,A,E11.4,A,E11.4,A,E11.4,A,E11.4)') &
            ' jt=',jt,', jt_t=',jt_total,', t=',tt, ', CFL=',cfl,', MKE=',mke,&
            ', u_star=',u_star_diag,', divu=',rmsdivu
    else
       !do nothing
    endif
    
    if(flag_add2buffer)then
       buffer_pos=buffer_pos+1
       buffer_jt(buffer_pos)=jt_total
       buffer_tt(buffer_pos)=tt
       buffer_cfl(buffer_pos)=cfl
       buffer_mke(buffer_pos)=mke
       buffer_u_star(buffer_pos)=u_star_diag
       buffer_divu(buffer_pos)=rmsdivu
    endif

    if(buffer_pos==buffer_size)then
       if((nrank==0).and.(jt_total==buffer_size))then
          open(98,file='./running_diagnostics.txt')
          do n=1,buffer_size
             write(98,'(I10,E13.6,E16.9,E16.9,E16.9,E16.9)')&
                  buffer_jt(n),buffer_tt(n),buffer_cfl(n),buffer_mke(n),buffer_u_star(n),buffer_divu(n)
          enddo
          close(98)
       elseif((nrank==0))then
          open(98,file='./running_diagnostics.txt',access='append')
          do n=1,buffer_size
             write(98,'(I10,E13.6,E16.9,E16.9,E16.9,E16.9)')&
                  buffer_jt(n),buffer_tt(n),buffer_cfl(n),buffer_mke(n),buffer_u_star(n),buffer_divu(n)
          enddo
          close(98)
       endif
       buffer_pos=0
    endif


  end subroutine compute_running_diagnistics

  !------------------------------------------------------------------------------
  !! subroutine: compute CFL numer 
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 04/07/2014
  !!
  !! - description:
  !!  - output: cfl 
  !!
  !------------------------------------------------------------------------------
  subroutine compute_CFL(cfl)
    implicit none

    real(rprec),intent(OUT) :: cfl
    real(rprec) :: u_max,u_max_l
    real(rprec) :: v_max,v_max_l
    real(rprec) :: w_max,w_max_l
    integer :: i,j,k
    integer :: ierror

    u_max_l = -1234567890._rprec
    v_max_l = -1234567890._rprec
    w_max_l = -1234567890._rprec 

    !$omp parallel do reduction(max:u_max_l,v_max_l,w_max_l)
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             if(abs(u(i,j,k))>u_max_l) u_max_l=abs(u(i,j,k))
             if(abs(v(i,j,k))>v_max_l) v_max_l=abs(v(i,j,k))
             if(abs(w(i,j,k))>w_max_l) w_max_l=abs(w(i,j,k))
          enddo
       enddo
    enddo
    !$omp end parallel do

    call MPI_allreduce(u_max_l,u_max,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
    call MPI_allreduce(v_max_l,v_max,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
    call MPI_allreduce(w_max_l,w_max,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)

    cfl=dt*(u_max/dx+v_max/dy+w_max/dz)

    return
  end subroutine compute_CFL

  !------------------------------------------------------------------------------
  !! subroutine: compute mean kinetic energy 
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 04/07/2014
  !!
  !! - description:
  !!  - output: mke,ke
  !!
  !------------------------------------------------------------------------------
  subroutine compute_MKE(mke)
    implicit none

    real(rprec),intent(OUT) :: mke
    real(rprec) :: u_sum,u_sum_l
    real(rprec) :: v_sum,v_sum_l
    real(rprec) :: w_sum,w_sum_l
    real(rprec) :: ke_l
    integer :: i,j,k,k_min,k_max
    integer :: ierror

    u_sum = 0._rprec;u_sum_l = 0._rprec;
    v_sum = 0._rprec;v_sum_l = 0._rprec;
    w_sum = 0._rprec;w_sum_l = 0._rprec;
    !ke = 0._rprec;ke_l = 0._rprec;

    !$omp parallel do reduction(+:u_sum_l,v_sum_l,w_sum_l)
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             u_sum_l=u_sum_l+u(i,j,k)
             v_sum_l=v_sum_l+v(i,j,k)
             w_sum_l=w_sum_l+0.5_rprec*(w(i,j,k)+w(i,j,k+1))
          enddo
       enddo
    enddo
    !$omp end parallel do

    u_sum_l=u_sum_l/real(nx*ny*nz,rprec)
    v_sum_l=v_sum_l/real(nx*ny*nz,rprec)
    w_sum_l=w_sum_l/real(nx*ny*nz,rprec)

    call MPI_allreduce(u_sum_l,u_sum,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    call MPI_allreduce(v_sum_l,v_sum,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    call MPI_allreduce(w_sum_l,w_sum,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    !call MPI_allreduce(ke_l,ke,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)

    mke=0.5_rprec*(u_sum*u_sum + v_sum*v_sum + w_sum*w_sum)
    !ke=0.5_rprec*ke/(nx*ny*nz) 

    return
  end subroutine compute_MKE

  !------------------------------------------------------------------------------
  !! subroutine: compute u_star average 
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 04/07/2014
  !!
  !! - description:
  !!  - output: u_star_diag
  !!
  !------------------------------------------------------------------------------
  subroutine compute_u_star_diag(u_star_diag)
    implicit none

    real(rprec),intent(OUT) :: u_star_diag
    real(rprec) :: u_sum,u_sum_l
    integer :: i,j
    integer :: ierror

    u_sum = 0._rprec;u_sum_l = 0._rprec;

    if(dcp%xst(3)==1)then
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             u_sum_l=u_sum_l+u_star_avg(i,j)
          enddo
       enddo
       call MPI_reduce(u_sum_l,u_sum,1,real_type,MPI_SUM,0,DECOMP_2D_LAYER_X,ierror)
    endif
    if(nrank==0)then
       u_star_diag=u_sum/(nx*ny)
    endif
    call MPI_Bcast(u_star_diag,1,real_type,0,MPI_COMM_WORLD, ierror);
    return
  end subroutine compute_u_star_diag

  !------------------------------------------------------------------------------
  !! subroutine: compute average div u
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 13/01/2017
  !!
  !! - description:
  !!  - output: rmsdivu
  !!
  !------------------------------------------------------------------------------
  subroutine compute_rms_divu(rms)
    implicit none

    real(rprec),intent(OUT) :: rms
    real(rprec) :: rms_l
    integer :: i,j,k
    integer :: ierror

    rms=0._rprec;rms_l=0._rprec;
    
    !$omp parallel do reduction(+:rms_l)
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             rms_l=rms_l+(dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))**2
          enddo
       enddo
    enddo
    !$omp end parallel do

    call MPI_allreduce(rms_l,rms,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    rms=sqrt(rms/real(nx*ny*nz,rprec))

    return
  end subroutine compute_rms_divu

end module running_diagnositcs

