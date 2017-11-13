!------------------------------------------------------------------------------
!! module: initial conditions
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------
module initial_conditions

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use decomp_2d_io
  use boundary_conditions

  !! - global variables
  use system_variables,only:u,v,w,lbc_u,lbc_v,lbc_w,ubc_u,ubc_v,ubc_w

  implicit none
  include 'mpif.h'

  private

  public :: set_initial_conditions

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: set initial conditions
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - ,u_star,pfx,pfy
  !!  - constatn velocity fields:
  !!      set_initial_conditions('constant',Re,0._rprec,0._rprec,0._rprec,dcp)
  !!  - Log law:     
  !!      set_initial_conditions('log_law',Re,u_star,pfx,pfy,dcp)
  !!  - Poseuille:   
  !!      set_initial_conditions('poiseuille',Re,0._rprec,0._rprec,0._rprec,dcp)
  !!
  !------------------------------------------------------------------------------
  subroutine set_initial_conditions(ic_type,Re,u_star,pfx,pfy,dcp)
    use compute_wall_law,only:lbc_z0
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp
    character(*),intent(in) :: ic_type
    real(rprec), intent(in) :: Re,u_star,pfx,pfy
    

    integer :: i,j,k
    integer :: reclen
    integer :: errorcode
    real(rprec) :: x,y,z
    real(rprec) :: dspl_h=0._rprec

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)) :: arg,arg2
    real(rprec),dimension(nx,ny,nz) :: vel_in
    real(rprec) :: dir_x,dir_y,norm

    ! to change when scalar are implemented
    real(rprec) :: w_star
    real(rprec),parameter :: T_scale=290._rprec,wt_s= 0.1_rprec,T_init=290._rprec

    !normalised direction
    norm=(pfx**2+pfy**2)**0.5_rprec
    if(norm .gt. 0.0_rprec)then
       dir_x=pfx/norm
       dir_y=pfy/norm
    else
       dir_x=1.0_rprec
       dir_y=0.0_rprec
    endif

    w_star=(9.81d0/T_init*abs(wt_s)*zi)**(1._rprec/3._rprec)/u_scale;

    if(ic_type.eq.'constant')then
       ! constant velocity ======================================
       u=u0
       v=v0
       w=w0
    elseif(coriolis_flag)then
       ! coriolis_forcing ======================================
       do k=1,dcp%xsz(3)
          !define height
          z=(dcp%xst(3)+k-1.5_rprec)*dz
          if(z.le.0.5)then
             u(:,:,k)=ug
             v(:,:,k)=vg
             w(:,:,k)=0.0_rprec
          else
             u(:,:,k)=ug
             v(:,:,k)=0.0_rprec
             w(:,:,k)=0.0_rprec
          endif
       enddo
    elseif(ic_type.eq.'log_law')then
       ! log law ============================================
       do k=1,dcp%xsz(3)
          !define height
          z=(dcp%xst(3)+k-1.5_rprec)*dz
          if(z.ge.dspl_h)then
             ! impose log law (vectors as lbc_z0 might vary in space)
             arg2(:,:)=(z-dspl_h)/lbc_z0(:,:)
             arg=(1._rprec/KvonK)*dlog(arg2)
             u(:,:,k)=dir_x*arg(:,:)*u_star
             v(:,:,k)=dir_y*arg(:,:)*u_star
             w(:,:,k)=0.0_rprec
          else
             u(:,:,k)=0.0_rprec
             v(:,:,k)=0.0_rprec
             w(:,:,k)=0.0_rprec
          endif
       enddo
    elseif(ic_type.eq.'poiseuille')then
       ! Poiseuille flow  ======================================
       do k=1,dcp%xsz(3)
          z=(dcp%xst(3)+k-1.5_rprec)*dz
          u(:,:,k)=Re*z*(1.0_rprec-0.5_rprec*z)  !non-dimensional
       enddo
       v=0.0_rprec
       w=0.0_rprec
    elseif(ic_type.eq.'do_nothing')then
       ! nothing to do ======================================
       return
    elseif(ic_type.eq.'files')then
       ! load form files ======================================
       inquire(iolength=reclen) vel_in

       open(98,file='vel_u.dat',form='unformatted',access='direct',recl=reclen)
       read(98,rec=1) vel_in
       do k=1,dcp%xsz(3)
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1)  
                u(i,j,k)=real(vel_in(dcp%xst(1)-1+i,dcp%xst(2)-1+j,dcp%xst(3)-1+k),rprec)
             enddo
          enddo
       enddo
       close(98)

       open(98,file='vel_v.dat',form='unformatted',access='direct',recl=reclen)
       read(98,rec=1) vel_in
       do k=1,dcp%xsz(3)
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1)  
                v(i,j,k)=real(vel_in(dcp%xst(1)-1+i,dcp%xst(2)-1+j,dcp%xst(3)-1+k),rprec)
             enddo
          enddo
       enddo
       close(98)

       open(98,file='vel_w.dat',form='unformatted',access='direct',recl=reclen)
       read(98,rec=1) vel_in
       do k=1,dcp%xsz(3)
          do j=1,dcp%xsz(2)
             do i=1,dcp%xsz(1)  
                w(i,j,k)=real(vel_in(dcp%xst(1)-1+i,dcp%xst(2)-1+j,dcp%xst(3)-1+k),rprec)
             enddo
          enddo
       enddo
       close(98)

    else
       errorcode = 1
       call decomp_2d_abort(errorcode,'Invalid initial flow condition type')
    endif

    if(ic_type.ne.'files')then
       ! add noise to initial conditions
       call add_noise1(u,u_star,z_turb*lz,dcp)
       call add_noise1(v,u_star,z_turb*lz,dcp)
       call add_noise1(w,u_star,z_turb*lz,dcp)

       !call add_noise2(u,u_star,w_star,z_turb*lz,dcp)
       !call add_noise2(v,u_star,w_star,z_turb*lz,dcp)
       !call add_noise2(w,u_star,w_star,z_turb*lz,dcp)
    endif


    ! update ghost cells after initial conditions
    call update_ghost(u,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(v,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)
    call update_ghost(w,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3),dcp)

    ! boundary conditions
    if(dcp%xst(3)==1)then
       call compute_bc('bot',lbc_u,lbc_v,lbc_w,dcp)

       call set_bc('bot',u,lbc%typ_u,lbc_u,.true.,dcp)
       call set_bc('bot',v,lbc%typ_v,lbc_v,.true.,dcp)
       call set_bc('bot',w,lbc%typ_w,lbc_w,.false.,dcp)
    endif
    if(dcp%xen(3)==nz)then
       call compute_bc('top',ubc_u,ubc_v,ubc_w,dcp)

       call set_bc('top',u,ubc%typ_u,ubc_u,.true.,dcp)
       call set_bc('top',v,ubc%typ_v,ubc_v,.true.,dcp)
       call set_bc('top',w,ubc%typ_w,ubc_w,.false.,dcp)
    endif

    ! set BOGUS to useless nodes
    if(dcp%xst(3)==1)then
       w(:,:,0)=BOGUS
    endif
    !if(dcp%xen(3)==nz)then
    !   u(:,:,dcp%xsz(3)+1)=BOGUS
    !   v(:,:,dcp%xsz(3)+1)=BOGUS   
    !   w(:,:,dcp%xsz(3)+1)=BOGUS
    !endif

    call output_initial_conditions(dcp)

    return
  end subroutine set_initial_conditions

  !------------------------------------------------------------------------------
  !! subroutine: add noise to initial conditions
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine output_initial_conditions(dcp)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp

    character(512) :: filename
    integer :: ierror, fh
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
  
    ! open file for IO system
    write(filename,'(A)') trim(out_path)//'initial_conditions.dat'

    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    
    call decomp_2d_write_var(fh,disp,1,u(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,v(:,:,1:dcp%xsz(3)),dcp)
    call decomp_2d_write_var(fh,disp,1,w(:,:,1:dcp%xsz(3)),dcp)
    
    call MPI_FILE_CLOSE(fh,ierror)

  end subroutine output_initial_conditions

  !------------------------------------------------------------------------------
  !! subroutine: add noise to initial conditions
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine add_noise1(x,ustar,zturb,dcp)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(inout) :: x
    real(rprec),intent(in) :: ustar,zturb

    real(rprec) :: noise,z,ran3,variance_fact,avg_noise
    integer :: i,i2,j,j2,k,k_glob
    integer :: seed,id
    integer,save :: rnd=0

    if(ustar.gt.0.0_rprec)then
       rnd=rnd+4
       id=0
       avg_noise=0
       do k=1,dcp%xsz(3)
          k_glob=dcp%xst(3)+k-1
          z=(k_glob-0.5_rprec)*dz
          if(z.le.zturb)then
             seed=-80-k_glob-rnd
             !noise=ran3(seed)*2.0_rprec-1.0_rprec
             !write (*,*) nrank, k, k_glob, seed
             do j=1,ny
                do i=1,nx
                   noise=ran3(seed)*2.0_rprec-1.0_rprec
                   if((i.ge.dcp%xst(1)).and.(i.le.dcp%xen(1)).and. &
                        (j.ge.dcp%xst(2)).and.(j.le.dcp%xen(2)))then
                      i2=i-dcp%xst(1)+1;j2=j-dcp%xst(2)+1;
                      variance_fact=3.0_rprec*vnf*u_star**2
                      x(i2,j2,k)=x(i2,j2,k)+variance_fact*(noise*(lz-z*0.9_rprec)/lz)
                   endif
                enddo
             enddo
          else
             ! no turbulence...
          endif
       enddo
    endif

    return
  end subroutine add_noise1
  !------------------------------------------------------------------------------
  !! subroutine: add noise to initial conditions
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine add_noise2(x,ustar,wstar,zturb,dcp)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(inout) :: x
    real(rprec),intent(in) :: ustar,wstar,zturb

    real(rprec) :: noise,z,ran3,variance_fact,avg_noise,rms
    integer :: i,j,k,k_glob
    integer :: seed,id
    integer,save :: rnd=0

    rms = 3._rprec
    rnd=rnd+4
    id=0
    avg_noise=0
    do k=1,dcp%xsz(3)
       k_glob=dcp%xst(3)+k-1
       z=(k_glob-0.5_rprec)*dz
       seed=-80-k_glob
       !noise=ran3(seed)*2.0_rprec-1.0_rprec
       !write (*,*) nrank, k, k_glob, seed
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             if(z.le.zturb)then
                noise=rms/.289_rprec*(ran3(seed)-.5_rprec)
                x(i,j,k)=x(i,j,k)+noise*(1._rprec-z)*wstar
             else
                noise=rms/.289_rprec*(ran3(seed)-0.5_rprec)
                x(i,j,k)=x(i,j,k)+noise*wstar*.01_rprec
             endif
          enddo
       enddo
    enddo
 
    return
  end subroutine add_noise2

end module initial_conditions
