!------------------------------------------------------------------------------
!! module: parameters_IO
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------
module parameters_IO

  implicit none
  include 'mpif.h'

  !! - code global parameters
  integer,parameter,public :: rprec=kind(1.d0)
  logical,parameter,public :: VERBOSE=.true.
  logical,parameter,public :: DEBUG_TIMELOOP=.false.
#ifdef DNS_MOD
  logical,parameter,public :: SGS=.false.
#else
  logical,parameter,public :: SGS=.true.
#endif

  ! define info type for momentum boundary conditions
  type, public :: MOM_BC_INFO
     character(25) :: model,special
     character(25) :: typ_u,typ_v,typ_w
     real(rprec) :: u,v,w
  end type MOM_BC_INFO

  ! define info type for scalar
  type, public :: SC_BC_INFO
     real(rprec) :: s_init,dsdz_init
     real(rprec) :: nf,z_turb,z_ent
     real(rprec) :: ent_s_init,ent_dsdz_init
     character(25) :: lbc_typ,ubc_typ
     real(rprec) :: lbc,ubc
  end type SC_BC_INFO

  !! - useful constants
  real(rprec),parameter,public :: pi=3.1415926535897932384626433_rprec
  complex(rprec),parameter,public :: eye =(0._rprec,1._rprec)
  real(rprec),parameter,public :: g=9.81_rprec
  real(rprec),parameter,public :: F_rho=1.0_rprec
  real(rprec),parameter,public :: KvonK=0.4_rprec 
  real(rprec),parameter,public :: fcoriol=1.067E-04 !for latitude = 47.02N - DN
  integer,parameter,public :: IBOGUS=-1234567890
  real(rprec),parameter,public :: BOGUS=-1.23456789E20_rprec

  !! - numerical parameters 
  integer,public :: p_on_y,p_on_z

  integer,public :: nx,ny,nz
  integer,public :: nt
  integer,public :: lhx,lhy,ldx,ldy
  integer,public :: nx2,ny2

  real(rprec),public :: dx,dy,dz
  real(rprec),public :: dt
  logical :: dim_dt_flag
  logical,public :: restart_flag  

  !! - physical parameters
  real(rprec),public :: Re
  real(rprec),public :: u_scale,zi
  real(rprec),public :: g_hat
  logical :: XY_pi
  real(rprec),public :: lx,ly,lz 
  real(rprec),public :: pfx,pfy
  real(rprec),public :: dt_dim

  !! -SGS parameters
  real(rprec),public :: Cs0
  character(16),public :: sgs_model="cst_Smag"
  integer,public :: cs_count,dyn_init
  logical,public :: Lagran_init

  ! Coriolis focring
  logical,public :: coriolis_flag,geostr_rot
  real(rprec),public :: ug_dim,vg_dim
  integer,public :: geostr_start
  real(rprec),public :: ug,vg
  real(rprec),public :: coriol

  !! - initial conditions parameters 
  character(16),public :: ic_type
  real(rprec),public :: u0,v0,w0,u_star,dudz0,vnf,z_turb

  !! - boundary conditions parameters 
  TYPE(MOM_BC_INFO),public :: lbc
  TYPE(MOM_BC_INFO),public :: ubc

  character(25) :: lbc_model,lbc_special
  character(25) :: ubc_model,ubc_special
  logical,public :: surface_files

  !! - roughness length and wall conditions
  character(16),public :: wall_model
  real(rprec),public :: lbc_dz0,dspl_h,gt_u,gt_v,wall_c

  !! - output and stats parameters
  character(256),public :: out_path  
  integer,public :: checkpoint,scheckpoint
  logical,public :: output
  integer,public :: step_out
  logical,public :: stat_flag 
  integer,public :: start_stats,running_period
  logical,public :: tke_flag
  integer,public :: start_tke,tke_period
  logical,public :: spectra_flag=.false.

  private
  public :: read_parameters
  public :: write_parameters_restart,read_parameters_restart

contains !=======================================================================

  !------------------------------------------------------------------------------
  !! subroutine: read parameters
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!  - parameters are read form an input files using namelist (see example)
  !------------------------------------------------------------------------------
  subroutine read_parameters(me,nproc,ierror)
    implicit none

    integer, intent(inout) :: me,nproc,ierror
    integer :: op_len
    integer :: open_status,close_status,nlst_status
    integer :: errorcode


    ! lower boundary condition
    character(25) :: lbc_typ_u,lbc_typ_v,lbc_typ_w
    real(rprec) :: lbc_u,lbc_v,lbc_w

    ! upper boundary condition
    character(25) :: ubc_typ_u,ubc_typ_v,ubc_typ_w
    real(rprec) :: ubc_u,ubc_v,ubc_w

    lbc_typ_u='drc';lbc_typ_v='drc';lbc_typ_w='drc'
    lbc_u=0.0_rprec;lbc_v=0.0_rprec;lbc_w=0.0_rprec

    ubc_typ_u='neu';ubc_typ_v='neu';ubc_typ_w='drc'
    ubc_u=0.0_rprec;ubc_v=0.0_rprec;ubc_w=0.0_rprec

    namelist / simulation_param / p_on_y,p_on_z,nx,ny,nz,nt,dt,dim_dt_flag,restart_flag

    namelist / output_param / out_path,&
         output,step_out,&
         stat_flag,start_stats,running_period,&
         tke_flag,start_tke,tke_period,&
         checkpoint,scheckpoint

    namelist / physical_param / XY_pi,lx,ly,lz,u_scale,zi,lbc_dz0,Re,&
         pfx,pfy,coriolis_flag,ug_dim,vg_dim,geostr_rot,geostr_start

    namelist / sgs_param / Cs0,sgs_model,cs_count,dyn_init

    namelist / momentum_bc / lbc_model,lbc_special,ubc_model,ubc_special,wall_model,surface_files

    namelist / momentum_ic / ic_type,u0,v0,w0,u_star,dudz0,z_turb,vnf

    if(me==0)then
       write(*,*) '## Reading parameters from parameters.nlst'

       open(88,file="parameters.nlst",status='OLD',recl=80,delim='APOSTROPHE',&
            iostat=open_status,action='READ')
       if(open_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Could not open parameters.nlst for reading. (ioerror=',open_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       read(88,simulation_param,iostat=nlst_status)
       if(nlst_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Error in namelist simulation_param. (ioerror=',nlst_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       read(88,output_param,iostat=nlst_status)
       if(nlst_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Error in namelist ouptut_param. (ioerror=',nlst_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       read(88,physical_param,iostat=nlst_status)
       if(nlst_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Error in namelist physical_param. (ioerror=',nlst_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       read(88,sgs_param,iostat=nlst_status)
       if(nlst_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Error in namelist sgs_param. (ioerror=',nlst_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       read(88,momentum_bc,iostat=nlst_status)
       if(nlst_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Error in namelist momentum_bc. (ioerror=',nlst_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       read(88,momentum_ic,iostat=nlst_status)
       if(nlst_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Error in namelist momentum_ic. (ioerror=',nlst_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       close(88,iostat=close_status)
       if ( open_status /= 0 ) then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Could not close parameters.nlst.(ioerror=',open_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

    endif
    call broadcast_parameters(me,nproc,ierror)
    call MPI_barrier(MPI_COMM_WORLD, ierror)

    out_path=trim(out_path)
    op_len=len(out_path)
    if(out_path(op_len:).ne.'/')then
       out_path=out_path//'/'
    endif

    if(XY_pi)then
       lx=lx*pi;ly=ly*pi
    endif

    !nondimensionalize domain size
    lx=lx/zi;ly=ly/zi;lz=lz/zi;

    !nondimensionalize physical constant
    g_hat=g/(u_scale*u_scale/zi)

    ! Coriolis focring
    if(coriolis_flag)then
       ug=ug_dim/u_scale
       vg=vg_dim/u_scale
       coriol=1.067E-04*zi/u_scale 
    else
       ug=0.0_rprec
       vg=0.0_rprec
       coriol=0.0_rprec
    endif

    if(dim_dt_flag)then
       dt_dim = dt
       dt = dt*u_scale/zi
    else
       dt_dim = dt/u_scale*zi
    endif

    lhx=nx/2+1;lhy=ny/2+1
    ldx=2*lhx;ldy=2*lhy
    nx2=3*nx/2;ny2=3*ny/2

    dx=lx/real(nx,rprec)
    dy=ly/real(ny,rprec)
    dz=lz/real(nz,rprec)

    lbc%model=lbc_model
    lbc%special=lbc_special
    lbc%typ_u=lbc_typ_u
    lbc%typ_v=lbc_typ_v
    lbc%typ_w=lbc_typ_w
    lbc%u=lbc_u
    lbc%v=lbc_v
    lbc%w=lbc_w

    ubc%model=ubc_model
    ubc%special=ubc_special
    ubc%typ_u=ubc_typ_u
    ubc%typ_v=ubc_typ_v
    ubc%typ_w=ubc_typ_w
    ubc%u=ubc_u
    ubc%v=ubc_v
    ubc%w=ubc_w

    if(me==0)then
       write(*,*) '## Parameters summary'
       write(*,*) '=============================================================================='
       write(*,simulation_param)
       write(*,*) '=============================================================================='
       write(*,output_param)
       write(*,*) '=============================================================================='
       write(*,physical_param)
       write(*,*) '=============================================================================='
       write(*,sgs_param)
       write(*,*) '=============================================================================='
       write(*,momentum_bc)
       write(*,*) '=============================================================================='
       write(*,momentum_ic)
       write(*,*) '=============================================================================='   
       write(*,*) ' '       
    endif

    ! needed so all param are set on all nodes
    call MPI_barrier(MPI_COMM_WORLD, ierror)

    return
  end subroutine read_parameters

  !------------------------------------------------------------------------------
  !! subroutine: write parameters for restart
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine write_parameters_restart(jt_total,me)
    implicit none

    integer,intent(IN) :: jt_total,me
    integer :: p_on_y_f,p_on_z_f,nx_f,ny_f,nz_f
    integer :: open_status,close_status,nlst_status,errorcode,ierror

    p_on_y_f=p_on_y;p_on_z_f=p_on_z;
    nx_f=nx;ny_f=ny;nz_f=nz;

    namelist / restart_param / p_on_y_f,p_on_z_f,nx_f,ny_f,nz_f,jt_total,dyn_init,Lagran_init,sgs_model

    if(me==0)then
       write(*,*) '## Parameters summary ...'

       open(99,file="parameters_restart.nlst",iostat=open_status,action='WRITE')
       if(open_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Could not open parameters_restart.nlst for reading. (ioerror=',open_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       write(99,restart_param)

       close(99,iostat=close_status)
       if ( open_status /= 0 ) then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Could not close parameters_restar.nlst.(ioerror=',open_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif
    endif

    call MPI_barrier(MPI_COMM_WORLD, ierror)

    return
  end subroutine write_parameters_restart

  !------------------------------------------------------------------------------
  !! subroutine: read parameters for restart
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !! - description:
  !!
  !!
  !------------------------------------------------------------------------------
  subroutine read_parameters_restart(jt_total,me)
    implicit none

    integer,intent(IN) :: me
    integer,intent(OUT) :: jt_total
    integer :: p_on_y_f,p_on_z_f,nx_f,ny_f,nz_f

    integer :: open_status,close_status,nlst_status
    integer :: errorcode,ierror
    namelist / restart_param / p_on_y_f,p_on_z_f,nx_f,ny_f,nz_f,jt_total,dyn_init,Lagran_init,sgs_model

    if(me==0)then
       open(88,file="parameters_restart.nlst",status='OLD',recl=80,delim='APOSTROPHE',&
            iostat=open_status,action='READ')
       if(open_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Could not open parameters_restart.nlst for reading. (ioerror=',open_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       read(88,restart_param,iostat=nlst_status)
       if(nlst_status/=0)then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Error in namelist restart_param. (ioerror=',nlst_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       close(88,iostat=close_status)
       if ( open_status /= 0 ) then
          write(*,'(A,I3.0,A)')&
               '[ParametresIO error] Could not close parameters_restar.nlst.(ioerror=',open_status,')'
          errorcode=22
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif

       ! checking parameters
       if((p_on_y.ne.p_on_y_f) .or.(p_on_y.ne.p_on_y_f) .or. &
            (nx.ne.nx_f) .or. (nz.ne.nz_f) .or. (nz.ne.nz_f))then
          write(*,'(A)')&
               '[Restart error] invalid restart parameters'
          errorcode=1
          call MPI_Abort(MPI_COMM_WORLD,errorcode)
       endif
    endif
 
    call MPI_BCAST(jt_total,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(dyn_init,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(Lagran_init,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(sgs_model, 16, MPI_CHARACTER,0,MPI_COMM_WORLD,ierror)
    call MPI_barrier(MPI_COMM_WORLD, ierror)

    return
  end subroutine read_parameters_restart

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !!
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !!
  !------------------------------------------------------------------------------
  subroutine broadcast_parameters(me,nproc,ierr)
    
    implicit none
    integer, intent (inout) :: ierr,me,nproc
    
    integer,parameter :: bufferlength=1500
    integer :: position
    character :: buffer(bufferlength)

    if(me==0)then
       ! on node 0
       position = 0
       
       !namelist / simulation_param
       call MPI_pack(p_on_y, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(p_on_z, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(nx, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ny, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(nz, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(nt, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(dt, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(dim_dt_flag, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(restart_flag, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)

       !namelist / physical_param
       call MPI_pack(XY_pi, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(lx, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ly, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(lz, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(u_scale, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(zi, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(lbc_dz0, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(Re, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(pfx, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(pfy, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(coriolis_flag, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ug_dim, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(vg_dim, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(geostr_rot, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(geostr_start, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)

       !namelist / sgs_param
       call MPI_pack(Cs0, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(sgs_model, 16, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(cs_count, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(dyn_init, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)

       !namelist / momentum_bc
       call MPI_pack(lbc_model, 25, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(lbc_special, 25, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ubc_model, 25, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(ubc_special, 25, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(wall_model, 16, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(surface_files, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)

       !namelist / momentum_ic
       call MPI_pack(ic_type, 16, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(u0, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(v0, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(w0, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(u_star, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(dudz0, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(z_turb, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(vnf, 1, MPI_DOUBLE_PRECISION, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
    else
       ! on the other nodes
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

       position = 0

       !namelist / simulation_param
       call MPI_unpack(buffer, bufferlength, position, p_on_y, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, p_on_z, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, nx, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ny, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, nz, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, nt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, dt, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, dim_dt_flag, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, restart_flag, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)

       !namelist /  physical_param
       call MPI_unpack(buffer, bufferlength, position, XY_pi, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, lx, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ly, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, lz, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, u_scale, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, zi, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, lbc_dz0, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, Re, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, pfx, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, pfy, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, coriolis_flag, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ug_dim, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, vg_dim, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, geostr_rot, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, geostr_start, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
                           
       !namelist / sgs_param 
       call MPI_unpack(buffer, bufferlength, position, Cs0, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, sgs_model, 16, MPI_CHARACTER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, cs_count, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, dyn_init, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       
       !namelist / momentum_bc
       call MPI_unpack(buffer, bufferlength, position, lbc_model, 25, MPI_CHARACTER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, lbc_special, 25, MPI_CHARACTER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ubc_model, 25, MPI_CHARACTER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, ubc_special, 25, MPI_CHARACTER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, wall_model, 16, MPI_CHARACTER, MPI_COMM_WORLD, ierr)       
       call MPI_unpack(buffer, bufferlength, position, surface_files, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)

       !namelist / momentum_ic
       call MPI_unpack(buffer, bufferlength, position, ic_type, 16, MPI_CHARACTER, MPI_COMM_WORLD, ierr)       
       call MPI_unpack(buffer, bufferlength, position, u0, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, v0, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)      
       call MPI_unpack(buffer, bufferlength, position, w0, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, u_star, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, dudz0, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)       
       call MPI_unpack(buffer, bufferlength, position, z_turb, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, vnf, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)       
    endif

    if(me==0)then
       !! on node 0
       position = 0

       !namelist / output_param 
       call MPI_pack(out_path, 256, MPI_CHARACTER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(output, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(step_out, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(stat_flag, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(start_stats, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(running_period, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(tke_flag, 1, MPI_LOGICAL, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(start_tke, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(tke_period, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(checkpoint, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)
       call MPI_pack(scheckpoint, 1, MPI_INTEGER, buffer, bufferlength, position, MPI_COMM_WORLD, ierr)

       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)
    else
       !! on the other nodes
       call MPI_bcast(buffer, bufferlength, MPI_PACKED, 0, MPI_COMM_WORLD, ierr)

       position = 0

       !namelist / output_param 
       call MPI_unpack(buffer, bufferlength, position, out_path, 256, MPI_CHARACTER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, output, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, step_out, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, stat_flag, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, start_stats, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, running_period, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, tke_flag, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, start_tke, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, tke_period, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, checkpoint, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_unpack(buffer, bufferlength, position, scheckpoint, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    endif

    

  end subroutine broadcast_parameters

end module parameters_IO
