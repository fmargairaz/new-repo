module mod_IO_statistics

  use decomp_2d
  use decomp_2d_custom
  use decomp_2d_io
  use utility_tools

  use system_variables, only : u,v,w,p,txx,txy,txz,tyy,tyz,tzz,dudz,dvdz
  use compute_wall_law, only : u_star_avg
  use compute_sgs, only : nu_t,cs_opt2

  use parameters_IO

  implicit none
  include 'mpif.h'

  private
  TYPE(DECOMP_INFO),save :: dcp
  integer, save :: iterations_av 
  
  ! variables for running averages
  real(rprec),dimension(:,:,:),allocatable :: Mav_u,Mav_v,Mav_w,Mav_p
  real(rprec),dimension(:,:,:),allocatable :: Mav_u2,Mav_v2,Mav_w2
  real(rprec),dimension(:,:,:),allocatable :: Mav_u3,Mav_v3,Mav_w3
  real(rprec),dimension(:,:,:),allocatable :: Mav_u4,Mav_v4,Mav_w4
  real(rprec),dimension(:,:,:),allocatable :: Mav_uv,Mav_uw,Mav_vw
  real(rprec),dimension(:,:,:),allocatable :: Mav_txx, Mav_tyy, Mav_tzz
  real(rprec),dimension(:,:,:),allocatable :: Mav_txy, Mav_txz, Mav_tyz
  real(rprec),dimension(:,:,:),allocatable :: Mav_dudz,Mav_dvdz
  real(rprec),dimension(:,:,:),allocatable :: Mav_nut,Mav_cs

  real(rprec),dimension(:,:),allocatable :: Mav_ustar
  
  public :: io_stat_init,io_stat_finalize
  public :: io_stat_averages
  public :: io_stat_spectra

contains !=======================================================================

  subroutine io_stat_init(dcp_main)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: dcp_main

    integer :: status, errorcode

    dcp = dcp_main

    ! u,v,w,p
    allocate(Mav_u(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       write(*,*) status
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variable: u')
    end if
    allocate(Mav_v(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variable: v')
    end if
    allocate(Mav_w(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variable: w')
    end if
    allocate(Mav_p(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variable: p')
    end if
    
    ! u2,v2,w2
    allocate(Mav_u2(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_v2(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_w2(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if

    ! u3,v3,w3
    allocate(Mav_u3(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_v3(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_w3(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if

    ! u4,v4,w4
    allocate(Mav_u4(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_v4(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_w4(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if

    ! uv,uw,vw
    allocate(Mav_uv(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_uw(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_vw(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    
    ! txx,tyy,tzz,txy,txz,tyz
    allocate(Mav_txx(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_tyy(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_tzz(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_txy(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_txz(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_tyz(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    
    ! dudz,dudz
    allocate(Mav_dudz(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_dvdz(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if

    ! dudz,dudz
    allocate(Mav_nut(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if
    allocate(Mav_cs(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variables')
    end if

    ! 2D ground array
    allocate(Mav_ustar(dcp%xsz(1),dcp%xsz(2)), STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising io module variable: ustar')
    end if

    ! initialisation of variable to 0.0
    Mav_u = 0.d0; Mav_v = 0.d0; Mav_w = 0.d0; Mav_p = 0.d0;
    Mav_u2 = 0.d0; Mav_v2 = 0.d0; Mav_w2 = 0.d0
    Mav_u3 = 0.d0; Mav_v3 = 0.d0; Mav_w3 = 0.d0
    Mav_u4 = 0.d0; Mav_v4 = 0.d0; Mav_w4 = 0.d0
    Mav_uv = 0.d0; Mav_uw = 0.d0; Mav_vw = 0.d0
    Mav_txx = 0.d0; Mav_tyy = 0.d0; Mav_tzz = 0.d0
    Mav_txy = 0.d0; Mav_txz = 0.d0; Mav_tyz = 0.d0
    Mav_dudz = 0.d0; Mav_dvdz = 0.d0;
    Mav_nut = 0.d0; Mav_cs = 0.d0;
    Mav_ustar = 0.d0;
    iterations_av = 0 
    
    return
  end subroutine io_stat_init
  
  subroutine io_stat_finalize()

    deallocate(Mav_u,Mav_v,Mav_w,Mav_p)
    deallocate(Mav_u2,Mav_v2,Mav_w2)
    deallocate(Mav_u3,Mav_v3,Mav_w3,Mav_u4,Mav_v4,Mav_w4)
    deallocate(Mav_uv,Mav_uw,Mav_vw)
    deallocate(Mav_txx,Mav_tyy,Mav_tzz,Mav_txy, Mav_txz, Mav_tyz)
    deallocate(Mav_dudz,Mav_dvdz)
    deallocate(Mav_nut,Mav_cs)

    deallocate(Mav_ustar)

    return
  end subroutine io_stat_finalize

  !------------------------------------------------------------------------------
  !! subroutine: io main running averages
  !------------------------------------------------------------------------------
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !------------------------------------------------------------------------------
  subroutine io_stat_averages(jt,jt_total)
    implicit none

    integer,intent(in) :: jt,jt_total
    
    character(512) :: filename 

    iterations_av = iterations_av + 1
    
    if((mod(jt_total,running_period)==0).and.(jt_total.gt.start_stats))then
       if(nrank==0)then
          write(*,*) '=============================================================================='
          write(*,'(A,I12)') '# writing momentum stats at jt_t=',jt_total
          write(*,*) '=============================================================================='

          write(filename,'(A,I9.9,A)') trim(out_path)//'RAV_',jt_total,'_info.txt'
          open(99,file=filename)
          write(99,'(I6,I6,I6,I10)') nx,ny,nz,iterations_av
          close(99)
       endif
       iterations_av = 0  
    endif

    call io_2D_stat_averages(jt,jt_total)
    call io_3D_stat_averages(jt,jt_total)
    
    return 
  end subroutine io_stat_averages
  
  !------------------------------------------------------------------------------
  !! subroutine: io 2D running averages
  !------------------------------------------------------------------------------
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !------------------------------------------------------------------------------
  subroutine io_2D_stat_averages(jt,jt_total)
    implicit none

    integer,intent(in) :: jt,jt_total
    
    integer :: i,j,k,k_min,k_max
    character(512) :: filename 

    integer :: ierror, fh
    integer, dimension(2) :: sizes, subsizes, starts
    integer :: newtype, data_type
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp

    sizes(1) = dcp%xsz(1)
    sizes(2) = dcp%ysz(2)
    subsizes(1) = dcp%xsz(1)
    subsizes(2) = dcp%xsz(2)
    starts(1) = 0
    starts(2) = dcp%xst(2)-1

    if(dcp%xst(3)==1)then
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1) 
             Mav_ustar(i,j) = Mav_ustar(i,j) + u_star_avg(i,j)
          enddo
       enddo

       if((mod(jt_total,running_period)==0).and.(jt_total.gt.start_stats))then
          ! open file for IO system
          write(filename,'(A,I9.9,A)') trim(out_path)//'RAV_',jt_total,'_2D.dat'
 
          call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts,  &
               MPI_ORDER_FORTRAN, real_type, newtype, ierror)
          call MPI_TYPE_COMMIT(newtype,ierror)         

          call MPI_FILE_OPEN(DECOMP_2D_LAYER_X, filename, &
               MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
               fh, ierror)
          
          filesize = 0_MPI_OFFSET_KIND
          call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
          disp = 0_MPI_OFFSET_KIND
          call MPI_FILE_SET_VIEW(fh,disp,real_type, &
               newtype,'native',MPI_INFO_NULL,ierror)
          
          call MPI_FILE_WRITE_ALL(fh, Mav_ustar, &
               subsizes(1)*subsizes(2),real_type, MPI_STATUS_IGNORE, ierror)
          
          call MPI_FILE_CLOSE(fh,ierror)
          call MPI_TYPE_FREE(newtype,ierror)
          
          Mav_ustar = 0.d0
       endif
    endif
    
    return
  end subroutine io_2D_stat_averages

  !------------------------------------------------------------------------------
  !! subroutine: io 3D running averages
  !------------------------------------------------------------------------------
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !------------------------------------------------------------------------------
  subroutine io_3D_stat_averages(jt,jt_total)
    implicit none

    integer,intent(in) :: jt,jt_total
    
    integer :: i,j,k,k_min,k_max
    character(512) :: filename 

    integer :: ierror, fh
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp

    if(dcp%xst(3)==1)then
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1) 
             Mav_uw(i,j,1) = 0.d0
             Mav_vw(i,j,1) = 0.d0
          enddo
       enddo
       k_min=2
    else
       k_min=1
    endif

    !$omp parallel
    !$omp do
    do k=k_min,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             !For computing the vertical Reynolds stresses we, compute "u" and "v"
             !at the "w" node. This way we can add, the reynolds stresses dirctly to the
             !SGSTxz stresses, which arec omputed at the 'w' nodes.
             !(At the first node, there is no interpolation needed, it is '0' anway.)
             Mav_uw(i,j,k) = Mav_uw(i,j,k) + ( w(i,j,k)*(0.5d0*(u(i,j,k-1)+u(i,j,k))) )
             Mav_vw(i,j,k) = Mav_vw(i,j,k) + ( w(i,j,k)*(0.5d0*(v(i,j,k-1)+v(i,j,k))) )
          enddo
       enddo
    enddo
    !$omp end do
    
    !$omp do
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)    
             !For computing Mean values, and variances:
             Mav_u(i,j,k) = Mav_u(i,j,k) + u(i,j,k)
             Mav_v(i,j,k) = Mav_v(i,j,k) + v(i,j,k)
             Mav_w(i,j,k) = Mav_w(i,j,k) + w(i,j,k)
             Mav_p(i,j,k) = Mav_p(i,j,k) + p(i,j,k)

             Mav_u2(i,j,k) = Mav_u2(i,j,k) + (u(i,j,k)**2)
             Mav_v2(i,j,k) = Mav_v2(i,j,k) + (v(i,j,k)**2)
             Mav_w2(i,j,k) = Mav_w2(i,j,k) + (w(i,j,k)**2)
             
             Mav_u3(i,j,k) = Mav_u3(i,j,k) + (u(i,j,k)**3)
             Mav_v3(i,j,k) = Mav_v3(i,j,k) + (v(i,j,k)**3)
             Mav_w3(i,j,k) = Mav_w3(i,j,k) + (w(i,j,k)**3)

             Mav_u4(i,j,k) = Mav_u4(i,j,k) + (u(i,j,k)**4)
             Mav_v4(i,j,k) = Mav_v4(i,j,k) + (v(i,j,k)**4)
             Mav_w4(i,j,k) = Mav_w4(i,j,k) + (w(i,j,k)**4)
             
             !For computing variances:
             Mav_uv(i,j,k) = Mav_uv(i,j,k) + (u(i,j,k)*v(i,j,k))
             
             !SGS stress tensor, important near the WALL. (On the first grid point = Wall Stress.)
             Mav_txx(i,j,k) = Mav_txx(i,j,k) + txx(i,j,k)
             Mav_tyy(i,j,k) = Mav_tyy(i,j,k) + tyy(i,j,k)
             Mav_tzz(i,j,k) = Mav_tzz(i,j,k) + tzz(i,j,k)
             Mav_txy(i,j,k) = Mav_txy(i,j,k) + txy(i,j,k)
             Mav_txz(i,j,k) = Mav_txz(i,j,k) + txz(i,j,k)
             Mav_tyz(i,j,k) = Mav_tyz(i,j,k) + tyz(i,j,k)
             
             Mav_dudz(i,j,k) = Mav_dudz(i,j,k) + dudz(i,j,k)
             Mav_dvdz(i,j,k) = Mav_dvdz(i,j,k) + dvdz(i,j,k)

             Mav_nut(i,j,k) = Mav_nut(i,j,k) + nu_t(i,j,k)
             Mav_cs(i,j,k) = Mav_cs(i,j,k) + sqrt(cs_opt2(i,j,k))
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    if((mod(jt_total,running_period)==0).and.(jt_total.gt.start_stats))then
       ! open file for IO system
       write(filename,'(A,I9.9,A)') trim(out_path)//'RAV_',jt_total,'_3D.dat'
       
       call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND

       ! u,v,w,p
       call decomp_2d_write_var(fh,disp,1,Mav_u,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_v,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_w,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_p,dcp)
       
       ! u2,v2,w2
       call decomp_2d_write_var(fh,disp,1,Mav_u2,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_v2,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_w2,dcp)

       ! u3,v3,w3
       call decomp_2d_write_var(fh,disp,1,Mav_u3,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_v3,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_w3,dcp)

       ! u4,v4,w3
       call decomp_2d_write_var(fh,disp,1,Mav_u4,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_v4,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_w4,dcp)

       ! uv,uw,vw
       call decomp_2d_write_var(fh,disp,1,Mav_uv,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_uw,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_vw,dcp)
       
       ! txx,tyy,tzz,txy,txz,tyz
       call decomp_2d_write_var(fh,disp,1,Mav_txx,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_tyy,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_tzz,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_txy,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_txz,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_tyz,dcp)

       call decomp_2d_write_var(fh,disp,1,Mav_dudz,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_dvdz,dcp)
       
       call decomp_2d_write_var(fh,disp,1,Mav_nut,dcp)
       call decomp_2d_write_var(fh,disp,1,Mav_cs,dcp)
       
       call MPI_FILE_CLOSE(fh,ierror)

       !All variables are cleaned, so next output starts from scrsatch again.
       Mav_u = 0.d0; Mav_v = 0.d0; Mav_w = 0.d0; Mav_p = 0.d0;
       Mav_u2 = 0.d0; Mav_v2 = 0.d0; Mav_w2 = 0.d0
       Mav_u3 = 0.d0; Mav_v3 = 0.d0; Mav_w3 = 0.d0
       Mav_u4 = 0.d0; Mav_v4 = 0.d0; Mav_w4 = 0.d0
       Mav_uv = 0.d0; Mav_uw = 0.d0; Mav_vw = 0.d0
       Mav_txx = 0.d0; Mav_tyy = 0.d0; Mav_tzz = 0.d0
       Mav_txy = 0.d0; Mav_txz = 0.d0; Mav_tyz = 0.d0
       Mav_dudz = 0.d0; Mav_dvdz = 0.d0;
       Mav_nut = 0.d0; Mav_cs = 0.d0;   
    endif
    
    return
  end subroutine io_3D_stat_averages
  
  !------------------------------------------------------------------------------
  !! subroutine: io spectra
  !------------------------------------------------------------------------------
  !! - fabien margairaz (fabien.margairaz@gmail.com)
  !------------------------------------------------------------------------------
  subroutine io_stat_spectra(jt,jt_total)
    implicit none

    integer,intent(in) :: jt,jt_total

    if((jt_total.ge.1500000).and.(jt_total.le.1502000).and.(mod(jt,10).eq.0))then 
       call io_stat_spectra_out(jt,jt_total)
    endif
    if((jt_total.ge.1550000).and.(jt_total.le.1552000).and.(mod(jt,10).eq.0))then 
       call io_stat_spectra_out(jt,jt_total)
    endif
    if((jt_total.ge.1600000).and.(jt_total.le.1602000).and.(mod(jt,10).eq.0))then 
       call io_stat_spectra_out(jt,jt_total)
    endif
    if((jt_total.ge.1650000).and.(jt_total.le.1652000).and.(mod(jt,10).eq.0))then 
       call io_stat_spectra_out(jt,jt_total)
    endif
  end subroutine io_stat_spectra
  
  subroutine io_stat_spectra_out(jt,jt_total)
    implicit none

    integer,intent(in) :: jt,jt_total

    !real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)) :: u_tmp
    integer :: i,j,k

    character(512) :: filename 
    integer :: ierror, fh
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp

    !if(nrank==0) write(*,*) 'spectra output'

    !!$omp parallel do 
    ! do k=1,dcp%xsz(3)
    !   do j=1,dcp%xsz(2)
    !      do i=1,dcp%xsz(1)
    !         u_tmp(i,j,k) = u(i,j,k)
    !      enddo
    !   enddo
    !enddo
    !!$omp end parallel do

    ! open file for IO system
    write(filename,'(A,I9.9,A)') trim(out_path)//'SP_',jt_total,'.dat'
    !write(*,*) nrank,dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)
    !call decomp_2d_write_one(1,u_tmp,filename,dcp)

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

  end subroutine io_stat_spectra_out

end module mod_IO_statistics
