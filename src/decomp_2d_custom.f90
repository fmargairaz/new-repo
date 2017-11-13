!------------------------------------------------------------------------------
!!  module: decomp_2d_custom
!------------------------------------------------------------------------------
!!
!! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
!!
!! - description: 
!!  - this module is based on 2Decomp library and ghost cell functionality to 
!!    the library
!!
!------------------------------------------------------------------------------

module decomp_2d_custom
  
  use decomp_2d  ! 2D decomposition module

  implicit none
  include 'mpif.h'

  private

  !number of ghost cells
  integer level

  ! define neighboring blocks (to be used for ghost cell)
  ! first dimension 1=X-pencil, 2=Y-pencil,
  ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom 
  integer, save, dimension(2,6) :: neighbour
  integer, save, dimension(2),public :: dims, coord
  integer, allocatable, dimension(:),public :: ranks_proc

  integer, save, public :: DECOMP_2D_LAYER_X
  integer, save, public :: nrank_2d_layer_x
  integer, allocatable, dimension(:,:), public :: ranks_2d_layer_x

  integer, save, public :: DECOMP_2D_VRTCL_X
  integer, save, public :: nrank_2d_vrtcl_x
  integer, allocatable, dimension(:,:), public :: ranks_2d_vrtcl_x

  integer,allocatable,dimension(:,:,:),public :: dcp_ph_sizes,dcp_sp_sizes

  public :: decomp_info_ghost_init, decomp_2d_partition_info
  public :: alloc_x_ghost, alloc_y_ghost, update_ghost, update_top_ghost, update_halo3 

  interface alloc_x_ghost
     module procedure alloc_x_ghost_real
     module procedure alloc_x_ghost_complex
  end interface alloc_x_ghost
  
  interface alloc_y_ghost
     module procedure alloc_y_ghost_real
     module procedure alloc_y_ghost_complex
  end interface alloc_y_ghost
  
  interface update_ghost
     module procedure update_ghost_real
     module procedure update_ghost_complex
  end interface update_ghost

  interface update_top_ghost
     module procedure update_ghost_real_top
     module procedure update_ghost_complex_top
  end interface update_top_ghost

contains !=======================================================================   
  
  !------------------------------------------------------------------------------
  !! subroutine: initialise
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------

  subroutine decomp_info_ghost_init(nx,ny,nz,ghost_cell,p_on_z,decomp,decomp_ghost)
    implicit none

    integer, intent(in) :: nx,ny,nz,ghost_cell,p_on_z
    TYPE(DECOMP_INFO), intent(OUT) :: decomp,decomp_ghost

    logical, dimension(2) :: dummy_periods,dummy_coords
    integer, allocatable, dimension(:) :: tmp_nodes
    integer :: MPI_size,i,icount,ierror
    integer :: new_rank, sendbuf, recvbuf
    integer :: m,n

    ! local MPI group
    integer base_grp, grp1, grp2
    
    ! get MPI size
    call MPI_COMM_SIZE(MPI_COMM_WORLD,MPI_size,ierror)
    
    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, ierror)
    call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
    
    call MPI_COMM_GROUP(MPI_COMM_WORLD,base_grp,ierror)

    ! liste of processor ranks
    allocate(ranks_proc(nproc))
    do m=1,nproc
       ranks_proc(m)=m-1
    enddo

    ! create a group commuicator per layer
    allocate(tmp_nodes(dims(1)))
    
    tmp_nodes(1) = coord(2)
    do i=2,dims(1)
       tmp_nodes(i) = tmp_nodes(i-1) + dims(2)
    enddo

    call MPI_GROUP_INCL(base_grp,dims(1),tmp_nodes,grp1,ierror)
    call MPI_COMM_CREATE(MPI_COMM_WORLD,grp1,DECOMP_2D_LAYER_X,ierror)
    call MPI_GROUP_RANK(grp1,nrank_2d_layer_x,ierror)
        
    deallocate(tmp_nodes)
    
    ! create 2D array with ranks info for 2D_layer_x groups
    allocate(ranks_2d_layer_x(dims(2),dims(1)))
    do m=1,dims(2)
       do n=1,dims(1)
          ranks_2d_layer_x(m,n)=(m-1) + (n-1)*dims(2)
       enddo
    enddo
       
    ! create a group commuicator vertical
    allocate(tmp_nodes(dims(2)))

    tmp_nodes(1) = coord(1)*p_on_z
    do i=2,dims(2)
       tmp_nodes(i) = tmp_nodes(i-1) + 1
    enddo

    call MPI_GROUP_INCL(base_grp,dims(2),tmp_nodes,grp2,ierror)
    call MPI_COMM_CREATE(MPI_COMM_WORLD,grp2,DECOMP_2D_VRTCL_X,ierror)
    call MPI_GROUP_RANK(grp2,nrank_2d_vrtcl_x,ierror)
    
    deallocate(tmp_nodes)

    ! create 2D array with ranks info for 2D_layer_x groups
    allocate(ranks_2d_vrtcl_x(dims(1),dims(2)))
    do m=1,dims(1)
       do n=1,dims(2)
          ranks_2d_vrtcl_x(m,n)=(n-1) + (m-1)*dims(2)
       enddo
    enddo

    !save the number of ghost cells
    level = ghost_cell
    
    ! create deomp info
    call decomp_info_init(nx, ny, nz, decomp)
    call decomp_info_init(nx, ny, nz+2*ghost_cell*p_on_z, decomp_ghost)
    
    ! inssure that ghostcell are ok
    decomp_ghost%xsz(3)=decomp%xsz(3)+2*ghost_cell
    decomp_ghost%ysz(3)=decomp%ysz(3)+2*ghost_cell
 
    !decomp_ghost%xst(3)=decomp%xst(3)-1
    !decomp_ghost%xen(3)=decomp%xen(3)+1
    !decomp_ghost%yst(3)=decomp%yst(3)-1
    !decomp_ghost%yen(3)=decomp%yen(3)+1
    
        
    ! get neighbour for ghost cell comm
    call init_neighbour()

    return
  end subroutine decomp_info_ghost_init

  !------------------------------------------------------------------------------
  !! subroutine: 
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 22/05/2015
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine decomp_2d_partition_info(ph)
    use decomp_2d_fft_2d, only : decomp_2d_fft_2d_get_sp_info
    implicit none

    TYPE(DECOMP_INFO),intent(IN) :: ph
    TYPE(DECOMP_INFO) :: sp
    integer, dimension(3,3) :: dcp_size_local
    integer :: i,ierror

    dcp_size_local(1,:)=ph%xsz(:)
    dcp_size_local(2,:)=ph%ysz(:)
    dcp_size_local(3,:)=ph%zsz(:)
    allocate(dcp_ph_sizes(3,3,nproc))
    call MPI_ALLGATHER(dcp_size_local,9,MPI_INT,dcp_ph_sizes,9,MPI_INT,MPI_COMM_WORLD,ierror)

    call decomp_2d_fft_2d_get_sp_info(sp)

    dcp_size_local(1,:)=sp%xsz(:)
    dcp_size_local(2,:)=sp%ysz(:)
    dcp_size_local(3,:)=sp%zsz(:)
    allocate(dcp_sp_sizes(3,3,nproc))
    call MPI_ALLGATHER(dcp_size_local,9,MPI_INT,dcp_sp_sizes,9,MPI_INT,MPI_COMM_WORLD,ierror)

    return
  end subroutine decomp_2d_partition_info

  !------------------------------------------------------------------------------
  !! subroutine: allocate real varialbe on X-pencil with ghost cells
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine alloc_x_ghost_real(var,decomp)
    implicit none

    TYPE(DECOMP_INFO), intent(IN) :: decomp
    real(mytype), allocatable, dimension(:,:,:),intent(out) :: var
    integer :: status, errorcode

    allocate(var(decomp%xsz(1),decomp%xsz(2),0:decomp%xsz(3)+1),STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising variables')
    endif

    return
  end subroutine alloc_x_ghost_real

  !------------------------------------------------------------------------------
  !! subroutine: allocate real varialbe on Y-pencil with ghost cells
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine alloc_y_ghost_real(var,decomp)
    implicit none

    TYPE(DECOMP_INFO), intent(IN) :: decomp
    real(mytype), allocatable, dimension(:,:,:),intent(out) :: var
    integer :: status, errorcode

    allocate(var(decomp%ysz(1),decomp%ysz(2),0:decomp%ysz(3)+1),STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising variables')
    endif

    return
  end subroutine alloc_y_ghost_real

  !------------------------------------------------------------------------------
  !! subroutine: allocate complex varialbe on X-pencil with ghost cells
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine alloc_x_ghost_complex(var,decomp)
    implicit none

    TYPE(DECOMP_INFO), intent(IN) :: decomp
    complex(mytype), allocatable, dimension(:,:,:),intent(out) :: var
    integer :: status, errorcode

    allocate(var(decomp%xsz(1),decomp%xsz(2),0:decomp%xsz(3)+1),STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising variables')
    endif

    return
  end subroutine alloc_x_ghost_complex

  !------------------------------------------------------------------------------
  !! subroutine: allocate complex varialbe on Y-pencil with ghost cells
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine alloc_y_ghost_complex(var,decomp)
    implicit none

    TYPE(DECOMP_INFO), intent(IN) :: decomp
    complex(mytype), allocatable, dimension(:,:,:),intent(out) :: var
    integer :: status, errorcode

    allocate(var(decomp%ysz(1),decomp%ysz(2),0:decomp%ysz(3)+0),STAT=status)
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising variables')
    endif

    return
  end subroutine alloc_y_ghost_complex

  !------------------------------------------------------------------------------
  !! subroutine: update ghost cells
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!  - this function is inpiered form 2decomp
  !!
  !------------------------------------------------------------------------------
  subroutine update_ghost_real(inout,s1,s2,s3,dcp)
    implicit none

    TYPE(DECOMP_INFO),intent(IN) :: dcp
    integer, intent(IN) :: s1, s2, s3 
    real(mytype), dimension(1:s1,1:s2,0:s3+1), intent(INOUT) :: inout   

    ! starting/ending index of array with ghost cells
    integer :: xs, ys, zs, xe, ye, ze
    integer :: i, j, k
    integer :: ierror

    integer :: icount              
    integer, dimension(4) :: requ_4
    integer, dimension(2) :: requ_2
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_t, tag_b

    integer :: data_type
    
    data_type = real_type
    
    if (s1==dcp%xsz(1)) then  ! X-pencil input
       xs = 1 
       xe = s1
       ys = 1
       ye = s2
       zs = 1 - level
       ze = s3 + level
    elseif (s2==dcp%ysz(2)) then  ! Y-pencil input
       xs = 1
       xe = s1
       ys = 1
       ye = s2
       zs = 1 - level
       ze = s3 + level
    else ! invalid input
       call decomp_2d_abort(10, &
            'Invalid data passed to update_ghost')
    endif

    if(dims(2)==1)then
       !nothing to do with ghost cell

    elseif(s1==dcp%xsz(1)) then ! X-pencil
       
       ! *** top/bottom ***
       tag_b = coord(2)
       if (coord(2)==dims(2)-1) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = s1*s2*level
  
       if(coord(2)==0)then
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_2(1), ierror)
          ! send to top
          call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       elseif(coord(2)==dims(2)-1)then
          ! receive from bottom
          call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_2(1), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       else
          ! receive from bottom
          call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_4(1), ierror)
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_4(2), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_4(3), ierror)
          ! send to top
          call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_4(4), ierror)
          call MPI_WAITALL(4, requ_4, status, ierror)
       endif

       
    else if (s2==dcp%ysz(2)) then ! Y-pencil
       ! *** top/bottom ***
       tag_b = coord(2)
       if (coord(2)==dims(2)-1) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = s2*s1*level
       ! receive from bottom
       if(coord(2)==0)then
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_2(1), ierror)
          ! send to top
          call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       elseif(coord(2)==dims(2)-1)then
          ! receive from bottom
          call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_2(1), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       else
          ! receive from bottom
          call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_4(1), ierror)
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_4(2), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_4(3), ierror)
          ! send to top
          call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_4(4), ierror)
          call MPI_WAITALL(4, requ_4, status, ierror)
       endif

    endif
    
    return
  end subroutine update_ghost_real
  
  !------------------------------------------------------------------------------
  !! subroutine: update ghost cells complex
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!  - this function is inpiered form 2decomp
  !!
  !------------------------------------------------------------------------------
  subroutine update_ghost_complex(inout,s1,s2,s3,dcp)
    implicit none

    TYPE(DECOMP_INFO),intent(IN) :: dcp
    integer, intent(IN) :: s1, s2, s3 
    complex(mytype), dimension(1:s1,1:s2,0:s3+1), intent(INOUT) :: inout   

    ! starting/ending index of array with ghost cells
    integer :: xs, ys, zs, xe, ye, ze
    integer :: i, j, k
    integer :: ierror

    integer :: icount              
    integer, dimension(4) :: requ_4
    integer, dimension(2) :: requ_2
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_t, tag_b

    integer :: data_type

    data_type = complex_type
    
    if (s1==dcp%xsz(1)) then  ! X-pencil input
       xs = 1 
       xe = s1
       ys = 1
       ye = s2
       zs = 1 - level
       ze = s3 + level
    elseif (s2==dcp%ysz(2)) then  ! Y-pencil input
       xs = 1
       xe = s1
       ys = 1
       ye = s2
       zs = 1 - level
       ze = s3 + level
    else ! invalid input
       call decomp_2d_abort(10, &
            'Invalid data passed to update_ghost')
    endif

    if(dims(2)==1)then
       !nothing to do with ghost cell
    elseif (s1==dcp%xsz(1)) then ! X-pencil
       
       ! *** top/bottom ***
       tag_b = coord(2)
       if (coord(2)==dims(2)-1) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = s1*s2*level
  
       if(coord(2)==0)then
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_2(1), ierror)
          ! send to top
          call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       elseif(coord(2)==dims(2)-1)then
          ! receive from bottom
          call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_2(1), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       else
          ! receive from bottom
          call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_4(1), ierror)
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_4(2), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_4(3), ierror)
          ! send to top
          call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_4(4), ierror)
          call MPI_WAITALL(4, requ_4, status, ierror)
       endif

       
    else if (s2==dcp%ysz(2)) then ! Y-pencil
       ! *** top/bottom ***
       tag_b = coord(2)
       if (coord(2)==dims(2)-1) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = s2*s1*level
       ! receive from bottom
       if(coord(2)==0)then
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_2(1), ierror)
          ! send to top
          call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       elseif(coord(2)==dims(2)-1)then
          ! receive from bottom
          call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_2(1), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       else
          ! receive from bottom
          call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_4(1), ierror)
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_4(2), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_4(3), ierror)
          ! send to top
          call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_4(4), ierror)
          call MPI_WAITALL(4, requ_4, status, ierror)
       endif
    endif

    return
  end subroutine update_ghost_complex
  
  !------------------------------------------------------------------------------
  !! subroutine: update top ghostcells
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!  - this function is inpiered form 2decomp
  !!
  !------------------------------------------------------------------------------
  subroutine update_ghost_real_top(inout,s1,s2,s3,dcp)
    implicit none

    TYPE(DECOMP_INFO),intent(IN) :: dcp
    integer, intent(IN) :: s1, s2, s3 
    real(mytype), dimension(1:s1,1:s2,1:s3+1), intent(INOUT) :: inout   

    ! starting/ending index of array with ghost cells
    integer :: xs, ys, zs, xe, ye, ze
    integer :: i, j, k
    integer :: ierror

    integer :: icount              
    integer, dimension(2) :: requ_2
    integer, dimension(1) :: requ_1
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_t, tag_b

    integer :: data_type
    
    data_type = real_type
    
    if (s1==dcp%xsz(1)) then  ! X-pencil input
       xs = 1 
       xe = s1
       ys = 1
       ye = s2
       zs = 1
       ze = s3 + 1
    elseif (s2==dcp%ysz(2)) then  ! Y-pencil input
       xs = 1
       xe = s1
       ys = 1
       ye = s2
       zs = 1
       ze = s3 + 1
    else ! invalid input
       call decomp_2d_abort(10, &
            'Invalid data passed to update_ghost')
    endif

    if(dims(2)==1)then
       !nothing to do with ghost cell

    elseif(s1==dcp%xsz(1)) then ! X-pencil
       
       ! *** top/bottom ***
       tag_b = coord(2)
       if (coord(2)==dims(2)-1) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = s1*s2
  
       if(coord(2)==0)then
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_1(1), ierror)
          call MPI_WAITALL(1, requ_1, status, ierror)
       elseif(coord(2)==dims(2)-1)then
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_1(1), ierror)
          call MPI_WAITALL(1, requ_1, status, ierror)
       else
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_2(1), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       endif

       
    else if (s2==dcp%ysz(2)) then ! Y-pencil
       ! *** top/bottom ***
       tag_b = coord(2)
       if (coord(2)==dims(2)-1) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = s2*s1
       ! receive from bottom
       if(coord(2)==0)then
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_1(1), ierror)
          call MPI_WAITALL(1, requ_1, status, ierror)
       elseif(coord(2)==dims(2)-1)then
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_1(1), ierror)
          call MPI_WAITALL(1, requ_1, status, ierror)
       else
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_2(1), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       endif

    endif
    
    return
  end subroutine update_ghost_real_top

  !------------------------------------------------------------------------------
  !! subroutine: update top ghostcells
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!  - this function is inpiered form 2decomp
  !!
  !------------------------------------------------------------------------------
  subroutine update_ghost_complex_top(inout,s1,s2,s3,dcp)
    implicit none

    TYPE(DECOMP_INFO),intent(IN) :: dcp
    integer, intent(IN) :: s1, s2, s3 
    complex(mytype), dimension(1:s1,1:s2,1:s3+1), intent(INOUT) :: inout   

    ! starting/ending index of array with ghost cells
    integer :: xs, ys, zs, xe, ye, ze
    integer :: i, j, k
    integer :: ierror

    integer :: icount              
    integer, dimension(2) :: requ_2
    integer, dimension(1) :: requ_1
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_t, tag_b

    integer :: data_type
    
    data_type = complex_type
    
    if (s1==dcp%xsz(1)) then  ! X-pencil input
       xs = 1 
       xe = s1
       ys = 1
       ye = s2
       zs = 1
       ze = s3 + 1
    elseif (s2==dcp%ysz(2)) then  ! Y-pencil input
       xs = 1
       xe = s1
       ys = 1
       ye = s2
       zs = 1
       ze = s3 + 1
    else ! invalid input
       call decomp_2d_abort(10, &
            'Invalid data passed to update_ghost')
    endif

    if(dims(2)==1)then
       !nothing to do with ghost cell

    elseif(s1==dcp%xsz(1)) then ! X-pencil
       
       ! *** top/bottom ***
       tag_b = coord(2)
       if (coord(2)==dims(2)-1) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = s1*s2
  
       if(coord(2)==0)then
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_1(1), ierror)
          call MPI_WAITALL(1, requ_1, status, ierror)
       elseif(coord(2)==dims(2)-1)then
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_1(1), ierror)
          call MPI_WAITALL(1, requ_1, status, ierror)
       else
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_2(1), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       endif

       
    else if (s2==dcp%ysz(2)) then ! Y-pencil
       ! *** top/bottom ***
       tag_b = coord(2)
       if (coord(2)==dims(2)-1) then
          tag_t = 0
       else
          tag_t = coord(2) + 1
       end if
       icount = s2*s1
       ! receive from bottom
       if(coord(2)==0)then
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_Y, &
               requ_1(1), ierror)
          call MPI_WAITALL(1, requ_1, status, ierror)
       elseif(coord(2)==dims(2)-1)then
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_Y, &
               requ_1(1), ierror)
          call MPI_WAITALL(1, requ_1, status, ierror)
       else
          ! receive from top
          call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
               neighbour(2,5), tag_t, DECOMP_2D_COMM_CART_X, &
               requ_2(1), ierror)
          ! send to bottom
          call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
               neighbour(2,6), tag_b, DECOMP_2D_COMM_CART_X, &
               requ_2(2), ierror)
          call MPI_WAITALL(2, requ_2, status, ierror)
       endif

    endif
    
    return
  end subroutine update_ghost_complex_top

  !------------------------------------------------------------------------------
  !! subroutine: update halo for Lagrangian interpolation
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 08/03/2016
  !!
  !! - description:
  !!  - this function is inpiered form 2decomp
  !!
  !------------------------------------------------------------------------------
  subroutine update_halo3(in,out,level,dcp)

    implicit none
    TYPE(DECOMP_INFO),intent(in) :: dcp
    integer, intent(IN) :: level      ! levels of halo cells required
    real(mytype), dimension(dcp%xsz(1),dcp%xsz(2),dcp%xsz(3)+1), intent(IN) :: in    
    real(mytype), allocatable, dimension(:,:,:), intent(OUT) :: out

    ! starting/ending index of array with halo cells
    integer :: xs, ys, zs, xe, ye, ze

    integer :: i, j, k, s1, s2, s3, ierror
    integer :: data_type

    integer :: icount, ilength, ijump 
    integer :: halo12, halo21, halo31, halo32
    
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

    data_type = real_type

    s1 = dcp%xsz(1)
    s2 = dcp%xsz(2)
    s3 = dcp%xsz(3)

    xs = 1 - level
    xe = s1 + level
    ys = 1 - level
    ye = s2 + level 
    zs = 1 - level
    ze = s3 + level
    
    allocate(out(xs:xe, ys:ye, zs:ze))
    
    do k=1,s3
       do j=1,s2
          do i=1,s1
             out(i,j,k) = in(i,j,k)
          enddo
       enddo
    enddo

    ! *** east/west ***
    do k=1,s3
       do j=1,s2
          do i=1-level,0
             out(i,j,k) = in(s1+i,j,k)
          enddo
       enddo
    enddo
    do k=1,s3
       do j=1,s2
          do i=1,level
             out(s1+i,j,k) = in(i,j,k)
          enddo
       enddo
    enddo

    ! *** north/south *** 
    tag_s = coord(1)
    if (coord(1)==dims(1)-1) then
       tag_n = 0
    else
       tag_n = coord(1) + 1
    end if
    icount = s3 + 2*level
    ilength = level * (s1+2*level)
    ijump = (s1+2*level)*(s2+2*level)
    call MPI_TYPE_VECTOR(icount, ilength, ijump, &
         data_type, halo12, ierror)
    call MPI_TYPE_COMMIT(halo12, ierror)
    ! receive from south
    call MPI_IRECV(out(xs,ys,zs), 1, halo12, &
         neighbour(1,4), tag_s, DECOMP_2D_COMM_CART_X, &
         requests(1), ierror)
    ! receive from north
    call MPI_IRECV(out(xs,ye-level+1,zs), 1, halo12, &
         neighbour(1,3), tag_n, DECOMP_2D_COMM_CART_X, &
         requests(2), ierror)
    ! send to south
    call MPI_ISSEND(out(xs,ys+level,zs), 1, halo12, &
         neighbour(1,4), tag_s, DECOMP_2D_COMM_CART_X, &
         requests(3), ierror)
    ! send to north
    call MPI_ISSEND(out(xs,ye-level-level+1,zs), 1, halo12, &
         neighbour(1,3), tag_n, DECOMP_2D_COMM_CART_X, &
         requests(4), ierror)
    call MPI_WAITALL(4, requests, status, ierror)
    call MPI_TYPE_FREE(halo12, ierror)

    ! *** top/bottom ***
    ! no need to define derived data type as data on xy-planes
    ! all contiguous in memory, which can be sent/received using
    ! MPI directly
    tag_b = coord(2)
    if (coord(2)==dims(2)-1) then
       tag_t = 0
    else
       tag_t = coord(2) + 1
    end if
    icount = ((s1+2*level) * (s2+2*level)) * level
    ! receive from bottom
    call MPI_IRECV(out(xs,ys,zs), icount, data_type, &
         neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
         requests(1), ierror)
    ! receive from top
    call MPI_IRECV(out(xs,ys,ze-level+1), icount, data_type, &
         neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
         requests(2), ierror)
    ! send to bottom
    call MPI_ISSEND(out(xs,ys,zs+level), icount, data_type, &
         neighbour(1,6), tag_b, DECOMP_2D_COMM_CART_X, &
         requests(3), ierror)
    ! send to top
    call MPI_ISSEND(out(xs,ys,ze-level-level+1), icount, data_type, &
         neighbour(1,5), tag_t, DECOMP_2D_COMM_CART_X, &
         requests(4), ierror)
    call MPI_WAITALL(4, requests, status, ierror)

    return
  end subroutine update_halo3

  !------------------------------------------------------------------------------
  !! subroutine: initialise neighbour for ghost cell
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 09/06/2014
  !!
  !! - description:
  !!  - this function is inpiered form 2decomp
  !!
  !------------------------------------------------------------------------------
  subroutine init_neighbour
    implicit none

    integer :: ierror
    
    ! For X-pencil
    neighbour(1,1) = MPI_PROC_NULL               ! east
    neighbour(1,2) = MPI_PROC_NULL               ! west
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
         neighbour(1,4), neighbour(1,3), ierror) ! north & south
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
         neighbour(1,6), neighbour(1,5), ierror) ! top & bottom

    ! For Y-pencil
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 0, 1, &
         neighbour(2,2), neighbour(2,1), ierror) ! east & west
    neighbour(2,3) = MPI_PROC_NULL               ! north
    neighbour(2,4) = MPI_PROC_NULL               ! south
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 1, 1, &
         neighbour(2,6), neighbour(2,5), ierror) ! top & bottom

    return
  end subroutine init_neighbour
  

  !------------------------------------------------------------------------------
  !! subroutine: distribute slice at given z
  !------------------------------------------------------------------------------
  !!
  !! - last modified: fabien margairaz ( fabien.margairaz@gmail.com ), 19/06/2014
  !!
  !! - description:
  !!
  !------------------------------------------------------------------------------

  subroutine distrib_slice_z()
    implicit none

    return
  end subroutine distrib_slice_z

end module decomp_2d_custom
