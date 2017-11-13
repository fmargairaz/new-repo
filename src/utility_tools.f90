!------------------------------------------------------------------------------
!!  module: utility tools
!------------------------------------------------------------------------------
!!
!! - fabien margairaz (fabien.margairaz@gmail.com)
!!
!! - description:
!!
!------------------------------------------------------------------------------

module utility_tools

  use parameters_IO
  use decomp_2d
  use decomp_2d_custom
  use fft_engine

  implicit none
  include 'mpif.h'
  include "fftw3.f"

  private
  public :: compute_planar_average,compute_planar_average_wghost
  public :: get_layer_x,get_layer_fft_x
  public :: get_layer_average,get_layer_filter
  public :: gather_layer_spz

  interface get_layer_average
     module procedure get_layer_average_3Din
     module procedure get_layer_average_2Din
  end interface get_layer_average

  interface get_layer_filter
     module procedure get_layer_filter_3Din
     module procedure get_layer_filter_2Din
  end interface get_layer_filter

contains !=======================================================================
  
  !------------------------------------------------------------------------------
  !! subroutine: initialise utility tools module
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine utility_tools_init()
    implicit none
    
    return
  end subroutine utility_tools_init
  
  !------------------------------------------------------------------------------
  !! subroutine: compute planar average 
  !------------------------------------------------------------------------------
  !! - description:
  !!  - var has to be passed without the ghost cell
  !!
  !------------------------------------------------------------------------------
  subroutine compute_planar_average(plavg,var,dcp,root)
    implicit none

    TYPE(DECOMP_INFO),intent(IN) :: dcp    
    real(rprec),dimension(:),intent(OUT) :: plavg
    real(rprec),dimension(:,:,:),intent(IN) :: var
    integer,intent(IN),optional :: root

    real(rprec),dimension(dcp%xsz(3)) :: var_sum_l
    integer :: i,j,k
    integer :: ierror

    plavg = 0._rprec;var_sum_l = 0._rprec;
    
    !$omp parallel do reduction(+:var_sum_l)
    do k=1,dcp%xsz(3)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             var_sum_l(k)=var_sum_l(k)+var(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    if(present(root))then
       call MPI_reduce(var_sum_l,plavg,dcp%xsz(3),real_type,MPI_SUM,root,DECOMP_2D_LAYER_X,ierror)   
    else
       call MPI_allreduce(var_sum_l,plavg,dcp%xsz(3),real_type,MPI_SUM,DECOMP_2D_LAYER_X,ierror)
    endif
    plavg=plavg/(nx*ny) 
    
    return
  end subroutine compute_planar_average

  !------------------------------------------------------------------------------
  !! subroutine: compute planar average with ghost cell
  !------------------------------------------------------------------------------
  !! - description:
  !!  - var has to be passed with the ghost cell
  !!
  !------------------------------------------------------------------------------
  subroutine compute_planar_average_wghost(plavg,var,dcp,root)
    implicit none

    TYPE(DECOMP_INFO),intent(IN) :: dcp
    real(rprec),dimension(0:dcp%xsz(3)+1),intent(OUT) :: plavg
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(IN) :: var
    integer,intent(IN),optional :: root

    real(rprec),dimension(0:dcp%xsz(3)+1) :: var_sum_l
    integer :: i,j,k
    integer :: ierror

    plavg = 0._rprec;var_sum_l = 0._rprec;
    
    !$omp parallel do reduction(+:var_sum_l)
    do k=0,dcp%xsz(3)+1
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             var_sum_l(k)=var_sum_l(k)+var(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    if(present(root))then
       call MPI_reduce(var_sum_l,plavg,dcp%xsz(3)+2,real_type,MPI_SUM,root,DECOMP_2D_LAYER_X,ierror)   
    else
       call MPI_allreduce(var_sum_l,plavg,dcp%xsz(3)+2,real_type,MPI_SUM,DECOMP_2D_LAYER_X,ierror)
    endif
    plavg=plavg/(nx*ny) 
    
    return
  end subroutine compute_planar_average_wghost

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine get_layer_average_3Din(var_out,var_in,layer,dcp)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(2)+1),intent(IN) :: var_in
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(OUT) :: var_out
    integer,intent(IN) :: layer

    real(rprec) :: var_s,var_sum,const
    integer :: i,j,ierror

    var_s=0._rprec;var_sum=0._rprec
       
    do j=1,dcp%xsz(2)
       do i=1,dcp%xsz(1)             
          var_s=var_s+var_in(i,j,layer)
       enddo
    enddo
    
    call MPI_Allreduce(var_s,var_sum,1,real_type,MPI_SUM,DECOMP_2D_LAYER_X,ierror)
    const=1._rprec/(nx*ny)
    var_sum=const*var_sum
    
    do j=1,dcp%xsz(2)
       do i=1,dcp%xsz(1)
          var_out(i,j)=var_sum
       enddo
    enddo
    
    return
  end subroutine get_layer_average_3Din

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine get_layer_average_2Din(var_out,var_in,dcp)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(IN) :: var_in
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(OUT) :: var_out

    real(rprec) :: var_s,var_sum,const
    integer :: i,j,ierror

    var_s=0._rprec;var_sum=0._rprec
       
    do j=1,dcp%xsz(2)
       do i=1,dcp%xsz(1)             
          var_s=var_s+var_in(i,j)
       enddo
    enddo
    
    call MPI_Allreduce(var_s,var_sum,1,real_type,MPI_SUM,DECOMP_2D_LAYER_X,ierror)
    const=1._rprec/(nx*ny)
    var_sum=const*var_sum
    
    do j=1,dcp%xsz(2)
       do i=1,dcp%xsz(1)
          var_out(i,j)=var_sum
       enddo
    enddo
    
    return
  end subroutine get_layer_average_2Din

  !------------------------------------------------------------------------------
  !! subroutine: extract slice at given z
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------
  subroutine get_layer_x(layer,layer_loc,dcp,root)
    implicit none
    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(IN) :: layer_loc
    real(rprec),dimension(nx,ny),intent(OUT) ::layer
    integer,intent(IN),optional :: root

    integer,dimension(:),allocatable :: rdispl,rcount
    integer :: k,m,n
    integer :: icount,ierror

    if(dims(1)>1)then
       ! gather on root node
       icount = dcp%xsz(1)*dcp%xsz(2) 
       allocate(rdispl(dims(1)));allocate(rcount(dims(1)))
       do m=1,dims(1)
          n=ranks_2d_layer_x(1,m)+1
          rcount(m) = dcp_ph_sizes(1,1,n)*dcp_ph_sizes(1,2,n); 
       enddo
       rdispl(1)=0
       do k=2,dims(1)
          rdispl(k) = rdispl(k-1)+rcount(k-1); 
       enddo
       if(present(root))then
          call MPI_GATHERV(layer_loc,icount,real_type,&
               layer,rcount,rdispl,MPI_DOUBLE_PRECISION,root,DECOMP_2D_LAYER_X,ierror)
       else
          call MPI_GATHERV(layer_loc,icount,real_type,&
               layer,rcount,rdispl,MPI_DOUBLE_PRECISION,0,DECOMP_2D_LAYER_X,ierror)
          icount = nx*ny
          call MPI_BCAST(layer,icount,MPI_DOUBLE_PRECISION,0,DECOMP_2D_LAYER_X,ierror)
       endif
       
       deallocate(rdispl);deallocate(rcount)
    else
       layer=layer_loc
    endif
    return
  end subroutine get_layer_x

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !! - description:                                                                                              
  !!
  !------------------------------------------------------------------------------
  subroutine get_layer_filter_3Din(var_out,var_in,layer,G,dcp)
    use testfilter
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2),0:dcp%xsz(3)+1),intent(IN) :: var_in
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(INOUT) :: var_out
    integer,intent(IN) :: layer
    real(rprec),dimension(lhx,ny),intent(IN) :: G 

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)) :: var_loc
    real(rprec),dimension(nx,ny) :: var_glob
    integer :: i,j,ierror

    var_glob=0.0_rprec
    
    if(dims(1)==1)then ! if dims(1)=1  then (dcp%xsz(1),dcp%xsz(2))=(nx,ny)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             var_glob(i,j)=var_in(i,j,layer)
          enddo
       enddo
    else
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             var_loc(i,j)=var_in(i,j,layer)
          enddo
       enddo
       !get var over the whole domaine at each node
       call get_layer_x(var_glob,var_loc,dcp)
    endif
    
    ! each node has the data for the whole domaine
    call test_filter_layer(var_glob,G)

    if(dims(1)==1)then
       var_out=var_glob
    else
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             var_out(i,j)=var_glob(dcp%xst(1)-1+i,dcp%xst(2)-1+j)
          enddo
       enddo
    endif

    return
  end subroutine get_layer_filter_3Din

  !------------------------------------------------------------------------------
  !! subroutine:
  !------------------------------------------------------------------------------
  !! - description:                                                                                              
  !!
  !------------------------------------------------------------------------------
  subroutine get_layer_filter_2Din(var_out,var_in,G,dcp)
    use testfilter
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(IN) :: var_in
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(INOUT) :: var_out
    real(rprec),dimension(lhx,ny),intent(IN) :: G 

    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)) :: var_loc
    real(rprec),dimension(nx,ny) :: var_glob
    integer :: i,j,ierror

    var_glob=0.0_rprec

    if(dims(1)==1)then ! if dims(1)=1  then (dcp%xsz(1),dcp%xsz(2))=(nx,ny)
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             var_glob(i,j)=var_in(i,j)
          enddo
       enddo
    else
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             var_loc(i,j)=var_in(i,j)
          enddo
       enddo
       !get var over the whole domaine at each node
       call get_layer_x(var_glob,var_loc,dcp)
    endif
    
    ! each node has the data for the whole domaine
    call test_filter_layer(var_glob,G)

    if(dims(1)==1)then
       var_out=var_glob
    else
       do j=1,dcp%xsz(2)
          do i=1,dcp%xsz(1)
             var_out(i,j)=var_glob(dcp%xst(1)-1+i,dcp%xst(2)-1+j)
          enddo
       enddo
    endif

    return
  end subroutine get_layer_filter_2Din

  !------------------------------------------------------------------------------
  !! subroutine: extract slice at given z
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------

  subroutine get_layer_fft_x(layer,layer_loc,dcp)
    implicit none
    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(IN) :: layer_loc
    complex(rprec),dimension(lhx,ny),intent(OUT) ::layer
    real(rprec),dimension(nx,ny) :: scr

    call get_layer_x(scr,layer_loc,dcp)
    call dfftw_execute_dft_r2c(plan2d_f,scr,layer)
    
    return
  end subroutine get_layer_fft_x

  !------------------------------------------------------------------------------
  !! subroutine: extract slice at given z
  !------------------------------------------------------------------------------
  !! - description:
  !!
  !------------------------------------------------------------------------------

  subroutine get_layer_y(layer,layer_loc,dcp,root)
    implicit none
    TYPE(DECOMP_INFO),intent(in) :: dcp
    real(rprec),dimension(dcp%xsz(1),dcp%xsz(2)),intent(IN) :: layer_loc
    real(rprec),dimension(nx,ny),intent(OUT) ::layer
    integer,intent(IN),optional :: root
    
    layer = (0._rprec,0._rprec)

!!$    complex(rprec),dimension(nx/2+1,ny) :: bot_layer,top_layer
!!$    complex(rprec),dimension(sp%ysz(1),sp%ysz(2)) :: bot_layer_loc,top_layer_loc
!!$    integer,dimension(:),allocatable :: rdispl,rcount 
!!$    integer :: icount,root,ierror
!!$    
!!$    if(sp%yst(3)==1)then
!!$       if(dims(1)==1)then
!!$          do j=1,sp%ysz(2)
!!$             do i=1,sp%ysz(1)
!!$                bot_layer(i,j)=wk_press(i,j,0)
!!$             enddo
!!$          enddo
!!$       else
!!$          do j=1,sp%ysz(2)
!!$             do i=1,sp%ysz(1)
!!$                bot_layer_loc(i,j)=wk_press(i,j,0)
!!$             enddo
!!$          enddo
!!$          ! gather on root node
!!$          icount = sp%ysz(1)*sp%ysz(2)
!!$          allocate(rdispl(dims(1)));allocate(rcount(dims(1)))
!!$          do m=1,dims(1)
!!$             n=ranks_2d_layer_x(1,m)+1
!!$             rcount(m) = dcp_sp_sizes(2,1,n)*dcp_sp_sizes(2,2,n); 
!!$          enddo
!!$          rdispl(1)=0
!!$          do k=2,dims(1)
!!$             rdispl(k) = rdispl(k-1)+rcount(k-1); 
!!$          enddo
!!$          call MPI_GATHERV(bot_layer_loc,icount,complex_type,&
!!$               bot_layer,rcount,rdispl,complex_type,0,DECOMP_2D_LAYER_X,ierror)
!!$          deallocate(rdispl);deallocate(rcount)
!!$       endif
!!$    endif
!!$    ! -- Z-pencil ---------------------------------
!!$    icount = (nx/2+1)*ny
!!$    root = 0 
!!$    call MPI_BCAST(bot_layer,icount,complex_type,root,MPI_COMM_WORLD,ierror)
!!$    do j=1,sp%zsz(2)
!!$       do i=1,sp%zsz(1)
!!$          RHS_col_z(i,j,1) = bot_layer(sp%zst(1)-1+i,sp%zst(2)-1+j)
!!$       enddo
!!$    enddo
    
    return
  end subroutine get_layer_y

  subroutine gather_layer_spz(olayer,ilayer,sp)
    implicit none
    TYPE(DECOMP_INFO),intent(in) :: sp
    complex(rprec),dimension(sp%zsz(1),sp%zsz(2)),intent(IN) :: ilayer
    complex(rprec),dimension(lhx,ny),intent(OUT) :: olayer
    complex(rprec),dimension(lhx,ny) :: wk_layer
    integer :: i,j,ierror,newtype1,newtype2
    !integer, dimension(2) :: sizes
    !integer, dimension(2) :: subsizes
    !integer, dimension(2) :: starts
    !integer, dimension(nproc) :: scount,rcount,displ

    olayer = (0._rprec,0._rprec)
    wk_layer = (0._rprec,0._rprec)

    do j=1,sp%zsz(2)
       do i=1,sp%zsz(1)
          wk_layer(sp%zst(1)-1+i,sp%zst(2)-1+j)=ilayer(i,j)
       enddo
    enddo
    call MPI_allreduce(wk_layer,olayer,lhx*ny,complex_type,MPI_SUM,MPI_COMM_WORLD,ierror)

!!$    write(*,*)
!!$    sizes(1)=lhx;sizes(2)=ny; !size of global array
!!$    subsizes(1)=sp%zsz(1);subsizes(2)=sp%zsz(2); !size of sub-region
!!$    starts(1)=sp%zst(1)-1;starts(2)=sp%zst(2)-1; !starts of each array
!!$    !starts(1)=1;starts(2)=1;
!!$
!!$    write(*,*) nrank,sizes,subsizes,starts
!!$
!!$    scount=subsizes(1)*subsizes(2)
!!$    rcount=sizes(1)*sizes(2)
!!$    !displ=
!!$
!!$    call MPI_Type_create_subarray(2,sizes,subsizes,starts,&
!!$         MPI_ORDER_FORTRAN,MPI_DOUBLE_COMPLEX,newtype2,ierror);
!!$     call MPI_Type_create_resized(newtype2, 0, 2*sizeof(rprec),newtype1,ierror);
!!$     call MPI_Type_commit(newtype1);
!!$
!!$    olayer = (0._rprec,0._rprec)
!!$    !call MPI_ALLGATHER()
!!$
!!$    call MPI_ALLGATHERV(ilayer,scount,newtype1,&
!!$         olayer,rcount,MPI_DOUBLE_COMPLEX,&
!!$         DECOMP_2D_COMM_CART_Z,ierror)
!!$
!!$    call MPI_TYPE_FREE(newtype1,ierror)

    return
  end subroutine gather_layer_spz
  
end module utility_tools

