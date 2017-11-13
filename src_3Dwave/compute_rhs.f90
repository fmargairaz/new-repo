module rhs_mod

  use derivatives_mod
  use fft_mod
  use decomp_2d

  implicit none
  include 'mpif.h'

  private
  
  public :: add2rhs_ddx, add2rhs_ddy, add2rhs_ddz
  public :: add2rhs_ddy_notrans
  
contains
  
  subroutine add2rhs_ddx(u,rhs,vector)
    implicit none
    
    TYPE(DECOMP_INFO),intent(in) :: vector
    !use local indexes
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: u
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: rhs

    real(rprec), allocatable, dimension(:,:,:) :: dudx
   
    integer :: j,k     

    real(rprec) :: test_time
    
    call alloc_x(dudx,vector)
    
    
    test_time = mpi_wtime()
    !compute FT in x to compute the derivative
    !$OMP PARALLEL DO 
    do k=1,vector%xsz(3)
       do j=1,vector%xsz(2)
          call ddx(dudx(:,j,k),u(:,j,k))
       enddo
    enddo
    !$OMP END PARALLEL DO
    test_time = mpi_wtime() - test_time
    
    !write(*,*) 'FFT on x',test_time
    
    rhs = rhs + dudx
    
    deallocate(dudx)

  end subroutine add2rhs_ddx

  subroutine add2rhs_ddy(u,rhs,vector)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: vector
    !use local indexes
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: u
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: rhs

    real(rprec), allocatable, dimension(:,:,:) :: u_4y
    real(rprec), allocatable, dimension(:,:,:) :: dudy_4x, dudy_4y    

    integer :: i,k 

    real(rprec) :: test_time

    call alloc_y(u_4y,vector)
    call alloc_x(dudy_4x,vector)
    call alloc_y(dudy_4y,vector)

    call transpose_x_to_y(u,u_4y)

    test_time = mpi_wtime()
    !compute FT in y to compute the derivative
    !$OMP PARALLEL DO 
    do k=1,vector%ysz(3)
       do i=1,vector%ysz(1)
          call ddy(dudy_4y(i,:,k),u_4y(i,:,k))
       enddo
    enddo
    !$OMP END PARALLEL DO
    test_time = mpi_wtime() - test_time

    !write(*,*) 'FFT on y',test_time

    call transpose_y_to_x(u_4y,u)
    call transpose_y_to_x(dudy_4y,dudy_4x)
    
    rhs = rhs + dudy_4x

    deallocate(u_4y)
    deallocate(dudy_4x)
    deallocate(dudy_4y)

  end subroutine add2rhs_ddy

  subroutine add2rhs_ddy_notrans(u,rhs,vector)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: vector
    !use local indexes
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: u
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: rhs

    real(rprec), allocatable, dimension(:,:,:) :: dudy

    integer :: i,k 

    real(rprec) :: test_time

    call alloc_x(dudy,vector)

    test_time = mpi_wtime()
    !compute FT in y to compute the derivative
    !$OMP PARALLEL DO 
    do k=1,vector%ysz(3)
       do i=1,vector%ysz(1)
          call ddy(dudy(i,:,k),u(i,:,k))
       enddo
    enddo
    !$OMP END PARALLEL DO
    test_time = mpi_wtime() - test_time

    !write(*,*) 'FFT on y',test_time

    rhs = rhs + dudy

    deallocate(dudy)
    
  end subroutine add2rhs_ddy_notrans


  subroutine add2rhs_ddz(u,rhs,vector)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: vector
    !use local indexes
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: u
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: rhs
    
    real(rprec), allocatable, dimension(:,:,:) :: dudz
    real(rprec), allocatable, dimension(:,:,:) :: dudz_halo

    integer :: k 

    real(rprec) :: test_time

    call alloc_x(dudz,vector)

    call ddz_uv_4x(dudz,u,vector)
    call update_halo(dudz,dudz_halo,level=1)

    test_time = mpi_wtime()
    !$OMP PARALLEL DO 
    do k=1,vector%xsz(3)
       rhs(:,:,k)=rhs(:,:,k)+0.5_rprec*(dudz_halo(1:vector%xsz(1),1:vector%xsz(2),k-1)+ &
            dudz_halo(1:vector%xsz(1),1:vector%xsz(2),k)) 
    end do
    !$OMP END PARALLEL DO
    test_time = mpi_wtime() - test_time

    !write(*,*) 'DF on z',test_time

    deallocate(dudz)
  end subroutine add2rhs_ddz
    
  subroutine add2rhs_ddz_nohalo(u,rhs,vector)
    implicit none

    TYPE(DECOMP_INFO),intent(in) :: vector
    !use local indexes
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: u
    real(rprec),dimension(vector%xsz(1),vector%xsz(2),vector%xsz(3)),intent(inout) :: rhs
    
    real(rprec), allocatable, dimension(:,:,:) :: dudz
    
    integer :: k 

    call alloc_x(dudz,vector)

    call ddz_uv_4x(dudz,u,vector)
    
    !$OMP PARALLEL DO
    do k=2,vector%xsz(3)
       rhs(:,:,k)=rhs(:,:,k)+0.5_rprec*(dudz(1:vector%xsz(1),1:vector%xsz(2),k-1)+ &
            dudz(1:vector%xsz(1),1:vector%xsz(2),k)) 
    end do
    !$OMP END PARALLEL DO

    deallocate(dudz)

  end subroutine add2rhs_ddz_nohalo
    
end module rhs_mod
