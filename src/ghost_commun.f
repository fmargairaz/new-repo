  ! This file contain common code to be included by subroutines                                                                                                                                            
  ! 'update_ghost_...' in decomp_2d_cunstom.f90 

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
  
  if (s1==dcp%xsz(1)) then ! X-pencil
     
     ! *** top/bottom ***
     tag_b = coord(2)
     if (coord(2)==dims(2)-1) then
        tag_t = 0
     else
        tag_t = coord(2) + 1
     end if
     icount = s1*s2*level
     
     ! receive from bottom
     call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
          neighbour(1,2), tag_b, DECOMP_2D_COMM_CART_X, &
          requests(1), ierror)
     ! receive from top
     call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
          neighbour(1,1), tag_t, DECOMP_2D_COMM_CART_X, &
          requests(2), ierror)
     ! send to bottom
     call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
          neighbour(1,2), tag_b, DECOMP_2D_COMM_CART_X, &
          requests(3), ierror)
     ! send to top
     call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
          neighbour(1,1), tag_t, DECOMP_2D_COMM_CART_X, &
          requests(4), ierror)
     call MPI_WAITALL(4, requests, status, ierror)
     
     
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
     call MPI_IRECV(inout(xs,ys,zs), icount, data_type, &
          neighbour(2,2), tag_b, DECOMP_2D_COMM_CART_Y, &
          requests(1), ierror)
     ! receive from top
     call MPI_IRECV(inout(xs,ys,ze-level+1), icount, data_type, &
          neighbour(2,1), tag_t, DECOMP_2D_COMM_CART_Y, &
          requests(2), ierror)
     ! send to bottom
     call MPI_ISSEND(inout(xs,ys,zs+level), icount, data_type, &
          neighbour(2,2), tag_b, DECOMP_2D_COMM_CART_Y, &
          requests(3), ierror)
     ! send to top
     call MPI_ISSEND(inout(xs,ys,ze-level-level+1), icount, data_type, &
          neighbour(2,1), tag_t, DECOMP_2D_COMM_CART_Y, &
          requests(4), ierror)
     call MPI_WAITALL(4, requests, status, ierror)
  endif
  
  
