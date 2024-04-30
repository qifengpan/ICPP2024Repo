
!*******************************************************************!
!             Update Bounds of subdomain with me process            !
!*******************************************************************!

subroutine updateBound(x, size_tot_x, size_tot_y, size_tot_z, &
                       neighBor, comm3d, matrix_type_oxz, matrix_type_oxy, matrix_type_oyz, &
                       me, xs, ys, zs, xe, ye, ze, nproc)

   include 'mpif.h'

   ! Input parameters
   integer size_tot_x, size_tot_y, size_tot_z, nproc, me
   double precision, dimension(0:size_tot_x-1,0:size_tot_y-1,0:size_tot_z-1) :: x
   !Type Matrix type
   integer matrix_type_oxz, matrix_type_oxy, matrix_type_oyz
   integer, parameter :: nbvi = 4
   integer, parameter :: S=1, E=2, N=3, W=4, Zd=5, Zu=6
   integer, dimension(6) :: neighBor
   integer infompi, comm3d
   integer flag
   integer, dimension(mpi_status_size) :: status
   integer xs(0:nproc-1), ys(0:nproc-1), zs(0:nproc-1), xe(0:nproc-1), ye(0:nproc-1), ze(0:nproc-1)

!********* North/South communication *********************************
   flag = 1
   !Send my boundary to North and receive from South
   call MPI_Sendrecv(x(xs(me), ys(me), zs(me)), 1, matrix_type_oxz, neighBor(N), flag, &
      x(xs(me), ye(me)+1, zs(me)), 1, matrix_type_oxz, neighBor(S), flag, comm3d, status, infompi)

   !Send my boundary to South and receive from North
   call MPI_Sendrecv(x(xs(me), ye(me), zs(me)), 1, matrix_type_oxz, neighBor(S), flag, &
      x(xs(me), ys(me)-1, zs(me)), 1, matrix_type_oxz, neighBor(N), flag, comm3d, status, infompi)

!********* Est/West communication ************************************
   flag = 2
   !Send my boundary to Est and receive from West
   call MPI_Sendrecv(x(xe(me), ys(me), zs(me)), 1, matrix_type_oyz, neighBor(E), flag, &
      x(xs(me)-1, ys(me), zs(me)), 1, matrix_type_oyz, neighBor(W), flag, comm3d, status, infompi)

   !Send my boundary to West and receive from Est
   call MPI_Sendrecv(x(xs(me), ys(me), zs(me)), 1, matrix_type_oyz, neighBor(W), flag, &
      x(xe(me)+1, ys(me), zs(me)), 1, matrix_type_oyz, neighBor(E), flag, comm3d, status, infompi)

!********* Zdown/Zup communication ***********************************
   flag = 3
   !Send my boundary to Zup and receive from Zdown
   call MPI_Sendrecv(x(xs(me), ys(me), ze(me)), 1, matrix_type_oxy, neighBor(Zu), flag, &
      x(xs(me), ys(me), zs(me)-1), 1, matrix_type_oxy, neighBor(Zd), flag, comm3d, status, infompi)

   !Send my boundary to Zdown and receive from Zup
   call MPI_Sendrecv(x(xs(me), ys(me), zs(me)), 1, matrix_type_oxy, neighBor(Zd), flag, &
      x(xs(me), ys(me), ze(me)+1), 1, matrix_type_oxy, neighBor(Zu), flag, comm3d, status, infompi)
end
