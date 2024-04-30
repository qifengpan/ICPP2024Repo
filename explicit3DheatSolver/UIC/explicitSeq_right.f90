program heat

   implicit none

   ! Sizes for the discretization
   integer size_x, size_y, size_z
   integer nsample

   ! Arrays
   double precision, allocatable :: x(:,:,:,:), x0(:,:,:,:)

   ! Space and time steps
   double precision dt, dt1, dt2, hx, hy, hz, epsilon

   ! Current error and limit convergence
   double precision, dimension(128) :: result

   ! Index variables
   integer i, j, k

   ! Time variable and step
   double precision t
   integer step

   ! Max step
   integer maxStep

   ! Variables for clock
   integer count_0, count_1
   integer count_rate, count_max
   double precision time_init, time_final, elapsed_time

   ! Number of initial borders layers to
   ! avoid artifacts problem on corners
   integer nb_layers

   ! Physical parameters
   double precision temp1_init, temp2_init
   double precision ,dimension(128) :: k0
   ! Initial number of borders layers
   nb_layers = 1

   ! temp1_init: temperature init on borders
   temp1_init = 10.0

   ! temp2_init: temperature init inside
   temp2_init = -10.0

   ! Diffusivity coefficient
   k0 = 1
   nSample=size(k0)
   ! Get input parameters
   ! print *,' Size x of the square '
   ! read(*,*) size_x
   ! print *,' Size y of the square '
   ! read(*,*) size_y
   ! print *,' Size z of the square '
   ! read(*,*) size_z
   ! print *, 'Max. number of steps '
   ! read(*,*) maxStep
   ! print *, 'Time step'
   ! read(*,*) dt1
   ! print *,'Convergence'
   ! read(*,*) epsilon
   size_x = 50
   size_y = 50
   size_z = 5
   maxStep= 100000
   dt1    = 0.001
   epsilon= 0.00001

   ! Compute space and time steps
   hx = 1.0d0/dble(size_x+2)
   hy = 1.0d0/dble(size_y+2)
   hz = 1.0d0/dble(size_z+2)
   dt2 = minval(0.125*(min(min(hx,hy),hz)**2)/k0)

   ! Take a right time step for convergence
   if(dt1.ge.dt2) then
    write(*,*)
    write(*,*) 'Time step too large in param file - ', &
               'Taking convergence criterion'
    dt=dt2
   else
    dt=dt1
   endif

   ! Allocation of 3D arrays
   allocate(x(0:size_x+1,0:size_y+1,0:size_z+1,1:nSample))
   allocate(x0(0:size_x+1,0:size_y+1,0:size_z+1,1:nSample))

   ! Initialize values
   call initValues(nb_layers, x0, size_x, size_y, size_z, temp1_init, temp2_init,nSample)

   ! Initialize step and time
   step = 0
   t = 0.0

   ! Starting time
   call system_clock(count_0, count_rate, count_max)
   time_init = count_0*1.0/count_rate

   ! Main loop
10 continue

   ! Increment step and time
   step = step + 1
   t = t + dt

   ! Perform one step of the explicit scheme
   call computeNext(x0, x, size_x, size_y, size_z, dt, hx, hy, hz, result, k0,nSample)

   ! Current error
   result = sqrt(result)

   ! Break conditions of main loop
   if(.not.((maxval(result).lt.epsilon).or.(step.gt.maxStep))) goto 10

   ! Ending time
   call system_clock(count_1, count_rate, count_max)
   time_final = count_1*1.0/count_rate
   ! Elapsed time
   elapsed_time = time_final - time_init

   ! Print results
   write(*,*)
   write(*,*) '  Time step = ',dt
   write(*,*)
   write(*,1001) epsilon,step
   write(*,*)
   write(*,*) '  Problem size = ',size_x*size_y*size_z
   write(*,*)
   write(*,1002) elapsed_time
   write(*,*)
   write(*,*) '  Computed solution in outputSeq.dat '
   write(*,*)

!   ! Store solution into output file
!   open(5,file='outputSeq.dat',action='write',status='replace')
!   do 40, k=0,size_z+1
!    do 30, j=0,size_y+1
!     write(5,1000,advance='no') (x0(:,i,j,k),i=0,size_x)
!     write(5,999,advance='no') x0(:,size_x+1,j,k)
!     write(5,*)
!30  continue
!    if(k.ne.(size_z+1)) then
!     write(5,*)
!    endif
!40 continue
!   close(5)

   ! Free arrays
   deallocate(x)
   deallocate(x0)

   ! Formats available to display the computed values on the grid
999 format(1000(f15.11))
1000 format(1000(f15.11,1x))
1001 format('   Convergence = ',f11.9,' after ',i9,' steps ')
1002 format('   Wall Clock = ',f15.6)

   stop
end

