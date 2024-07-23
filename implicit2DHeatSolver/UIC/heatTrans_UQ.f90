program heat2d
implicit none
integer::i,j,k,nx,ny,nt,ns,nf,dim_a
character(len=32) ::arg
real*8 ::dx,dy,dt,x0,xL,y0,yL,t,Tmax
real*8 ::alpha(64)
real*8,allocatable ::u(:,:,:),x(:),y(:)
real*8 :: start, finish
! call getarg(1,arg)
! read(arg,*)alpha

!Domain
x0 =-1.0d0 !left
xL = 1.0d0 !right

y0 =-1.0d0 !bottom
yL = 1.0d0 !up

!number of points
nx = 1000
ny = 1000

!grid spacing (spatial)
dx = (xL-x0)/dfloat(nx)
dy = (yL-y0)/dfloat(ny)

!spatial coordinates 
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

allocate(y(0:ny))
do j=0,ny
y(j) = y0 + dfloat(j)*dy
end do


!maximum time desired
Tmax = 0.1d0
dim_a = size(alpha)
!dim_a = 18000
!call generateGauss(alpha,dim_a)
alpha = 1.0d0
!diffusion constant:
!time step
dt = 0.05d0

!number of points in time
nt = nint(Tmax/dt)

!number of snapshots to plot
ns = nt

!frequency for plotting
nf = nint(dfloat(nt)/dfloat(ns))

!u: convected variable 
allocate(u(1:dim_a,0:nx,0:ny))

!initial condition
t = 0.0d0
do j=0,ny
do i=0,nx
u(:,i,j) = 0.0d0
end do
end do

!Plot initial condition
! open(18,file='u.plt')
! write(18,*) 'variables ="x","y","u"'
! write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
! do j=0,ny
! do i=0,nx
! write(18,*) x(i),y(j),u(i,j)
! end do
! end do

!time integration
do k=1,nt
    call CN(nx,ny,dx,dy,dt,t,alpha,x,y,u,dim_a)
      
    !update t
    t = t+dt 

    !plot field
 !    if (mod(k,nf).eq.0) then
	! write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
	! do j=0,ny
	! do i=0,nx
	! write(18,*) x(i),y(j),u(i,j)
	! end do
	! end do
 !    end if

   


end do
call cpu_time(start)
open(10,file = 'output_UQ.dat')
write(10,*)u(:,:,10)
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
!do i =0,nx
!do j =1,100
!print *,u(j,i,10)
!end do
!end do
!close(18)


! 100 format(a16,i8,a4,i8,a10,f10.4,a3)
end



!-----------------------------------------------------------------------------!
!source term (nonhomogenous forcing term in heat equation)
!-----------------------------------------------------------------------------!
real*8 function f(t,x,y)
implicit none
real*8::t,x,y
f = 2.0d0*(2.0d0-x*x-y*y)
end

subroutine generateGauss(k0,nSample)
   implicit none
   integer nSample
   double precision k0(nSample)
   double precision :: x
   double precision, parameter :: pi=3.14159265
   double precision :: u1,u2
   integer :: i
   do i = 1 ,nSample
   call random_stduniform(u1)
   call random_stduniform(u2)
   x = sqrt(-2*log(u1))*cos(2*pi*u2)+ 1
   k0(i) = abs(x)
   end do
end subroutine generateGauss
subroutine random_stduniform(u)
   implicit none
   double precision,intent(out) :: u
   double precision :: r
   call random_number(r)
   u = 1 - r
end subroutine random_stduniform

