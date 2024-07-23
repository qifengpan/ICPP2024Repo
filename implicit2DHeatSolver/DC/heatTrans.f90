!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> ADI Method for solving unsteady heat equation in 2D
!     du/dt = a*(d2u/dx2 + d2u/dy2) + q(x,y) 
!     Approximate factorization + Crank-Nicolson
!     Drichlet b.c.
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) 
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 13, 2015
!-----------------------------------------------------------------------------!

program heat2d
implicit none
integer::i,j,k,nx,ny,nt,ns,nf
real*8 :: acu_time,start,finish
real*8 ::dx,dy,dt,x0,xL,y0,yL,t,Tmax,alpha
character(len = 32) :: arg
real*8,allocatable ::u(:,:),x(:),y(:)
! get argv
alpha = 1
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
Tmax = 1.0d0

!diffusion constant:
!read(*,*) alpha
!alpha = 1.0d0

!time step
dt = 0.05d0

!number of points in time
nt = nint(Tmax/dt)

!number of snapshots to plot
ns = nt

!frequency for plotting
nf = nint(dfloat(nt)/dfloat(ns))

!u: convected variable 
allocate(u(0:nx,0:ny))

!initial condition
t = 0.0d0
do j=0,ny
do i=0,nx
u(i,j) = 0.0d0
end do
end do

!Plot initial condition
!open(18,file='u.plt')
!write(18,*) 'variables ="x","y","u"'
!write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!do j=0,ny
!do i=0,nx
!write(18,*) x(i),y(j),u(i,j)
!end do
!end do
acu_time = 0
!time integration
do k=1,nt
    call CPU_time(start)
    call CN(nx,ny,dx,dy,dt,t,alpha,x,y,u)
    call CPU_time(finish)
    acu_time = acu_time + (finish-start)  
    !update t
    t = t+dt 

    !plot field
!    if (mod(k,nf).eq.0) then
!	write(18,100)'zone f=point i=',nx+1,',j=',ny+1,',t="time',t,'"'
!	do j=0,ny
!	do i=0,nx
!	write(18,*) x(i),y(j),u(i,j)
!	end do
!	end do
!    end if

    
	!print*,k,t,maxval(u)


end do
!print*,u(:,10)  
!close(18)
print *,acu_time
open(10,file = 'output.dat')
write(10,*)u(:,10)
close(10)
!100 format(a16,i8,a4,i8,a10,f10.4,a3)
end

!-----------------------------------------------------------------------------!
!source term (nonhomogenous forcing term in heat equation)
!-----------------------------------------------------------------------------!
real*8 function f(t,x,y)
implicit none
real*8::t,x,y
f = 2.0d0*(2.0d0-x*x-y*y)
end
