!-----------------------------------------------------------------------------!
!Crank-Nicolson(CN) scheme with ADI
!-----------------------------------------------------------------------------!
subroutine CN(nx,ny,dx,dy,dt,t,alpha,x,y,u,dim_a)
implicit none
integer::nx,ny,i,j,dim_a
real*8 ::dx,dy,dt,t,f
real*8,dimension(dim_a) :: alpha,bx,by,ua,ub
!real*8 ::u(1:dim_a,0:nx,0:ny),x(0:nx),y(0:ny)
real*8 ::x(0:nx),y(0:ny)
real*8 ,intent(inout) :: u(1:dim_a,0:nx,0:ny)
!real*8,dimension(1:dim_a,1:ny-1) ::a,b,c,r,q
!real*8 ::s(1:dim_a,0:nx,0:ny),p(1:dim_a,0:nx,0:ny),z(1:dim_a,0:nx,0:ny)
real*8,allocatable,dimension(:,:,:) :: s,p,z
real*8,allocatable,dimension(:,:) :: a,b,c,r,q
if (.not. allocated(s)) then
allocate(s(1:dim_a,0:nx,0:ny),p(1:dim_a,0:nx,0:ny),z(1:dim_a,0:nx,0:ny))
endif
do j=0,ny
do i=0,nx
s(:,i,j) = 0.0d0
p(:,i,j) = 0.0d0
z(:,i,j) = 0.0d0
end do
end do

bx = 0.5d0*alpha*dt/(dx*dx)
by = 0.5d0*alpha*dt/(dy*dy)
!compute source term in known step
do j=1,ny-1
do i=0,nx
s(:,i,j) = u(:,i,j) + by*(u(:,i,j+1)-2.0d0*u(:,i,j)+u(:,i,j-1))
end do
end do

do j=1,ny-1
do i=1,nx-1
p(:,i,j) = s(:,i,j) + bx*(s(:,i+1,j)-2.0d0*s(:,i,j)+s(:,i-1,j)) &
       + 0.5d0*dt*(f(t,x(i),y(j))+f(t+dt,x(i),y(j))) 
end do
end do

if (.not. allocated(a)) then
allocate(a(1:dim_a,1:nx-1),b(1:dim_a,1:nx-1),c(1:dim_a,1:nx-1),r(1:dim_a,1:nx-1),q(1:dim_a,1:nx-1))
endif

!x-sweep to compute intermediate values:
do j=1,ny-1

  
	!Build coefficient matrix:
!	allocate(a(1:dim_a,1:nx-1),b(1:dim_a,1:nx-1),c(1:dim_a,1:nx-1),r(1:dim_a,1:nx-1),q(1:dim_a,1:nx-1))

	do i=1,nx-1
	a(:,i) = -bx
	b(:,i) = (1.0d0+2.0d0*bx)
	c(:,i) = -bx
	r(:,i) = p(:,i,j) 
	end do
    
	!apply boundary conditions
    ua = u(:,0,j)  - by*(u(:,0,j+1)-2.0d0*u(:,0,j)+u(:,0,j-1))
    ub = u(:,nx,j) - by*(u(:,nx,j+1)-2.0d0*u(:,nx,j)+u(:,nx,j-1))
    
	r(:,1)   = r(:,1) - a(:,1)*ua        !b.c.
	r(:,nx-1) = r(:,nx-1) - c(:,nx-1)*ub !b.c.
	call tdma(a,b,c,r,q,1,nx-1,dim_a)

	!assign solutions for as z
	do i=1,nx-1
	z(:,i,j)=q(:,i)
	end do
    z(:,0,j) =ua
    z(:,nx,j)=ub

!    deallocate(a,b,c,r,q)

end do
!deallocate(a,b,c,r,q)

	!allocate(a(1:dim_a,1:ny-1),b(1:dim_a,1:ny-1),c(1:dim_a,1:ny-1),r(1:dim_a,1:ny-1),q(1:dim_a,1:ny-1))
!y-sweep to compute final solution:
do i=1,nx-1

  
	!Build coefficient matrix:
	!allocate(a(1:dim_a,1:ny-1),b(1:dim_a,1:ny-1),c(1:dim_a,1:ny-1),r(1:dim_a,1:ny-1),q(1:dim_a,1:ny-1))

	do j=1,ny-1
	a(:,j) = -by
	b(:,j) = (1.0d0+2.0d0*by)
	c(:,j) = -by
	r(:,j) = z(:,i,j) 
	end do
    
	!apply boundary conditions
    ua = u(:,i,0)  
    ub = u(:,i,ny) 
    
	r(:,1)   = r(:,1) - a(:,1)*ua        !b.c.
	r(:,ny-1) = r(:,ny-1) - c(:,ny-1)*ub !b.c.
    
	call tdma(a,b,c,r,q,1,ny-1,dim_a)

	!assign solutions for as z
	do j=1,ny-1
	u(:,i,j)=q(:,j)
	end do

    !deallocate(a,b,c,r,q)
end do

!   deallocate(a,b,c,r,q)


end 


