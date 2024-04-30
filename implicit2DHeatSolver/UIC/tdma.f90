subroutine tdma(a,b,c,r,x,s,e,dim_a)
implicit none
integer s,e,i,dim_a
real*8, dimension(1:dim_a,s:e) ::a,b,c,r,x

! forward elimination phase
do i=s+1,e
b(:,i) = b(:,i) - a(:,i)/b(:,i-1)*c(:,i-1)
r(:,i) = r(:,i) - a(:,i)/b(:,i-1)*r(:,i-1)
end do
! backward substitution phase
x(:,e) = r(:,e)/b(:,e)
do i=e-1,s,-1
x(:,i) = (r(:,i)-c(:,i)*x(:,i+1))/b(:,i)
end do

return
end

