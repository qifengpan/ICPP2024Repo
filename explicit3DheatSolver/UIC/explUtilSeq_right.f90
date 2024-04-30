
!*******************************************************|
!                                                       |
! This subroutine computes next values on global domain |
!                                                       |
!*******************************************************|

subroutine computeNext(x0, x, size_x, size_y, size_z, dt, hx, hy, hz, r, k0,nSample)

   ! Input parameters
   integer nSample
   integer size_x, size_y, size_z
   double precision x(0:size_x+1,0:size_y+1,0:size_z+1,1:nSample)
   double precision dt, hx, hy, hz
   double precision,dimension(1:nSample) :: k0
   ! Output array
   double precision x0(0:size_x+1,0:size_y+1,0:size_z+1,1:nSample)

   ! Index variables
   integer i, j, k,l

   ! Output parameter for error
   ! Local variable for computing error
   double precision,dimension(1:nSample):: r,rk

   ! Factors for the stencil
   double precision,dimension(1:nSample):: diagx, diagy, diagz, weightx, weighty, weightz

   ! The stencil of the explicit operator for the heat equation
   ! on a regular rectangular grid using a seven point finite difference
   ! scheme in space is :
   !
   ! |                                       wx * x[i-1][j][k]     wz * x[i][j][k-1]               |
   ! |                                                                  /                          |
   ! | wy * x[i][j-1][k]   (diagx * wx + diagy * wy + diagz * wz) * x[i][j][k]   wx * x[i][j+1][k] |
   ! |                            /                                                                |
   ! |                 wz * x[i][j][k+1]     wy * x[i+1][j][k]                                     |

   diagx = -2.0 + hx*hx/(3*k0*dt)
   diagy = -2.0 + hy*hy/(3*k0*dt)
   diagz = -2.0 + hz*hz/(3*k0*dt)
   weightx = k0*dt/(hx*hx)
   weighty = k0*dt/(hy*hy)
   weightz = k0*dt/(hz*hz)

   ! Perform an explicit update on the points within the domain
   do l = 1,nSample
   do 31, k=1,size_z
    do 20, j=1,size_y
     do 30, i=1,size_x
      
      x(i,j,k,l) = weightx(l)*(x0(i-1,j,k,l) + x0(i+1,j,k,l) + x0(i,j,k,l)*diagx(l)) + &
                   weighty(l)*(x0(i,j-1,k,l) + x0(i,j+1,k,l) + x0(i,j,k,l)*diagy(l)) + &
                   weightz(l)*(x0(i,j,k-1,l) + x0(i,j,k+1,l) + x0(i,j,k,l)*diagz(l))

30   continue
20  continue
31 continue
end do

   ! Copy back the computed value : x0(n) <-- x(n)
   ! and compute the 2_norm of the 'residual'
   r = 0.0
   do l = 1,nSample
   do 51, k=1,size_z
    do 40, j=1,size_y
     do 50, i=1,size_x
      rk(l) = x0(i,j,k,l) - x(i,j,k,l)
      r(l) = r(l) + rk(l)*rk(l)
      x0(i,j,k,l) = x(i,j,k,l)
50   continue
40  continue
51 continue
end do
   return
end

!***********************************************************************!
!                                                                       !
! This subroutine sets up the initial temperature on borders and inside !
!                                                                       !
!***********************************************************************!

subroutine initValues(nb_layers, x0, x_dim, y_dim, z_dim, temp1_init, temp2_init,nSample)

   ! Input parameters
   integer nSample
   integer x_dim, y_dim, z_dim, nb_layers
   ! Init temperatures
   double precision temp1_init, temp2_init

   ! Output array
   double precision x0(0:x_dim+1,0:y_dim+1,0:z_dim+1,1:nSample)

   ! Index variables
   integer i, j, k, l,np

   ! Setup temp1_init on borders
   do 788, l=1,nb_layers
    do np = 1,nSample
    do 10, j=(l-1),y_dim+1-(l-1)
     do 12, k=(l-1),z_dim+1-(l-1)

      x0(l-1,j,k,np) = temp1_init
      x0(x_dim+1-(l-1),j,k,np) = temp1_init

12   continue
10  continue
    end do
    do np = 1,nSample
    do 11, i=(l-1),x_dim+1-(l-1)
     do 13, k=(l-1),z_dim+1-(l-1)

      x0(i,l-1,k,np) = temp1_init
      x0(i,y_dim+1-(l-1),k,np) = temp1_init

13   continue
11  continue
    end do
    do np = 1,nSample
    do 20, i=(l-1),y_dim+1-(l-1)
     do 21, k=(l-1),x_dim+1-(l-1)

      x0(k,i,l-1,np) = temp1_init
      x0(k,i,z_dim+1-(l-1),np) = temp1_init

21   continue
20  continue
     end do
788 continue

   ! Setup temp2_init inside
   do np = 1,nSample
   do 60, i=nb_layers,x_dim-(nb_layers-1)
    do 70, j=nb_layers,y_dim-(nb_layers-1)
     do 71, k=nb_layers,z_dim-(nb_layers-1)

      x0(i,j,k,np) = temp2_init

71   continue
70  continue
60 continue
    end do
   return
end

