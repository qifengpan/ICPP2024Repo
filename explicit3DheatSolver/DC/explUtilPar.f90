
!************************************************************************!
!                                                                        !
! This subroutin computes next values in subdomain of current process me !
!                                                                        !
!************************************************************************!

subroutine computeNext(x0, x, lda1, lda2, lda3, dt, hx, hy, hz, r, me, &
                       xs, ys, zs, xe, ye, ze, nproc, k0)

   ! Input parameters
   integer lda1, lda2, lda3, me, nproc
   integer xs(0:nproc-1), ys(0:nproc-1), xe(0:nproc-1), ye(0:nproc-1), ze(0:nproc-1), &
           zs(0:nproc-1)
   double precision x(0:lda1-1,0:lda2-1,0:lda3-1)
   double precision dt, hx, hy, hz, k0

   ! Output array
   double precision x0(0:lda1-1,0:lda2-1,0:lda3-1)

   ! Output error
   double precision r
   ! Local variable for computing error
   double precision rk

   ! Index variables
   integer i, j, k

   ! Factors for the stencil
   double precision diagx, diagy, diagz, weightx, weighty, weightz

   !   The stencil of the explicit operator for the heat equation
   !   on a regular rectangular grid using a seven point finite difference
   !   scheme in space is :
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
   do 31, k=zs(me),ze(me)
    do 20, j=ys(me),ye(me)
     do 30, i=xs(me),xe(me)

      x(i,j,k) = weightx*(x0(i-1,j,k) + x0(i+1,j,k) + x0(i,j,k)*diagx) + &
                 weighty*(x0(i,j-1,k) + x0(i,j+1,k) + x0(i,j,k)*diagy) + &
                 weightz*(x0(i,j,k-1) + x0(i,j,k+1) + x0(i,j,k)*diagz)

30   continue
20  continue
31 continue

   ! Copy back the computed value : x0(n) <-- x(n)
   ! and compute the 2_norm of the 'residual'
   r = 0.0
   do 51, k=zs(me),ze(me)
    do 40, j=ys(me),ye(me)
     do 39, i=xs(me),xe(me)

      rk = x0(i,j,k) - x(i,j,k)
      r = r + rk*rk
      x0(i,j,k) = x(i,j,k)

39   continue
40  continue
51 continue

   return
end

!************************************************************************!
!                                                                        !
! This subroutine sets up the initial temperatures on borders and inside !
!                                                                        !
!************************************************************************!

subroutine initValues(nb_layers, x0, size_tot_x, size_tot_y, size_tot_z, temp1_init, temp2_init)

   ! Input parameters
   integer size_tot_x, size_tot_y, size_tot_z
   ! Init Temperatures
   double precision temp1_init, temp2_init

   ! Output array
   double precision x0(0:size_tot_x-1,0:size_tot_y-1,0:size_tot_z-1)

   ! Index variables
   integer i, j, k, l

   ! Number of initial borders layers
   integer nb_layers

   ! Setup temp1_init on borders
   do 788, l=1,nb_layers+1

    do 10, j=(l-1),size_tot_y-l
     do 12, k=(l-1),size_tot_z-l

      x0(l-1,j,k) = temp1_init
      x0(size_tot_x-l,j,k) = temp1_init

12   continue
10  continue

    do 11, i=(l-1),size_tot_x-l
     do 13, k=(l-1),size_tot_z-l

      x0(i,l-1,k) = temp1_init
      x0(i,size_tot_y-l,k) = temp1_init

13   continue
11  continue

    do 15, i=(l-1),size_tot_x-l
     do 14, j=(l-1),size_tot_y-l

      x0(i,j,l-1) = temp1_init
      x0(i,j,size_tot_z-l) = temp1_init

14   continue
15  continue

788 continue

   ! Setup temp2_init inside
   do 60, i=nb_layers+1,size_tot_x-(nb_layers+2)
    do 70, j=nb_layers+1,size_tot_y-(nb_layers+2)
     do 71, k=nb_layers+1,size_tot_z-(nb_layers+2)

            x0(i,j,k) = temp2_init

71   continue
70  continue
60 continue

   return
end

!**************************************************************!
!                                                              !
! This subroutine computes the coordinates xs, xe, ys, ye, zs, !
! ze for each cell on the grid                                 !
!                                                              !
!**************************************************************!

subroutine processToMap(xs, xe, ys, ye, zs, ze, xcell, ycell, zcell, &
                        x_domains, y_domains, z_domains, nproc)

   ! Input parameters
   integer nproc
   integer xcell, ycell, zcell, x_domains, y_domains, z_domains

   ! Output arrays
   integer xs(0:nproc-1), xe(0:nproc-1), ys(0:nproc-1), ye(0:nproc-1), &
           zs(0:nproc-1), ze(0:nproc-1)

   ! Index variables
   integer i, j, k, l, m, v, p

   ! Computation of xs and xe with processes topology
   xs(0:(z_domains*y_domains)-1) = 2
   xe(0:(z_domains*y_domains)-1) = xs(0:(z_domains*y_domains)-1)+xcell-1

   do 774, j=1,x_domains-1

    xs(j*(z_domains*y_domains):(j+1)*(z_domains*y_domains)-1) = xs((j-1)*(z_domains*y_domains))+xcell+2
    xe(j*(z_domains*y_domains):(j+1)*(z_domains*y_domains)-1) = xs(j*(z_domains*y_domains))+xcell-1

   ! Computation of ys and ye with processes topology
   do 700, i=1,y_domains

    ys((i-1)*z_domains) = y_domains*(ycell+2)-ycell*i-2*(i-1)
    ye((i-1)*z_domains) = ys((i-1)*z_domains)+ycell-1

    do 780, l=1,z_domains-1

     ys((i-1)*z_domains+l) = ys((i-1)*z_domains)
     ye((i-1)*z_domains+l) = ys((i-1)*z_domains+l)+ycell-1

780 continue
700 continue
774 continue

   ! Prolongation along y_domain
   do 703, m=1,y_domains

    ys((m-1)*z_domains) = y_domains*(ycell+2)-ycell*m-2*(m-1)
    ye((m-1)*z_domains) = ys((m-1)*z_domains)+ycell-1

    do 701, i=1,x_domains-1

     ys(i*(y_domains*z_domains)+(m-1)*z_domains) = ys((m-1)*z_domains)
     ye(i*(y_domains*z_domains)+(m-1)*z_domains) = &
           ys(i*(y_domains*z_domains)+(m-1)*z_domains)+ycell-1

     do 702, l=1,z_domains-1

      ys(i*(y_domains*z_domains)+(m-1)*z_domains+l) = &
      ys(i*(y_domains*z_domains)+(m-1)*z_domains)

      ye(i*(y_domains*z_domains)+(m-1)*z_domains+l) = &
      ys(i*(y_domains*z_domains)+(m-1)*z_domains+l)+ycell-1

702  continue
701 continue
703 continue

   ! Computation of zs and ze with processes topology
   do 704, k=0,y_domains-1

    v = k*z_domains
    zs(v) = 2
    ze(v) = 2+zcell-1

    do 705, p=1,x_domains-1

     zs(v+p*(y_domains*z_domains)) = zs(v)
     ze(v+p*(y_domains*z_domains)) = ze(v)

705 continue
704 continue

   ! Prolongation along z_domain
   do 706, m=1,z_domains-1
    do 707, i=0,y_domains-1

     l = m+i*z_domains
     zs(l) = zs(l-1)+zcell+2
     ze(l) = zs(l)+zcell-1

     do 708, v=1,x_domains-1

      zs(l+v*(y_domains*z_domains)) = zs(l)
      ze(l+v*(y_domains*z_domains)) = zs(l+v*(y_domains*z_domains))+zcell-1

708  continue
707 continue
706 continue

end
