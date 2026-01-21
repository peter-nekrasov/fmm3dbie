subroutine modified_flex2d_g(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  real *8 :: dr 
  integer ipars(ndi)
  real *8 :: dx, dy, r2
  complex *16 :: val
  complex *16 :: zpars(2), zk1, zk2 
  complex *16 :: zt1, zt2 
  complex *16 :: g1, g2, h1 


  zk1 = zpars(1)
  zk2 = zpars(2)

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)



  call helmdiffgreen_g(zk1,dx,dy,g1)
  call helmdiffgreen_g(zk2,dx,dy,g2)

  val = (g1-g2)/(zk1*zk1-zk2*zk2)



end subroutine modified_flex2d_g
