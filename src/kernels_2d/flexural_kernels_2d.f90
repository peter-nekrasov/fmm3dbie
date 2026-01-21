subroutine modified_flex2d_g(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  real *8 :: dr 
  integer ipars(ndi)
  real *8 :: dx, dy, r2
  real *8 :: over4pi
  complex *16 :: val
  complex *16 :: zpars(2), zk1, zk2 
  complex *16 :: zt1, zt2 
  complex *16 :: g1, g2, h1 
  complex *16 :: ima 

  ima = (0,1)

  zk1 = zpars(1)
  zk2 = zpars(2)

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)

  dr = sqrt(dx**2+dy**2)
  zt1 = zk1*dr
  zt2 = zk2*dr

  call hank101(zt1,g1,h1)
  call hank101(zt2,g2,h1)

  val = ima*(g1-g2)/(4*(zk1*zk1-zk2*zk2))



end subroutine modified_flex2d_g
