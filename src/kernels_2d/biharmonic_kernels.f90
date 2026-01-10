
! New begins here.

subroutine bh2d_g(src, ndt,targ, ndd,dpars,ndz,zk,ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: src(*), targ(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2
  real *8 :: over4pi
  real *8 :: val
  data over4pi/0.07957747154594767d0/
  !
  ! returns the biharmonic volumetric potential kernel
  !

  dx=targ(1)-src(1)
  dy=targ(2)-src(2)

  r2=dx**2+dy**2

  !val = -over4pi*log(r2)
  val = r2 * log(r2) * over4pi / 4  

  return
end subroutine bh2d_g 

subroutine bh2d_gdn(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: val
  real *8 :: over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the normal derivative of the biharmonic volumetric kernel
  !

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  rdotn = dx*targinfo(10) + dy*targinfo(11)
  r2 = dx**2 + dy**2

  !val =  -2*rdotn/r2*over4pi
  val = rdotn*(log(r2)+1)*over4pi/2 

  return
end subroutine bh2d_gdn
