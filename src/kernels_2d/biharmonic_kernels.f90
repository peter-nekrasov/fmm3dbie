
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

subroutine bh2d_gsupp2(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: val, nu, gsxx, gsxy, gsyy
  real *8 :: over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the second supported plate condition of the 
  ! biharmonic volumetric kernel
  !

  call bh2d_green_hess(dx,dy,gxx,gxy,gyy)

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  taux = targinfo(4)
  tauy = targinfo(5)

  nx = targinfo(10)
  ny = targinfo(11)

  nu = dpars(1)  

  val = gxx*nx*nx + 2*gxy*nx*ny + gyy*ny*ny 
     1  nu*(gxx*taux*taux + 2*gxy*taux*tauy + gyy*tauy*tauy)
            
  return
end subroutine bh2d_gsupp2

subroutine bh2d_gfree1(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(ndd)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: val, nu, gsxx, gsxy ,gsyy
  real *8 :: over4pi
  data over4pi/0.07957747154594767d0/
  !
  ! returns the first free plate condition of the 
  ! biharmonic volumetric kernel
  !

  call bh2d_green_hess(dx,dy,gxx,gxy,gyy)

  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)

  taux = targinfo(4)
  tauy = targinfo(5)

  nx = targinfo(10)
  ny = targinfo(11)

  nu = dpars(1)

  val = gxx*nx*nx + 2*gxy*nx*ny + gyy*ny*ny
     1  nu*(gxx*taux*taux + 2*gxy*taux*tauy + gyy*tauy*tauy)

  return
end subroutine bh2d_gfree1

