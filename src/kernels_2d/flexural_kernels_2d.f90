subroutine flex2d_g(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  real *8 :: dr 
  integer ipars(ndi)
  real *8 :: dx, dy, r2
  complex *16 :: val
  complex *16 :: zpars(2), zk1, zk2 
  complex *16 :: zt1, zt2 
  complex *16 :: g1, g2, h1 
  complex *16 :: h01,h0x1,h0y1
  complex *16 :: h02,h0x2,h0y2
  complex *16 :: h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy
  complex *16 :: gradx, grady


  zk1 = zpars(1)
  zk2 = zpars(2)

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)

  call helmdiffgreen(zk1,dx,dy,g1,h0x1,h0y1,h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)
  call helmdiffgreen(zk2,dx,dy,g2,h0x2,h0y2,h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)

  val = (g1-g2)/(zk1*zk1-zk2*zk2)


end subroutine flex2d_g
!
!
!
!
!
subroutine flex2d_gdn(src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars,val)
  real *8 :: src(*), targ(ndt), dpars(ndd)
  real *8 :: dr 
  integer ipars(ndi)
  real *8 :: dx, dy, r2
  real *8 :: nx, ny 
  complex *16 :: val
  complex *16 :: zpars(2), zk1, zk2 
  complex *16 :: zt1, zt2 
  complex *16 :: h01,h0x1,h0y1
  complex *16 :: h02,h0x2,h0y2
  complex *16 :: h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy
  complex *16 :: gradx, grady
  

  zk1 = zpars(1)
  zk2 = zpars(2)

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)

  nx = targ(10)
  ny = targ(11)


  call helmdiffgreen(zk1,dx,dy,h01,h0x1,h0y1,h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h0x2,h0y2,h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)

  gradx = h0x1-h0x2
  grady = h0y1-h0y2

  val = (nx*gradx + ny*grady)/(zk1*zk1-zk2*zk2)

end subroutine flex2d_gdn
!
!
!
!
!
subroutine flex2d_gsupp2(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk, &
   ndi,ipars,val)
  implicit real *8 (a-h,o-z)
  real *8 :: srcinfo(*),targinfo(ndt),dpars(1)
  integer ipars(ndi)
  real *8 :: dx, dy, r2, rdotn
  real *8 :: nu, nx, ny
  complex *16 :: zk(2), zk1, zk2 
  complex *16 :: val, gsxx, gsxy, gsyy
  complex *16 :: h01,h01x,h01y,h02,h02x,h02y
  complex *16 :: h01xx,h01xy,h01yy,h01xxx,h01xxy,h01xyy,h01yyy
  complex *16 :: h02xx,h02xy,h02yy,h02xxx,h02xxy,h02xyy,h02yyy


  dx = targinfo(1) - srcinfo(1)
  dy = targinfo(2) - srcinfo(2)
  r2 = dx**2 + dy**2

  zk1 = zk(1)
  zk2 = zk(2)

!  rdotn = dx*targinfo(10) + dy*targinfo(11)
  
  call helmdiffgreen(zk1,dx,dy,h01,h01x,h01y,h01xx,h01xy,h01yy,&
    h01xxx,h01xxy,h01xyy,h01yyy)
  call helmdiffgreen(zk2,dx,dy,h02,h02x,h02y,h02xx,h02xy,h02yy,&
    h02xxx,h02xxy,h02xyy,h02yyy)

  gsxx = (h01xx-h02xx)/(zk1*zk1-zk2*zk2)
  gsxy = (h01xy-h02xy)/(zk1*zk1-zk2*zk2)
  gsyy = (h01yy-h02yy)/(zk1*zk1-zk2*zk2)


  nx = targinfo(10)
  ny = targinfo(11)

  nu = dpars(1)  


  val = nu*(gsxx + gsyy) + &
    (1.0d0 - nu)*(nx*nx*gsxx + 2*nx*ny*gsxy + ny*ny*gsyy)

            
  return
end subroutine flex2d_gsupp2