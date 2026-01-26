
program test_vpp
  implicit real *8 (a-h,o-z)

  integer ipars(5)
  real *8 :: dpars0
  complex * 16 :: zpars, zpars0, ima

  real *8, allocatable :: rtest(:), dpars(:)
  complex *16, allocatable :: vtest(:,:), vtestt(:,:)

  real *8 :: src(3), targ(3)
  
  external fun1, h3d_slp
  integer :: ipars0
  data ima /(0d0,1d0)/


  ldpars = 1000000
  allocate(dpars(ldpars))
  
  done = 1
  pi = atan(done)*4

  nf = 2
  a = 1d-12
  b = 2*pi
  n = 16
  tol = 1d-12
  maxsub = 1000
  maxdepth = 50
  ietype = 1
  ier = 0

  ! test build for simple function handle 
  call cpu_time(t1)
  call vpp_build(fun1,nf,a,b,n,tol,maxsub,maxdepth,ietype, &
       ipars,ndd,ldpars,dpars,ier)
  call cpu_time(t2)

  write(*,*) t2-t1
  

  write(*,*) ier
  do j = 1,5
     write(*,*) ipars(j)
  enddo

  ntest = 10000000

  allocate(rtest(ntest),vtest(nf,ntest),vtestt(nf,ntest))
  
  do j = 1,ntest
     h = (b-a)*(cos(j*1d0)+1d0)/2d0
     rtest(j) = a + h
  enddo

  src(:) = 0
  targ(:) = 0

  n3 = 3
  ndz = 0
  ndi = 5

  call cpu_time(t1)
  do j = 1,ntest
     targ(1) = rtest(j)
     call vpp_kern(src,n3,targ,ndd,dpars,ndz,zpars,ndi,ipars,vtest(1,j))
  enddo
  call cpu_time(t2)

  write(*,*) ntest/(t2-t1), 'vpp throughput '
  
  call cpu_time(t1)
  do j = 1,ntest
     call fun1(rtest(j),vtestt(1,j))
  enddo
  call cpu_time(t2)

  write(*,*) ntest/(t2-t1), 'fun1 throughput '

  do j = 1,min(ntest,10)
     write(*,*) abs(vtest(1,j)-vtestt(1,j))
  enddo

  ! test build for kernel routine
  write(*,*) 'kernel version ....'

  ndi0 = 0
  ndd0 = 0
  ndz0 = 1
  zpars0 = 6d0 + ima*1d0
  nf = 2
  call cpu_time(t1)
  call vpp_buildkern(h3d_slp,ndd0,dpars0,ndz0,zpars0,ndi0,ipars0,&
       nf,a,b,n,tol,maxsub,maxdepth,ietype, &
       ipars,ndd,ldpars,dpars,ier)

  call cpu_time(t2)

  write(*,*) t2-t1
  

  write(*,*) ier
  do j = 1,5
     write(*,*) ipars(j)
  enddo

  do j = 1,ntest
     h = (b-a)*(cos(j*1d0)+1d0)/2d0
     rtest(j) = a + h
  enddo

  src(:) = 0
  targ(:) = 0

  n3 = 3
  ndz = 0
  ndi = 5

  call cpu_time(t1)
  do j = 1,ntest
     call vpp_eval(rtest(j),ndd,dpars,ndi,ipars,vtest(1,j))
  enddo
  call cpu_time(t2)

  write(*,*) ntest/(t2-t1), 'vpp throughput '
  
  call cpu_time(t1)
  do j = 1,ntest
     targ(1) = rtest(j)
     call h3d_slp(src,n3,targ,ndd0,dpars0,ndz0,zpars0,ndi0,ipars0, &
          vtestt(1,j))
  enddo
  call cpu_time(t2)

  write(*,*) ntest/(t2-t1), 'e^(ikr)/r throughput'

  do j = 1,min(ntest,10)
     write(*,*) abs(vtest(1,j)-vtestt(1,j))
  enddo
  
  stop
end program test_vpp

subroutine fun1(r,val)
  implicit real *8 (a-h,o-z)
  complex *16 :: val, zim
  data zim /(0d0,1d0)/

  val = log(r)*cos(10*r)*(2d0+ zim*sin(3d0*r))

  return
end subroutine fun1

