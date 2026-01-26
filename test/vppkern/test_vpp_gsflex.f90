
program test_vpp_gsflex
  implicit real *8 (a-h,o-z)

  integer ipars(5)
  real *8 :: dpars0
  complex * 16 :: zpars, zpars0(14), ima

  real *8, allocatable :: rtest(:), dpars(:)  
  complex *16, allocatable :: vtest(:,:), vtestt(:,:)

  real *8 :: src(3), targ(3)
  
  external fun1, gphiflexkern
  integer :: ipars0
  data ima /(0d0,1d0)/


  ldpars = 10000000
  allocate(dpars(ldpars))

  done = 1
  pi = atan(done)*4

  ipars0 = 0
  nf = 2
  a = 1d-12
  b = 2*pi
  n = 16
  tol = 1d-12
  maxsub = 1000
  maxdepth = 50
  ietype = 1
  ier = 0

  ! test build for kernel routine
  write(*,*) 'kernel version ....'

  zpars0(1) =  1.423605848552332d0 + 0.000000000000000d0*ima
  zpars0(2) =  0.246729256910564d0 + 1.320816347450247d0*ima
  zpars0(3) =  0.246729256910564d0 - 1.320816347450247d0*ima
  zpars0(4) = -0.958532181186730d0 + 0.498427779031846d0*ima
  zpars0(5) = -0.958532181186730d0 - 0.498427779031846d0*ima

  zpars0(6) =  0.539472550644810d0 + 0.000000000000000d0*ima
  zpars0(7) =  0.453879061238694d0 + 0.495600250144307d0*ima
  zpars0(8) =  0.453879061238694d0 - 0.495600250144307d0*ima
  zpars0(9) = -0.723615336561098d0 + 1.073365249353908d0*ima
  zpars0(10) = -0.723615336561098d0 - 1.073365249353908d0*ima

  ndi0 = 0
  ndd0 = 0
  ndz0 = 14
  nf = 2

        write(*,*) ndd0
        write(*,*) dpars0
        write(*,*) ndz0
        write(*,*) zpars0(1)
        write(*,*) ndi0
        write(*,*) ipars0
        write(*,*) nf
        write(*,*) a
        write(*,*) b
        write(*,*) n
        write(*,*) tol
        write(*,*) maxsub
        write(*,*) maxdepth
        write(*,*) ietype
        write(*,*) ipars(1)
        write(*,*) ndd
        write(*,*) ldpars
        write(*,*) dpars(1)
        write(*,*) ier

  ntest = 1000000  
  allocate(rtest(ntest),vtest(nf,ntest),vtestt(nf,ntest))

  call cpu_time(t1)
  call vpp_buildkern(gphiflexkern,ndd0,dpars0,ndz0,zpars0,ndi0,ipars0,&
       nf,a,b,n,tol,maxsub,maxdepth,ietype, &
       ipars,ndd,ldpars,dpars,ier)

  call cpu_time(t2)

  
        write(*,*) ndd0
        write(*,*) dpars0
        write(*,*) ndz0
        write(*,*) zpars0(1)
        write(*,*) ndi0
        write(*,*) ipars0
        write(*,*) nf
        write(*,*) a
        write(*,*) b
        write(*,*) n
        write(*,*) tol
        write(*,*) maxsub
        write(*,*) maxdepth
        write(*,*) ietype
        write(*,*) ipars(1)
        write(*,*) ndd
        write(*,*) ldpars
        write(*,*) dpars(1)
        write(*,*) ier


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
     targ(1) = rtest(j)
     call vpp_kern(src,n3,targ,ndd,dpars,ndz,zpars,ndi,ipars,vtest(1,j))
  enddo
  call cpu_time(t2)

  write(*,*) ntest/(t2-t1), 'vpp throughput '
  
  call cpu_time(t1)
  do j = 1,ntest
     targ(1) = rtest(j)
     call gphiflexkern(src,n3,targ,ndd0,dpars0,ndz0,zpars0,ndi0,ipars0, &
          vtestt(1,j))
  enddo
  call cpu_time(t2)

  write(*,*) ntest/(t2-t1), 'G_S flex throughput'

  do j = 1,min(ntest,10)
     write(*,*) abs(vtest(1,j) - vtestt(1,j) )
  enddo


  stop
end program test_vpp_gsflex

subroutine fun1(r,val)
  implicit real *8 (a-h,o-z)
  complex *16 :: val, zim
  data zim /(0d0,1d0)/

  val = log(r)*cos(10*r)*(2d0+ zim*sin(3d0*r))

  return
end subroutine fun1

