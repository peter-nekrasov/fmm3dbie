      implicit real *8 (a-h,o-z)
      real *8 src(3), targ(12)
      real *8 val1, val2
      data ima /(0,1)/

      src(1) = 0
      src(2) = 0
      src(3) = 0

      targ(1) = 0.5
      targ(2) = -0.25
      targ(3) = 0

      targ(10) = 0.4
      targ(11) = 0.6

      call l2d_g(src,12,targ,0,0,0,0,0,0,val1)
      call l2d_gdn(src,12,targ,0,0,0,0,0,0,val2)

      open(unit=33,file='print_test_lap2d.txt')
      write(33,'(2x,e22.16)'),val1,val2
      close(33)

      stop
      end
