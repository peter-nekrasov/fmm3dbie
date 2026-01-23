      implicit real *8 (a-h,o-z)
      complex *16 zk, zpars(2)
      complex *16 h0,h0x,h0y
      complex *16 h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy
      real *8 dx,dy 
      complex *16 val, z, zt, ima, cvals(4,1)
      dimension src(2), targ(12) 
        
      call prini(6,13)
      
      ima = (0,1)
      
      


      src(1) = 1.2d0 
      src(2) = 2.3d0 

      targ(1) = 2d0 
      targ(2) = 3.5d0  


      dx = targ(1)-src(1)
      dy = targ(2)-src(2)

      dr = sqrt(dx**2+dy**2)

      zk = 0.01d0

      call helmdiffgreen(zk,dx,dy,h0,h0x,h0y,h0xx,
     1      h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)


      call prin2_long('h0=*',h0,2)

      end 
