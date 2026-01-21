      implicit real *8 (a-h,o-z)
      complex *16 zk, zpars(2)
      complex *16 h0,h0x,h0y
      complex *16 h0xx,h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy
      real *8 dx,dy 
      complex *16 val ,z 
      dimension src(2), targ(2) 
        
      call prini(6,13)
      src(1) = 1.2d0 
      src(2) = 2.3d0 

      targ(1) = 2d0 
      targ(2) = 3.5d0 
      


      dx = targ(1)-src(1)
      dy = targ(2)-src(2)
      
      

      if (1.eq.1) then 
        zk = (10d0,0)

c        call helmdiffgreen(zk,dx,dy,h0,h0x,h0y,h0xx,
c     1      h0xy,h0yy,h0xxx,h0xxy,h0xyy,h0yyy)
        call helmdiffgreen_g(zk,dx,dy,h0)

        call prin2_long('h0=*',h0,2)

        if (1.eq.0) then 
        call prin2_long('h0=*',h0,2)
        call prin2_long('h0x=*',h0x,2)
        call prin2_long('h0y=*',h0y,2)
        call prin2_long('h0xx=*',h0xx,2)
        call prin2_long('h0xy=*',h0xy,2)
        call prin2_long('h0yy=*',h0yy,2)
        call prin2_long('h0xxx=*',h0xxx,2)
        call prin2_long('h0xxy=*',h0xxy,2)
        call prin2_long('h0xyy=*',h0xyy,2)
        call prin2_long('h0yyy=*',h0yyy,2)
        endif 


      endif 


ccc   test green's function of modified flex2d problem
      if (1.eq.0) then 

        src(1) = 1.2d0 
        src(2) = 2.3d0 

        targ(1) = 2.5d0 
        targ(2) = 3.5d0 

        zpars(1) = 1
        zpars(2) = 2
        call modified_flex2d_g(src,2,targ,0,0,0,zpars,0,0,val)

        z = (0.024751196814465d0,  -0.061915622771211d0)
        call prin2_long('gfunc of modified flex (fortran) = *',val,2)
        call prin2_long('gfunc of modified flex (matlab) = *',z,2)

        err = abs(z-val)/abs(z)
        call prin2_long('relative error is *',err,1)
      endif  


      end 
