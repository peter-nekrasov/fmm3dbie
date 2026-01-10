      implicit real *8 (a-h,o-z)
      dimension src(12), targ(12) 
        
      call prini(6,13)
      src(1) = 1.2d0 
      src(2) = 2.3d0 

      targ(1) = 2.5d0 
      targ(2) = 3.5d0 
      
      targ(10) = 0.3d0 
      targ(11) = -0.2d0
          


      ndt = 12 
      ndd = 1 
      dpars = 1 
      ndz = 1 
      
      zk = 1 
      ndi = 1 
      ipars = 1 
      

      val = 0 
      call bh2d_g(src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipars,val)
      call prin2_long('val=*',val,1)
      
      val = 0 
      call bh2d_gdn(src,ndt,targ,ndd,dpars,ndz,zk,ndi,ipals,val)
      call prin2_long('val=*',val,1)
      end 
