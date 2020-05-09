c       This is the spectral ccde for three-D isotropic turbulence
c       using the message-passing in SP machines
c       Oct. 18, 1994 - December, 1994
c       lx = ly = lz = 2**n ( n = 4,5,...)
c       lx = x-diretion mesh size in each node
c       ly = y-diretion mesh size in each node
c       lz = z-direction mesh size in each node
c       nproc = number of process (nodes)
c       

      program three_spectral 
      complex*16,allocatable,dimension(:,:,:)::vx,vy,vz
      complex*16,allocatable,dimension(:,:,:)::wx,wy,wz
      complex*16,allocatable,dimension(:,:,:)::ox,oy,oz
      real*8,allocatable,dimension(:,:,:)::kx,ky,kz,k2
      real*8,allocatable,dimension(:,:,:)::tmp,tmp1,k2_e
      real*8,allocatable,dimension(:)::seed 
      integer*4 isign,icontxt,ip(40),nstep,lx,ly,lz,nnodes
      integer*4 lxp,lzp,lxp2,lzp2,nproc
      integer*4 ieout,itout,nek
      real*8 t_ex,t_ey,t_ez,t_wx,t_wy,t_wz
      real*8 to_ex,to_ey,to_ez,to_wx,to_wy,to_wz
      real*8 scale,scale_r,scale_s
      real*8 tout,eout,dt,rnu,u0,time,sd,rr,ek,e_t
      real*8 force(2),pi,cc1
      real*8 time0,time1,time2,time3,time4,time5,time6

      character*40 fdata(10),filename
      logical*4 new,timing

      data   fdata /'/database/syc/E1.sp',
     1 '/database/syc/E2.sp','/database/syc/E3.sp',
     2 '/database/syc/E4.sp','/database/syc/E5.sp',
     3 '/database/syc/E6.sp','/database/syc/E7.sp',
     4 '/database/syc/E8.sp','/database/syc/E9.sp',
     5 '/database/syc/E10.sp'/

      data istep, is, ak0 / -1, -1, 4.75683 /
      
      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2
       
      external d_vadd

c  setup MPP environment

        call blacs_pinfo(m,nnodes)
        call blacs_gridinit(icontxt,1,nnodes)

        nallgrp = 0
        call mp_sync(nallgrp)

c   read parameters from node m = 0

      pi = 4.D0*datan(1.D0)
 
      if(m.eq.0)then
        open(90,file='input_sp',status='unknown')
        write(*,*)'enter lx,ly,lz'
        read(90,*)lx,ly,lz
        write(*,*)'enter number of processors'
        read(90,*)nproc
        close(90)
        write(*,*)'lx = ',lx,' ly =',ly,' lz =',lz
        write(*,*)'nproc = ',nproc
      endif 

      call mp_bcast(lx,4,0,nallgrp) 
      call mp_bcast(ly,4,0,nallgrp) 
      call mp_bcast(lz,4,0,nallgrp) 
      call mp_bcast(nproc,4,0,nallgrp)
      
      lxp = lx/nproc
      lzp = lz/nproc
      lxp2 = lxp/2
      lzp2 = lzp/2

      if(m.eq.0)then
       write(*,*)'lx/nproc = ',lxp
       write(*,*)'lz/nproc = ',lzp
       write(*,*)'lxp2 = ',lxp2
       write(*,*)'lzp2 = ',lzp2
      end if 
  
      scale = 1.D0/(lx*ly*lz)
      scale_s = dsqrt(scale)
      scale_r = 1.D0
    
c allocate memory 
      
      allocate( vx(lz,ly,lxp2) )
      allocate( vy(lz,ly,lxp2) )
      allocate( vz(lz,ly,lxp2) )
      allocate( wx(lz,ly,lxp2) )
      allocate( wy(lz,ly,lxp2) )
      allocate( wz(lz,ly,lxp2) )
      allocate( ox(lz,ly,lxp2) )
      allocate( oy(lz,ly,lxp2) )
      allocate( oz(lz,ly,lxp2) )

      allocate( kx(lz,ly,lxp2) )
      allocate( ky(lz,ly,lxp2) )
      allocate( kz(lz,ly,lxp2) )
      allocate( k2(lz,ly,lxp2) )
      allocate( k2_e(lz,ly,lxp2) )
      allocate( tmp(lz,ly,lxp2) )
      allocate( tmp1(lz,ly,lxp2) )
      allocate( seed(nproc) )

      if(m.eq.0)print*,' memory allocated'

        is = 0
        iii = 1
        jjj = 1
        kkk = 1

      if(m.eq.0)then

       open(1,file='parameter.d',status='old')
       read(1,*) nstep,tout,eout,dt,rnu,u0,time
     1    ,new,timing,force,sd
       write(*,*)nstep,tout,eout,dt,rnu,u0,time,
     1    new,timing,force,sd
       close(1)

       open(8,file='input.file',status='old')
       write(*,*)'enter filename'
       read(8,1005)filename
       write(*,*)filename
       close(8)
      end if

1005   format(a)
 
        call mp_bcast(nstep,4,0,nallgrp)
        call mp_bcast(tout,8,0,nallgrp)
        call mp_bcast(eout,8,0,nallgrp)
        call mp_bcast(dt,8,0,nallgrp)
        call mp_bcast(rnu,8,0,nallgrp)
        call mp_bcast(u0,8,0,nallgrp)
        call mp_bcast(time,8,0,nallgrp)
        call mp_bcast(new,4,0,nallgrp)
        call mp_bcast(timing,4,0,nallgrp)
        call mp_bcast(force,16,0,nallgrp)
        call mp_bcast(filename,40,0,nallgrp)
        write(*,*)'m=',m,filename
        call mp_bcast(sd,8,0,nallgrp)

        ieout = int(eout/dt)
        itout = int(tout/dt)
        nek=0.9*sqrt(lx*lx+.25*ly*ly+.25*lz*lz)
        
        if(m.eq.0)then
         write(*,*)ieout,itout,nek
         open(11,file='energy.d',status='unknown')
         open(20,file='spectrum.d',status='unknown') 
         open(22,file='transfer.d',status='unknown') 
         open(70,file='initial.sp',status='unknown')
        end if
        
c---Generate initial condition -------

      if(new) then

c---Generate random number seed for each process ---

        if(m.eq.0)then

         call srand(sd)
         do i=1,nproc
           seed(i) = rand()
         end do

        end if 

       call mp_bcast(seed,nproc*8,0,nallgrp)

       rr = seed(m+1)
       call srand(rr)

c--------------------------------------------
         
         call gaussian(vx,pi)
         call gaussian(vy,pi)
         call gaussian(vz,pi)

        call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,m)

         tmp = (kx*real(vx) + ky*real(vy) + kz*real(vz))/k2
         vx = cmplx(real(vx) - kx*tmp, aimag(vx))
         vy = cmplx(real(vy) - ky*tmp, aimag(vy))
         vz = cmplx(real(vz) - kz*tmp, aimag(vz))
         tmp = (kx*aimag(vx) + ky*aimag(vy) + kz*aimag(vz))/k2
         vx = cmplx(real(vx), aimag(vx) - kx*tmp)
         vy = cmplx(real(vy), aimag(vy) - ky*tmp)
         vz = cmplx(real(vz), aimag(vz) - kz*tmp)
      
         cc1 = u0 * sqrt(8.*sqrt(2./pi)/(3.*pi*ak0**5))
         tmp = cc1 * dsqrt(k2) * dexp (-k2/(ak0*ak0))

         vx = vx * tmp
         vy = vy * tmp
         vz = vz * tmp

         call symmetrize(vx,m)
         call symmetrize(vy,m)
         call symmetrize(vz,m)

c         call dealiasing(vx,k2)
c         call dealiasing(vy,k2)
c         call dealiasing(vz,k2)
        
c     calculating initial spectrum

        tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
        do i=1,nek
           ek = 0.D0
           e_t = 0.D0
           ek=sum(tmp,mask=(abs(dsqrt(k2)-i-0.49999).lt.0.5))
           call mp_reduce(ek,e_t,8,0,d_vadd,nallgrp)
         if(m.eq.0)then
            write(70,*)i,e_t
         end if
        end do

        if(m.eq.0)close(70)

      else

c ------ Rreading from DISK ------ 

      open(2,file=filename,status='old',
     1    form='unformatted')
         read(2)vx
         read(2)vy
         read(2)vz

        nallgrp = 0
        call mp_sync(nallgrp)

        close(2)

        if(m.eq.0)then
          write(*,*)'after reading'
        end if
c
      endif

      call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,m)

      ox=vx
      oy=vy
      oz=vz

      dt_h=.5*dt
      ip = 0

 1    istep = istep + 1

c   write out the total energy in K space

       time1 = gettime() 
       tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
       ek = sum(tmp)
       call mp_reduce(ek,e_t,8,0,d_vadd,nallgrp)
       if(m.eq.0)then
       write(*,*)istep,time,e_t
       write(76,*)time,e_t
       end if
c
      if(mod(istep,itout).eq.0.and.istep.ne.0) then

      open(2,file=fdata(iii),status='unknown',
     1   form='unformatted')
          write(2)vx
          write(2)vy
          write(2)vz
        nallgrp = 0
        call mp_sync(nallgrp)

        close(2)

         iii = iii + 1

      endif

      wx = (0.,1.) * (ky*vz - kz*vy)
      wy = (0.,1.) * (kz*vx - kx*vz)
      wz = (0.,1.) * (kx*vy - ky*vx)

c     energy spectrum calculation 
     
      if(mod(istep,ieout).eq.0) then

       tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
         
       if(m.eq.0)then   
         write(20,*)(kkk-1)
       end if

        do i=1,nek
           ek = 0.D0
           e_t = 0.D0
           ek=sum(tmp,mask=(abs(dsqrt(k2)-i-0.49999).lt.0.5))
           call mp_reduce(ek,e_t,8,0,d_vadd,nallgrp)

         if(m.eq.0)then
            write(20,*)i,e_t
         end if

       end do
         kkk = kkk + 1
      end if
     
      kx = real(vx)
      ky = aimag(vx)
      kz = real(vy)
      k2 = aimag(vy)
      k2_e = real(vz)
      tmp1 = aimag(vz)
c
      ip(1) = 0
      isign = -1
c      
      call pdcrft3(vx,vx,lx,ly,lz,isign,scale_r,icontxt,ip) 
      call pdcrft3(vy,vy,lx,ly,lz,isign,scale_r,icontxt,ip)
      call pdcrft3(vz,vz,lx,ly,lz,isign,scale_r,icontxt,ip)

      call pdcrft3(wx,wx,lx,ly,lz,isign,scale_r,icontxt,ip)
      call pdcrft3(wy,wy,lx,ly,lz,isign,scale_r,icontxt,ip)
      call pdcrft3(wz,wz,lx,ly,lz,isign,scale_r,icontxt,ip)

c  ---Check energy and vorticity in physical space -----
 
c      t_ex = 0.
c      t_ey = 0.
c      t_ez = 0.
      
c      t_wx = 0.
c      t_wy = 0.
c      t_wz = 0.

c        t_ex=sum(vx*vx)
c        t_ey=sum(vy*vy)
c        t_ez=sum(vz*vz)

c        t_wx=sum(wx*wx)
c        t_wy=sum(wy*wy)
c        t_wz=sum(wz*wz)

c       call mp_reduce(t_ex,to_ex,8,0,d_vadd,nallgrp)
c       call mp_reduce(t_ey,to_ey,8,0,d_vadd,nallgrp)
c       call mp_reduce(t_ez,to_ez,8,0,d_vadd,nallgrp)
 
c       call mp_reduce(t_wx,to_wx,8,0,d_vadd,nallgrp)
c       call mp_reduce(t_wy,to_wy,8,0,d_vadd,nallgrp)
c       call mp_reduce(t_wz,to_wz,8,0,d_vadd,nallgrp)

c      if(m.eq.0)then
c
c         to_ex = to_ex*scale
c        to_ey = to_ey*scale
c         to_ez = to_ez*scale
c
c         to_wx = to_wx*scale
c         to_wy = to_wy*scale
c         to_wz = to_wz*scale
c
c         to_energy = to_ex + to_ey + to_ez
c         to_omega  = to_wx + to_wy + to_wz

c         write(11,*)time,to_ex,to_ey,to_ez,0.5*to_energy
c         write(11,*)time,to_wx,to_wy,to_wz,0.5*to_omega

c      end if
       
      tmp = real(vy)*real(wz) - real(vz)*real(wy)
      wz = cmplx(real(vz)*real(wx) - real(vx)*real(wz), aimag(wz))
      wx = cmplx(real(vx)*real(wy) - real(vy)*real(wx), aimag(wx))
      wy = cmplx(tmp, aimag(wy))
      tmp = aimag(vy)*aimag(wz) - aimag(vz)*aimag(wy)
      wz = cmplx(real(wz), aimag(vz)*aimag(wx) - aimag(vx)*aimag(wz))
      wx = cmplx(real(wx), aimag(vx)*aimag(wy) - aimag(vy)*aimag(wx))
      wy = cmplx(real(wy), tmp)

      vx = cmplx(kx,ky)
      vy = cmplx(kz,k2)
      vz = cmplx(k2_e,tmp1)

      call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,m)

      isign = 1

      ip(1) = 0

c      wx = wx*scale_s
c      wy = wy*scale_s
c      wz = wz*scale_s

      call pdrcft3(wx,wx,lx,ly,lz,isign,scale_r,icontxt,ip)
      call pdrcft3(wy,wy,lx,ly,lz,isign,scale_r,icontxt,ip)
      call pdrcft3(wz,wz,lx,ly,lz,isign,scale_r,icontxt,ip)

        wx = wx*scale
        wy = wy*scale
        wz = wz*scale

       call symmetrize(wx,m)
       call symmetrize(wy,m)
       call symmetrize(wz,m)

c       call dealiasing(vx,k2)
c       call dealiasing(vy,k2)
c       call dealiasing(vz,k2)

      tmp = (kx*real(wy) + ky*real(wz) + kz*real(wx))/k2
      wy = cmplx(real(wy) - kx*tmp, aimag(wy))
      wz = cmplx(real(wz) - ky*tmp, aimag(wz))
      wx = cmplx(real(wx) - kz*tmp, aimag(wx))
      tmp = (kx*aimag(wy) + ky*aimag(wz) + kz*aimag(wx))/k2
      wy = cmplx(real(wy), aimag(wy) - kx*tmp)
      wz = cmplx(real(wz), aimag(wz) - ky*tmp)
      wx = cmplx(real(wx), aimag(wx) - kz*tmp)

       if(mod(istep,ieout).eq.0) then
        tmp = real(vx)*real(wy) + real(vy)*real(wz) + 
     &        real(vz)*real(wx) +
     &        aimag(vx)*aimag(wy) + aimag(vy)*aimag(wz) + 
     &        aimag(vz)*aimag(wx)
 
         if(m.eq.0)then
          write(22,*)(jjj-1)
         end if
 
         do i=1,nek

           ek = 0.
           e_t = 0. 
           ek=sum(tmp,mask=(abs(dsqrt(k2)-i-0.49999).lt.0.5))
           call mp_reduce(ek,e_t,8,0,d_vadd,nallgrp)
             
           if(m.eq.0)then
           write(22,*) i,ek
           end if
 
         end do

         jjj = jjj + 1

       endif

      if(istep.eq.0) then

         vx = vx * (1.-rnu*dt_h*k2) + dt_h*wy
         vy = vy * (1.-rnu*dt_h*k2) + dt_h*wz
         vz = vz * (1.-rnu*dt_h*k2) + dt_h*wx

       call supfor(vx,vy,vz,tmp,k2,force,m)

        time2 = gettime()

       if(m.eq.0)then
          write(6,fmt='(5X,A,F20.4,A)')
     &   'elapsed time =',1.0D3*(time2-time1),' msec'
        endif


       open(2,file='/database/syc/scratch',status='unknown',
     1   form='unformatted')
            write(2)wy
            write(2)wz
            write(2)wx
       close(2)

         time = time + dt_h

      elseif(istep.eq.1) then
 
         vx = ox + dt*wy - rnu*dt*k2*vx
         vy = oy + dt*wz - rnu*dt*k2*vy
         vz = oz + dt*wx - rnu*dt*k2*vz

         call supfor(vx,vy,vz,tmp,k2,force,m)

         time2 = gettime()

       if(m.eq.0)then
          write(6,fmt='(5X,A,F20.4,A)')
     &   'elapsed time =',1.0D3*(time2-time1),' msec'
        endif

       open(2,file='/database/syc/scratch',status='old',
     1   form='unformatted')
           read(2)ox
           read(2)oy
           read(2)oz
        close(2)

         time = time + dt_h

      else

         vx = vx + dt_h*(3.*wy - k2_e*ox)
         vy = vy + dt_h*(3.*wz - k2_e*oy)
         vz = vz + dt_h*(3.*wx - k2_e*oz)

         vx = vx*k2_e
         vy = vy*k2_e
         vz = vz*k2_e
 
        call supfor(vx,vy,vz,tmp,k2,force,m)

         ox = wy
         oy = wz
         oz = wx

        time = time + dt

        time2 = gettime()

       if(m.eq.0)then
          write(6,fmt='(5X,A,F20.4,A)') 
     &   'elapsed time =',1.0D3*(time2-time1),' msec'
        endif

      endif

      if(istep.le.nstep) goto 1

      if(m.eq.0)then
        close(11)
        close(20)
        close(22)
      end if

      stop

 100  format(' step',i4,':  t = ',e14.6,' E = ',e14.6,' W = ',e14.6)
 200  format(8e13.5)
 102  format(2(1x,E16.7))

      end
c
      subroutine supfor(vx,vy,vz,tmp,k2,force,m)

      complex*16,dimension(lz,ly,lxp2)::vx,vy,vz
      real*8,dimension(lz,ly,lxp2)::k2,tmp
      real*8 force(2),ek,ff
      integer*4 lx,ly,lz,lxp,lzp,lxp2,lzp2
      real*8 t0,t1,t2
  
      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2
      external d_vadd

       nallgrp = 0
       call mp_sync(nallgrp)

      tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      
      do is = 1, 2

         ek = 0.D0
         ff = 0.D0
         ek=sum(tmp,mask=(abs(dsqrt(k2)-is-0.49999).lt.0.5))

          call mp_reduce(ek,ff,8,0,d_vadd,nallgrp)
          call mp_bcast(ff,8,0,nallgrp)

         ff = dsqrt(force(is)/ff)

        where(abs(dsqrt(k2)-is-0.499999).lt.0.5)
          vx=vx*ff
          vy=vy*ff
          vz=vz*ff
         end where

      end do

      return
      end

c
      subroutine wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,m)

      real*8,dimension(lz,ly,lxp2)::kx,ky,kz,k2,k2_e
      integer*4 lx,ly,lz,lxp,lzp,lxp2,lzp2
      real*8 dt,rnu

      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2
 
       do i=1,ly
       ky(:,i,:) = mod(i-1+ly/2,ly)-ly/2
       end do
   
       do i=1,lz
       kz(i,:,:) = mod(i-1+lz/2,lz)-lz/2
       end do
 
       do i=1,lxp2
       kx(:,:,i) = m*lxp2+i-1
       end do

      k2=kx*kx+ky*ky+kz*kz
      k2_e = dexp(-k2*dt*rnu)
 
      if(m.eq.0)then
       k2(1,1,1)=1.e-5
       k2_e(1,1,1) = 1.D0
      end if

      return
      end
c

      subroutine gaussian(u,pi)
      complex*16,dimension(lz,ly,lxp2) :: u
      integer*4 lx,ly,lz,lxp,lzp,lxp2,lzp2
      real*8 pi,t1,t2
      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2

      u  = 0.

      do i=1,lz
       do j=1,ly
        do k=1,lxp2
         t1 = rand()
         t2 = rand()
         if(t1.le.1.e-20) t1 = 1.e-20
         if(t2.le.1.e-20) t2 = 1.e-20
         t2 = 2*pi*t2
         u(i,j,k)=dsqrt(-2.0D0*dlog(t1))*cmplx(dcos(t2),dsin(t2))
        end do
       end do
      end do
  
      return
      end


      subroutine symmetrize(c,m)
      complex*16,dimension(lz,ly,lxp2)::c
      integer*4 lx,ly,lz,lxp,lzp,lxp2,lzp2

      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2

      c(:,ly/2+1,:) = (0.D0,0.D0)
      c(lz/2+1,:,:) = (0.D0,0.D0)

      if(m.eq.0) then
       c(1,1,1) = (0.D0,0.D0)
       
        do j = ly/2+1,ly
         c(:,j,1) = (0.D0,0.D0)
        end do
 
       do i = lz/2+1,lz
        c(i,1,1) = (0.D0,0.D0)
       end do

      end if

      return
      end

      subroutine dealiasing(vx,k2)
         
c   2/3 rule     

      complex*16,dimension(lz,ly,lxp2)::vx
      real*8,dimension(lz,ly,lxp2)::k2
      integer*4 lx,ly,lz,lxp,lzp,lxp2,lzp2
      real*8 ktr
  
      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2

      ktr = lx/3.D0 

      where(dsqrt(k2).ge.ktr)vx = (0.D0,0.D0)
      
      return
      end


      subroutine esvemonp
      return
      end








------- End of Forwarded Message
