c       This is the spectral ccde for three-D isotropic turbulence
c       using the message-passing in SP machines. A random
c       forcing is used, and a scalar field is advected.
c       March, 25th, 1997.
c       lx = ly/2 ; ly = lz = 2**n ( n = 4,5,...)
c       lx = x-diretion mesh size in each node
c       ly = y-diretion mesh size in each node
c       lz = z-direction mesh size in each node
c       nproc = number of process (nodes)
c       

      program pass_evlt 
      complex*16,allocatable,dimension(:,:,:)::vx,vy,vz
      complex*16,allocatable,dimension(:,:,:)::wx,wy,wz
      complex*16,allocatable,dimension(:,:,:)::ox,oy,oz
      complex*16,allocatable,dimension(:,:,:)::p_x,p_y,p_z
      complex*16,allocatable,dimension(:,:,:)::pass,pass0
      real*8,allocatable,dimension(:,:,:)::kx,ky,kz,k2
      real*8,allocatable,dimension(:,:,:)::tmp,tmp1,k2_e
      real*8,allocatable,dimension(:)::seed 
      integer*4 isign,icontxt,ip(40),nstep,lx,ly,lz,nnodes
      integer*4 lxp,lzp,lxp2,lzp2,nproc,ifijk(200,3,2)
      integer*4 ieout,itout,nek,nvfor,npfor,ifcnt(2)
      real*8 t_ex,t_ey,t_ez,t_wx,t_wy,t_wz,kd
      real*8 distp(100),distv(100),hdu(100),hdp(100)
      real*8 to_ex,to_ey,to_ez,to_wx,to_wy,to_wz
      real*8 scale,scale_r,scale_s,eru(10),erp(10)
      real*8 tout,eout,dt,rnu,u0,time,sd,aa,rr,eku,ekp,e_tu,e_tp
      real*8 force_p(2),force_u(2),pi,cc1,akapa
      real*8 time0,time1,time2,time3,time4,time5,time6

c     ic would be to input cases: pass-new,both-new, old & so on.
      character*40 forcing,fdata(5),filename
c      character fdatan(10)
      logical*4 new,timing

c      data   fdata /'/database/syc/E1.sp',
c     1 '/database/syc/E2.sp','/database/syc/E3.sp',
c     2 '/database/syc/E4.sp','/database/syc/E5.sp',
c     3 '/database/syc/E6.sp','/database/syc/E7.sp',
c     4 '/database/syc/E8.sp','/database/syc/E9.sp',
c     5 '/database/syc/E10.sp'/

c      watson

      data   fdata /'/tmp/H1.tmpsp',
     1 '/tmp/H2.tmpsp','/tmp/H3.tmpsp',
     2 '/tmp/H4.tmpsp','/tmp/H5.tmpsp'/

c      cornell

c      data   fdata /'/tmp/scratch/H1.sp',
c     1 '/tmp/scratch/H2.sp','/tmp/scratch/H3.sp',
c     2 '/tmp/scratch/H4.sp','/tmp/scratch/H5.sp',
c     3 '/tmp/scratch/H6.sp','/tmp/scratch/H7.sp',
c     4 '/tmp/scratch/H8.sp','/tmp/scratch/H9.sp',
c     5 '/tmp/scratch/H10.sp'/
 

      data istep, is, ak0, akp / -1, -1, 4.75683, 4.75683 /
 
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
 

c==========================
c     piofs. 
c
c           `parallel I/O file system' at watson
c
c     bound to number of procs < 100. trivial extension.
    
c       do i=1,9
 
c         if (m.lt.10) then
         
c           fdatan(i)='/piofs/syc/H'//char(i+48)//'.0'
c     1//char(m+48)//'.spd'   
c           write(*,2) fdatan(i)  
c         else

c           m1=int(m/10)
c           m2=mod(m,10)
c           fdatan(i)='/piofs/syc/H'//char(i+48)//'.'
c     1//char(m1+48)//char(m2+48)//'.spd'
c            write (*,2)fdatan(i)
c         endif

c       enddo
c 2     format(a40)

c       if (m.lt.10) then

c           fdatan(10)='/piofs/syc/H10.0'
c     1//char(m+48)//'.spd'
c           write(*,2) fdatan(10)
c       else

c           m1=int(m/10)
c           m2=mod(m,10)
c           fdatan(10)='/piofs/syc/H10.'
c     1//char(m1+48)//char(m2+48)//'.spd'
c            write (*,2)fdatan(10)
c       endif

c    
c=============================== 
      
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
      allocate( p_z(lz,ly,lxp2) )
      allocate( p_y(lz,ly,lxp2) )
      allocate( p_x(lz,ly,lxp2) )
      allocate( pass(lz,ly,lxp2) )
      allocate( pass0(lz,ly,lxp2) )

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
c
       open(1,file='parameter.d',status='old')
       read(1,*) nstep,tout,eout,dt,rnu,akapa,u0,ap0,time,nvfor
     1    ,npfor,new,timing,sd,force_u,force_p  
        forcing='random'

       write(*,*)'nste= ',nstep,' to= ',tout,' eo= ',eout,' dt= ',dt
       write(*,*)'nu= ',rnu,' akp= ',akapa,' u0= ',u0,' ap0= ',ap0
       write(*,*)'ti= ',time,' nv= ',nvfor,' np= ',npfor,' new ',new
        write(*,*)'timing= ',timing,' forcin= ',forcing
        write(*,*)'f_p= ',force_p,' f_u= ',force_u,' sd= ',sd
       close(1)
c
c      Commented IF using piofs Parallel I/O File System.
c      if not )/tmp, or whatever, use this feature
c
c      next will contain filename where all data is enclosed
c      (if old)
       open(8,file='input.file',status='old')
       write(*,*)'enter filename'
       read(8,1005)filename
       write(*,*)filename
       close(8)
      end if

1005   format(a)

c       bcast param
c
        call mp_bcast(nstep,4,0,nallgrp)
        call mp_bcast(tout,8,0,nallgrp)
        call mp_bcast(eout,8,0,nallgrp)
        call mp_bcast(dt,8,0,nallgrp)
        call mp_bcast(rnu,8,0,nallgrp)
        call mp_bcast(akapa,8,0,nallgrp)
        call mp_bcast(u0,8,0,nallgrp)
        call mp_bcast(ap0,8,0,nallgrp)
        call mp_bcast(time,8,0,nallgrp)
        call mp_bcast(new,4,0,nallgrp)
        call mp_bcast(timing,4,0,nallgrp)
        call mp_bcast(force_p,16,0,nallgrp)
        call mp_bcast(force_u,16,0,nallgrp)
        call mp_bcast(akapa,8,0,nallgrp)
c NOTc        now using piofs.
        call mp_bcast(filename,40,0,nallgrp)
        call mp_bcast(forcing,40,0,nallgrp)
        call mp_bcast(nvfor,4,0,nallgrp)
        call mp_bcast(npfor,4,0,nallgrp)
        write(*,*)'m=',m,filename
        call mp_bcast(sd,8,0,nallgrp)

        ieout = int(eout/dt)
        itout = int(tout/dt)
c
c       will contain up to which modes are 'useful'
c
        nek=0.9*sqrt(lx*lx+.25*ly*ly+.25*lz*lz)
c
c       files for data processing        
        if(m.eq.0)then
         write(*,*)'ieo= ',ieout,' itou= ',itout,' nek ',nek
c         open(11,file='energy.d',status='unknown')
         open(20,file='spectrum_u.d',status='unknown') 
         open(21,file='spectrum_p.d',status='unknown')
         open(22,file='transfer.d',status='unknown') 
         open(70,file='initialu.sp',status='unknown')
         open(80,file='initialp.sp',status='unknown')
        end if
        
c---  Generate initial condition -------

      if(new) then

c---     Generate random number seed for each process ---

         if (m.eq.0)then
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
         call gaussian(pass,pi)

         call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,m)
c
c        in this version, i apply an homogeneous forcing.
c        supset is overriden.
c   
c         call supset(k2,ifijk,ifcnt,lx,ly,lxp2)
c        proyector.

         tmp = (kx*real(vx) + ky*real(vy) + kz*real(vz))/k2
         vx = dcmplx(real(vx) - kx*tmp, aimag(vx))
         vy = dcmplx(real(vy) - ky*tmp, aimag(vy))
         vz = dcmplx(real(vz) - kz*tmp, aimag(vz))
         tmp = (kx*aimag(vx) + ky*aimag(vy) + kz*aimag(vz))/k2
         vx = dcmplx(real(vx), aimag(vx) - kx*tmp)
         vy = dcmplx(real(vy), aimag(vy) - ky*tmp)
         vz = dcmplx(real(vz), aimag(vz) - kz*tmp)
c        initial spectra. later, include some 'if' for 'ic'
c        cases. now "both - new"

c          velocity:

          cc1 = u0 * dsqrt(8.*sqrt(2./pi)/(3.*pi*ak0**5))
          tmp = cc1 * dsqrt(k2) * dexp (-k2/(ak0*ak0))

          vx = vx * tmp
          vy = vy * tmp
          vz = vz * tmp
c
c          passive scalar:

          cc1 = ap0 * dsqrt(8.*dsqrt(2./pi)/(3.*pi*akp**5))
          k2_e = cc1 * dsqrt(k2) * dexp(-k2/(akp*akp))

          pass = pass * k2_e
c
c           symmetryze all
         call symmetrize(vx,m)
         call symmetrize(vy,m)
         call symmetrize(vz,m)
         call symmetrize(pass,m)

c         call dealiasing(vx,k2)
c         call dealiasing(vy,k2)
c         call dealiasing(vz,k2)
        
c        calculating initial spectrum
c
c        velocity:

         tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
         do i=1,nek
           eku = 0.D0
           e_tu = 0.D0
           eku=sum(tmp,mask=(abs(dsqrt(k2)-i-0.49999).lt.0.5))
           call mp_reduce(eku,e_tu,8,0,d_vadd,nallgrp)
           if(m.eq.0)then
              write(70,*)i,e_tu
           end if
         end do
c        
c         passive scalar:

         tmp = pass*conjg(pass)
         do i=1,nek
           ekp = 0.D0
           e_tp = 0.D0
           ekp=sum(tmp,mask=(abs(dsqrt(k2)-i-0.49999).lt.0.5))
           call mp_reduce(ekp,e_tp,8,0,d_vadd,nallgrp)
           if(m.eq.0)then
              write(80,*)i,e_tp
           end if
         end do

         if(m.eq.0)close(70)
         if(m.eq.0)close(80)
      else               
c                
c          not new
c ------ reading from DISK ------ 
c
c             usin piofs. read last file index.
c           if (m.eq.0) then
c               write(*,*) 'what was the index of the piofs file'
c               open(3,file='input.indexpiofs',status='unknown')
c               read (3,*) indf
c               write(*,*) indf,' filenamemaster= ',fdatan(indf)
c               close(3)
c           endif
c           call mp_bcast(indf,4,0,nallgrp)
c         open(2,file=fdatan(indf),status='old',

         open(2,file=filename,status='old',
     1    form='unformatted')
         read(2)vx
         read(2)vy
         read(2)vz
         read(2)pass
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
      pass0=pass

      dt_h=.5*dt
      ip = 0

 1    istep = istep + 1

c     write out the total energy in K space (fort.76)

c      add fort.77 for passive sc. spect      
      time1 = gettime() 
      tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      eku = sum(tmp)
      call mp_reduce(eku,e_tu,8,0,d_vadd,nallgrp)
      if(m.eq.0)then
         write(*,*)istep,time,' e_u= ',e_tu
         write(76,*)time,e_tu
         if (e_tu.gt.1000) then
             write(*,*) 'blow!'
             istep=nstep  
             call mp_bcast(istep,4,0,nallgrp)
         endif
      end if
c
         tmp=pass*conjg(pass)
        ekp=sum(tmp)
        call mp_reduce(ekp,e_tp,8,0,d_vadd,nallgrp)
        if(m.eq.0)then
          write(*,*)istep,time,' e_p= ',e_tp
          write(77,*)time,e_tp
          if (e_tp.gt.1000) then
              write(*,*) 'blow ep !'
              istep=nstep  
              call mp_bcast(istep,4,0,nallgrp)
          end if

        end if
c        
      
      if(mod(istep,itout).eq.0.and.istep.ne.0) then
c
c 
          open(2,file=fdata(iii),status='unknown',
     1    form='unformatted')
          write(2)vx
          write(2)vy
          write(2)vz
          write(2)pass
          nallgrp = 0
          call mp_sync(nallgrp)

          close(2)

          iii = iii + 1

      endif
c
c      rot u =>  ik x v
c
      wx = (0.,1.) * (ky*vz - kz*vy)
      wy = (0.,1.) * (kz*vx - kx*vz)
      wz = (0.,1.) * (kx*vy - ky*vx)

c     energy spectrum calculation 
c      velocity:     
c     
      if(mod(istep,ieout).eq.0) then

        tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
        if(m.eq.0)then   
          write(20,*)'#k ',(kkk-1)
        end if
        do i=1,nek
            eku = 0.D0
            e_tu = 0.D0
            eku=sum(tmp,mask=(abs(dsqrt(k2)-i-0.49999).lt.0.5))
            call mp_reduce(eku,e_tu,8,0,d_vadd,nallgrp)
            if(m.eq.0)then
              write(20,*)i,e_tu
            end if
        end do
c
c       passive scalar:
c
        tmp = pass*conjg(pass)
        if(m.eq.0)then
           write(21,*) '#k ',(kkk-1)
        end if
        do i=1,nek
            ekp = 0.D0
            e_tp = 0.D0
            ekp=sum(tmp,mask=(abs(dsqrt(k2)-i-0.49999).lt.0.5))
            call mp_reduce(ekp,e_tp,8,0,d_vadd,nallgrp)
            if(m.eq.0)then
               write(21,*)i,e_tp
            end if
        end do
   
        kkk = kkk + 1
      end if
c     
c     save values for v      
     
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

      call pdcrft3(pass,pass,lx,ly,lz,isign,scale_r,icontxt,ip)

c     here is where all statistics should be entered.
c
c    for now, i'll leave things as they are, once the 
c    code runs, i'll get the statistics part coded.
c       
      tmp = real(vy)*real(wz) - real(vz)*real(wy)
      wz = dcmplx(real(vz)*real(wx) - real(vx)*real(wz), aimag(wz))
      wx = dcmplx(real(vx)*real(wy) - real(vy)*real(wx), aimag(wx))
      wy = dcmplx(tmp, aimag(wy))
      tmp = aimag(vy)*aimag(wz) - aimag(vz)*aimag(wy)
      wz = dcmplx(real(wz), aimag(vz)*aimag(wx) - aimag(vx)*aimag(wz))
      wx = dcmplx(real(wx), aimag(vx)*aimag(wy) - aimag(vy)*aimag(wx))
      wy = dcmplx(real(wy), tmp)

      p_x = dcmplx(real(pass)*real(vx),aimag(pass)*aimag(vx))
      p_y = dcmplx(real(pass)*real(vy),aimag(pass)*aimag(vy))
      p_z = dcmplx(real(pass)*real(vz),aimag(pass)*aimag(vz))

      vx = dcmplx(kx,ky)
      vy = dcmplx(kz,k2)
      vz = dcmplx(k2_e,tmp1)
      call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,m)

      isign = 1

      ip(1) = 0

c      wx = wx*scale_s
c      wy = wy*scale_s
c      wz = wz*scale_s
       call pdrcft3(wx,wx,lx,ly,lz,isign,scale_r,icontxt,ip)
       call pdrcft3(wy,wy,lx,ly,lz,isign,scale_r,icontxt,ip)
       call pdrcft3(wz,wz,lx,ly,lz,isign,scale_r,icontxt,ip)
       call pdrcft3(pass,pass,lx,ly,lz,isign,scale_r,icontxt,ip)
       call pdrcft3(p_x,p_x,lx,ly,lz,isign,scale_r,icontxt,ip)
       call pdrcft3(p_y,p_y,lx,ly,lz,isign,scale_r,icontxt,ip)
       call pdrcft3(p_z,p_z,lx,ly,lz,isign,scale_r,icontxt,ip)

       wx = wx*scale
       wy = wy*scale
       wz = wz*scale
       pass=pass*scale
       p_z=p_z*scale
       p_y=p_y*scale
       p_x=p_x*scale

       call symmetrize(vx,m)
       call symmetrize(vy,m)
       call symmetrize(vz,m)
       call symmetrize(wx,m)
       call symmetrize(wy,m)
       call symmetrize(wz,m)
       call symmetrize(p_x,m)
       call symmetrize(p_y,m)
       call symmetrize(p_z,m)
       call symmetrize(pass,m)

c       call dealiasing(vx,k2)
c       call dealiasing(vy,k2)
c       call dealiasing(vz,k2)

       p_x = (0.,-1.)*(kx*p_x + ky*p_y + kz*p_z)
c
      tmp = (kx*real(wy) + ky*real(wz) + kz*real(wx))/k2
      wy = dcmplx(real(wy) - kx*tmp, aimag(wy))
      wz = dcmplx(real(wz) - ky*tmp, aimag(wz))
      wx = dcmplx(real(wx) - kz*tmp, aimag(wx))
      tmp = (kx*aimag(wy) + ky*aimag(wz) + kz*aimag(wx))/k2
      wy = dcmplx(real(wy), aimag(wy) - kx*tmp)
      wz = dcmplx(real(wz), aimag(wz) - ky*tmp)
      wx = dcmplx(real(wx), aimag(wx) - kz*tmp)
c
c     energy transfer..
c
       if(mod(istep,ieout).eq.0) then

c           p-proyection

            tmp = (kx*real(vx) + ky*real(vy) + kz*real(vz))/k2
            vx = dcmplx(real(vx) - kx*tmp, aimag(vx))
            vy = dcmplx(real(vy) - ky*tmp, aimag(vy))
            vz = dcmplx(real(vz) - kz*tmp, aimag(vz))
            tmp = (kx*aimag(vx) + ky*aimag(vy) + kz*aimag(vz))/k2
            vx = dcmplx(real(vx), aimag(vx) - kx*tmp)
            vy = dcmplx(real(vy), aimag(vy) - ky*tmp)
            vz = dcmplx(real(vz), aimag(vz) - kz*tmp)

            tmp = real(vx)*real(wy) + real(vy)*real(wz) + 
     &      real(vz)*real(wx) +
     &      aimag(vx)*aimag(wy) + aimag(vy)*aimag(wz) + 
     &      aimag(vz)*aimag(wx)
             
            if(m.eq.0)then
               write(22,*)(jjj-1)
            end if
 
            do i=1,nek
              eku = 0.
              e_tu = 0. 
              eku=sum(tmp,mask=(abs(dsqrt(k2)-i-0.49999).lt.0.5))
              call mp_reduce(eku,e_tu,8,0,d_vadd,nallgrp)
              if(m.eq.0)then
                 write(22,*) i,eku
              end if
            end do

            jjj = jjj + 1

         endif
c
c     time step.

      if(istep.eq.0) then

         vx = vx * (1.-rnu*dt_h*k2) + dt_h*wy
         vy = vy * (1.-rnu*dt_h*k2) + dt_h*wz
         vz = vz * (1.-rnu*dt_h*k2) + dt_h*wx
         pass = pass*(1.-akapa*dt_h*k2) + dt_h*p_x

          if(m.eq.0) write(*,*) 'b calling forcing ' 
        call supfor(vx,vy,vz,tmp,k2,force_u,m)
        call psvfor(pass,tmp,k2,force_p,m)

c       remember we are trying the homoge. forcing.
 
c         call supfor(vx,vy,vz,tmp,k2,force_u,m,ifcnt,ifijk,nvfor)

c         call psvfor(pass,tmp,k2,force_p,m,ifcnt,ifijk,forcing,npfor)
         time2 = gettime()

         if(m.eq.0)then
            write(6,fmt='(5X,A,F20.4,A)')
     &     'elapsed time =',1.0D3*(time2-time1),' msec'
         endif
c
c        watson

        open(2,file='/tmp/scratmp',status='unknown',
     1   form='unformatted')

c          cornell
c

c         open(2,file='/tmp/scratch/scra',status='unknown',
c     1   form='unformatted')
            write(2)wy
            write(2)wz
            write(2)wx
            write(2)p_x
         close(2)

         time = time + dt_h

      elseif(istep.eq.1) then
 
         vx = ox + dt*wy - rnu*dt*k2*vx
         vy = oy + dt*wz - rnu*dt*k2*vy
         vz = oz + dt*wx - rnu*dt*k2*vz
         pass = pass0 +dt*p_x - akapa*dt*k2*pass

         call supfor(vx,vy,vz,tmp,k2,force_u,m)
         call psvfor(pass,tmp,k2,force_p,m)
         
c         call supfor(vx,vy,vz,tmp,k2,force_u,m,ifcnt,ifijk,nvfor)
c         call psvfor(pass,tmp,k2,force_p,m,ifcnt,ifijk,forcing,npfor)

         time2 = gettime()

         if(m.eq.0)then
            write(6,fmt='(5X,A,F20.4,A)')
     &      'elapsed time =',1.0D3*(time2-time1),' msec'
         endif

c        open(2,file='/database/syc/scratch',status='old',
c     1   form='unformatted')
c    
c        watson

         open(2,file='/tmp/scratmp',status='unknown',
     1   form='unformatted')

c         cornell

c         open(2,file='/tmp/scratch/scra',status='unknown',
c     1   form='unformatted')
           read(2)ox
           read(2)oy
           read(2)oz
           read(2)pass0
         close(2)

         time = time + dt_h

      else

         vx = vx + dt_h*(3.*wy - k2_e*ox)
         vy = vy + dt_h*(3.*wz - k2_e*oy)
         vz = vz + dt_h*(3.*wx - k2_e*oz)

         vx = vx*k2_e
         vy = vy*k2_e
         vz = vz*k2_e
 
        call supfor(vx,vy,vz,tmp,k2,force_u,m)
        call psvfor(pass,tmp,k2,force_p,m)
        
c         call supfor(vx,vy,vz,tmp,k2,force_u,m,ifcnt,ifijk,nvfor)
c         call psvfor(pass,tmp,k2,force_p,m,ifcnt,ifijk,forcing,npfor)

         call wavenumber(kx,ky,kz,k2,k2_e,dt,rnu,m)

         pass = pass + dt_h*(3*p_x - k2_e*pass0)
         pass = pass * k2_e

         ox = wy
         oy = wz
         oz = wx
         pass0 = p_x  
        
        time = time + dt

        time2 = gettime()

        if(m.eq.0)then
          write(6,fmt='(5X,A,F20.4,A)') 
     &   'elapsed time =',1.0D3*(time2-time1),' msec'
        endif

      endif

      if(istep.le.nstep) goto 1

      if(m.eq.0)then
c        close(11)
        close(20)
        close(22)
      end if

      stop

 100  format(' step',i4,':  t = ',e14.6,' E = ',e14.6,' W = ',e14.6)
 200  format(8e13.5)
 102  format(2(1x,E16.7))

      end

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
      subroutine psvfor(pass,tmp,k2,force,m)

      complex*16,dimension(lz,ly,lxp2)::pass
      real*8,dimension(lz,ly,lxp2)::k2,tmp
      real*8 force(2),ek,ff
      integer*4 lx,ly,lz,lxp,lzp,lxp2,lzp2
      real*8 t0,t1,t2
  
      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2
      external d_vadd

       nallgrp = 0
       call mp_sync(nallgrp)

      tmp = pass*conjg(pass)
      
      do is = 1, 2

         ek = 0.D0
         ff = 0.D0
         ek=sum(tmp,mask=(abs(dsqrt(k2)-is-0.49999).lt.0.5))

          call mp_reduce(ek,ff,8,0,d_vadd,nallgrp)
          call mp_bcast(ff,8,0,nallgrp)

         ff = dsqrt(force(is)/ff)

        where(abs(dsqrt(k2)-is-0.499999).lt.0.5)
          pass=pass*ff
         end where

      end do

      return
      end


c
c      subroutine supset(k2,ifijk,ifcnt,nx,ny,nz)
c      integer  ifcnt(2), ifijk(200,3,2)
c      real*8,dimension(nx,ny,nz)::k2
c
c      eps = 1.0e-6
c      ifcnt(1) = 0
c      ifcnt(2) = 0
c      do i=1,3
c         do k=1,nz
c            if(k.gt.3.and.k.lt.nz-1)goto 2
c            do j=1,ny
c               if(j.gt.3.and.j.lt.ny-1)goto 1
c               irk = int(sqrt(k2(i,j,k))+eps)
c               if ((irk.ge.3).or.(irk.le.0)) goto 1
c               ifcnt(irk) = ifcnt(irk) + 1
c              ifijk(ifcnt(irk),1,irk) = i
c              ifijk(ifcnt(irk),2,irk) = j
c               ifijk(ifcnt(irk),3,irk) = k
c 1             continue
c            enddo
c 2          continue
c         enddo
c      enddo

c      return
c      end

cc    modified by Inigo San Gil, march 3rd 97
cc    parameters added: ifcnt, ifijk,nvfor
c      subroutine supfor(vx,vy,vz,tmp,k2,force,m,ifcnt,ifijk,nvfor)

c      complex*16,dimension(lz,ly,lxp2)::vx,vy,vz
c      real*8,dimension(lz,ly,lxp2)::k2,tmp
c      real*8 force(2),ek,ff
c      integer*4 lx,ly,lz,lxp,lzp,lxp2,lzp2
c      integer*4 nvfor,ifcnt(2),ifijk(200,3,2),nproc,m
c      real*8 t0,t1,t2
c 
c      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2
c      external d_vaddc
c
c       nallgrp = 0
c       call mp_sync(nallgrp)
c      
c      do is =1, nvforccc
c
c        tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
c      
c         ek = 0.D0
c         ff = 0.D0
c        ek=sum(tmp,mask=(abs(dsqrt(k2)-is-0.49999).lt.0.5))
c
c         call mp_reduce(ek,ff,8,0,d_vadd,nallgrp)
c         call mp_bcast(ff,8,0,nallgrp)
c         ff = dsqrt(force(is)/ff)
c 
c         do i = 1, ifcnt(is)
c           i1=ifijk(i,1,is)
c           i2=ifijk(i,2,is)
c           i3=ifijk(i,3,is)
c           where(abs(dsqrt(k2)-is-0.499999).lt.0.5)
c             vx=vx*ff
c             vy=vy*ff
c             vz=vz*ff
c           end where
c           vx(i1,i2,i1)=vx(i1,i2,i3)*ff
c           vy(i1,i2,i1)=vy(i1,i2,i3)*ff
c           vz(i1,i2,i1)=vz(i1,i2,i3)*ff
c         end do
c      enddo
c      return
c      end

c      subroutine added
c     still check validity of ifijk.

c      subroutine psvfor(pass,tmp,k2,force,m,ifcnt,ifijk,
c     . forcing,npfor,nproc,seed)cc

c      complex*16,dimension(lz,ly,lxp2)::pass
c      complex*16,dimension(5,10,10):: for_r
c      complex*16 fff,ppp 
c      real*8,dimension(nproc)::seed
c      real*8,dimension(lz,ly,lxp2)::k2,tmp
c      real*8 force(2),ek,ff
c      integer*4 lx,ly,lz,lxp,lzp,lxp2,lzp2,m
c      integer*4 npfor,nproc,ifcnt(2),ifijk(200,3,2)
c      real*8 t0,t1,t2
c      character*40 forcing
c      real*8 a(5,10,10)
c
c     it makes more sense (10,10,5)
c      common / dim / lx,ly,lz,lxp,lzp,lxp2,lzp2
c      external d_vadd
c
c      nallgrp = 0
c       call mp_sync(nallgrp)
      
c      if(forcing.eq.'random') then

c        aa = seed(m+1)
c        call srand(aa)
c        aa = 6.2831853*aa
c        for_r = dcmplx(cos(aa),sin(aa))
c      endif
 
c      tmp= pass*conjg(pass)

c      do is =1, npfor

c         ek = 0.D0
c         ff = 0.D0
c         ek=sum(tmp,mask=(abs(dsqrt(k2)-is-0.49999).lt.0.5))

c         call mp_reduce(ek,ff,8,0,d_vadd,nallgrp)
c         call mp_bcast(ff,8,0,nallgrp)

c         ff = dsqrt(force(i)/ff)
c         do i = 1,ifcnt(is)
c             i1=ifijk(i,1,is)
c             i2=ifijk(i,2,is)
c             i3=ifijk(i,3,is)
c             ppp=pass(i1,i2,i3)
c             if (forcing.eq.'ker') then
c                 fff=ppp/abs(ppp)
c             elseif (forcing.eq.'random') then
c                 fff=ppp*for_r(i1,i2,i3)*ff
c             elseif (forcing.eq.'shell') then
c                 fff=ppp*ff
c             else

c            if(m.eq.0) write(*,*) 'forcing:',forcing,' not recognized!'
c                 stop

c             endif 
c             pass(i1,i2,i3) = fff 
           
c       where(abs(dsqrt(k2)-is-0.499999).lt.0.5)
c            pass=pass*ff
c       end where   && before new routine, and for vel, only
         
c         enddo
c      enddo
      
c      if (forcing.eq.'ker') then
c            tmp=pass*conjg(pass)
c            do is=1,npfor 
c                ek=sum(tmp,mask=(abs(dsqrt(k2)-is-0.49999).lt.0.5))
c                call mp_reduce(ek,ff,8,0,d_vadd,nallgrp)
c                call mp_bcast(ff,8,0,nallgrp)
c                ff = dsqrt(force(i)/ff)
c                where(abs(dsqrt(k2)-is-0.499999).lt.0.5)
c                    pass=pass*ff
c                end where
c            enddo
c      endif
c      return
c      end

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
         u(i,j,k)=dsqrt(-2.0D0*dlog(t1))*dcmplx(dcos(t2),dsin(t2))
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






