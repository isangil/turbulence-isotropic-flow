        include 'mpif.h'

        integer nx,ny,nz,lx,ly,lz,lx2,lx1,lly
        integer nallgrp,ierror,nproc,nprocs
        integer isign,idp,id,lz1,ly1
        real ek,e_t,scale,pi
        common /dim/nx,ny,nz,lx,ly,lz,lx2,lx1,lly,nproc

