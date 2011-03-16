******************************************
      program moment
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'mpif.h'

      include 'param.escan_real'
      integer status(MPI_STATUS_SIZE)
******************************************
       real*8 AL(3,3)
********************************************
       character*20 fwg_in,fkpt_in
*************************************************


       complex*16,allocatable,dimension(:) :: ugtemp
       real*8,allocatable,dimension(:,:) :: pmatrix

       complex*16 cc,ccx,ccy,ccz

**************************************************
c      atime00=mclock()/100.d0
c
c initialize mpi and number each node
c
       call mpi_init(ierr)
       call mpi_comm_rank(MPI_COMM_WORLD,inode,ierr)
       call mpi_comm_size(MPI_COMM_WORLD,nnodes,ierr)
       inode=inode+1

       pi=4.0d0*datan(1.0d0)

       open(10,file="moment.input",status='old',
     &  action='read',iostat=ierr)
       if(ierr.ne.0) then
       if(inode.eq.1) 
     &  write(6,*) "file ***moment.input*** does not exist, stop"
       stop
       endif
      
       rewind(10)
       read(10,*) fwg_in
       read(10,*) fkpt_in
       close(10)

*******************************************************

       if(inode==nnodes) then
       open(11,file=fwg_in,form='unformatted',status='old',
     &  action='read',iostat=ierr)
       if(ierr.ne.0) then
       write(6,*) "file ***",fwg_in,"*** does not exist, stop"
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif
       
       rewind(11)
       read(11)n1,n2,n3,mx
       read(11)Ecut
       read(11)AL
       read(11)nnodes_o

       if(nnodes_o.ne.nnodes) then
         write(6,*) "nnodes_o changed, stop", nnodes_o
         call mpi_abort(MPI_COMM_WORLD,ierr)
       endif
       endif   ! out inode.eq.nnodes

       call mpi_bcast(n1,1,MPI_INTEGER,nnodes-1,MPI_COMM_WORLD,ierr)
       call mpi_bcast(n2,1,MPI_INTEGER,nnodes-1,MPI_COMM_WORLD,ierr)
       call mpi_bcast(n3,1,MPI_INTEGER,nnodes-1,MPI_COMM_WORLD,ierr)
       call mpi_bcast(mx,1,MPI_INTEGER,nnodes-1,MPI_COMM_WORLD,ierr)
       call mpi_bcast(Ecut,1,MPI_REAL8,nnodes-1,MPI_COMM_WORLD,ierr)
       call mpi_bcast(AL,9,MPI_REAL8,nnodes-1,MPI_COMM_WORLD,ierr)

       call mpi_barrier(MPI_COMM_WORLD,ierr)

       allocate(pmatrix(mx,mx))

**********************************
       nr1x=n1
       nr2x=n2
       nr3x=n3+2
       nr=n1*n2*n3

       mr=n1*n2*(n3+2)

       nr_n = nr/nnodes

       mr_n=mr/nnodes

       mr=1
       mr_n=1
       nr_n=1  
       nr=1

c     Calculate the approx. no. of g points, so we can dimension arrays
      vol=al(3,1)*(al(1,2)*al(2,3)-al(1,3)*al(2,2))
     &     +al(3,2)*(al(1,3)*al(2,1)-al(1,1)*al(2,3))
     &     +al(3,3)*(al(1,1)*al(2,2)-al(1,2)*al(2,1))
      vol=dabs(vol)
      delta_k=(2*pi)**3/vol
      totg=(4.0d0/3.0d0*pi*(dsqrt(2.0d0*Ecut))**3)/delta_k
      mg_nx=int(1.1d0*totg/nnodes)+100

       call get_ALI(AL,ALI)


       open(12,file=fkpt_in,status='old',action='read',iostat=ierr)
       if(ierr.ne.0) then
       if(inode.eq.1) 
     & write(6,*) "file ***",fkpt_in,"*** does not exist, stop"
       stop
       endif

       rewind(12)
       read(12,*) nkpt

       call data_allocate_akx(nkpt)
       call ngtot_allocate(nnodes,nkpt)

       read(12,*) iflag, ALxyz
       if(iflag.eq.1)  then
         if(inode.eq.1) write(6,*) "input kpts in Cartesian Coord"
         do kpt=1,nkpt
         read(12,*) akx(kpt),aky(kpt),akz(kpt),weighkpt(kpt)
         akx(kpt)=akx(kpt)*2*pi/ALxyz
         aky(kpt)=aky(kpt)*2*pi/ALxyz
         akz(kpt)=akz(kpt)*2*pi/ALxyz
         enddo
      endif

      if(iflag.eq.2) then
         if(inode.eq.1) write(6,*) "input kpts in primary cell unit"

         do kpt=1,nkpt
         read(12,*) ak1_t,ak2_t,ak3_t,weighkpt(kpt)
         akx(kpt)=2*pi*(ALI(1,1)*ak1_t+ALI(1,2)*ak2_t+
     &           ALI(1,3)*ak3_t)
         aky(kpt)=2*pi*(ALI(2,1)*ak1_t+ALI(2,2)*ak2_t+
     &           ALI(2,3)*ak3_t)
         akz(kpt)=2*pi*(ALI(3,1)*ak1_t+ALI(3,2)*ak2_t+
     &           ALI(3,3)*ak3_t)
         enddo
       endif

       close(12)
************************************************
       ilocal=1
       call fft_allocate(n1,n2,n3,nnodes)
       call load_allocate(ncolx,nnodes,mg_nx,mr_n)
       islda=1
       call data_allocate(mg_nx,mx,mr_n,ilocal,inode,nkpt,
     &   islda,1)
************************************************

        if(inode.eq.1) then
        open(15,file="pmatrix")
        rewind(15)
        endif

       do 100 kpt=1,nkpt
       call gen_G_comp(kpt,0) 

       ng_n=ngtotnod(inode,kpt)

       call get_ug()
**************************************************

         do 200 is1=1,mx
         do 200 is2=1,is1

         ccx=dcmplx(0.d0,0.d0)
         ccy=dcmplx(0.d0,0.d0)
         ccz=dcmplx(0.d0,0.d0)

         do i=1,ng_n
         cc=ug_n(i,is1)*dconjg(ug_n(i,is2))
         ccx=ccx+cc*gkx_n(i,kpt)
         ccy=ccy+cc*gky_n(i,kpt)
         ccz=ccz+cc*gkz_n(i,kpt)
         enddo

         call mpi_barrier(MPI_COMM_WORLD,ierr)

         call global_sumc(ccx)
         call global_sumc(ccy)
         call global_sumc(ccz)
        
         pmatrix(is1,is2)=cdabs(ccx*vol)**2+cdabs(ccy*vol)**2+
     &                    cdabs(ccz*vol)**2
         pmatrix(is2,is1)=pmatrix(is1,is2)

200      continue

         if(inode.eq.1) then
         write(15,*) "kpt,akx,y,z (Cartesian Coord,a.u),mx "
         write(15,330) kpt,akx(kpt),aky(kpt),akz(kpt),mx
330      format(i4,2x,3(f12.6,2x),i4) 
         write(15,*) "**** <|P_x|>^2+<|P_y|>^2+<|P_z|>^2 (a.u) ***"
         do is1=1,mx
         write(15,331) (pmatrix(is1,is2),is2=1,mx)
         enddo
331      format(10(E8.3,1x))
         write(15,*) "********************************************"
         endif

100    continue      ! do kpt
    


       if(inode==nnodes) close(11)
       if(inode==1) close(15)
*********************************************************
      
      call data_deallocate()
      call load_deallocate()
c      call fft_deallocate()   ! some deallocated arrays are allocated in fftprep

      call mpi_finalize(ierr)

      stop

      contains
      
**************************************************
      subroutine get_ug()
      implicit double precision (a-h,o-z)

       if(inode==nnodes) then
            allocate(ugtemp(mg_nx))

            do i=0,nnodes-2
               do iwavefun=1,mx
                  read(11)ugtemp    ! assuming the mg_nx is the same as before
                  do j=1,mg_nx
                  ug_n(j,iwavefun)=ugtemp(j)
                  enddo
               end do


          call mpi_send(ug_n,mg_nx*mx,MPI_DOUBLE_COMPLEX,i,
     &            102,MPI_COMM_WORLD,ierr)
            end do
c     Now read in the data for node nnodes
          do iwavefun=1,mx
               read(11)ugtemp
               do j=1,mg_nx
               ug_n(j,iwavefun)=ugtemp(j)
               enddo
            end do
            deallocate(ugtemp)

         else   ! for other nodes

         call mpi_recv(ug_n,mg_nx*mx,MPI_DOUBLE_COMPLEX,nnodes-1,
     &       102,MPI_COMM_WORLD,status,ierr)

         end if
**************************************************
         call mpi_barrier(MPI_COMM_WORLD,ierr)
         return
         end subroutine get_ug

         end
