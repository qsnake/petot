      program anaG_PEtot
****************************************
****  written by Lin-Wang Wang
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

****  this is a G space analysition program
******************************************
      use fft_data
      use load_data
      use data
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

      integer status(MPI_STATUS_SIZE)
******************************************
       real*8 AL(3,3),ALt(3,3),AL0(3,3),ALI(3,3),AG0(3,3)
***********************************************
      integer mumG
      parameter (mumG=2000)

c      integer,parameter :: mkk=50
      integer,parameter :: mkk=20

      real*4,allocatable,dimension(:,:,:,:) :: dd
      real*4,allocatable,dimension(:,:,:) :: ddtemp
      real*4,allocatable,dimension(:,:,:) :: akx_st,aky_st,akz_st
      complex*16,allocatable,dimension(:,:) :: ugt_n

      integer status(MPI_STATUS_SIZE)
      real*8 akxt(mumG),akyt(mumG),akzt(mumG)

      character*20 filepot,filewg_out,filewg_in,kpt_file
      character*20 f_xatom
      character*20 outfile

******************************************
      call mpi_init(ierr)
      inode = my_pe()+1
      nnodes=n$pes

      pi=4.0d0*datan(1.0d0)

***********************************************************
***** read in the wavefunction
***** First, just read in the header
***********************************************************
cccccccc all processors read together
            open(10,file='input.AL0')
            read(10,*) AL0(1,1),AL0(2,1),AL0(3,1)
            read(10,*) AL0(1,2),AL0(2,2),AL0(3,2)
            read(10,*) AL0(1,3),AL0(2,3),AL0(3,3)
            read(10,*) if_so
            read(10,*) filewg_in
            read(10,*) kpt_file 
            close(10)
           call get_ALI(AL0,AG0)
           AG0=AG0*2*pi

            open(11,file='akxyz_u0.input')
cccccccc input the k-point for which the graphu.G0.kp will be output
cccccccc They must be in xyz direction, and in a.u.
            read(11,*) akx_u0,aky_u0,akz_u0
            close(11)
            
    

           open(12,file=kpt_file,status='old',action='read',iostat=ierr)
           rewind(12)
           read(12,*) nkpt

           

           call data_allocate_akx(nkpt)
           call ngtot_allocate(nnodes,nkpt)

           read(12,*) iflag, ALxyz
           do kpt=1,nkpt
           read(12,*) akx(kpt),aky(kpt),akz(kpt),weighkpt(kpt)
           enddo
           close(12)

        if(inode==nnodes) then
           open(11,file=filewg_in,form='unformatted',status='old')
  	   rewind(11)
           read(11)n1,n2,n3,mx
           read(11)Ecut
           read(11)AL
           read(11)nnodes_o

           if(nnodes_o.ne.nnodes) then
           write(6,*) "nnodes_o changed, stop", nnodes_o
           call mpi_abort(MPI_COMM_WORLD,ierr)
           endif
         endif     ! take the information from header, keep (11) open for inode=nnodes
    

	  call mpi_barrier(mpi_comm_world,ierr)

	  call mpi_bcast(n1,1,MPI_INTEGER,nnodes-1,MPI_COMM_WORLD,ierr)
          call mpi_bcast(n2,1,MPI_INTEGER,nnodes-1,MPI_COMM_WORLD,ierr)
          call mpi_bcast(n3,1,MPI_INTEGER,nnodes-1,MPI_COMM_WORLD,ierr)
          call mpi_bcast(mx,1,MPI_INTEGER,nnodes-1,MPI_COMM_WORLD,ierr)
          call mpi_bcast(AL,9,MPI_REAL8,nnodes-1,MPI_COMM_WORLD,ierr)
          call mpi_bcast(Ecut,1,MPI_REAL8,nnodes-1,MPI_COMM_WORLD,ierr)
          call mpi_bcast(AL0,1,MPI_REAL8,nnodes-1,MPI_COMM_WORLD,ierr)
	  call mpi_barrier(mpi_comm_world,ierr)

         call get_ALI(AL,ALI)
         do kpt=1,nkpt
         if(iflag.eq.1) then
         akx(kpt)=akx(kpt)*2*pi/ALxyz
         aky(kpt)=aky(kpt)*2*pi/ALxyz
         akz(kpt)=akz(kpt)*2*pi/ALxyz
         else
         ak1_t=akx(kpt)
         ak2_t=aky(kpt)
         ak3_t=akz(kpt)
         akx(kpt)=2*pi*(ALI(1,1)*ak1_t+ALI(1,2)*ak2_t+
     &           ALI(1,3)*ak3_t)
         aky(kpt)=2*pi*(ALI(2,1)*ak1_t+ALI(2,2)*ak2_t+
     &           ALI(2,3)*ak3_t)
         akz(kpt)=2*pi*(ALI(3,1)*ak1_t+ALI(3,2)*ak2_t+
     &           ALI(3,3)*ak3_t)
         endif
         enddo
cccccccccccccccccccccccccccccccccc


       nwg_in=mx

************************************************************
****** doing this initialization after the node nnodes get n1,n2,n3,AL
************************************************************
      nr1x=n1
      nr2x=n2
      nr3x=n3+2
      nr=n1*n2*n3
      mr=n1*n2*(n3+2)
      mr_n=mr/nnodes


c     Calculate the approx. no. of g points, so we can dimension arrays
      volume=al(3,1)*(al(1,2)*al(2,3)-al(1,3)*al(2,2))
     &     +al(3,2)*(al(1,3)*al(2,1)-al(1,1)*al(2,3))
     &     +al(3,3)*(al(1,1)*al(2,2)-al(1,2)*al(2,1))
      volume=dabs(volume)
      delta_k=(2*pi)**3/volume
      totg=(4.0d0/3.0d0*pi*(dsqrt(2.0d0*Ecut))**3)/delta_k
      mg_nx=int(1.1*totg/nnodes)+100

c     Allocate all the arrays that depend on n1,n2,n3,nnodes

      ilocal=1
      islda=1
      call fft_allocate(n1,n2,n3,nnodes)
      call load_allocate(ncolx,nnodes,mg_nx,mr_n)
      call data_allocate(mg_nx,mx,mr_n,ilocal,inode,nkpt,
     & islda,1)
      
       allocate(dd(-mkk:mkk,-mkk:mkk,-mkk:mkk,1:mx))
       allocate(akx_st(-mkk:mkk,-mkk:mkk,-mkk:mkk))
       allocate(aky_st(-mkk:mkk,-mkk:mkk,-mkk:mkk))
       allocate(akz_st(-mkk:mkk,-mkk:mkk,-mkk:mkk))
       allocate(ugt_n(mumG,mx))

      do 2000 kpt=1,nkpt

      dd=0.d0
      akx_st=0.d0
      aky_st=0.d0
      akz_st=0.d0
      ugt_n=dcmplx(0.d0,0.d0)
      akxt=0.d0
      akyt=0.d0
      akzt=0.d0

      call gen_G_comp(kpt,0)
      call get_ug()
**********************************************************
      if (inode==1) write(6,*)'Initialisation complete!'
**********************************************************
      vol=volume

       AG01=dsqrt(AG0(1,1)**2+AG0(2,1)**2+AG0(3,1)**2)
       AG02=dsqrt(AG0(1,2)**2+AG0(2,2)**2+AG0(3,2)**2)
       AG03=dsqrt(AG0(1,3)**2+AG0(2,3)**2+AG0(3,3)**2)
       

       num_warn=0

       numG=0
       do 1000 ig=1,ngtotnod(inode,kpt)

	  gkx=gkx_n(ig,kpt)
	  gky=gky_n(ig,kpt)
	  gkz=gkz_n(ig,kpt)
  
        s1=(gkx*AL0(1,1)+gky*AL0(2,1)+gkz*AL0(3,1))/(2*pi)+300
        s2=(gkx*AL0(1,2)+gky*AL0(2,2)+gkz*AL0(3,2))/(2*pi)+300
        s3=(gkx*AL0(1,3)+gky*AL0(2,3)+gkz*AL0(3,3))/(2*pi)+300

      is1=s1
      s1=s1-is1
      is2=s2
      s2=s2-is2
      is3=s3
      s3=s3-is3

      num_test=0
ccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc put the k' point inside the BZ, fine tune i1,j1,k1

      do 303 i1=-1,1
      do 303 j1=-1,1
      do 303 k1=-1,1
      s11=s1+i1+0.0001
      s22=s2+j1+0.0002
      s33=s3+k1+0.00013

      gkxt=s11*AG0(1,1)+s22*AG0(1,2)+s33*AG0(1,3)
      gkyt=s11*AG0(2,1)+s22*AG0(2,2)+s33*AG0(2,3)
      gkzt=s11*AG0(3,1)+s22*AG0(3,2)+s33*AG0(3,3)


      do 301 ii1=-1,1
      do 301 jj1=-1,1
      do 301 kk1=-1,1

      if(ii1.eq.0.and.jj1.eq.0.and.kk1.eq.0) goto 301

      AGx=AG0(1,1)*ii1+AG0(1,2)*jj1+AG0(1,3)*kk1
      AGy=AG0(2,1)*ii1+AG0(2,2)*jj1+AG0(2,3)*kk1
      AGz=AG0(3,1)*ii1+AG0(3,2)*jj1+AG0(3,3)*kk1

      ss=(gkxt*AGx+gkyt*AGy+gkzt*AGz)/
     &  (AGx**2+AGy**2+AGz**2)

      if(ss.lt.-0.5d0.or.ss.gt.0.5d0) goto 303 

301   continue
      num_test=num_test+1
      i1sh=i1
      j1sh=j1
      k1sh=k1
303   continue

      if(num_test.ne.1) then
      write(6,*) "num_test.ne.1, stop", inode,num_test
      call abort()
      endif
      
      s11=s1+i1sh
      s22=s2+j1sh
      s33=s3+k1sh
cccccccc primary cell reciprocal vector
      akx_t=s11*AG0(1,1)+s22*AG0(1,2)+s33*AG0(1,3)
      aky_t=s11*AG0(2,1)+s22*AG0(2,2)+s33*AG0(2,3)
      akz_t=s11*AG0(3,1)+s22*AG0(3,2)+s33*AG0(3,3)
      akx_t2=akx_t+akx(kpt)
      aky_t2=aky_t+aky(kpt)
      akz_t2=akz_t+akz(kpt)

ccc akx_t2,.. and i11,j11,.. has the one to one correspondence

      i11=(akx_t2*AL(1,1)+aky_t2*AL(2,1)+
     &                    akz_t2*AL(3,1))/(2*pi)*1.002
      j11=(akx_t2*AL(1,2)+aky_t2*AL(2,2)+
     &                    akz_t2*AL(3,2))/(2*pi)*1.002
      k11=(akx_t2*AL(1,3)+aky_t2*AL(2,3)+
     &                    akz_t2*AL(3,3))/(2*pi)*1.002

      akx_st(i11,j11,k11)=akx_t
      aky_st(i11,j11,k11)=aky_t
      akz_st(i11,j11,k11)=akz_t


      if(i11.gt.mkk.or.i11.lt.-mkk.or.j11.gt.mkk.or.j11.lt.-mkk.
     &  or.k11.gt.mkk.or.k11.lt.-mkk) then
      write(6,*) "i11,j11,k11, out range, stop", i11,j11,k11
      call abort()
      endif


       do iwave=1,mx
      dd(i11,j11,k11,iwave)=dd(i11,j11,k11,iwave)+
     & cdabs(ug_n(ig,iwave))**2*vol
      end do


***********************************************************
****** store the wavefunction of k=0, for later analysis
cc      if(i11.eq.0.and.j11.eq.0.and.k11.eq.0) then
      if(dabs(akx_t-akx_u0)+dabs(aky_t-aky_u0)+
     &   dabs(akz_t-akz_u0).lt.0.00001) then

      numG=numG+1

      if(numG.gt.mumG) then
      write(6,*) "numG.gt.mumG, stop", numG, mumG
      stop
      endif

      do iwave=1,mx
      ugt_n(numG,iwave)=ug_n(ig,iwave)
      enddo
      akxt(numG)=gkx
      akyt(numG)=gky
      akzt(numG)=gkz
      endif
***********************************************************

1000  continue
      if(num_warn.gt.0) then
      write(6,*) "num_warning > 0", num_warn
      endif


****************************************************************
***** add up the contribution from different node
****************************************************************
       allocate(ddtemp(-mkk:mkk,-mkk:mkk,-mkk:mkk))
       ncount=(2*mkk+1)**3

       do 200 iwave=1,mx
       call mpi_reduce(dd(-mkk,-mkk,-mkk,iwave),ddtemp,ncount,
     &  MPI_REAL4,MPI_SUM,0,MPI_COMM_WORLD,ierr)

	if(inode.eq.1) then
	do k=-mkk,mkk
	do j=-mkk,mkk
	do i=-mkk,mkk
	dd(i,j,k,iwave)=ddtemp(i,j,k)
	enddo
	enddo
	enddo
	endif

	call mpi_barrier(mpi_comm_world,ierr)
200     continue
        deallocate(ddtemp)

*********** final output
**********************************************
	if(inode.eq.1) then
	open(12,file='graph.G.kp'//char(kpt+48))
	rewind(12)

        numk=0
	do 300 k=-mkk,mkk
	do 300 j=-mkk,mkk
	do 300 i=-mkk,mkk
	s=0.d0
	do iwave=1,mx
	s=s+abs(dd(i,j,k,iwave))
	enddo

	if(s.lt.1.d-20)  goto 300

	 akx_t=akx_st(i,j,k)
	 aky_t=aky_st(i,j,k)
	 akz_t=akz_st(i,j,k)

	 write(12,305) akx_t,aky_t,akz_t,(dd(i,j,k,iwave),iwave=1,mx)
	 numk=numk+1

300     continue
305    format(3(f6.4,1x),2x,10(E8.3,1x))
        write(6,*) "num of line in graph.G=", numk

	close(12)
	endif
       
**********************************************************
****** output the k=0 wavefunction in graph.G0
**********************************************************
        call mpi_allreduce(numG,numG_tot,1, 
     &  MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD,ierr)

	if(inode.eq.1) then
	open(13,file='graphu.G0.kp'//char(kpt+48))
	rewind(13)
        write(13,*) numG_tot,mx
	endif

        do 700 j=1,nnodes

	if(j.eq.inode.and.j.ne.1) then
	call mpi_send(numG,1,MPI_INTEGER,0,200,MPI_COMM_WORLD,ierr)
	call mpi_send(ugt_n,mumG*mx,MPI_DOUBLE_COMPLEX,0,
     &             201,MPI_COMM_WORLD,ierr)
	call mpi_send(akxt,mumG,MPI_REAL8,0,202,MPI_COMM_WORLD,ierr)
	call mpi_send(akyt,mumG,MPI_REAL8,0,203,MPI_COMM_WORLD,ierr)
	call mpi_send(akzt,mumG,MPI_REAL8,0,204,MPI_COMM_WORLD,ierr)
	endif

	if(inode.eq.1.and.j.ne.1) then
	call mpi_recv(numG,1,MPI_INTEGER,j-1,200,MPI_COMM_WORLD,
     &  	status,ierr)
	call mpi_recv(ugt_n,mumG*mx,MPI_DOUBLE_COMPLEX,j-1,
     &              201,MPI_COMM_WORLD,status,ierr)
	call mpi_recv(akxt,mumG,MPI_REAL8,j-1,202,MPI_COMM_WORLD,
     &  	status,ierr)
	call mpi_recv(akyt,mumG,MPI_REAL8,j-1,203,MPI_COMM_WORLD,
     &  	status,ierr)
	call mpi_recv(akzt,mumG,MPI_REAL8,j-1,204,MPI_COMM_WORLD,
     &  	status,ierr)
	endif

        call mpi_barrier(mpi_comm_world,ierr)

	if(inode.eq.1.and.numG.gt.0) then

        write(6,*) "numG=", numG

	 do i=1,numG
	 write(13,302) akxt(i),akyt(i),akzt(i),
     &   (ugt_n(i,iwave),iwave=1,mx)
	 enddo

	endif

700     continue
302    format(3(f7.4,1x),2x,10(E10.4,1x,E10.4,3x))

	if(inode.eq.1) close(13)

2000    continue


      call data_deallocate()
      call load_deallocate()
      call mpi_finalize(ierr)

      stop


      contains
*************************************************

****************************************
cccccc this is just for one kpt. 
****************************************

       subroutine get_ug()

       complex*16,allocatable,dimension(:) :: ugtemp

        call mpi_barrier(MPI_COMM_WORLD,ierr)

cccccc  file (11) should have been opened before the do loop for kpt=1,nkpt
ccccccccccccccccccccccccccccccccccccccccc

         call mpi_barrier(MPI_COMM_WORLD,ierr)

         if(inode==nnodes) then
            allocate(ugtemp(mg_nx))

            do i=0,nnodes-2
               do iwavefun=1,mx
                  read(11)ugtemp    ! assuming the mg_nx is the same as before
                  ug_n(:,iwavefun)=ugtemp
               end do


          call mpi_send(ug_n,mg_nx*mx,MPI_DOUBLE_COMPLEX,i,
     &            102,MPI_COMM_WORLD,ierr)
            end do
c     Now read in the data for node nnodes
            do iwavefun=1,mx
               read(11)ugtemp
               ug_n(:,iwavefun)=ugtemp
            end do
            deallocate(ugtemp)

         else   ! for other nodes

         call mpi_recv(ug_n,mg_nx*mx,MPI_DOUBLE_COMPLEX,nnodes-1,
     &       102,MPI_COMM_WORLD,status,ierr)

         end if

        return
        end subroutine get_ug

        end

