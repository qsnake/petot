      subroutine fftprep_comp(nr1,nr2,nr3)

c Written by A. Canning (CRAY-EPFL) 25th July 94  
c
c Small and large sphere wavefunction prep prog for FFT.
c (Small only for escan code )
c this subroutine does the preparation work for the 
c 3d fft adapted to the CP code. 
c The first pass of the FFT ie x dir (nr1) is done on psi
c in the format x column load balanced format, the full x columns
c of g vectors are present and padded to zero. It will 
c only do the FFT on the non-zero g columns in the 
c four (two for gamma) corners of the cube giving quarter cylinders. 
c Only these cylinders are transpoed into the y dir (nr2)
c which gives two (one gamma) slabs after the FFT in the y dir. 
c the width of the two slabs is mgz ie max radius of sphere
c in the z direction including the zero point. One slab
c (the one with the zero component) will actually be one point
c wider than the other but we just padd the other one so they
c are both of width mgz (not used for gamma).
c The two (one gamma) slices are then transposed.
c This is done by a chunking local transpose of which the 
c correct chunk of the chunk is sent transposed to its target processor.
c the z column index will tell us the target node.
c This is done via a strided put so the transpose is only
c explicitly done when the z piece of column remains on the 
c node. Messy but runs fast on T3D/T3D as transposes are 
c done locally and not written to memory before sending. 
c The FFT is then done 
c in the Z dir (nr3) and left as z columns in real space. 
c
c In this routine we choose the distributions for the 
c y-dir and z-dir FFT's and then calculate the target 
c nodes for the shmem puts for the two transposes. 
c and also the chunk addresses and widths for the chunking
c transpose. 
c in the y and z direction we just treat the system as a 
c series of y and z dir columns and hand out the same 
c number of columns to each PE (load 
c-----------------------------------------------------
c revised ADV 12/8/94: call to nz_find included to
c                           determine nonzero columns
c
c last revised ADV 7/7/95 to parametrise the no. of 
c                         chunks on each PE
c AMC mpi version for SP  7/9/00
c
      use fft_data
      use load_data
      
      implicit none

      include 'mpif.h'
c
c
c scalars used 
c

      integer i,j,ic,ico,idum,inodec,itaradd,itnode,ix,iy,iz,
     c        jcol,ngy,ngyadd,ngz,igdum,nr1,nr2,nr3,inode,nr3u,ierr,
     c        ilocadd,izb,iw,ib,itar,isc,ii,mgztemp
      integer nnodes
      integer inode_tot,nnodes_tot,icolor,num_group,MPI_COMM_K,
     &       MPI_COMM_N

      integer ifirst
      data ifirst/0/
      save ifirst

c
c mpi version arrays 
c
      integer mpistatus(mpi_status_size)
      integer ivdum(mnrx),ivpacn_cum(nnodes)
      integer ireq(nnodes)

      common /mpi_data/inode,nnodes,inode_tot,nnodes_tot,
     &  icolor,num_group,MPI_COMM_K,MPI_COMM_N


      call mpi_barrier(MPI_COMM_K,ierr)

c     nr3u = nr3/2+1
c
c set ups for different direction FFT's
c
      call system_ccfft(1,1,nr1,1.0d0,0,0,tabnr1fw,0,0,ntabnr1,0)
      call system_ccfft(1,1,nr2,1.0d0,0,0,tabnr2fw,0,0,ntabnr2,0)
      call system_ccfft(1,1,nr3,1.0d0,0,0,tabnr3fw,0,0,ntabnr3,0)

      call system_ccfft(1,-1,nr1,1.0d0,0,0,tabnr1in,0,0,ntabnr1,0)
      call system_ccfft(1,-1,nr2,1.0d0,0,0,tabnr2in,0,0,ntabnr2,0)
      call system_ccfft(1,-1,nr3,1.0d0,0,0,tabnr3in,0,0,ntabnr3,0)
c
c first make list of non-zero columns of psi wavefunctions
c and get max value radius of ngz ie. mgz 
c only treat non-zero cols in escan 
c 
      mgz = 0
c
      do jcol = 1,ncol(inode)
         ngz = izcol(jcol,inode)
ccccccc
ccccc         if(ngz.gt.mgz.and.ngz.lt.nr3/2) mgz = ngz
ccccccc modified by Lin-Wang, for odd number nr3
        if(ngz.gt.mgz.and.ngz.lt.(nr3+1)/2) mgz = ngz
      enddo 
c
      call mpi_allreduce(mgz,mgztemp,1,mpi_integer,mpi_max,
     &                    MPI_COMM_K, ierr)
      mgz = mgztemp
      call mpi_barrier(MPI_COMM_K,ierr)

      if(inode_tot.eq.1.and.ifirst.eq.0) then
         write(*,*)
         write(*,*) 'minimum possible value of mgz = ',mgz
      endif
 
 100  if(mod(nr1*mgz,nnodes).ne.0) then
         mgz=mgz+1
         goto 100
      endif
 
      if(inode_tot.eq.1.and.ifirst.eq.0) then
         write(*,*) 'optimal (and used) value of mgz = ',mgz
cc         if(mgz.gt.(nr3/2)) then 
         if(mgz.gt.(nr3+1)/2) then 
           write(*,*) 'load balancing failed, change n3 ' 
         endif
      endif
 
c     ncolnz = ico

      call mpi_barrier(MPI_COMM_K,ierr)
c
c now decide on blocking structure for the y dir FFT
c We have two slabs of mgz*nr1 y columns 
c We will do the x dir FFT in place in psi 
c then transpose to psiy then do FFT in y dir then 
c transpose back into psi and do z dir FFT
c
c For each FFT we assume each PE has the same number of
c columns  Check this is true !!
c
       if(mod(nr1*mgz,nnodes).ne.0) then
        write(*,*) 'nr1*mgz not multiple of nnodes, FFT stopped'
        stop
       endif

       call mpi_barrier(MPI_COMM_K,ierr)
c 
c
       if(mod(nr1*nr2,nnodes).ne.0) then
        write(*,*) 'nr1*nr2 not multiple of nnodes, FFT stopped'
        stop
       endif

       call mpi_barrier(MPI_COMM_K,ierr)
c 
c  choose data layout for  slab  for y dir FFT
c  each PE will take (nr1*mgz)/nnodes consecutive 
c  columns 
c
c  first slab goes from 1 to mgz in z direction 
c ASSUME mgz < nnodes so that each PE contains at
c max two parts of slab.
c 

       ncoly = (nr1*mgz)/nnodes
       ico = 0
       inodec = 1

       do ix = 1,nr1
        do iz = 1,mgz
         ico = ico + 1
         izcol_y(ico,inodec) = iz 
         ixcol_y(ico,inodec) = ix 
         ipe_y(ix,iz) = inodec
         ibca_y(ix,iz) = (ico-1)*nr2
         if(ico.eq.ncoly) then 
          inodec = inodec + 1
          ico = 0
         endif !(ico.eq.ncoly)
        enddo !iz
       enddo !ix
c
c put second slab on top of first with same type
c of layout in memory
c
c
       ico = ncoly
       inodec = 1

       do ix = 1,nr1
        do iz = 1+(nr3-mgz),nr3
         ico = ico + 1
         izcol_y(ico,inodec) = iz
         ixcol_y(ico,inodec) = ix
         ipe_y(ix,iz) = inodec
         ibca_y(ix,iz) = (ico-1)*nr2
         if(ico.eq.2*ncoly) then
          inodec = inodec + 1
          ico = ncoly
         endif !(ico.eq.2*ncoly)
        enddo !ix
       enddo !iz
c
c calculate starting addresses for chunks and also 
c their width in the z direction 
c each proc has no more than   nofchks   chunks 
c these chunks will be transposed locally 
c before being sent in z strips to the appropriate
c processor and memory location  before z dir fft
      
      if(ifirst.eq.1) then    ! ifirst.eq.1, not the first call of fftprep_comp
      deallocate(izchb)
      deallocate(ixch)
      deallocate(ichw)
      endif


      allocate(izchb(ncoly))
      allocate(ixch(ncoly))
      allocate(ichw(ncoly))

      izchb(1) = izcol_y(1,inode)
       ixch(1) = ixcol_y(1,inode)
       ix = ixcol_y(1,inode)
       ichunk = 1
       ichw=0
       do ic = 1,ncoly
         iz = izcol_y(ic,inode)
         if(ix.ne.ixcol_y(ic,inode)) then
          ix = ixcol_y(ic,inode) 
          ichunk = ichunk + 1
          izchb(ichunk) = izcol_y(ic,inode)
          ixch(ichunk) = ix
         endif
         ichw(ichunk) = ichw(ichunk)+1
       enddo


       izchb(ichunk+1) = izcol_y(1+ncoly,inode)
       ixch(ichunk+1) = ixcol_y(1+ncoly,inode)
       ix = ixcol_y(1+ncoly,inode)
       ichunk = ichunk+1
       do ic = 1+ncoly,ncoly*2
         iz = izcol_y(ic,inode)
         if(ix.ne.ixcol_y(ic,inode)) then
          ix = ixcol_y(ic,inode)
          ichunk = ichunk + 1
          izchb(ichunk) = izcol_y(ic,inode)
          ixch(ichunk) = ix
         endif
         ichw(ichunk) = ichw(ichunk)+1
       enddo

c
c check no more than nofchks chunks on each PE
c
       if(ichunk.gt.ncoly) then
       write(*,*) 'More than ncoly chunks on PEs '
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif


c choose layout for final FFT in z direction 
c this will be on the whole system ie nr1*nr2 columns
c divided by the number of processors 

       ncolz = nr1*nr2/nnodes

       ico = 0
       inodec = 1
 
       do ix = 1,nr1
        do iy = 1,nr2
         ico = ico + 1
         iycol_z(ico,inodec) = iy
         ixcol_z(ico,inodec) = ix
         ipe_z(iy,ix) = inodec
         ibca_z(iy,ix) = (ico-1)*nr3
         if(ico.eq.ncolz) then 
          inodec = inodec + 1
          ico = 0
         endif !(ico.eq.ncolz)
        enddo !iy
       enddo !ix

c
c calculate ivecadd for the first scatterred put
c in the fft's
c
       icount = 0
 
       do jcol = 1,ncol(inode)
        ngy = iycol(jcol,inode)
        ngz = izcol(jcol,inode)
        ngyadd = 1+(ngy-1)
        do ix = 1,nr1
          itnode = ipe_y(ix,ngz)
          itaradd = ngyadd + ibca_y(ix,ngz) - 1
          idum = icount(itnode)+1

          ivecadd(idum,itnode) = itaradd
 
          icount(itnode) = icount(itnode)+1
 
        enddo !ix
       enddo !i
c
c   check that we have allocated enough memory for vec,ivecadd
c
       idum = mnrx/(nnodes)
       do i = 1,nnodes
        if(icount(i).gt.idum) then
         write(*,*) 'ERROR not enough memory allocated for vec,ivecadd'
         write(*,*) ' array dim = ',idum,' min dim = ',icount(i)
         call mpi_abort(MPI_COMM_WORLD,ierr)
        endif
       enddo

       call mpi_barrier(MPI_COMM_K,ierr)
c
c put ivecadd into 1 dim array  ivdum
c
      idum = 1
      do i = 1,nnodes
        ivpacn1(i) = icount(i)
        do j = 1,icount(i)
          ivdum(idum) = ivecadd(j,i)
          idum = idum + 1
        enddo
      enddo

      idum = idum -1

      if(idum.gt.mnrx) then
        write(*,*) " ivpac,ivdum arrays not large enough "
        call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

c
c get cumulative address for ivpac1
c
      ivpacn_cum(1) = 1
      do i= 1,nnodes-1
        ivpacn_cum(i+1)= ivpacn_cum(i)+ivpacn1(i)
      enddo

      idum = idum + ivpacn1(nnodes) -1


      if(idum.gt.mnrx) then
        write(*,*) " ivpac,ivdum arrays not large enough "
        call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

c   
c set up gather indexes for mpi first send
c   
      ivpacn1  = 0

      do jcol = 1,ncol(inode)
         ilocadd = 1+(jcol-1)*nr1
         ngy = iycol(jcol,inode)
         ngz = izcol(jcol,inode)
         ngyadd = 1+(ngy-1)
         do ix = 1,nr1
            itnode = ipe_y(ix,ngz)
            idum = ivpacn_cum(itnode)+ivpacn1(itnode)

            ivpac1(idum) = ilocadd

            ivpacn1(itnode) = ivpacn1(itnode)+1

            ilocadd = ilocadd + 1

         enddo                  !ix
      enddo                     !i
c
c communicate scatter indexes for mpi first send
c   ie transformation of cylinder
c
       do i = 1,nnodes
        call mpi_isend(ivpacn1(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
       enddo

       do i = 1,nnodes
        call mpi_recv(ivunpn1(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       idum = 1
       do i = 1,nnodes
        call mpi_isend(ivdum(idum),ivpacn1(i),mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
        idum = idum + ivpacn1(i)
       enddo

       idum = 1
       do i = 1,nnodes
        call mpi_recv(ivunp1(idum),ivunpn1(i),mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
        idum = idum + ivunpn1(i)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       ivunp1 = ivunp1 + 1

c
c set up packing index for second transformation of slice
c
      ivpacn2 = 0

      do ic = 1,ichunk
         izb = izchb(ic)
         ix = ixch(ic)
         iw = ichw(ic)
         ib = ibca_y(ix,izb)
         do iy = 1,nr2
            itnode = ipe_z(iy,ix)
            itaradd = ibca_z(iy,ix)+(izb-1)+1
            itar = itaradd
            isc = ib+1
            do ii = 1,iw
              ivpacn2(itnode) = ivpacn2(itnode)+1
              itar = itar + 1
              isc = isc + nr2
            enddo
            ib = ib + 1
         enddo
      enddo
c
c get cumulative address for ivpac2
c
      ivpacn_cum(1) = 1
      do i= 1,nnodes-1
        ivpacn_cum(i+1)= ivpacn_cum(i)+ivpacn2(i)
      enddo

      ivpacn2 = 0

      do ic = 1,ichunk
         izb = izchb(ic)
         ix = ixch(ic)
         iw = ichw(ic)
         ib = ibca_y(ix,izb)
         do iy = 1,nr2
            itnode = ipe_z(iy,ix)
            itaradd = ibca_z(iy,ix)+(izb-1)+1
            itar = itaradd
            isc = ib+1
            do ii = 1,iw
              idum = ivpacn_cum(itnode)+ivpacn2(itnode)
              ivpac2(idum) = isc
              ivdum(idum) = itar
              ivpacn2(itnode) = ivpacn2(itnode)+1
              itar = itar + 1
              isc = isc + nr2
            enddo
            ib = ib + 1
         enddo
      enddo

c
c communicate scatter indexes for mpi second send
c
       do i = 1,nnodes
        call mpi_isend(ivpacn2(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
       enddo

       do i = 1,nnodes
        call mpi_recv(ivunpn2(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       idum = 1
       do i = 1,nnodes
        call mpi_isend(ivdum(idum),ivpacn2(i),mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
        idum = idum + ivpacn2(i)
       enddo

       idum = 1
       do i = 1,nnodes
        call mpi_recv(ivunp2(idum),ivunpn2(i),mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
        idum = idum + ivunpn2(i)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       if(ifirst.eq.0) ifirst=1

       call mpi_barrier(MPI_COMM_K,ierr)

       return
       end

      
