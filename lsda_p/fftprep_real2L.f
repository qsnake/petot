      subroutine fftprep_real2L(nr1L,nr2L,nr3L,mr_nL)

ccccccccccccccccccccccccccccccccccccccccccccccc
cccc written by Andew Canning, 2001
ccccccccccccccccccccccccccccccccccccccccccccccc


c Large sphere potential  prep prog for FFT starts here
c This FFT is only used in vofrho
c this part does the preparation work for the 
c 3d fft adapted to the CP code. 
c The first pass of the FFT ie x dir (nr1) is done on psi
c in the format of vofrho. 
c This prep prog is similar to fftprep in that 
c we go from prog format to y,z,x and then to z,y,x
c
c In this routine we choose the distributions for the 
c y-dir and z-dir FFT's and then calculate the target 
c nodes for the shmem puts for the two transposes. 
c and also the chunk addresses and widths for the chunking
c transpose which takes y,z,x to z,y,x 
c
c first make list of non-zero columns of psi wavefunctions
c and get max value radius of ngz ie. mgz
c only treat non-zero cols in escan
c
      use fft_data
      use load_data

      implicit none

      include 'param.escan_real'
      include 'mpif.h'
c
c scalars used
c

      integer i,j,ic,ico,idum,inodec,itaradd,itnode,ix,iy,iz,
     c        jcol,ngy,ngyadd,ngz,igdum,nr1L,nr2L,nr3L,nr3uL,ierr,
     c        ilocadd,izb,iw,ib,itar,isc,ii,mgz2Ltmp,
     c        n1_inv,n2_inv,indepg_d,jjnode_dum,indepg,ig,
     c        k0npac_max,k0nunp_max,isum,imax,mr_nL

c
c mpi version arrays
c
      integer mpistatus(mpi_status_size)
      integer ivdum(mnr2xL),ivpacn_cum(nnodes)
ccccc      integer k0idum(nr1*nr2,nnodes)    ! old stupid way, too much memory
      integer ireq(nnodes)

      integer, allocatable,dimension(:,:) ::  k0idum


      call mpi_barrier(MPI_COMM_K,ierr)

      nr3uL = nr3L/2+1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call system_ccfft(1,1,nr1L,1.0d0,0,0,tabnr1lfwL,
     &     0,0,ntabnr1L,0)
      call system_ccfft(1,1,nr2L,1.0d0,0,0,tabnr2lfwL,
     &     0,0,ntabnr2L,0)
      call system_scfft(1,1,nr3L,1.0d0,0,0,tabnr3lrcL,
     &     0,0,ntabnr3L,0)

      call system_ccfft(1,-1,nr1L,1.0d0,0,0,tabnr1linL,
     &     0,0,ntabnr1L,0)
      call system_ccfft(1,-1,nr2L,1.0d0,0,0,tabnr2linL,
     &     0,0,ntabnr2L,0)
      call system_csfft(1,-1,nr3L,1.0d0,0,0,tabnr3lcrL,
     &     0,0,ntabnr3L,0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mgz2L = 0
c
       do jcol = 1,ncol2L(inode)
         ngz = izcol2L(jcol,inode)
         if(ngz.gt.mgz2L) mgz2L = ngz
       enddo
c
       call global_maxi(mgz2L)
c
       call mpi_barrier(MPI_COMM_K,ierr)
 
       call mpi_allreduce(mgz2L,mgz2Ltmp,1,mpi_integer,mpi_max,
     &                    MPI_COMM_K, ierr)
       mgz2L = mgz2Ltmp
       call mpi_barrier(MPI_COMM_K,ierr)

       if(inode_tot.eq.1) then
         write(*,*)
         write(*,*) 'minimum possible value of (pot)mgz2L=',mgz2L
       endif

 101   if(mod(nr1L*mgz2L,nnodes).ne.0) then
         mgz2L=mgz2L+1
         goto 101
       endif

       if(inode_tot.eq.1) then
         write(*,*) 'optimal(and used)value of(pot) mgz2L=',mgz2L
         if(mgz2L.gt.(1+nr3L/2)) then
           write(*,*) 'load balancing failed, change n3L '
cccccc PROBLEMS ????  WARNING
         endif
       endif

       call mpi_barrier(MPI_COMM_K,ierr)

c      ncolnz = ico

c
c now decide on blocking structure for the y dir FFT
c We have one slabs of mgz2*nr1 y columns
c We will do the x dir FFT in place in psi
c then transpose to psiy then do FFT in y dir then
c transpose back into psi and do z dir FFT
c
c For each FFT we assume each PE has the same number of
c columns  Check this is true !!
c
       if(mod(nr1L*mgz2L,nnodes).ne.0) then
        write(*,*) 'nr1L*mgz2L not multiple of(pot)nnodes, FFT stop'
        call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       call mpi_barrier(MPI_COMM_K,ierr)

c  choose data layout for  slab  for y dir FFT
c  each PE will take (nr1*mgz2)/nnodes consecutive
c  columns
c
c  first slab goes from 1 to mgz2 in z direction
c

       ncoly2L = (nr1L*mgz2L)/nnodes
       ico = 0
       inodec = 1

       do ix = 1,nr1L
        do iz = 1,mgz2L
         ico = ico + 1
         izcol_y2L(ico,inodec) = iz
         ixcol_y2L(ico,inodec) = ix
         ipe_y2L(ix,iz) = inodec
         ibca_y2L(ix,iz) = (ico-1)*2*nr2L
         if(ico.eq.ncoly2L) then
          inodec = inodec + 1
          ico = 0
         endif !(ico.eq.ncoly2L)
        enddo !iz
       enddo !ix

c calculate starting addresses for chunks and also
c their width in the z direction
c each proc has no more than   nofchks   chunks
c these chunks will be transposed locally
c before being sent in z strips to the appropriate
c processor and memory location  before z dir fft

      allocate(izchb2L(ncoly2L))
      allocate(ixch2L(ncoly2L))
      allocate(ichw2L(ncoly2L))

      izchb2L(1) = izcol_y2L(1,inode)
       ixch2L(1) = ixcol_y2L(1,inode)
       ix = ixcol_y2L(1,inode)
       ichunk2L = 1
       ichw2L=0
       do ic = 1,ncoly2L
         iz = izcol_y2L(ic,inode)
         if(ix.ne.ixcol_y2L(ic,inode)) then
          ix = ixcol_y2L(ic,inode)
          ichunk2L = ichunk2L + 1
          izchb2L(ichunk2L) = izcol_y2L(ic,inode)
          ixch2L(ichunk2L) = ix
         endif
         ichw2L(ichunk2L) = ichw2L(ichunk2L)+1
       enddo
c
c check no more than nofchks chunks on each PE
c
       if(ichunk2L.gt.ncoly2L) then
       write(*,*) 'More than ncoly2L chunks on PEs '
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

c choose layout for final FFT in z direction
c this will be on the whole system ie nr1*nr2 columns
c divided by the number of processors

       ncolz2L = nr1L*nr2L/nnodes

       ico = 0
       inodec = 1

       do ix = 1,nr1L
        do iy = 1,nr2L
         ico = ico + 1
         iycol_z2L(ico,inodec) = iy
         ixcol_z2L(ico,inodec) = ix
         ipe_z2L(iy,ix) = inodec
         ibca_z2L(iy,ix) = (ico-1)*2*nr3uL
         if(ico.eq.ncolz2L) then
          inodec = inodec + 1
          ico = 0
         endif !(ico.eq.ncolz2L)
        enddo !iy
       enddo !ix
c
c calculate ivecadd for the first scatterred put
c in the fft's
c
       icount2L = 0

       do jcol = 1,ncol2L(inode)
        ngy = iycol2L(jcol,inode)
        ngz = izcol2L(jcol,inode)
        ngyadd = 1+(ngy-1)*2
        do ix = 1,nr1L
          itnode = ipe_y2L(ix,ngz)
          itaradd = ngyadd + ibca_y2L(ix,ngz) - 1
          idum = icount2L(itnode)+1

          ivecadd2L(idum,itnode) = itaradd
          ivecadd2L(idum+1,itnode) = itaradd +1

          icount2L(itnode) = icount2L(itnode)+2

        enddo !ix
       enddo !i
c
c   check that we have allocated enough memory for vec,ivecadd
c
       idum = mnr2xL/(nnodes)
       do i = 1,nnodes
        if(icount2L(i).gt.idum) then
         write(*,*) 'not enough memory allocated for vec2L,ivecadd2L'
         write(*,*) ' array dim = ',idum,' min dim = ',icount2L(i)
         call mpi_abort(MPI_COMM_WORLD,ierr)
        endif
       enddo

       call mpi_barrier(MPI_COMM_K,ierr)
c
c put ivecadd into 1 dim array  ivdum
c
      idum = 1
      do i = 1,nnodes
        ivpacn1lL(i) = icount2L(i)
        do j = 1,icount2L(i)
          ivdum(idum) = ivecadd2L(j,i)
          idum = idum + 1
        enddo
      enddo

      idum = idum -1

      if(idum.gt.mnr2xL) then
        write(*,*) " ivpac1lL,ivdum arrays not large enough "
        call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

c
c get cumulative address for ivpac1
c
      ivpacn_cum(1) = 1
      do i= 1,nnodes-1
        ivpacn_cum(i+1)= ivpacn_cum(i)+ivpacn1lL(i)
      enddo

      idum = idum + ivpacn1lL(nnodes) -1


      if(idum.gt.mnr2xL) then
        write(*,*) " ivpac1lL,ivdum arrays not large enough "
        call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
c  
c set up gather indexes for mpi first send
c
      ivpacn1lL  = 0

      do jcol = 1,ncol2L(inode)
         ilocadd = 1+(jcol-1)*(2*nr1L)
         ngy = iycol2L(jcol,inode)
         ngz = izcol2L(jcol,inode)
         ngyadd = 1+(ngy-1)*2
         do ix = 1,nr1L
            itnode = ipe_y2L(ix,ngz)
            idum = ivpacn_cum(itnode)+ivpacn1lL(itnode)

            ivpac1lL(idum) = ilocadd
            ivpac1lL(idum+1) = ilocadd + 1

            ivpacn1lL(itnode) = ivpacn1lL(itnode)+2

            ilocadd = ilocadd + 2

         enddo                  !ix
      enddo                     !i
c
c communicate scatter indexes for mpi first send
c   ie transformation of cylinder
c

       do i = 1,nnodes
        call mpi_isend(ivpacn1lL(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
       enddo

       do i = 1,nnodes
        call mpi_recv(ivunpn1lL(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       idum = 1
       do i = 1,nnodes
        call mpi_isend(ivdum(idum),ivpacn1lL(i),mpi_integer,i-1,
     &               inode,MPI_COMM_K,ireq(i),ierr)
        idum = idum + ivpacn1lL(i)
       enddo

       idum = 1
       do i = 1,nnodes
        call mpi_recv(ivunp1lL(idum),ivunpn1lL(i),mpi_integer,
     &              i-1,i,MPI_COMM_K,mpistatus,ierr)
        idum = idum + ivunpn1lL(i)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       ivunp1lL = ivunp1lL + 1
c
c set up packing index for second transformation of slice
c
      ivpacn2lL = 0

      do ic = 1,ichunk2L
         izb = izchb2L(ic)
         ix = ixch2L(ic)
         iw = ichw2L(ic)
         ib = ibca_y2L(ix,izb)
         do iy = 1,nr2L
            itnode = ipe_z2L(iy,ix)
            itaradd = ibca_z2L(iy,ix)+(izb-1)+1
            itar = itaradd
            isc = ib+1
            do ii = 1,iw
              ivpacn2lL(itnode) = ivpacn2lL(itnode)+2
              itar = itar + 2
              isc = isc + nr2L*2
            enddo
            ib = ib + 2
         enddo
      enddo
c
c get cumulative address for ivpac2
c
      ivpacn_cum(1) = 1
      do i= 1,nnodes-1
        ivpacn_cum(i+1)= ivpacn_cum(i)+ivpacn2lL(i)
      enddo

      ivpacn2lL = 0

      do ic = 1,ichunk2L
         izb = izchb2L(ic)
         ix = ixch2L(ic)
         iw = ichw2L(ic)
         ib = ibca_y2L(ix,izb)
         do iy = 1,nr2L
            itnode = ipe_z2L(iy,ix)
            itaradd = ibca_z2L(iy,ix)+(izb-1)*2+1
            itar = itaradd
            isc = ib+1
            do ii = 1,iw
              idum = ivpacn_cum(itnode)+ivpacn2lL(itnode)
              ivpac2lL(idum) = isc
              ivpac2lL(idum+1) = isc+1
              ivdum(idum) = itar
              ivdum(idum+1) = itar+1
              ivpacn2lL(itnode) = ivpacn2lL(itnode)+2
              itar = itar + 2
              isc = isc + nr2L*2
            enddo
            ib = ib + 2
         enddo
      enddo

c
c communicate scatter indexes for mpi second send
c
       do i = 1,nnodes
        call mpi_isend(ivpacn2lL(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
       enddo

       do i = 1,nnodes
        call mpi_recv(ivunpn2lL(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       idum = 1
       do i = 1,nnodes
        call mpi_isend(ivdum(idum),ivpacn2lL(i),mpi_integer,
     &              i-1,inode,MPI_COMM_K,ireq(i),ierr)
        idum = idum + ivpacn2lL(i)
       enddo

       idum = 1
       do i = 1,nnodes
        call mpi_recv(ivunp2lL(idum),ivunpn2lL(i),mpi_integer,
     &             i-1,i,MPI_COMM_K,mpistatus,ierr)
        idum = idum + ivunpn2lL(i)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

      call mpi_barrier(MPI_COMM_K,ierr)
c
c mpi indexes for filling in of half sphere k.eq.0 for ffts
c

      k0npacL = 0

      do ig = 1, ngtotnod2L(inode)

         indepg=(jjcol2L(n2p2_nL(ig),n3p2_nL(ig))-1)*nr1L+n1p2_nL(ig)

         if(n3p2_nL(ig).eq.1) then
            n1_inv = mod(n1L-n1p2_nL(ig)+1,nr1L)+1
            n2_inv = mod(n2L-n2p2_nL(ig)+1,nr2L)+1
            indepg_d = (jjcol2L(n2_inv,1)-1)*n1L+n1_inv
            jjnode_dum = jjnode2L(n2_inv,1)
            k0npacL(jjnode_dum)=k0npacL(jjnode_dum)+1
ccccc            k0idum(k0npacL(jjnode_dum),jjnode_dum)=2*indepg_d-1
            k0npacL(jjnode_dum)=k0npacL(jjnode_dum)+1
cccccc           k0idum(k0npacL(jjnode_dum),jjnode_dum)=2*indepg_d
         endif

      enddo


      do i = 1,nnodes
        call mpi_isend(k0npacL(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
      enddo

      do i = 1,nnodes
        call mpi_recv(k0nunpL(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      call mpi_barrier(MPI_COMM_K,ierr)
*******************************************************
      k0npac_max=1
      k0nunp_max=1
      do i=1,nnodes
      if(k0npacL(i).gt.k0npac_max) k0npac_max=k0npacL(i)
      if(k0nunpL(i).gt.k0nunp_max) k0nunp_max=k0nunpL(i)
      enddo

      allocate(k0idum(k0npac_max,nnodes))
      call fft_allocate2L(k0nunp_max,nnodes)  ! allocate k0iunp(k0nunp_max,nnodes)
*******************************************************
ccccc recalculate k0npac, and this time: k0idum 

      k0npacL = 0

      do ig = 1, ngtotnod2L(inode)

         indepg=(jjcol2L(n2p2_nL(ig),n3p2_nL(ig))-1)*nr1L+n1p2_nL(ig)

         if(n3p2_nL(ig).eq.1) then
            n1_inv = mod(n1L-n1p2_nL(ig)+1,nr1L)+1
            n2_inv = mod(n2L-n2p2_nL(ig)+1,nr2L)+1
            indepg_d = (jjcol2L(n2_inv,1)-1)*n1L+n1_inv
            jjnode_dum = jjnode2L(n2_inv,1)
            k0npacL(jjnode_dum)=k0npacL(jjnode_dum)+1
            k0idum(k0npacL(jjnode_dum),jjnode_dum)=2*indepg_d-1
            k0npacL(jjnode_dum)=k0npacL(jjnode_dum)+1
            k0idum(k0npacL(jjnode_dum),jjnode_dum)=2*indepg_d
         endif

      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_K,ierr)


      do i = 1,nnodes
        call mpi_isend(k0idum(1,i),k0npacL(i),mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
      enddo

      do i = 1,nnodes
        call mpi_recv(k0iunpL(1,i),k0nunpL(i),mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      call mpi_barrier(MPI_COMM_K,ierr)

cccccc  safety check
      isum=0
      imax=0
      do i=1,nnodes
      isum=isum+ivpacn2lL(i)
        do j=1,k0nunpL(i)
        if(k0iunpL(j,i).gt.imax) imax=k0iunpL(j,i)
        enddo
      enddo
        if(isum.gt.mr_nL.or.imax.gt.mr_nL) then
       write(6,*) "increase mr_nL, due to imbalance, stop",
     & inode,isum,imax,mr_nL
cccccc This is probably nature, due to the load imbalance, some
cccccc node might have slightly more FFT columns than other nodes,
cccccc  then isum,imax > (n1L*n2L*(n3L+2)/nnodes)
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif


      return
      end
