      subroutine fftprep_real2L2(nr1L2,nr2L2,nr3L2,mr_nL2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc written by Andew Canning, 2001
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc This subroutine need to be checked more carefully, especially 
cccc for why it might need mr_nL2 larger than n1L2*n2L2*(n3L2+2)/nnodes
cccc Is it really due to imbalance? Perhaps. 

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
     c     jcol,ngy,ngyadd,ngz,igdum,nr1L2,nr2L2,nr3L2,nr3uL2,ierr,
     c        ilocadd,izb,iw,ib,itar,isc,ii,mgz2L2tmp,
     c        n1_inv,n2_inv,indepg_d,jjnode_dum,indepg,ig,
     c        k0npac_max,k0nunp_max,isum,imax,mr_nL2

c
c mpi version arrays
c
      integer mpistatus(mpi_status_size)
      integer ivdum(mnr2xL2),ivpacn_cum(nnodes)
ccccc      integer k0idum(nr1*nr2,nnodes)    ! old stupid way, too much memory
      integer ireq(nnodes)

      integer, allocatable,dimension(:,:) ::  k0idum


      call mpi_barrier(MPI_COMM_K,ierr)

      nr3uL2 = nr3L2/2+1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call system_ccfft(1,1,nr1L2,1.0d0,0,0,tabnr1lfwL2,
     &     0,0,ntabnr1L2,0)
      call system_ccfft(1,1,nr2L2,1.0d0,0,0,tabnr2lfwL2,
     &     0,0,ntabnr2L2,0)
      call system_scfft(1,1,nr3L2,1.0d0,0,0,tabnr3lrcL2,
     &     0,0,ntabnr3L2,0)

      call system_ccfft(1,-1,nr1L2,1.0d0,0,0,tabnr1linL2,
     &     0,0,ntabnr1L2,0)
      call system_ccfft(1,-1,nr2L2,1.0d0,0,0,tabnr2linL2,
     &     0,0,ntabnr2L2,0)
      call system_csfft(1,-1,nr3L2,1.0d0,0,0,tabnr3lcrL2,
     &     0,0,ntabnr3L2,0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mgz2L2 = 0
c
       do jcol = 1,ncol2L2(inode)
         ngz = izcol2L2(jcol,inode)
         if(ngz.gt.mgz2L2) mgz2L2 = ngz
       enddo
c
       call global_maxi(mgz2L2)
c
       call mpi_barrier(MPI_COMM_K,ierr)
 
       call mpi_allreduce(mgz2L2,mgz2L2tmp,1,mpi_integer,mpi_max,
     &                    MPI_COMM_K, ierr)
       mgz2L2 = mgz2L2tmp
       call mpi_barrier(MPI_COMM_K,ierr)

       if(inode_tot.eq.1) then
         write(*,*)
         write(*,*) 'minimum possible value of (pot)mgz2L2=',mgz2L2
       endif

 101   if(mod(nr1L2*mgz2L2,nnodes).ne.0) then
         mgz2L2=mgz2L2+1
         goto 101
       endif

       if(inode_tot.eq.1) then
         write(*,*) 'optimal(and used)value of(pot) mgz2L2=',mgz2L2
         if(mgz2L2.gt.(1+nr3L2/2)) then
           write(*,*) 'load balancing failed, change n3L2 '
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
       if(mod(nr1L2*mgz2L2,nnodes).ne.0) then
        write(*,*) 'nr1L2*mgz2L not multiple of(pot)nnodes, FFT stop'
        call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       call mpi_barrier(MPI_COMM_K,ierr)

c  choose data layout for  slab  for y dir FFT
c  each PE will take (nr1*mgz2)/nnodes consecutive
c  columns
c
c  first slab goes from 1 to mgz2 in z direction
c

       ncoly2L2 = (nr1L2*mgz2L2)/nnodes
       ico = 0
       inodec = 1

       do ix = 1,nr1L2
        do iz = 1,mgz2L2
         ico = ico + 1
         izcol_y2L2(ico,inodec) = iz
         ixcol_y2L2(ico,inodec) = ix
         ipe_y2L2(ix,iz) = inodec
         ibca_y2L2(ix,iz) = (ico-1)*2*nr2L2
         if(ico.eq.ncoly2L2) then
          inodec = inodec + 1
          ico = 0
         endif !(ico.eq.ncoly2L2)
        enddo !iz
       enddo !ix

c calculate starting addresses for chunks and also
c their width in the z direction
c each proc has no more than   nofchks   chunks
c these chunks will be transposed locally
c before being sent in z strips to the appropriate
c processor and memory location  before z dir fft


      allocate(izchb2L2(ncoly2L2))
      allocate(ixch2L2(ncoly2L2))
      allocate(ichw2L2(ncoly2L2))

      izchb2L2(1) = izcol_y2L2(1,inode)
       ixch2L2(1) = ixcol_y2L2(1,inode)
       ix = ixcol_y2L2(1,inode)
       ichunk2L2 = 1
       ichw2L2=0
       do ic = 1,ncoly2L2
         iz = izcol_y2L2(ic,inode)
         if(ix.ne.ixcol_y2L2(ic,inode)) then
          ix = ixcol_y2L2(ic,inode)
          ichunk2L2 = ichunk2L2 + 1
          izchb2L2(ichunk2L2) = izcol_y2L2(ic,inode)
          ixch2L2(ichunk2L2) = ix
         endif
         ichw2L2(ichunk2L2) = ichw2L2(ichunk2L2)+1
       enddo
c
c check no more than nofchks chunks on each PE
c
       if(ichunk2L2.gt.ncoly2L2) then
       write(*,*) 'More than ncoly2L2 chunks on PEs '
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

c choose layout for final FFT in z direction
c this will be on the whole system ie nr1*nr2 columns
c divided by the number of processors

       ncolz2L2 = nr1L2*nr2L2/nnodes

       ico = 0
       inodec = 1

       do ix = 1,nr1L2
        do iy = 1,nr2L2
         ico = ico + 1
         iycol_z2L2(ico,inodec) = iy
         ixcol_z2L2(ico,inodec) = ix
         ipe_z2L2(iy,ix) = inodec
         ibca_z2L2(iy,ix) = (ico-1)*2*nr3uL2
         if(ico.eq.ncolz2L2) then
          inodec = inodec + 1
          ico = 0
         endif !(ico.eq.ncolz2L2)
        enddo !iy
       enddo !ix
c
c calculate ivecadd for the first scatterred put
c in the fft's
c
       icount2L2 = 0

       do jcol = 1,ncol2L2(inode)
        ngy = iycol2L2(jcol,inode)
        ngz = izcol2L2(jcol,inode)
        ngyadd = 1+(ngy-1)*2
        do ix = 1,nr1L2
          itnode = ipe_y2L2(ix,ngz)
          itaradd = ngyadd + ibca_y2L2(ix,ngz) - 1
          idum = icount2L2(itnode)+1

          ivecadd2L2(idum,itnode) = itaradd
          ivecadd2L2(idum+1,itnode) = itaradd +1

          icount2L2(itnode) = icount2L2(itnode)+2

        enddo !ix
       enddo !i
c
c   check that we have allocated enough memory for vec,ivecadd
c
       idum = mnr2xL2/(nnodes)
       do i = 1,nnodes
        if(icount2L2(i).gt.idum) then
         write(*,*) 'not enough memory allocated for vec2L2,ivecadd2L2'
         write(*,*) ' array dim = ',idum,' min dim = ',icount2L2(i)
         call mpi_abort(MPI_COMM_WORLD,ierr)
        endif
       enddo

       call mpi_barrier(MPI_COMM_K,ierr)
c
c put ivecadd into 1 dim array  ivdum
c
      idum = 1
      do i = 1,nnodes
        ivpacn1lL2(i) = icount2L2(i)
        do j = 1,icount2L2(i)
          ivdum(idum) = ivecadd2L2(j,i)
          idum = idum + 1
        enddo
      enddo

      idum = idum -1

      if(idum.gt.mnr2xL2) then
        write(*,*) " ivpac1lL,ivdum arrays not large enough "
        call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

c
c get cumulative address for ivpac1
c
      ivpacn_cum(1) = 1
      do i= 1,nnodes-1
        ivpacn_cum(i+1)= ivpacn_cum(i)+ivpacn1lL2(i)
      enddo

      idum = idum + ivpacn1lL2(nnodes) -1


      if(idum.gt.mnr2xL2) then
        write(*,*) " ivpac1lL2,ivdum arrays not large enough "
        call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
c  
c set up gather indexes for mpi first send
c
      ivpacn1lL2  = 0

      do jcol = 1,ncol2L2(inode)
         ilocadd = 1+(jcol-1)*(2*nr1L2)
         ngy = iycol2L2(jcol,inode)
         ngz = izcol2L2(jcol,inode)
         ngyadd = 1+(ngy-1)*2
         do ix = 1,nr1L2
            itnode = ipe_y2L2(ix,ngz)
            idum = ivpacn_cum(itnode)+ivpacn1lL2(itnode)

            ivpac1lL2(idum) = ilocadd
            ivpac1lL2(idum+1) = ilocadd + 1

            ivpacn1lL2(itnode) = ivpacn1lL2(itnode)+2

            ilocadd = ilocadd + 2

         enddo                  !ix
      enddo                     !i
c
c communicate scatter indexes for mpi first send
c   ie transformation of cylinder
c

       do i = 1,nnodes
        call mpi_isend(ivpacn1lL2(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
       enddo

       do i = 1,nnodes
        call mpi_recv(ivunpn1lL2(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       idum = 1
       do i = 1,nnodes
        call mpi_isend(ivdum(idum),ivpacn1lL2(i),mpi_integer,
     &             i-1,inode,MPI_COMM_K,ireq(i),ierr)
        idum = idum + ivpacn1lL2(i)
       enddo

       idum = 1
       do i = 1,nnodes
        call mpi_recv(ivunp1lL2(idum),ivunpn1lL2(i),mpi_integer,
     &              i-1,i,MPI_COMM_K,mpistatus,ierr)
        idum = idum + ivunpn1lL2(i)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       ivunp1lL2 = ivunp1lL2 + 1
c
c set up packing index for second transformation of slice
c
      ivpacn2lL2 = 0

      do ic = 1,ichunk2L2
         izb = izchb2L2(ic)
         ix = ixch2L2(ic)
         iw = ichw2L2(ic)
         ib = ibca_y2L2(ix,izb)
         do iy = 1,nr2L2
            itnode = ipe_z2L2(iy,ix)
            itaradd = ibca_z2L2(iy,ix)+(izb-1)+1
            itar = itaradd
            isc = ib+1
            do ii = 1,iw
              ivpacn2lL2(itnode) = ivpacn2lL2(itnode)+2
              itar = itar + 2
              isc = isc + nr2L2*2
            enddo
            ib = ib + 2
         enddo
      enddo
c
c get cumulative address for ivpac2
c
      ivpacn_cum(1) = 1
      do i= 1,nnodes-1
        ivpacn_cum(i+1)= ivpacn_cum(i)+ivpacn2lL2(i)
      enddo

      ivpacn2lL2 = 0

      do ic = 1,ichunk2L2
         izb = izchb2L2(ic)
         ix = ixch2L2(ic)
         iw = ichw2L2(ic)
         ib = ibca_y2L2(ix,izb)
         do iy = 1,nr2L2
            itnode = ipe_z2L2(iy,ix)
            itaradd = ibca_z2L2(iy,ix)+(izb-1)*2+1
            itar = itaradd
            isc = ib+1
            do ii = 1,iw
              idum = ivpacn_cum(itnode)+ivpacn2lL2(itnode)
              ivpac2lL2(idum) = isc
              ivpac2lL2(idum+1) = isc+1
              ivdum(idum) = itar
              ivdum(idum+1) = itar+1
              ivpacn2lL2(itnode) = ivpacn2lL2(itnode)+2
              itar = itar + 2
              isc = isc + nr2L2*2
            enddo
            ib = ib + 2
         enddo
      enddo

c
c communicate scatter indexes for mpi second send
c
       do i = 1,nnodes
        call mpi_isend(ivpacn2lL2(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
       enddo

       do i = 1,nnodes
        call mpi_recv(ivunpn2lL2(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       idum = 1
       do i = 1,nnodes
        call mpi_isend(ivdum(idum),ivpacn2lL2(i),mpi_integer,
     &              i-1,inode,MPI_COMM_K,ireq(i),ierr)
        idum = idum + ivpacn2lL2(i)
       enddo

       idum = 1
       do i = 1,nnodes
        call mpi_recv(ivunp2lL2(idum),ivunpn2lL2(i),mpi_integer,
     &             i-1,i,MPI_COMM_K,mpistatus,ierr)
        idum = idum + ivunpn2lL2(i)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

      call mpi_barrier(MPI_COMM_K,ierr)
c
c mpi indexes for filling in of half sphere k.eq.0 for ffts
c

      k0npacL2 = 0

      do ig = 1, ngtotnod2L2(inode)

         indepg=(jjcol2L2(n2p2_nL2(ig),n3p2_nL2(ig))-1)*nr1L2+
     &        n1p2_nL2(ig)

         if(n3p2_nL2(ig).eq.1) then
            n1_inv = mod(n1L2-n1p2_nL2(ig)+1,nr1L2)+1
            n2_inv = mod(n2L2-n2p2_nL2(ig)+1,nr2L2)+1
            indepg_d = (jjcol2L2(n2_inv,1)-1)*n1L2+n1_inv
            jjnode_dum = jjnode2L2(n2_inv,1)
            k0npacL2(jjnode_dum)=k0npacL2(jjnode_dum)+1
ccccc            k0idum(k0npacL2(jjnode_dum),jjnode_dum)=2*indepg_d-1
            k0npacL2(jjnode_dum)=k0npacL2(jjnode_dum)+1
cccccc           k0idum(k0npacL2(jjnode_dum),jjnode_dum)=2*indepg_d
         endif

      enddo


      do i = 1,nnodes
        call mpi_isend(k0npacL2(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
      enddo

      do i = 1,nnodes
        call mpi_recv(k0nunpL2(i),1,mpi_integer,i-1,i,
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
      if(k0npacL2(i).gt.k0npac_max) k0npac_max=k0npacL2(i)
      if(k0nunpL2(i).gt.k0nunp_max) k0nunp_max=k0nunpL2(i)
      enddo


      allocate(k0idum(k0npac_max,nnodes))
      call fft_allocate2L2(k0nunp_max,nnodes)  ! allocate k0iunp(k0nunp_max,nnodes)
*******************************************************
ccccc recalculate k0npac, and this time: k0idum 

      k0npacL2 = 0

      do ig = 1, ngtotnod2L2(inode)

         indepg=(jjcol2L2(n2p2_nL2(ig),n3p2_nL2(ig))-1)*nr1L2+
     &         n1p2_nL2(ig)

         if(n3p2_nL2(ig).eq.1) then
            n1_inv = mod(n1L2-n1p2_nL2(ig)+1,nr1L2)+1
            n2_inv = mod(n2L2-n2p2_nL2(ig)+1,nr2L2)+1
            indepg_d = (jjcol2L2(n2_inv,1)-1)*n1L2+n1_inv
            jjnode_dum = jjnode2L2(n2_inv,1)
            k0npacL2(jjnode_dum)=k0npacL2(jjnode_dum)+1
            k0idum(k0npacL2(jjnode_dum),jjnode_dum)=2*indepg_d-1
            k0npacL2(jjnode_dum)=k0npacL2(jjnode_dum)+1
            k0idum(k0npacL2(jjnode_dum),jjnode_dum)=2*indepg_d
         endif

      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_K,ierr)


      do i = 1,nnodes
        call mpi_isend(k0idum(1,i),k0npacL2(i),mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
      enddo

      do i = 1,nnodes
        call mpi_recv(k0iunpL2(1,i),k0nunpL2(i),mpi_integer,i-1,i,
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
      isum=isum+ivpacn2lL2(i)
        do j=1,k0nunpL2(i)
        if(k0iunpL2(j,i).gt.imax) imax=k0iunpL2(j,i)
        enddo
      enddo

      if(isum.gt.mr_nL2.or.imax.gt.mr_nL2) then
       write(6,*) "increase mrL2, due to imbalance, stop", 
     & inode,isum,imax,mr_nL2
cccccc This is probably nature, due to the load imbalance, some 
cccccc node might have slightly more FFT columns than other nodes, 
cccccc  then isum,imax > (n1L2*n2L2*(n3L2+2)/nnodes)

      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif



      return
      end
