      subroutine fftprep_real2(nr1,nr2,nr3,mr_n)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc written by Andew Canning, 2001
cccccccccccccccccccccccccccccccccccccccccccccccccc

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
     c        jcol,ngy,ngyadd,ngz,igdum,nr1,nr2,nr3,nr3u,ierr,
     c        ilocadd,izb,iw,ib,itar,isc,ii,mgz2tmp,isum,imax,
     c        mr_n,n1_inv,n2_inv,indepg_d,jjnode_dum,indepg,ig,
     c        k0npac_max,k0nunp_max

c
c mpi version arrays
c
      integer mpistatus(mpi_status_size)
      integer ivdum(mnr2x),ivpacn_cum(nnodes)
ccccc      integer k0idum(nr1*nr2,nnodes)    ! old stupid way, too much memory
      integer ireq(nnodes)

      integer, allocatable,dimension(:,:) ::  k0idum


      call mpi_barrier(MPI_COMM_K,ierr)

      nr3u = nr3/2+1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call system_ccfft(1,1,nr1,1.0d0,0,0,tabnr1lfw,0,0,ntabnr1,0)
      call system_ccfft(1,1,nr2,1.0d0,0,0,tabnr2lfw,0,0,ntabnr2,0)
      call system_scfft(1,1,nr3,1.0d0,0,0,tabnr3lrc,0,0,ntabnr3,0)

      call system_ccfft(1,-1,nr1,1.0d0,0,0,tabnr1lin,0,0,ntabnr1,0)
      call system_ccfft(1,-1,nr2,1.0d0,0,0,tabnr2lin,0,0,ntabnr2,0)
      call system_csfft(1,-1,nr3,1.0d0,0,0,tabnr3lcr,0,0,ntabnr3,0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mgz2 = 0
c
       do jcol = 1,ncol2(inode)
         ngz = izcol2(jcol,inode)
         if(ngz.gt.mgz2) mgz2 = ngz
       enddo
c
       call global_maxi(mgz2)
c
       call mpi_barrier(MPI_COMM_K,ierr)
 
       call mpi_allreduce(mgz2,mgz2tmp,1,mpi_integer,mpi_max,
     &                    MPI_COMM_K, ierr)
       mgz2 = mgz2tmp

       call mpi_barrier(MPI_COMM_K,ierr)

       if(inode_tot.eq.1) then
         write(*,*)
         write(*,*) 'minimum possible value of (pot) mgz2 = ',mgz2
       endif

 101   if(mod(nr1*mgz2,nnodes).ne.0) then
         mgz2=mgz2+1
         goto 101
       endif

       if(inode_tot.eq.1) then
         write(*,*) 'optimal (and used) value of (pot) mgz2 = ',mgz2
         if(mgz2.gt.(1+nr3/2)) then
           write(*,*) 'load balancing failed, change n3 '
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
       if(mod(nr1*mgz2,nnodes).ne.0) then
        write(*,*) 'nr1*mgz2 not multiple of (pot) nnodes, FFT stopped'
        call mpi_abort(MPI_COMM_WORLD,ierr)
       endif

       call mpi_barrier(MPI_COMM_K,ierr)

c  choose data layout for  slab  for y dir FFT
c  each PE will take (nr1*mgz2)/nnodes consecutive
c  columns
c
c  first slab goes from 1 to mgz2 in z direction
c

       ncoly2 = (nr1*mgz2)/nnodes
       ico = 0
       inodec = 1

       do ix = 1,nr1
        do iz = 1,mgz2
         ico = ico + 1
         izcol_y2(ico,inodec) = iz
         ixcol_y2(ico,inodec) = ix
         ipe_y2(ix,iz) = inodec
         ibca_y2(ix,iz) = (ico-1)*2*nr2
         if(ico.eq.ncoly2) then
          inodec = inodec + 1
          ico = 0
         endif !(ico.eq.ncoly2)
        enddo !iz
       enddo !ix

c calculate starting addresses for chunks and also
c their width in the z direction
c each proc has no more than   nofchks   chunks
c these chunks will be transposed locally
c before being sent in z strips to the appropriate
c processor and memory location  before z dir fft

      allocate(izchb2(ncoly2))
      allocate(ixch2(ncoly2))
      allocate(ichw2(ncoly2))

      izchb2(1) = izcol_y2(1,inode)
       ixch2(1) = ixcol_y2(1,inode)
       ix = ixcol_y2(1,inode)
       ichunk2 = 1
       ichw2=0
       do ic = 1,ncoly2
         iz = izcol_y2(ic,inode)
         if(ix.ne.ixcol_y2(ic,inode)) then
          ix = ixcol_y2(ic,inode)
          ichunk2 = ichunk2 + 1
          izchb2(ichunk2) = izcol_y2(ic,inode)
          ixch2(ichunk2) = ix
         endif
         ichw2(ichunk2) = ichw2(ichunk2)+1
       enddo
c
c check no more than nofchks chunks on each PE
c
       if(ichunk2.gt.ncoly2) then
       write(*,*) 'More than ncoly2 chunks on PEs '
       call mpi_abort(mpi_comm_world,ierr)
       endif

c choose layout for final FFT in z direction
c this will be on the whole system ie nr1*nr2 columns
c divided by the number of processors

       ncolz2 = nr1*nr2/nnodes

       ico = 0
       inodec = 1

       do ix = 1,nr1
        do iy = 1,nr2
         ico = ico + 1
         iycol_z2(ico,inodec) = iy
         ixcol_z2(ico,inodec) = ix
         ipe_z2(iy,ix) = inodec
         ibca_z2(iy,ix) = (ico-1)*2*nr3u
         if(ico.eq.ncolz2) then
          inodec = inodec + 1
          ico = 0
         endif !(ico.eq.ncolz2)
        enddo !iy
       enddo !ix
c
c calculate ivecadd for the first scatterred put
c in the fft's
c
       icount2 = 0

       do jcol = 1,ncol2(inode)
        ngy = iycol2(jcol,inode)
        ngz = izcol2(jcol,inode)
        ngyadd = 1+(ngy-1)*2
        do ix = 1,nr1
          itnode = ipe_y2(ix,ngz)
          itaradd = ngyadd + ibca_y2(ix,ngz) - 1
          idum = icount2(itnode)+1

          ivecadd2(idum,itnode) = itaradd
          ivecadd2(idum+1,itnode) = itaradd +1

          icount2(itnode) = icount2(itnode)+2

        enddo !ix
       enddo !i
c
c   check that we have allocated enough memory for vec,ivecadd
c
       idum = mnr2x/(nnodes)
       do i = 1,nnodes
        if(icount2(i).gt.idum) then
         write(*,*) 'ERROR notenough memory allocated for vec2,ivecadd2'
         write(*,*) ' array dim = ',idum,' min dim = ',icount2(i)
         call mpi_abort(mpi_comm_world,ierr)
        endif
       enddo

       call mpi_barrier(MPI_COMM_K,ierr)
c
c put ivecadd into 1 dim array  ivdum
c
      idum = 1
      do i = 1,nnodes
        ivpacn1l(i) = icount2(i)
        do j = 1,icount2(i)
          ivdum(idum) = ivecadd2(j,i)
          idum = idum + 1
        enddo
      enddo

      idum = idum -1

      if(idum.gt.mnr2x) then
        write(*,*) " ivpac1l,ivdum arrays not large enough "
        call mpi_abort(mpi_comm_world,ierr)
      endif

c
c get cumulative address for ivpac1
c
      ivpacn_cum(1) = 1
      do i= 1,nnodes-1
        ivpacn_cum(i+1)= ivpacn_cum(i)+ivpacn1l(i)
      enddo

      idum = idum + ivpacn1l(nnodes) -1


      if(idum.gt.mnr2x) then
        write(*,*) " ivpac1l,ivdum arrays not large enough "
        call mpi_abort(mpi_comm_world,ierr)
      endif
c  
c set up gather indexes for mpi first send
c
      ivpacn1l  = 0

      do jcol = 1,ncol2(inode)
         ilocadd = 1+(jcol-1)*(2*nr1)
         ngy = iycol2(jcol,inode)
         ngz = izcol2(jcol,inode)
         ngyadd = 1+(ngy-1)*2
         do ix = 1,nr1
            itnode = ipe_y2(ix,ngz)
            idum = ivpacn_cum(itnode)+ivpacn1l(itnode)

            ivpac1l(idum) = ilocadd
            ivpac1l(idum+1) = ilocadd + 1

            ivpacn1l(itnode) = ivpacn1l(itnode)+2

            ilocadd = ilocadd + 2

         enddo                  !ix
      enddo                     !i
c
c communicate scatter indexes for mpi first send
c   ie transformation of cylinder
c
       do i = 1,nnodes
        call mpi_isend(ivpacn1l(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
       enddo

       do i = 1,nnodes
        call mpi_recv(ivunpn1l(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       idum = 1
       do i = 1,nnodes
        call mpi_isend(ivdum(idum),ivpacn1l(i),mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
        idum = idum + ivpacn1l(i)
       enddo

       idum = 1
       do i = 1,nnodes
        call mpi_recv(ivunp1l(idum),ivunpn1l(i),mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
        idum = idum + ivunpn1l(i)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       ivunp1l = ivunp1l + 1
c
c set up packing index for second transformation of slice
c
      ivpacn2l = 0

      do ic = 1,ichunk2
         izb = izchb2(ic)
         ix = ixch2(ic)
         iw = ichw2(ic)
         ib = ibca_y2(ix,izb)
         do iy = 1,nr2
            itnode = ipe_z2(iy,ix)
            itaradd = ibca_z2(iy,ix)+(izb-1)+1
            itar = itaradd
            isc = ib+1
            do ii = 1,iw
              ivpacn2l(itnode) = ivpacn2l(itnode)+2
              itar = itar + 2
              isc = isc + nr2*2
            enddo
            ib = ib + 2
         enddo
      enddo
c
c get cumulative address for ivpac2
c
      ivpacn_cum(1) = 1
      do i= 1,nnodes-1
        ivpacn_cum(i+1)= ivpacn_cum(i)+ivpacn2l(i)
      enddo

      ivpacn2l = 0

      do ic = 1,ichunk2
         izb = izchb2(ic)
         ix = ixch2(ic)
         iw = ichw2(ic)
         ib = ibca_y2(ix,izb)
         do iy = 1,nr2
            itnode = ipe_z2(iy,ix)
            itaradd = ibca_z2(iy,ix)+(izb-1)*2+1
            itar = itaradd
            isc = ib+1
            do ii = 1,iw
              idum = ivpacn_cum(itnode)+ivpacn2l(itnode)
              ivpac2l(idum) = isc
              ivpac2l(idum+1) = isc+1
              ivdum(idum) = itar
              ivdum(idum+1) = itar+1
              ivpacn2l(itnode) = ivpacn2l(itnode)+2
              itar = itar + 2
              isc = isc + nr2*2
            enddo
            ib = ib + 2
         enddo
      enddo

c
c communicate scatter indexes for mpi second send
c
       do i = 1,nnodes
        call mpi_isend(ivpacn2l(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
       enddo

       do i = 1,nnodes
        call mpi_recv(ivunpn2l(i),1,mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

       call mpi_barrier(MPI_COMM_K,ierr)

       idum = 1
       do i = 1,nnodes
        call mpi_isend(ivdum(idum),ivpacn2l(i),mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
        idum = idum + ivpacn2l(i)
       enddo

       idum = 1
       do i = 1,nnodes
        call mpi_recv(ivunp2l(idum),ivunpn2l(i),mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
        idum = idum + ivunpn2l(i)
       enddo

       do i = 1, nnodes
          call mpi_wait(ireq(i), mpistatus, ierr)
       end do

      call mpi_barrier(MPI_COMM_K,ierr)
c
c mpi indexes for filling in of half sphere k.eq.0 for ffts
c

      k0npac = 0

      do ig = 1, ngtotnod2(inode)

         indepg=(jjcol2(n2p2_n(ig),n3p2_n(ig))-1)*nr1+n1p2_n(ig)

         if(n3p2_n(ig).eq.1) then
            n1_inv = mod(n1-n1p2_n(ig)+1,nr1)+1
            n2_inv = mod(n2-n2p2_n(ig)+1,nr2)+1
            indepg_d = (jjcol2(n2_inv,1)-1)*n1+n1_inv
            jjnode_dum = jjnode2(n2_inv,1)
            k0npac(jjnode_dum)=k0npac(jjnode_dum)+1
ccccc            k0idum(k0npac(jjnode_dum),jjnode_dum)=2*indepg_d-1
            k0npac(jjnode_dum)=k0npac(jjnode_dum)+1
cccccc            k0idum(k0npac(jjnode_dum),jjnode_dum)=2*indepg_d
         endif

      enddo


      do i = 1,nnodes
        call mpi_isend(k0npac(i),1,mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
      enddo

      do i = 1,nnodes
        call mpi_recv(k0nunp(i),1,mpi_integer,i-1,i,
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
      if(k0npac(i).gt.k0npac_max) k0npac_max=k0npac(i)
      if(k0nunp(i).gt.k0nunp_max) k0nunp_max=k0nunp(i)
      enddo

      allocate(k0idum(k0npac_max,nnodes))
      call fft_allocate2(k0nunp_max,nnodes)  ! allocate k0iunp(k0nunp_max,nnodes)
*******************************************************
ccccc recalculate k0npac, and this time: k0idum 

      k0npac = 0

      do ig = 1, ngtotnod2(inode)

         indepg=(jjcol2(n2p2_n(ig),n3p2_n(ig))-1)*nr1+n1p2_n(ig)

         if(n3p2_n(ig).eq.1) then
            n1_inv = mod(n1-n1p2_n(ig)+1,nr1)+1
            n2_inv = mod(n2-n2p2_n(ig)+1,nr2)+1
            indepg_d = (jjcol2(n2_inv,1)-1)*n1+n1_inv
            jjnode_dum = jjnode2(n2_inv,1)
            k0npac(jjnode_dum)=k0npac(jjnode_dum)+1
            k0idum(k0npac(jjnode_dum),jjnode_dum)=2*indepg_d-1
            k0npac(jjnode_dum)=k0npac(jjnode_dum)+1
            k0idum(k0npac(jjnode_dum),jjnode_dum)=2*indepg_d
         endif

      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_barrier(MPI_COMM_K,ierr)


      do i = 1,nnodes
        call mpi_isend(k0idum(1,i),k0npac(i),mpi_integer,i-1,inode,
     &                 MPI_COMM_K,ireq(i),ierr)
      enddo

      do i = 1,nnodes
        call mpi_recv(k0iunp(1,i),k0nunp(i),mpi_integer,i-1,i,
     &                MPI_COMM_K,mpistatus,ierr)
      enddo

      do i = 1, nnodes
         call mpi_wait(ireq(i), mpistatus, ierr)
      end do

      call mpi_barrier(MPI_COMM_K,ierr)

ccccccccccccccccccccccccccccccccc
cccccc  safety check
      isum=0
      imax=0
      do i=1,nnodes
      isum=isum+ivpacn2l(i)
        do j=1,k0nunp(i)
        if(k0iunp(j,i).gt.imax) imax=k0iunp(j,i)
        enddo
      enddo

      if(isum.gt.mr_n.or.imax.gt.mr_n) then
       write(6,*) "increase mr_n, due to imbalance, stop",
     & inode,isum,imax,mr_n
cccccc This is probably nature, due to the load imbalance, some
cccccc node might have slightly more FFT columns than other nodes,
cccccc  then isum,imax > (n1*n2*(n3+2)/nnodes)
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif


      return
      end
