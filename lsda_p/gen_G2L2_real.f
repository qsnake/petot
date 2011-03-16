      subroutine gen_G2L2_real()

************************************************
****  Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************


***   this subroutine is for real wavefunction program, escan_real.f only
***   calculates gkk_n, wg_n, length squared and weight factor for 
***   g vectors on my node also calculates index matrices  
***   wg(g): smooth cutoff truncation used in d3fft_real.f
      
      use data
      use load_data
      use fft_data
      
      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'

      real*8 tmp(3)

      integer fmin

*****************************************
*     *  \sum_i AL(i,j1)*ALI(i,j2)= \delta_j1,j2
*     *  2*pi*ALI(i,j) is the jth G vector unit in (i)x,y,z components
*****************************************

***********************************
*     select the G points within the Ecut2
***********************************
      pi=4*datan(1.d0)
      nh1L2=(n1L2+1)/2+1

      ngcountL2=0

*********************************************
***   ngcountL2=1 correspond to G=0
***   This special point is used in other subroutines
*********************************************

      iorg2L2 = 0
      jcolL2 = 0
      ncol2L2 = 0
      ngcol2L2=0
      ngcol02L2=0
      ngtotnod2L2=0
      jjcol2L2 = 0
      jjnode2L2 = 0

      do 10 k3=1,n3L2
         do 10 j3=1,n2L2
            jcolL2 = jcolL2 + 1
            iccL2 = 0
            do 5 i3=1,n1L2

c     
c     parallel version does not allow for remapping inside 
c     half sphere  
c     keep half sphere cut by plane k1.eq.0 
c     neccessary for parallel real to complex fft's 
c     
c     

               i1=i3-1
               if(i1.gt.n1L2/2) i1=i1-n1L2
               j1=j3-1
               if(j1.gt.n2L2/2) j1=j1-n2L2
               k1=k3-1
               if(k1.gt.n3L2/2) k1=k1-n3L2


               akk=0.5d0*(2*pi)**2*(
     &            (ALI2(1,1)*i1+ALI2(1,2)*j1+ALI2(1,3)*k1)**2
     &           +(ALI2(2,1)*i1+ALI2(2,2)*j1+ALI2(2,3)*k1)**2
     &           +(ALI2(3,1)*i1+ALI2(3,2)*j1+ALI2(3,3)*k1)**2)
               
cccccccc special full FFT, must use full FFT for accurate convolution
cc               if(akk.gt.Ecut2L) goto 5

               if(k1.lt.0) goto 5
c     
c     list all columns including ones that contain no g vecs
c     on k1.eq.0 plane, when we do ffts' these will be 
c     filled
c     
               ngcol02L2(jcolL2) = ngcol02L2(jcolL2)+1
               if(k1.eq.0.and.j1.lt.0) goto 5
               if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 5

*****************************************
***   keep only half the G sphere for real wavefunction
*****************************************

               ngcol2L2(jcolL2)=ngcol2L2(jcolL2)+1

               ngcountL2=ngcountL2+1

 5          continue            ! i3 
            
c     throw away cols with no g vecs in half sphere
            if(ngcol02L2(jcolL2).eq.0) then  
               jcolL2 = jcolL2 -1  
            else
c     list y and z position of cols
               ngycol2L2(jcolL2) = j3
               ngzcol2L2(jcolL2) = k3
            endif

 10      continue               ! j3,k3

c     total number of g vec
         ng2L2=ngcountL2

c     
c     check  Ecut is contained in box n1*n2*n3
c     by counting g vectors in box twice the size
c     
c         ngcountL2 = 0

c         do k3= 1,n3L2*2
c            do j3= 1,n2L2*2
c               do i3= 1,n1L2*2
 
c                  i1=i3-1
c                  if(i1.gt.n1L2) i1=i1-2*n1L2
c                  j1=j3-1
c                  if(j1.gt.n2L2) j1=j1-2*n2L2
c                  k1=k3-1
c                  if(k1.gt.n3L2) k1=k1-2*n3L2
c
c                  akk=0.5d0*(2*pi)**2*(
c     &                 (ALI2(1,1)*i1+ALI2(1,2)*j1+ALI2(1,3)*k1)**2
c     &                +(ALI2(2,1)*i1+ALI2(2,2)*j1+ALI2(2,3)*k1)**2
c     &                +(ALI2(3,1)*i1+ALI2(3,2)*j1+ALI2(3,3)*k1)**2)
c
c                  if(akk.gt.Ecut2L) goto 15
c                  if(k1.lt.0) goto 15
c                  if(k1.eq.0.and.j1.lt.0) goto 15
c                  if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 15
c
c                  ngcountL2 = ngcountL2 + 1
c
c 15               continue

c               enddo
c            enddo
c         enddo

c         if(ngcountL2.ne.ng2L2.and.inode.eq.1) then
c            write(*,*) " Ecut2L not contained in grid n1L2*n2L2*n3L2" 
c           call mpi_abort(MPI_COMM_WORLD,ierr)
c         endif


c     
c     sort columns using ngcol2 in ascending order of length 
c     
         call heapsort(ngcol2L2,jcolL2,ngycol2L2,ngzcol2L2,ngcol02L2)
c     
c     Distribute cols. between the nodes starting with longest
c     new column given to node with least g vectors 
c     and construct info. for column distr. for ffts etc.
c     

           ig=0
           icolL2=0


           do jp=jcolL2,1,-1

            jnodeL2=fmin(ngtotnod2L2,nnodes)
            if(ngcol2L2(jp).eq.0) jnodeL2 = mod(jp,nnodes)+1
            ngyL2=ngycol2L2(jp)
            ngzL2=ngzcol2L2(jp)
            jjnode2L2(ngyL2,ngzL2)=jnodeL2
            ncol2L2(jnodeL2) = ncol2L2(jnodeL2) + 1
            if(ncol2L2(jnodeL2).gt.ncolxL2) then 
              write(*,*) "ncolxL2  too small for potential sphere "
              stop
            endif
            iycol2L2(ncol2L2(jnodeL2),jnodeL2) = ngyL2
            izcol2L2(ncol2L2(jnodeL2),jnodeL2) = ngzL2
            jjcol2L2(ngyL2,ngzL2) = ncol2L2(jnodeL2)
            ngtotnod2L2(jnodeL2) = ngtotnod2L2(jnodeL2) + ngcol2L2(jp)

c     
c     set up local info for each node 
c     
            if(jnodeL2.eq.inode) then 
               icolL2 = icolL2 + 1

               igstar_jjcol2L2(ncol2L2(jnodeL2))=ig+1

               do ii = 1,n1L2

                  i1=ii-1
                  if(i1.gt.n1L2/2) i1=i1-n1L2
                  j1=ngyL2-1
                  if(j1.gt.n2L2/2) j1=j1-n2L2
                  k1=ngzL2-1
                  if(k1.gt.n3L2/2) k1=k1-n3L2

                  akxt=(2*pi)*(ALI2(1,1)*i1+ALI2(1,2)*j1+ALI2(1,3)*k1)
                  akyt=(2*pi)*(ALI2(2,1)*i1+ALI2(2,2)*j1+ALI2(2,3)*k1)
                  akzt=(2*pi)*(ALI2(3,1)*i1+ALI2(3,2)*j1+ALI2(3,3)*k1)

                  akkt=0.5d0*(akxt**2+akyt**2+akzt**2)

cccc special, a full FFT, must use full FFT for convolution
ccccc                  if(akkt.gt.Ecut2L) goto 55

                  if(k1.lt.0) goto 55
                  if(k1.eq.0.and.j1.lt.0) goto 55
                  if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 55

                  ig = 1 + ig 

                  if(ig.gt.mr_nL2/2) then 
                    write(*,*) " mr_nL2/2 too small, must be > ",ig
                  endif


                  n2p2_nL2(ig) = ngyL2
                  n3p2_nL2(ig) = ngzL2 

                  n1p2_nL2(ig) = ii  

                  gkk2_nL2(ig) = akkt
                  gkx2_nL2(ig) = akxt
                  gky2_nL2(ig) = akyt
                  gkz2_nL2(ig) = akzt


c     store node number and location of origin

                  if(akkt.eq.0.0) then
                     iorg2L2(1) = inode
                     iorg2L2(2) = ig 
                  endif

 55               continue

               enddo            !ii

               igfin_jjcol2L2(ncol2L2(jnodeL2))=ig
               
            endif               ! jnodeL.eq.inode

         enddo           ! do jp, the column

c     number of g vectors on my processor = ng_n

         ng2_nL2 = ig

         if(ng2_nL2.ne.ngtotnod2L2(inode)) then
         write(6,*) "ng2_nL2.ne.ngtotnod2L2, stop", 
     &           ng2_nL2, ngtotnod2L2(inode)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

ccccccccccccccccccccccccccccccccc

c     
c     write out load balancing info 
c     
         if(inode_tot.eq.1) then
            write(*,*) " load balancing info for Large potential sphere"
            write(*,*) "pe.  g's. cols. "
            do i = 1,nnodes
               write(*,*) i,ngtotnod2L2(i),ncol2L2(i)
            enddo
         endif
         

         return
         end
      
