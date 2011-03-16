      subroutine gen_G2_real()

*******************************************************
***  Written by Lin-Wang Wang, 2001
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
      nh1=(n1+1)/2+1

      ngcount=0

*********************************************
***   ngcount=1 correspond to G=0
***   This special point is used in other subroutines
*********************************************

      iorg2 = 0
      jcol = 0
      ncol2 = 0
      ngcol2=0
      ngcol02=0
      ngtotnod2=0
      jjcol2 = 0
      jjnode2 = 0

      do 10 k3=1,n3
         do 10 j3=1,n2
            jcol = jcol + 1
            icc = 0
            do 5 i3=1,n1

c     
c     parallel version does not allow for remapping inside 
c     half sphere  
c     keep half sphere cut by plane k1.eq.0 
c     neccessary for parallel real to complex fft's 
c     
c     

               i1=i3-1
               if(i1.gt.n1/2) i1=i1-n1
               j1=j3-1
               if(j1.gt.n2/2) j1=j1-n2
               k1=k3-1
               if(k1.gt.n3/2) k1=k1-n3


               akk=0.5d0*(2*pi)**2*(
     &              (ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)**2
     &              +(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)**2
     &              +(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)**2)
               
               if(akk.gt.Ecut2) goto 5
               if(k1.lt.0) goto 5
c     
c     list all columns including ones that contain no g vecs
c     on k1.eq.0 plane, when we do ffts' these will be 
c     filled
c     
               ngcol02(jcol) = ngcol02(jcol)+1
               if(k1.eq.0.and.j1.lt.0) goto 5
               if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 5

*****************************************
***   keep only half the G sphere for real wavefunction
*****************************************

               ngcol2(jcol)=ngcol2(jcol)+1

               ngcount=ngcount+1

 5          continue            ! i3 
            
c     throw away cols with no g vecs in half sphere
            if(ngcol02(jcol).eq.0) then  
               jcol = jcol -1  
            else
c     list y and z position of cols
               ngycol2(jcol) = j3
               ngzcol2(jcol) = k3
            endif

 10      continue               ! j3,k3

c     total number of g vec
         ng2=ngcount

         ngcount = 0
c     
c     check  Ecut is contained in box n1*n2*n3
c     by counting g vectors in box twice the size
c     
         do k3= 1,n3*2
            do j3= 1,n2*2
               do i3= 1,n1*2

                  i1=i3-1
                  if(i1.gt.n1) i1=i1-2*n1
                  j1=j3-1
                  if(j1.gt.n2) j1=j1-2*n2
                  k1=k3-1
                  if(k1.gt.n3) k1=k1-2*n3

                  akk=0.5d0*(2*pi)**2*(
     &                 (ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)**2
     &                 +(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)**2
     &                 +(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)**2)

                  if(akk.gt.Ecut2) goto 15
                  if(k1.lt.0) goto 15
                  if(k1.eq.0.and.j1.lt.0) goto 15
                  if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 15

                  ngcount = ngcount + 1

 15               continue

               enddo
            enddo
         enddo

         if(ngcount.ne.ng2.and.inode.eq.1) then
            write(*,*) " Ecut2 not contained in grid n1*n2*n3" 
           call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

c     
c     sort columns using ngcol2 in ascending order of length 
c     
         call heapsort(ngcol2,jcol,ngycol2,ngzcol2,ngcol02)
c     
c     Distribute cols. between the nodes starting with longest
c     new column given to node with least g vectors 
c     and construct info. for column distr. for ffts etc.
c     

         ig = 0
         icol = 0

         do jp=jcol,1,-1

            jnode=fmin(ngtotnod2,nnodes)
            if(ngcol2(jp).eq.0) jnode = mod(jp,nnodes)+1
            ngy=ngycol2(jp)
            ngz=ngzcol2(jp)
            jjnode2(ngy,ngz)=jnode
            ncol2(jnode) = ncol2(jnode) + 1
            if(ncol2(jnode).gt.ncolx) then 
              write(*,*) "ncolx  too small for potential sphere "
              stop
            endif
            iycol2(ncol2(jnode),jnode) = ngy
            izcol2(ncol2(jnode),jnode) = ngz
            jjcol2(ngy,ngz) = ncol2(jnode)
            ngtotnod2(jnode) = ngtotnod2(jnode) + ngcol2(jp)

c     
c     set up local info for each node 
c     
            if(jnode.eq.inode) then 
               icol = icol + 1

               igstar_jjcol2(ncol2(jnode))=ig+1

               do ii = 1,n1

                  i1=ii-1
                  if(i1.gt.n1/2) i1=i1-n1
                  j1=ngy-1
                  if(j1.gt.n2/2) j1=j1-n2
                  k1=ngz-1
                  if(k1.gt.n3/2) k1=k1-n3

                  akxt=(2*pi)*(ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)
                  akyt=(2*pi)*(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)
                  akzt=(2*pi)*(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)

                  akkt=0.5d0*(akxt**2+akyt**2+akzt**2)

                  if(akkt.gt.Ecut2) goto 55
                  if(k1.lt.0) goto 55
                  if(k1.eq.0.and.j1.lt.0) goto 55
                  if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 55

                  ig = 1 + ig 

                  if(ig.gt.mr_n/2) then 
                    write(*,*) " mr_n/2 too small, must be > ",ig
                  endif

                  n2p2_n(ig) = ngy
                  n3p2_n(ig) = ngz 

                  n1p2_n(ig) = ii  

                  gkk2_n(ig) = akkt
                  gkx2_n(ig) = akxt
                  gky2_n(ig) = akyt
                  gkz2_n(ig) = akzt


c     store node number and location of origin

                  if(akkt.eq.0.0) then
                     iorg2(1) = inode
                     iorg2(2) = ig 
                  endif

 55               continue

               enddo            !ii

               igfin_jjcol2(ncol2(jnode))=ig
               
            endif               ! jnode.eq.inode

         enddo

c     number of g vectors on my processor = ng_n

         ng2_n = ig

         if(ng2_n.ne.ngtotnod2(inode)) then
         write(6,*) "ng2_n.ne.ngtotnod2, stop", ng2_n, ngtotnod2(inode)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif
 

c     
c     write out load balancing info 
c     
         if(inode_tot.eq.1) then
            write(*,*) " load balancing info for potential sphere"
            write(*,*) "pe.  g's. cols. "
            do i = 1,nnodes
               write(*,*) i,ngtotnod2(i),ncol2(i)
            enddo
         endif
         

         return
         end
      
