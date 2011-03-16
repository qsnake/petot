      subroutine gen_G2L_real()

********************************************
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
      nh1L=(n1L+1)/2+1

      ngcountL=0

*********************************************
***   ngcountL=1 correspond to G=0
***   This special point is used in other subroutines
*********************************************

      iorg2L = 0
      jcolL = 0
      ncol2L = 0
      ngcol2L=0
      ngcol02L=0
      ngtotnod2L=0
      jjcol2L = 0
      jjnode2L = 0

      do 10 k3=1,n3L
         do 10 j3=1,n2L
            jcolL = jcolL + 1
            iccL = 0
            do 5 i3=1,n1L

c     
c     parallel version does not allow for remapping inside 
c     half sphere  
c     keep half sphere cut by plane k1.eq.0 
c     neccessary for parallel real to complex fft's 
c     
c     

               i1=i3-1
               if(i1.gt.n1L/2) i1=i1-n1L
               j1=j3-1
               if(j1.gt.n2L/2) j1=j1-n2L
               k1=k3-1
               if(k1.gt.n3L/2) k1=k1-n3L


               akk=0.5d0*(2*pi)**2*(
     &              (ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)**2
     &              +(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)**2
     &              +(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)**2)
               
               if(akk.gt.Ecut2L) goto 5
               if(k1.lt.0) goto 5
c     
c     list all columns including ones that contain no g vecs
c     on k1.eq.0 plane, when we do ffts' these will be 
c     filled
c     
               ngcol02L(jcolL) = ngcol02L(jcolL)+1
               if(k1.eq.0.and.j1.lt.0) goto 5
               if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 5

*****************************************
***   keep only half the G sphere for real wavefunction
*****************************************

               ngcol2L(jcolL)=ngcol2L(jcolL)+1

               ngcountL=ngcountL+1

 5          continue            ! i3 
            
c     throw away cols with no g vecs in half sphere
            if(ngcol02L(jcolL).eq.0) then  
               jcolL = jcolL -1  
            else
c     list y and z position of cols
               ngycol2L(jcolL) = j3
               ngzcol2L(jcolL) = k3
            endif

 10      continue               ! j3,k3

c     total number of g vec
         ng2L=ngcountL

         ngcountL = 0
c     
c     check  Ecut is contained in box n1*n2*n3
c     by counting g vectors in box twice the size
c     
         do k3= 1,n3L*2
            do j3= 1,n2L*2
               do i3= 1,n1L*2

                  i1=i3-1
                  if(i1.gt.n1L) i1=i1-2*n1L
                  j1=j3-1
                  if(j1.gt.n2L) j1=j1-2*n2L
                  k1=k3-1
                  if(k1.gt.n3L) k1=k1-2*n3L

                  akk=0.5d0*(2*pi)**2*(
     &                 (ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)**2
     &                 +(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)**2
     &                 +(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)**2)

                  if(akk.gt.Ecut2L) goto 15
                  if(k1.lt.0) goto 15
                  if(k1.eq.0.and.j1.lt.0) goto 15
                  if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 15

                  ngcountL = ngcountL + 1

 15               continue

               enddo
            enddo
         enddo

         if(ngcountL.ne.ng2L.and.inode.eq.1) then
            write(*,*) " Ecut2 not contained in grid n1*n2*n3" 
           call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

c     
c     sort columns using ngcol2 in ascending order of length 
c     
         call heapsort(ngcol2L,jcolL,ngycol2L,ngzcol2L,ngcol02L)
c     
c     Distribute cols. between the nodes starting with longest
c     new column given to node with least g vectors 
c     and construct info. for column distr. for ffts etc.
c     

ccccccccccc if a column is defined in smaller fft, gen_G2_real, 
ccccccccccc Then place this column on the same node, and make 
ccccccccccc a maping connection between them 
ccccccccccc jjnode2(ngy,ngz) is from the small fft, gen_G2_real. 

         ig = 0
         icolL = 0
         nh1=n1/2+1
         nh2=n2/2+1
         nh3=n3/2+1
         map_StoL=0
         map_LtoS=0

         do jp=jcolL,1,-1

            ngyL=ngycol2L(jp)
            ngzL=ngzcol2L(jp)

         jnodeL=0
         if((ngyL.le.nh2.or.ngyL.ge.(n2L-n2+nh2+1)).and.
     &      (ngzL.le.nh3.or.ngzL.ge.(n3L-n3+nh3+1)))  then     ! column within the small fft box
            if(ngyL.le.nh2) ngy=ngyL
            if(ngyL.ge.(n2L-n2+nh2+1)) ngy=ngyL+n2-n2L
            if(ngzL.le.nh3) ngz=ngzL
            if(ngzL.ge.(n3L-n3+nh3+1)) ngz=ngzL+n3-n3L
            jnodeL=jjnode2(ngy,ngz)
          endif
ccccc at this stage, if jnodeL.ne.0, jjnode2 has been assigned in small fft box
           
            imap=1
            if(jnodeL.eq.0) then     ! this column is not assigned yet in small fft box
            imap=0
            jnodeL=fmin(ngtotnod2L,nnodes)
            if(ngcol2L(jp).eq.0) jnodeL = mod(jp,nnodes)+1
            endif

ccccccccccc test
           if(jnodeL.gt.nnodes) write(6,*) "test1", inode,jnodeL,nnodes


            jjnode2L(ngyL,ngzL)=jnodeL
            ncol2L(jnodeL) = ncol2L(jnodeL) + 1
            if(ncol2L(jnodeL).gt.ncolxL) then 
              write(*,*) "ncolxL  too small for potential sphere "
              stop
            endif
            iycol2L(ncol2L(jnodeL),jnodeL) = ngyL
            izcol2L(ncol2L(jnodeL),jnodeL) = ngzL
            jjcol2L(ngyL,ngzL) = ncol2L(jnodeL)
            ngtotnod2L(jnodeL) = ngtotnod2L(jnodeL) + ngcol2L(jp)

c     
c     set up local info for each node 
c     
            if(jnodeL.eq.inode) then 
               icolL = icolL + 1

               igstar_jjcol2L(ncol2L(jnodeL))=ig+1

               do ii = 1,n1L

                  i1=ii-1
                  if(i1.gt.n1L/2) i1=i1-n1L
                  j1=ngyL-1
                  if(j1.gt.n2L/2) j1=j1-n2L
                  k1=ngzL-1
                  if(k1.gt.n3L/2) k1=k1-n3L

                  akxt=(2*pi)*(ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)
                  akyt=(2*pi)*(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)
                  akzt=(2*pi)*(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)

                  akkt=0.5d0*(akxt**2+akyt**2+akzt**2)

                  if(akkt.gt.Ecut2L) goto 55
                  if(k1.lt.0) goto 55
                  if(k1.eq.0.and.j1.lt.0) goto 55
                  if(k1.eq.0.and.j1.eq.0.and.i1.lt.0) goto 55

                  ig = 1 + ig 

                  if(ig.gt.mr_nL/2) then 
                    write(*,*) " mr_nL/2 too small, must be > ",ig
                  endif


                  n2p2_nL(ig) = ngyL
                  n3p2_nL(ig) = ngzL 

                  n1p2_nL(ig) = ii  

                  gkk2_nL(ig) = akkt
                  gkx2_nL(ig) = akxt
                  gky2_nL(ig) = akyt
                  gkz2_nL(ig) = akzt

                  map_LtoS(ig)=0      

                  if(imap.eq.1) then   ! find the possible mapping point
                   jjcol_S=jjcol2(ngy,ngz)   ! ngy,ngz has been defined above for imap.eq.1
                   do ig_s=igstar_jjcol2(jjcol_S),
     &                                igfin_jjcol2(jjcol_S)
                
                   ii_s=ii
                   if(ii_s.le.nh1.or.ii_s.ge.(n1L-n1+nh1+1)) then
                   if(ii_s.le.nh1) ii_s=ii
                   if(ii_s.ge.(n1L-n1+nh1+1)) ii_s=ii+n1-n1L
                   if(ii_s.eq.n1p2_n(ig_s))  then
                   map_LtoS(ig)=ig_s
                        if(map_StoL(ig_s).ne.0) then
                 write(6,*) "double mapping in gen_G2L_real, stop"
                 write(6,*) inode, ii,ii_s, ig, ig_s, map_StoL(ig_s),
     &       jjcol_S, n1, nh1, n1L,n1p2_n(ig_s) 
                        call mpi_abort(MPI_COMM_WORLD,ierr)
                        endif
                   map_StoL(ig_s)=ig
                   endif
                   endif
                   enddo
                  endif



c     store node number and location of origin

                  if(akkt.eq.0.0) then
                     iorg2L(1) = inode
                     iorg2L(2) = ig 
                  endif

 55               continue

               enddo            !ii

               igfin_jjcol2L(ncol2L(jnodeL))=ig
               
            endif               ! jnodeL.eq.inode

         enddo           ! do jp, the column

c     number of g vectors on my processor = ng_n

         ng2_nL = ig

         if(ng2_nL.ne.ngtotnod2L(inode)) then
         write(6,*) "ng2_nL.ne.ngtotnod2L, stop", 
     &           ng2_nL, ngtotnod2L(inode)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

cccccccc test the mapping
         do ig=1,ngtotnod2(inode)
         if(map_StoL(ig).eq.0) then
         write(6,*) "map_StoL.eq.0, stop", ig
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif
         enddo
ccccccccccccccccccccccccccccccccc

c     
c     write out load balancing info 
c     
         if(inode_tot.eq.1) then
            write(*,*) " load balancing info for Large potential sphere"
            write(*,*) "pe.  g's. cols. "
            do i = 1,nnodes
               write(*,*) i,ngtotnod2L(i),ncol2L(i)
            enddo
         endif
         

         return
         end
      
