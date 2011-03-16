      subroutine gen_G_comp(kpt,iflag)

***   calculates gkk_n, wg_n, length squared and weight factor for 
***   g vectors on my node also calculates index matrices  
***   wg(g): smooth cutoff truncation used in d3fft_comp.f
      
      use data
      use load_data
      use fft_data
      
      implicit double precision (a-h,o-z)

      include 'mpif.h'
      include 'param.escan_real'

c      real*8 ALI(3,3)     ! defined in param.escan_real

      integer fmin
      integer jcol_hf
      integer kpt
      data ifirst/0/
      save ifirst


*****************************************
*     *  \sum_i AL(i,j1)*ALI(i,j2)= \delta_j1,j2
*     *  2*pi*ALI(i,j) is the jth G vector unit in (i)x,y,z components
*****************************************

***********************************
*     select the G points within the Ecut
***********************************
      pi=4*datan(1.d0)
      nh1=(n1+1)/2+1

      ngcount=0
      ngc_in=0

*********************************************
***   ngcount=1 correspond to G=0
***   This special point is used in other subroutines
*********************************************

      iorg = 0
      jcol = 0
      ncol = 0
      ngcol=0
      ngcol0=0
      ngtotnod(:,kpt)=0
      jjcol = 0
      jjnode = 0

      do 10 k3=1,n3
	  k1=k3-1
          if(k1.gt.n3/2) k1=k1-n3
         do 10 j3=1,n2
	    j1=j3-1
	    if(j1.gt.n2/2) j1=j1-n2

	    jcol=jcol+1
            icc = 0
            do 5 i3=1,n1

	    i1=i3-1
	    if(i1.gt.n1/2) i1=i1-n1


       akkx=2*pi*(ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)-akx(kpt)
       akky=2*pi*(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)-aky(kpt)
       akkz=2*pi*(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)-akz(kpt)

       akk=0.5d0*(akkx**2+akky**2+akkz**2)
               
               if(akk.gt.Ecut) goto 5

               ngcol(jcol)=ngcol(jcol)+1
	       ngcount=ngcount+1

 5          continue            ! i3 
            
c     throw away cols with no g vecs 
            if(ngcol(jcol).eq.0) then  
	       jcol=jcol-1
            else
c     list y and z position of cols
               ngycol(jcol) = j3
               ngzcol(jcol) = k3
           endif

 10      continue               ! j3,k3

c     total number of g vec
         ng(kpt)=ngcount
cccc jcol=1, corresponds to ngycol(1)=1,ngzcol(1)=1,
cccc i.e, the zero points.
cccccccccccccccccccccccccccccccccccccccccccccc
c     
c     check  Ecut is contained in box n1*n2*n3
c     by counting g vectors in box twice the size
c     
      if(iflag.eq.0) then

         ngcount = 0
         do k3= 1,n3*2
            do j3= 1,n2*2
               do i3= 1,n1*2

                  i1=i3-1
                  if(i1.gt.n1) i1=i1-2*n1
                  j1=j3-1
                  if(j1.gt.n2) j1=j1-2*n2
                  k1=k3-1
                  if(k1.gt.n3) k1=k1-2*n3

       akkx=2*pi*(ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)-akx(kpt)
       akky=2*pi*(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)-aky(kpt)
       akkz=2*pi*(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)-akz(kpt)

       akk=0.5d0*(akkx**2+akky**2+akkz**2)

                  if(akk.gt.Ecut) goto 15

                  ngcount = ngcount + 1

 15               continue

               enddo
            enddo
         enddo

         if(ngcount.ne.ng(kpt).and.inode.eq.1) then
            write(*,*) " Ecut not contained in grid n1*n2*n3",
     &    ngcount,ng(kpt)
            stop
         endif
       endif   ! iflag=0
cccccccccccccccccccccccccccccccccccccccccccccc

c     
c     sort columns using ngcol in ascending order of length 
c     the first coloumn is not sorted (not changed, so use jcol-1)
c     
         call heapsort(ngcol(2),
     & 	 jcol-1,ngycol(2),ngzcol(2),ngcol0(2))
c     
ccccccccc ngcol0 not defined and not used
c     Distribute cols. between the nodes starting with longest
c     new column given to node with least g vectors 
c     and construct info. for column distr. for ffts etc.
c     

         ig = 0
         icol = 0

         do jp=jcol,1,-1

            jnode=fmin(ngtotnod(1,kpt),nnodes)
            if(ngcol(jp).eq.0) jnode = mod(jp,nnodes)+1
            ngy=ngycol(jp)
            ngz=ngzcol(jp)
            jjnode(ngy,ngz)=jnode
            ncol(jnode) = ncol(jnode) + 1
            iycol(ncol(jnode),jnode) = ngy
            izcol(ncol(jnode),jnode) = ngz
            jjcol(ngy,ngz) = ncol(jnode)
            ngtotnod(jnode,kpt) = ngtotnod(jnode,kpt) + ngcol(jp)

c     
c     set up local info for each node 
c     
            if(jnode.eq.inode) then 

               icol = icol + 1

               do 200 ii = 1,n1 

                  i1=ii-1
                  if(i1.gt.n1/2) i1=i1-n1
                  j1=ngy-1
                  if(j1.gt.n2/2) j1=j1-n2
                  k1=ngz-1
                  if(k1.gt.n3/2) k1=k1-n3

       akkx=2*pi*(ALI(1,1)*i1+ALI(1,2)*j1+ALI(1,3)*k1)-akx(kpt)
       akky=2*pi*(ALI(2,1)*i1+ALI(2,2)*j1+ALI(2,3)*k1)-aky(kpt)
       akkz=2*pi*(ALI(3,1)*i1+ALI(3,2)*j1+ALI(3,3)*k1)-akz(kpt)

       akk=0.5d0*(akkx**2+akky**2+akkz**2)

                  if(akk.gt.Ecut) goto 200

                  ig = 1 + ig 
                  n2p_n(ig) = ngy
                  n3p_n(ig) = ngz 
                  n1p_n(ig) = ii  
                  gkk_n(ig,kpt) = akk 
                  gkx_n(ig,kpt) = akkx
                  gky_n(ig,kpt) = akky
                  gkz_n(ig,kpt) = akkz


                  if(iflag.eq.0) then      ! wg_n need not be calc. again for iflag>0.
                  if(akk.lt.Ecut*Smth) then
                     wg_n(ig,kpt)=1.d0
                     ngc_in=ngc_in+1
                  else
                     x=(akk-Ecut*Smth)/((1-Smth)*Ecut)*pi
                     wg_n(ig,kpt)=(dcos(x)+1.d0)*0.5d0
                  endif
                  endif

c     store node number and location of origin

200           continue
               
            endif               ! jnode.eq.inode

        enddo      ! do jp  

c     number of g vectors on my processor = ng_n

         ng_n = ig

	 if(ng_n.ne.ngtotnod(inode,kpt)) then
	 write(6,*) "ng_n.ne.ngtotnod, stop", ng_n, ngtotnod(inode,kpt)
	 call mpi_abort(MPI_COMM_WORLD,ierr)
	 endif

c     
c     write out load balancing info 
c     
         if(inode.eq.1.and.iflag.eq.0.and.ifirst.eq.0) then
            write(*,*) " load balancing of the first kpt"
            write(*,*) "kpt, pe.  g's. cols. "
            do i = 1,nnodes
               write(*,*) i,ngtotnod(i,kpt),ncol(i)
            enddo
         endif

         if(ifirst.eq.0) ifirst=1
         
         return
         end
      
