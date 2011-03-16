      subroutine symmop(rhot)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************
ccccc  Now, this is in n1L,n2L,n3L

ccccc  symmetry operation of rho

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'mpif.h'
      include 'param.escan_real'

      real*8 rhot(mr_nL)
      integer iloc_st_c(500),iloc_c(500)
      complex*16 csum(500),csumtmp(500)

      real*8,allocatable,dimension(:)  :: rhotmp,workr_n

  

      allocate(rhotmp(mr_nL))
      allocate(workr_n(mr_nL))

ccccccccccccccccccccccccccccccccccccccccccc
ccc  num_gstar: the total number of gstar from all nodes (same for all nodes, global)
ccc  ng_ig_local: the total number of g in this node  =< 2*ngtotnod2L(inode)
ccc  ig_star(ng_ig_local): the index of gstar this g belongs to.
ccc  ig_lenstar(ng_ig_local): the num of equivalent g-point of the gstar=ig_star(ng_ig_local)
ccc  ig_local(ng_ig_local): the local ig index (within ngtotnod2L) of this g. 


      do i=1,nr_nL
      workr_n(i) = rhot(i)
      enddo

      call d3fft_real2L(rhot,workr_n,1,0)

      iloc=0

      ig_min_flag=1

      igstar_1=0
200   continue
      igstar_0=igstar_1+1
      if(igstar_0.gt.num_gstar) goto 201
      igstar_1=igstar_0+499
      if(igstar_1.gt.num_gstar) igstar_1=num_gstar
      
cccccc   break up the do loop, into chuncks of array(500), faster MPI

      do 100 igstar=igstar_0,igstar_1  

      igs_c=igstar-igstar_0+1

ccccc when ig_star_stmp=0, this ig_star is in this processor. 
ccccc when ig_star_stmp>1, no ig in this processor has ig_star

      if(ig_min_flag.eq.1) ig_star_stmp=ig_star_stop(iloc+1)

      if(ig_star_stmp.eq.0.and.iloc.lt.ng_ig_local) then
      ig_min_flag=1               !    this igstar is iloc+1, involved in summation
      ig_min=ig_star(iloc+1)
      else
      ig_min_flag=0               !    this igstar is not iloc+1, wait
      ig_star_stmp=ig_star_stmp-1
      endif


      csum(igs_c)=dcmplx(0.d0,0.d0)

      iloc_st_c(igs_c)=iloc

      if(ig_min_flag.eq.1) then
      do while (iloc.lt.ng_ig_local.and.ig_star(iloc+1).eq.ig_min)
       iloc=iloc+1
       ig=ig_local(iloc)
       if(ig.gt.0) then
       csum(igs_c)=csum(igs_c)+dcmplx(rhot(2*ig-1),rhot(2*ig))  
       else
       csum(igs_c)=csum(igs_c)+dcmplx(rhot(-2*ig-1),-rhot(-2*ig))  
       endif
      enddo
      endif

       iloc_c(igs_c)=iloc

100    continue

      call mpi_allreduce(csum,csumtmp,igstar_1-igstar_0+1,
     & MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_K,ierr)
      csum = csumtmp

       
      do 101 igstar=igstar_0,igstar_1   

      igs_c=igstar-igstar_0+1

      if(iloc_c(igs_c).gt.0) 
     &   csum(igs_c)=csum(igs_c)/ig_lenstar(iloc_c(igs_c)) 

cccccc put back the symmetrized rho, only for ig>0 (include the origin)

       do iloc_t=iloc_st_c(igs_c)+1,iloc_c(igs_c)
        ig=ig_local(iloc_t)
        if(ig.gt.0) then
        rhotmp(2*ig-1)=dreal(csum(igs_c))  
        rhotmp(2*ig)=dimag(csum(igs_c))   
        endif
       enddo

101    continue

       goto 200
201    continue


       if(iloc.ne.ng_ig_local) then
       write(6,*) "iloc.ne.ng_ig_local,stop,symmop",iloc,ng_ig_local
       stop
       endif

       call d3fft_real2L(rhotmp,workr_n,-1,0)

       do i=1,nr_nL
       rhot(i)=workr_n(i)
       enddo 

       deallocate(rhotmp)
       deallocate(workr_n)

       return
       end



      

      



       
