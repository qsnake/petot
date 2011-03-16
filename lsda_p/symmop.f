      subroutine symmop(rhot,workr_n)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

ccccc  symmetry operation of rho

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'mpif.h'
      include 'param.escan_real'

      real*8 rhot(mr_n),workr_n(mr_n)
      complex*16 csum

      real*8,allocatable,dimension(:)  :: rhotmp

      nlarge=100000000

      allocate(rhotmp(mr_n))

      do i=1,nr_n
      workr_n(i) = rhot(i)
      enddo

      call d3fft_real2(rhot,workr_n,1,0)

      iloc=0
    
      do 100 igstar=1,num_gstar    ! one g_star at a time, cross all the processors

      if(iloc.lt.ng_ig_local) then
      ig_min=ig_star(iloc+1)
      else
      ig_min=nlarge
      endif

      call mpi_allreduce(ig_min,ig_min,1,MPI_INTEGER,
     &  MPI_MIN,MPI_COMM_WORLD,ierr)

      csum=dcmplx(0.d0,0.d0)

      iloc_st=iloc
      
      do while (iloc.lt.ng_ig_local.and.ig_star(iloc+1).eq.ig_min)
       iloc=iloc+1
       ig=ig_local(iloc)
       if(ig.gt.0) then
       csum=csum+dcmplx(rhot(2*ig-1),rhot(2*ig))
       else
       csum=csum+dcmplx(rhot(-2*ig-1),-rhot(-2*ig))
       endif
      enddo

      call mpi_allreduce(csum,csum,1,MPI_DOUBLE_COMPLEX,
     & MPI_SUM,MPI_COMM_WORLD,ierr)

       
      if(iloc.gt.0)  csum=csum/ig_lenstar(iloc)   ! fin, even if csum is not used

cccccc put back the symmetrized rho, only for ig>0 (include the origin)

       do iloc_t=iloc_st+1,iloc
        ig=ig_local(iloc_t)
        if(ig.gt.0) then
        rhotmp(2*ig-1)=dreal(csum)     ! have to use rhotmp because of the -k values
        rhotmp(2*ig)=dimag(csum)
        endif
       enddo

100    continue


       if(iloc.ne.ng_ig_local) then
       write(6,*) "iloc.ne.ng_ig_local,stop,symmop",iloc,ng_ig_local
       stop
       endif

       call d3fft_real2(rhotmp,workr_n,-1,0)

       do i=1,nr_n
       rhot(i)=workr_n(i)
       enddo 

       deallocate(rhotmp)

       return
       end



      

      



       
