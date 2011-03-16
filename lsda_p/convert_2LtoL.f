      subroutine convert_2LtoL(v_pL,v_pL2,isign)

*************************************************
****   Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************


cccccc  This program works for all icoul=1,11,12,13

ccccc  isign=1, from v_pL to v_pL2
ccccc  isign=2, from v_pL2 to v_pL

      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)

      include 'param.escan_real'
      include 'mpif.h'

      real*8 v_pL(mr_nL),v_pL2(mr_nL2)
      real*8,allocatable,dimension(:) :: workr_nL,workr_nL2 
      integer icoul
      real*8 xcoul(3)
      
      common /comcoul/icoul,xcoul

cccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc

      n1c=xcoul(1)*n1L+1
      n2c=xcoul(2)*n2L+1
      n3c=xcoul(3)*n3L+1
      if(n1c.gt.n1L) n1c=n1c-n1L
      if(n2c.gt.n2L) n2c=n2c-n2L
      if(n3c.gt.n3L) n3c=n3c-n3L

cccccc from small L to large 2L

      if(isign.eq.1) then
      allocate(workr_nL(mr_nL))
      
      v_pL2=0.d0

      do 100 j=1,nnodes

      if(j.eq.inode) then 
      workr_nL=v_pL
      endif

      call mpi_bcast(workr_nL,mr_nL,MPI_REAL8,j-1,
     &    MPI_COMM_K,ierr)

          do ii=1,nr_nL

             iisw=ii+(j-1)*nr_nL
             is=(iisw-1)/(n2L*n3L)+1
             js=((iisw-1)-(is-1)*n2L*n3L)/n3L+1
             ks=iisw-(is-1)*n2L*n3L-(js-1)*n3L

            if(is.gt.n1c) is=is+n1L2-n1L
            if(js.gt.n2c) js=js+n2L2-n2L
            if(ks.gt.n3c) ks=ks+n3L2-n3L

             iiLw  = (is-1)*n2L2*n3L2 + (js-1)*n3L2 + ks

                 iLnode= (iiLw-1)/nr_nL2
                 iiL = iiLw - iLnode*nr_nL2

                 if(iLnode.eq.inode-1) then
                 v_pL2(iiL)=workr_nL(ii)
                 endif
            enddo


       call mpi_barrier(mpi_comm_k,ierr)

100     continue          ! Loop over nodes, j

        deallocate(workr_nL)
        return
        endif      ! for isign=1

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc from large 2L to small L

      if(isign.eq.2) then
      allocate(workr_nL2(mr_nL2))
      
      v_pL=0.d0


      do 200 j=1,nnodes

      if(j.eq.inode) then 
      workr_nL2=v_pL2
      endif

      call mpi_bcast(workr_nL2,mr_nL2,MPI_REAL8,j-1,
     &    MPI_COMM_K,ierr)

          do 150  ii=1,nr_nL2

             iiLw=ii+(j-1)*nr_nL2
             is=(iiLw-1)/(n2L2*n3L2)+1
             js=((iiLw-1)-(is-1)*n2L2*n3L2)/n3L2+1
             ks=iiLw-(is-1)*n2L2*n3L2-(js-1)*n3L2

            if(is.gt.n1c.and.is.le.n1c+n1L2-n1L)  goto 150
            if(is.gt.n1c+n1L2-n1L) is=is-n1L2+n1L

            if(js.gt.n2c.and.js.le.n2c+n2L2-n2L)  goto 150
            if(js.gt.n2c+n2L2-n2L) js=js-n2L2+n2L

            if(ks.gt.n3c.and.ks.le.n3c+n3L2-n3L)  goto 150
            if(ks.gt.n3c+n3L2-n3L) ks=ks-n3L2+n3L

             iisw  = (is-1)*n2L*n3L + (js-1)*n3L + ks
                 isnode= (iisw-1)/nr_nL

                 iis = iisw - isnode*nr_nL

                 if(isnode.eq.inode-1) then
                 v_pL(iis)=workr_nL2(ii)
                 endif
150        continue


       call mpi_barrier(mpi_comm_k,ierr)

200     continue          ! Loop over nodes, j

        deallocate(workr_nL2)
        return
        endif      ! for isign=1

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccc

        return
        end


      



       



      



       


