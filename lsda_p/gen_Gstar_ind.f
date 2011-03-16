      subroutine gen_Gstar_ind(smatr,nrot,AL)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
cc     Copyright 2001 The Regents of the University of California
cc     The United States government retains a royalty free license in this work
******************************************

ccccccc This subroutine generates mynode index arrays: ig_star(i),ig_local(i), and
ccccccc ig_lenstar(i), to be used for symmetry operation. 
ccccccc Here, i=1,ng2_n. ig_local(i) is the local index of g vector. 
ccccccc  e.g, rho_g(ig_local(i)). ig_star(i) is the
ccccccc index of global g_star for which the ig_local(i) belong.
ccccccc ig_star(i) is in ascending order. ig_lenstar(i) is the total number of g-vec
ccccccc in this g-star

ccccccc ALI is passed in through param.escan_real


      use fft_data
      use load_data
      use data

      implicit double precision (a-h,o-z)
      include 'mpif.h'
      include 'param.escan_real'

      integer smatr(3,3,48)
      real*8 smatrC(3,3,48),tmp(3,3)
      real*8 AL(3,3)
      real*8 akxyz_sym(3,48),akxyz(3)
      real*8 one2pi
      logical,allocatable,dimension (:) :: itreat
      logical inew,inew_gstar

**************************************************

cccccccc generate the sym op for Cartesian Coord.
       do irot=1,nrot

       do k=1,3
       do j=1,3
       s=0.d0
       do i=1,3
       s=s+AL(k,i)*smatr(j,i,irot)
       enddo
       tmp(k,j)=s
       enddo
       enddo

       do k=1,3
       do k2=1,3
       s=0.d0
       do j=1,3
       s=s+tmp(k,j)*ALI(k2,j)
       enddo
       smatrC(k2,k,irot)=s
       enddo
       enddo

       enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      one2pi=1.0002d0/(2*4*datan(1.d0))
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ng2_n=ngtotnod2(inode)
      ig_treat_loc=0

      allocate(itreat(mr_n))
      itreat=.false.
      num_gstar=0
      ig_loc_count=0

      do 200 iproc=1,nnodes

      if(iproc.eq.inode) ng2_tmp=ng2_n

      call mpi_bcast(ng2_tmp,1,MPI_INTEGER,iproc-1,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      do 100 i=1,ng2_tmp   ! all the processor doing the same loop as iproc

      if(iproc.eq.inode) then
      inew_gstar =.not.itreat(i)
      endif

      call mpi_bcast(inew_gstar,1,MPI_LOGICAL,iproc-1,
     &        MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if(inew_gstar) then    ! doing the new gstar
      num_gstar=num_gstar+1

      if(iproc.eq.inode) then
      akxyz(1)=gkx2_n(i)
      akxyz(2)=gky2_n(i)
      akxyz(3)=gkz2_n(i)
      endif

      call mpi_bcast(akxyz,3,MPI_REAL8,iproc-1,MPI_COMM_WORLD,ierr)
          
cccccccc now, all the processor are doing the same thing: to generate gstar

           len_gstar=0
           do 50 irot=1,nrot
           akxn=smatrC(1,1,irot)*akxyz(1)+smatrC(2,1,irot)*akxyz(2)+
     &      smatrC(3,1,irot)*akxyz(3)
           akyn=smatrC(1,2,irot)*akxyz(1)+smatrC(2,2,irot)*akxyz(2)+
     &      smatrC(3,2,irot)*akxyz(3)
           akzn=smatrC(1,3,irot)*akxyz(1)+smatrC(2,3,irot)*akxyz(2)+
     &      smatrC(3,3,irot)*akxyz(3)
           
ccccccc a small double loop, to find out whether this new ak is already in the
ccccccc n_gstar group
              inew=.true.
              do j=1,len_gstar     
              if(dabs(akxyz_sym(1,j)-akxn)+dabs(akxyz_sym(2,j)-akyn)+
     &           dabs(akxyz_sym(3,j)-akzn).lt.0.00001d0) inew=.false.
              enddo

              if(inew) then
              len_gstar=len_gstar+1
              akxyz_sym(1,len_gstar)=akxn
              akxyz_sym(2,len_gstar)=akyn
              akxyz_sym(3,len_gstar)=akzn
              endif
50         continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc  note that, all the processors have the same information for
ccccccc  num_gstar,len_gstar, akxyz_sym. 

ccccccccc now, each processor to find whether akxyz belongs to one of their g-vect, 
ccccccccc if yes, write it into ig_star, ig_local, and mark itreat.

            ncount_tmp=0

            do 60 j=1,len_gstar 
            call find_gindex(akxyz_sym(1,j),iflag,ig)

ccccccccccc note, ig could be both >0 (k vector), and <0 (-k vector). 
ccccccccccc but they (k and -k) are not necessarily in the same gstar
ccccccccccc but for akxyz=0, only one k (positive one) is returned)

            if(iflag.eq.1) then     ! yes, this akxyz_sym equals ig
            ncount_tmp=ncount_tmp+1
            ig_loc_count=ig_loc_count+1
            ig_star(ig_loc_count)=num_gstar
            ig_local(ig_loc_count)=ig      ! ig could be negative for -k
            ig_lenstar(ig_loc_count)=len_gstar
            if(ig.gt.0) then
              if(itreat(ig)) then
          write(6,*) "shouldn't be here, stop, gen_gstar_ind,itreat"
              stop
              endif
            itreat(ig)=.true.
            ig_treat_loc=ig_treat_loc+1
            endif
            endif
60          continue
          
cccccccccccc  check whether all g-star has found its and only its g-vect

          call mpi_allreduce(ncount_tmp,ncount_tmp,1,
     &        MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

          if(ncount_tmp.ne.len_gstar) then    ! not all gstar has found its g-vect
          write(6,*) "not all gstar has found its g-vect, stop",
     &               ncount_tmp,len_gstar

          stop
          endif
         
          endif    ! for inew_gstar

100      continue
200      continue

ccccccc  note that ig_loc_count can be larger than ngtotnod2(inode), 
ccccccc  because -k can also be in ig_loc_count, (however, not all -k are 
ccccccc  are necessarily in ig_loc_count)

ccccccc  ig_loc_count will be retained as ng_ig_local (data.f) in each processor
ccccccc  num_gstar will also be retained in data.f

         ng_ig_local=ig_loc_count

         if(ig_treat_loc.ne.ngtotnod2(inode)) then
         write(6,*) "ig_treat_loc.ne.ngtotnod2,stop",
     &      ig_treat_loc, ngtotnod2(inode)
         call mpi_abort(MPI_COMM_WORLD,ierr)
         endif

         if(inode.eq.1) write(6,*) "num_gstar=",num_gstar
         deallocate(itreat)
         
      return

      contains

       subroutine find_gindex(akxyzt,iflag2,ig2)

ccccccc this is reversion gen_G2_real, i.e, from akx,y,z, to find ig,jnode

       implicit double precision (a-h,o-z)
       real*8 akxyzt(3)

       iflag2=0

       imink=1
10     continue

       nxt=one2pi*(AL(1,1)*akxyzt(1)+AL(2,1)*akxyzt(2)+
     &             AL(3,1)*akxyzt(3))
       nyt=one2pi*(AL(1,2)*akxyzt(1)+AL(2,2)*akxyzt(2)+
     &             AL(3,2)*akxyzt(3))
       nzt=one2pi*(AL(1,3)*akxyzt(1)+AL(2,3)*akxyzt(2)+
     &             AL(3,3)*akxyzt(3))

       if(nxt.lt.0) nxt=nxt+n1
       if(nyt.lt.0) nyt=nyt+n2
       if(nzt.lt.0) nzt=nzt+n3
       nxt=nxt+1
       nyt=nyt+1
       nzt=nzt+1
       if(nxt.gt.n1.or.nyt.gt.n2.or.nzt.gt.n3) then
       write(6,*) "nxt,nyt,nzt out of range, stop, find_gindex",
     &          nxt,nyt,nzt
       stop
       endif

       if(jjnode2(nyt,nzt).eq.inode) then
        jjcol2_t=jjcol2(nyt,nzt)
        do ig_t=igstar_jjcol2(jjcol2_t),igfin_jjcol2(jjcol2_t)     ! short loop
        if(n1p2_n(ig_t).eq.nxt) then
        iflag2=1
        ig2=ig_t
        if(imink.eq.2) ig2=-ig_t
        return
        endif
        enddo
       endif

      
       if(imink.eq.2) return

       if(dabs(akxyzt(1))+
     &       dabs(akxyzt(2))+dabs(akxyzt(3)).gt.0.00001d0) then
       imink=2
       akxyzt=-akxyzt       ! try -k
       goto 10
       endif

       return

       end subroutine find_gindex

       end
      
