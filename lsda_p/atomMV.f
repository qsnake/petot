      subroutine atomMV(mov,line_step,xatom,fatom,Etot,
     & E_pred,dt,istop_line,dd_max,dd_limit,AL,imov_at,
     & convergE,message,dt_fact_nextlinemin,xatom_old,imv_cont)
******************************************
cc     Written by Lin-Wang Wang, March 30, 2001.  
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

******************************************

***********************************************
*** line_step=0, the begining of a new line minization direction
******************************************
***** The Hessian matrix has a 3*natom x 3*natom dimension, 
***** even if some of the degree's are not moved
*******************************************************
***** This subroutine contains many units, the values of mov and line_step
***** determine whether to perform each of these (logic) units.
******************************************
*****   $$$$$$$$$
***** fatom is the derivative of total energy, thus is the minus actual force
******************************************
       use data

       implicit double precision (a-h,o-z)

       include 'mpif.h'
       include 'param.escan_real'
******************************************
       real*8 px(3*matom),px1(3*matom),px2(3*matom)
       real*8 px_u(3*matom)
       real*8 dsx(3*matom)
       real*8 E_tot_st(0:200),dt_st(0:200)
       real*8 f_max,dd_max_prev
********************************************
       real*8 AL(3,3)
********************************************
       real*8 xatom(3,matom),fatom(3,matom)
       real*8 xatom_old(3,matom),fatom_old(3,matom)
       integer iatom(matom),ityatom(matom)
       integer smatr(3,3,48),nrot
       real*8 Ealpha(mtype)
*************************************************
cc       real*8 HessI(3*matom,3*matom)
       real*8,allocatable,dimension(:,:) :: HessI

       integer imov_at(3,matom)
       real*8 dE_dt
       character*60 message
**************************************************
       save HessI,fatom_old,E_tot_st,dt_st,
     & px,px_u,dE_dt,f_max,dd_max_prev

       common /comMainEtot/iatom,totNel,
     &     smatr,nrot,ilocal,Ealpha,Ealphat

*******************************************************
**** for line_step=0, all the saved quantities will be updated. 
**** for line_step=1, the saved quantities will not be changed, except E_tot_st, 
****  and dt_st, they will be kept updated for all line_step
*******************************************************


      if(line_step.gt.0) goto 3000


      if(mov.eq.1.and.line_step.eq.0) then

      allocate(HessI(3*natom,3*natom))

      if(imv_cont.eq.0) then
      HessI=0.d0
      do i1=1,3*natom
      HessI(i1,i1)=1.d0
      enddo
      else
        if(inode_tot.eq.1) then
        open(10,file="Hessian_cont",form="unformatted")
        rewind(10)
        read(10) natomt
        if(natomt.ne.natom) then
        write(6,*) "natomt.ne.natom,stop"
        stop
        endif
        read(10) HessI
        read(10) ((xatom_old(i,j),i=1,3),j=1,natom)
        read(10) ((fatom_old(i,j),i=1,3),j=1,natom)
        close(10)
        endif

        call mpi_bcast(HessI,9*natom**2,
     $   MPI_REAL8, 0, MPI_COMM_WORLD,ierr)     ! all the proc (including diff icolor) are doing the same
        call mpi_bcast(xatom_old,3*natom,
     $   MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
        call mpi_bcast(fatom_old,3*natom,
     $   MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
       endif

      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc determine the search direction, CG or BFGS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  move.gt.1 updating the Hess matrix using BFGS

      if((mov.gt.1.or.imv_cont.eq.1).and.line_step.eq.0) then
      do ia=1,natom
      ian=(ia-1)*3
      px(1+ian)=(fatom(1,ia)-fatom_old(1,ia))*imov_at(1,ia)
      px(2+ian)=(fatom(2,ia)-fatom_old(2,ia))*imov_at(2,ia)
      px(3+ian)=(fatom(3,ia)-fatom_old(3,ia))*imov_at(3,ia)
      enddo

      do ia=1,natom
      ian=(ia-1)*3
      dsx(1+ian)=AL(1,1)*(xatom(1,ia)-xatom_old(1,ia))+
     &           AL(1,2)*(xatom(2,ia)-xatom_old(2,ia))+
     &           AL(1,3)*(xatom(3,ia)-xatom_old(3,ia))
      dsx(2+ian)=AL(2,1)*(xatom(1,ia)-xatom_old(1,ia))+
     &           AL(2,2)*(xatom(2,ia)-xatom_old(2,ia))+
     &           AL(2,3)*(xatom(3,ia)-xatom_old(3,ia))
      dsx(3+ian)=AL(3,1)*(xatom(1,ia)-xatom_old(1,ia))+
     &           AL(3,2)*(xatom(2,ia)-xatom_old(2,ia))+
     &           AL(3,3)*(xatom(3,ia)-xatom_old(3,ia))
      enddo

      ss=0.d0
      do i1=1,3*natom
      ss=ss+dsx(i1)*px(i1)     ! This formula is correct, fatom is +dE/dx, not -phy.force
      enddo

      if(abs(ss).lt.1.D-25) goto 341

      ssI=1.d0/ss

      do i1=1,3*natom
      s1=0.d0
      do i2=1,3*natom
      s1=s1+HessI(i1,i2)*px(i2)
      enddo
      px1(i1)=s1
      enddo

      ss2=0.d0
      do i1=1,3*natom
      ss2=ss2+px1(i1)*px(i1)
      enddo

      if(abs(ss2).lt.1.D-25) goto 341

      ss2I=1.d0/ss2

      do i1=1,3*natom
      px2(i1)=dsx(i1)*ssI-px1(i1)*ss2I
      enddo

cccccccccc  BFGS formula
      do i1=1,3*natom
      do i2=1,3*natom
      HessI(i1,i2)=HessI(i1,i2)+dsx(i1)*dsx(i2)*ssI-
     &  px1(i1)*px1(i2)*ss2I+px2(i1)*px2(i2)*ss2
      enddo
      enddo

      endif

341   continue
ccccccccccc end of updating HessI  ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 

ccccccccc line_step.eq.0  generate the searching direction and calc. dE_dt

      if(line_step.eq.0) then
      do ia=1,natom
      xatom_old(1,ia)=xatom(1,ia)
      xatom_old(2,ia)=xatom(2,ia)
      xatom_old(3,ia)=xatom(3,ia)
      fatom_old(1,ia)=fatom(1,ia)
      fatom_old(2,ia)=fatom(2,ia)
      fatom_old(3,ia)=fatom(3,ia)
      enddo

        if(inode_tot.eq.1) then
        open(10,file="Hessian_cont",form="unformatted")
        rewind(10)
        write(10) natom
        write(10) HessI
        write(10) ((xatom_old(i,j),i=1,3),j=1,natom)
        write(10) ((fatom_old(i,j),i=1,3),j=1,natom)
        close(10)
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do ia=1,natom
      ian=(ia-1)*3
      px(1+ian)=fatom(1,ia)*imov_at(1,ia)
      px(2+ian)=fatom(2,ia)*imov_at(2,ia)
      px(3+ian)=fatom(3,ia)*imov_at(3,ia)
      enddo


      do i1=1,3*natom
      s=0.d0
      do i2=1,3*natom
      s=s+HessI(i1,i2)*px(i2)
      enddo
      px1(i1)=-s      ! px1 is the search dir. Note, it is -fatom, and fatom=+dE/dx
      enddo

      s1=0.d0
      s2=0.d0
      do i1=1,3*natom
      s1=s1+px1(i1)**2
      s2=s2+px(i1)**2
      enddo

      ss=dsqrt(s2)/dsqrt(s1)


      do i1=1,3*natom
      px(i1)=ss*px1(i1)
      enddo


********************************************************
***** px(i1) is the new direction for this mov, at line_step=0
********************************************************

      dE_dt=0.d0
      do ia=1,natom
      ian=(ia-1)*3
      dE_dt=dE_dt+px(1+ian)*fatom(1,ia)+px(2+ian)*fatom(2,ia)+
     & px(3+ian)*fatom(3,ia)
      enddo

      if(dE_dt.gt.0.d0) then     ! usually, this will not happen. 
      dE_dt=-dE_dt
      do i1=1,3*natom
      px(i1)=-px(i1)
      enddo
      endif


cccccccc change the unit

      do ia=1,natom
      ian=(ia-1)*3
      s1=px(1+ian)*ALI(1,1)+px(2+ian)*ALI(2,1)+px(3+ian)*ALI(3,1)
      s2=px(1+ian)*ALI(1,2)+px(2+ian)*ALI(2,2)+px(3+ian)*ALI(3,2)
      s3=px(1+ian)*ALI(1,3)+px(2+ian)*ALI(2,3)+px(3+ian)*ALI(3,3)
      px_u(1+ian)=s1
      px_u(2+ian)=s2
      px_u(3+ian)=s3
      enddo

      endif

cccccccccccccccc

       f_max=0.d0
       do ia=1,natom
       ian=(ia-1)*3
       f_tmp=px(1+ian)**2+px(2+ian)**2+px(3+ian)**2
       if(f_tmp.gt.f_max) f_max=f_tmp
       enddo
       f_max=dsqrt(f_max)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc end of search direction 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
3000   continue          !  the returning point of re-entry with line_step > 0


       E_tot_st(line_step)=Etot     

ccccccccc E_tot_st(0) is the energy at the begining of this line_search

ccccccccc istop_line=1, stop this line miniz at this line_step
ccccccccc istop_line=0, continue line miniz along this direction.

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       if(line_step.eq.0) then      ! determine dt for a initial probing run
       istop_line=0
       dt=dt*dt_fact_nextlinemin         ! using the dt of previous mov
        if(mov.gt.1.and.dabs(dd_max_prev-dd_limit).lt.0.1*dd_limit) then
        dt=dd_limit/f_max
        endif
       endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc E=E_0+ a*t + b*t^2

       if(line_step.eq.1) then       ! determine the update dt, and E_pred
       a=dE_dt
       b=(E_tot_st(1)-E_tot_st(0)-a*dt_st(1))/dt_st(1)**2
       dt=-a/(2*b)


       E_pred=E_tot_st(0)+a*dt+b*dt**2


         if(E_pred.le.E_tot_st(0)+convergE) then
            istop_line=1                              ! normal ending of the line_min
            message="normal ending of line_min"
            dt_fact_nextlinemin=1.d0
            if(E_pred.gt.E_tot_st(1)) then
            dt=dt_st(1)
            message="E_pred.lt.E_0, gt.E_1, use E_1 to end"   ! warning
            dt_fact_nextlinemin=3.d0
            endif
         else
            if(E_tot_st(1).le.E_tot_st(0)) then
            istop_line=1      ! stop by the first line_step
            dt=dt_st(1)
         message="E_pred.gt.E_0, E_1.lt.E_0,+curv.use E_1"    ! warning
            dt_fact_nextlinemin=3.d0
	    else
            istop_line=0        ! have to try it again, at half dt_st(1)
         message="E_pred.gt.E_0, E_1.gt.E_0, try again dt/2"       ! warning
            dt_fact_nextlinemin=1.d0
            dt=dt_st(1)/2
            endif
         endif
      endif       !  line_step=1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc for line_step.gt.1, as long as there is one calc with its Etot less than the
cccccccc original, the line miniz will be stopped.

       if(line_step.gt.1) then

	 E_min=1000.d0
	 do ii=1,line_step
	 if(E_tot_st(ii).lt.E_min) then
	 E_min=E_tot_st(ii)
	 line_min=ii
	 endif
	 enddo

	 if(E_min.le.E_tot_st(0)+20*convergE) then
         istop_line=1
	 dt=dt_st(line_min)
	 else
         istop_line=0
	 d_min=1.D+10
	 do ii=1,line_step
	 if(dt_st(ii).lt.d_min) d_min=dt_st(ii)
	 enddo
	 dt=d_min/2
	 endif

       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if(line_step.gt.7.and.istop_line.eq.0.and.inode.eq.1) then
       write(6,*) "line_step.gt.7, stop", mov,line_step
       call  mpi_abort(MPI_COMM_WORLD,ierr)
       endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if(f_max*dt.gt.dd_limit) then
       dt=dd_limit/f_max
         message="dd_max.gt.dd_limit,use dd_limit"
       if(line_step.eq.0) then
       E_pred=0.d0
       else
       E_pred=E_tot_st(0)+a*dt+b*dt**2
       endif
       endif
       dd_max=dt*f_max
       dd_max_prev=dd_max
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       line_step=line_step+1
       dt_st(line_step)=dt
ccccccccccc update the xatom in search direction and do a 
ccccccccccc Etot calculation

      do ia=1,natom
      ian=(ia-1)*3
      xatom(1,ia)=xatom_old(1,ia)+dt*px_u(1+ian)
      xatom(2,ia)=xatom_old(2,ia)+dt*px_u(2+ian)
      xatom(3,ia)=xatom_old(3,ia)+dt*px_u(3+ian)
      enddo

      ave_dd=0.d0
      num=0.d0
      do ia=1,natom
      ian=(ia-1)*3
      sum=(AL(1,1)*px_u(1+ian)+AL(1,2)*px_u(2+ian)+
     &     AL(1,3)*px_u(2+ian))**2
     &   +(AL(2,1)*px_u(1+ian)+AL(2,2)*px_u(2+ian)+
     &     AL(2,3)*px_u(2+ian))**2
     &   +(AL(3,1)*px_u(1+ian)+AL(3,2)*px_u(2+ian)+
     &     AL(3,3)*px_u(2+ian))**2
      ave_dd=ave_dd+dsqrt(sum)
      if(imov_at(1,ia).eq.1.or.imov_at(2,ia).eq.1.or.
     &     imov_at(3,ia).eq.1) then
      num=num+1
      endif
      enddo
      ave_dd=dt*ave_dd/num
ccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc write out the xatom information in xatom.log
ccccc file 23 should be opened now
      if(inode_tot.eq.1) then
      write(23,*) "**********************************"
      write(23,*) "***  mov, line_step=", mov, line_step
      write(23,*) natom
      write(23,203) AL(1,1),AL(2,1),AL(3,1)
      write(23,203) AL(1,2),AL(2,2),AL(3,2)
      write(23,203) AL(1,3),AL(2,3),AL(3,3)
      do i=1,natom
      write(23,204) iatom(i),xatom(1,i),xatom(2,i),xatom(3,i),
     &   imov_at(1,i),imov_at(2,i),imov_at(3,i)
      enddo
      write(23,*) "**********************************"
      endif
203   format(3(f12.6,1x))
204   format(i3,2x,3(f12.6,1x),2x,3(i2,1x))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc making sure the xatom for all the nodes are the same

        call mpi_bcast(xatom,3*natom,
     $   MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
      

      return       ! return to calculate Etotcalc at the new xatom
       
      end

