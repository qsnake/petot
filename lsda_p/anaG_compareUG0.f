        program anaG_compareUG0

***********************************************
***    written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

cccccccccc

cccccccc This program does a projection analysis
cccccccc based on the graphu.G0.kp files generated
cccccccc from the program anaG_PEtot for the 
cccccccc large system (move graphu.G0.kp to graphu.large)
cccccccc and the primary cell system (move graphu.G0.kp to graphu.primary)
cccccccc Lin-Wang Wang, April 2, 2002

        parameter (mg=4000,mx=20)
        implicit double precision (a-h,o-z)

        complex*16 ug0(mg,mx),ug1(mg,mx)
        real*8 ak0(3,mg),ak1(3,mg)
        real*8 tmp(2*mx)
        real*8 pp(mx,mx),weight(mx)
        complex*16 cpp(mx,mx)

        open(10,file="graphu.large")
        rewind(10)
        read(10,*) ng1,mx1
        if(ng1.gt.mg.or.mx1.gt.mx) then
        write(6,*) "ng1.gt.mg.or.mx1.gt.mx, stop", ng1,mg,mx1,mx
        stop
        endif
        do i=1,ng1
        read(10,*) ak1(1,i),ak1(2,i),ak1(3,i),
     &  (tmp(j),j=1,2*mx1)
         do j=1,mx1
         ug1(i,j)=dcmplx(tmp(j*2-1),tmp(j*2))
         enddo
        enddo
        close(10)

        open(11,file="graphu.primary")
        rewind(11)
        read(11,*) ng0,mx0
        if(ng0.gt.mg.or.mx0.gt.mx) then
        write(6,*) "ng0.gt.mg.or.mx0.gt.mx, stop", ng0,mg,mx0,mx
        stop
        endif
        if(ng0.ne.ng1) then
        write(6,*) "ng0.ne.ng1, stop", ng0,ng1
        stop
        endif
        do i=1,ng0
        read(11,*) ak0(1,i),ak0(2,i),ak0(3,i),
     &  (tmp(j),j=1,2*mx0)
         do j=1,mx0
         ug0(i,j)=dcmplx(tmp(j*2-1),tmp(j*2))
         enddo
        enddo
        close(11)

        do j=1,mx1
        sum=0.d0
        do i=1,ng1
        sum=sum+cdabs(ug1(i,j))**2 
        enddo
        weight(j)=sum
        sum=1.d0/dsqrt(sum)
        do i=1,ng1
        ug1(i,j)=ug1(i,j)*sum
        enddo
        enddo

        do j=1,mx0
        sum=0.d0
        do i=1,ng0
        sum=sum+cdabs(ug0(i,j))**2 
        enddo
        sum=1.d0/dsqrt(sum)
        do i=1,ng0
        ug0(i,j)=ug0(i,j)*sum
        enddo
        enddo
cccccccccccccccccccccccccccccccccccccccccccccc
        cpp=dcmplx(0.d0,0.d0)

        do 200 i1=1,ng1
        i10=0
        do i0=1,ng0
        if(ak0(1,i0).eq.ak1(1,i1).
     & and.ak0(2,i0).eq.ak1(2,i1).
     & and.ak0(3,i0).eq.ak1(3,i1)) then
        i10=i0
        goto 101
        endif
        enddo
101     continue
        if(i10.eq.0) then
        write(6,*) "ak1(i1).not found in ak0, stop", i1
        stop
        endif

        do j1=1,mx1
        do j0=1,mx0
        cpp(j0,j1)=cpp(j0,j1)+dconjg(ug0(i10,j0))*ug1(i1,j1)
        enddo
        enddo
200     continue
       
        do j1=1,mx1
        do j0=1,mx0
        pp(j0,j1)=cdabs(cpp(j0,j1))**2
        enddo
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        write(6,*) "each row is one state in graphu.large"
        write(6,*) "Weights are normalized for each state"
        write(6,*) "***************************"
        write(6,*) "The projections"
        do j1=1,mx1
        write(6,301) j1,(pp(j0,j1),j0=1,mx0)
        enddo
        write(6,*) "***************************"
301     format(i3,2x,10(E8.3,1x))
        stop
        end

        
       


        
        

