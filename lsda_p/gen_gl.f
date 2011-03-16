       program gen_gl

*************************************************************************
*** Written by Lin-Wang Wang, 2001
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************

       implicit double precision(a-h,o-z)
       real*8 y(1000),aj(1000,0:7),daj(1000,0:7)
       real*8 yz(2,0:7)
       real*8 gl(1001,0:7)
       real*8 aj_t1(0:7),aj_t2(0:7)
       real*8 y_cut(0:7),aj_cut(0:7)

       do i=1,1000
       y(i)=(i-1)*20.d0/1000
       if(i.eq.1) y(i)=1.D-3
       x=y(i)
       aj(i,0)=dsin(x)/x
       aj(i,1)=dsin(x)/x**2-dcos(x)/x
       do j=1,6
       aj(i,j+1)=(2*j+1)*aj(i,j)/x-aj(i,j-1)
       enddo
       enddo
ccccccc deal with the numerical error problem
ccccccc we didn't use Taylor expansion
       y(1)=0.d0
       aj(1,0)=1.d0
       do j=1,7
       aj(1,j)=0.d0
       enddo

       do j=1,7
       aj_min=100.d0
       do i=2,300
       if(y(i).lt.0.5d0) then
       if(aj(i,j).gt.0.d0.and.aj(i,j).lt.aj_min) then
       aj_min=aj(i,j)
       i_min=i
       endif
       endif
       enddo
       y_cut(j)=y(i_min)
       aj_cut(j)=aj_min
       do i=2,i_min
       aj(i,j)=aj_cut(j)*(y(i)/y_cut(j))**j
       if(aj(i,j).lt.1.D-30) aj(i,j)=0.d0
       enddo
       enddo
ccccccccccccccccccccccccccccccccccccccccccc
       
       do i=1,1000
       daj(i,0)=-aj(i,1) 
       do j=1,6
       daj(i,j)=(j*aj(i,j-1)-(j+1)*aj(i,j+1))/(2*j+1)
       enddo
       enddo

       do j=0,6
       ll=0
       do i=6,1000
       if(daj(i-1,j)*daj(i,j).lt.0.d0) then
       x=(daj(i,j)*y(i-1)-daj(i-1,j)*y(i))/
     &       (daj(i,j)-daj(i-1,j))
       if(ll.le.1) then
       ll=ll+1
       yz(ll,j)=x
       endif
       endif
       enddo
       enddo

       yz(2,0)=yz(1,0)
       yz(1,0)=0.d0

c       do j=0,6
c       write(6,*) j, yz(1,j),yz(2,j)
c       enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccc
       do j=0,6

       call get_aj(aj_t2,yz(2,j),aj_cut,y_cut)
       call get_aj(aj_t1,yz(1,j),aj_cut,y_cut)

       fac=aj_t2(j)/aj_t1(j)
      

       do i=1,1001
       x=(i-1)*1.d0/1000
       x1=x*yz(1,j)
       x2=x*yz(2,j)

       call get_aj(aj_t1,x1,aj_cut,y_cut)
       call get_aj(aj_t2,x2,aj_cut,y_cut)
       gl(i,j)=aj_t2(j)-aj_t1(j)*fac

ccccccccc  there is one minors for very small x and gl(i,5)
       gl(i,j)=dabs(gl(i,j))

       enddo

       enddo

ccccccccccccccccccccccccccccccccccccccccccc
       do j=0,6
       sum=0.d0
       do i=1,1001 
       x=(i-1)*1.d0/1000
       sum=sum+gl(i,j)*x**(j+2)
       enddo
       sum=sum*1.d0/1000
       sum=1.d0/sum
       do i=1,1001
       gl(i,j)=sum*gl(i,j)
       enddo
       enddo
ccccccccccccccccccccccccccccccccccccccccccc

       open(10,file="graph.j")
       rewind(10)
       do i=1,1001
       x=(i-1)*1.d0/1000
       write(10,200) x,(gl(i,j),j=0,6)
       enddo
       close(10)
200    format(9(E13.7,1x))

       stop
       end

       
       subroutine get_aj(aj_t,x1,aj_cut,y_cut)
       implicit double precision (a-h,o-z)
       real*8 x, aj_t(0:6), aj_cut(0:6), y_cut(0:6)

       x=x1
       if(x.lt.1.D-5) x=1.D-5

       aj_t(0)=dsin(x)/x
       aj_t(1)=dsin(x)/x**2-dcos(x)/x
       do j=1,6
       aj_t(j+1)=(2*j+1)*aj_t(j)/x-aj_t(j-1)
       enddo

       do j=1,6
       if(x.lt.y_cut(j)) then
       aj_t(j)=aj_cut(j)*(x1/y_cut(j))**j
       if(aj_t(j).lt.1.D-30) aj_t(j)=0.d0
       endif
       enddo

       if(x1.lt.1.D-5) aj_t(0)=1.d0

       return
       end
