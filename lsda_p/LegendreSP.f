        subroutine LegendreSP(Plm,x)

*************************************************************************
**  Written by Lin-Wang Wang, 2003
*************************************************************************
**  copyright (c) 2003, The Regents of the University of California,
**  through Lawrence Berkeley National Laboratory (subject to receipt of any
**  required approvals from the U.S. Dept. of Energy).  All rights reserved.
*************************************************************************


cccc Plm(m,l)[x]=sqrt[(2*l+1)/4pi*(l-m)!/(l+m)!] P_l^m(x)
cccc  here, P_l^m(x) is the Lengendre function

cccc Note the normalized spherical harmonic function is:
cccc Y_lm(th,phi)=Plm(m,l)[costh]*exp(i*m*phi)

        implicit double precision(a-h,o-z)
        parameter (pi=3.14159265358979312d0)
        real*8 plm(-6:6,0:6),F(-6:6,0:6)
        

ccccccc   F(m,l)=d^(l+m)/dx^(l+m) (x^2-1)^l
        x2=x**2
        x3=x2*x
        x4=x3*x
        x5=x4*x
        x6=x5*x
        x7=x6*x
        x8=x7*x
        x9=x8*x
        x10=x9*x
        x11=x10*x
        x12=x11*x

        F(0,0)=1.d0

        F(-1,1)=x2-1
        F(0,1)=2*x
        F(1,1)=2

        F(-2,2)=x4-2*x2+1
        F(-1,2)=4*x3-4*x
        F(0,2)=12*x2-4
        F(1,2)=24*x
        F(2,2)=24

        F(-3,3)=x6-3*x4+3*x2-1
        F(-2,3)=6*x5-12*x3+6*x   
        F(-1,3)=30*x4-36*x2+6
        F(0,3)=120*x3-72*x
        F(1,3)=360*x2-72
        F(2,3)=720*x
        F(3,3)=720

        F(-4,4)=x8-4*x6+6*x4-4*x2+1
        F(-3,4)=8*x7-24*x5+24*x3-8*x
        F(-2,4)=56*x6-120*x4+72*x2-8
        F(-1,4)=336*x5-480*x3+144*x
        F(0,4)=1680*x4-1440*x2+144
        F(1,4)=6720*x3-2880*x
        F(2,4)=20160*x2-2880
        F(3,4)=40320*x
        F(4,4)=40320

        F(-5,5)=x10-5*x8+10*x6-10*x4+5*x2-1
        F(-4,5)=10*x9-40*x7+60*x5-40*x3+10*x
        F(-3,5)=90*x8-280*x6+300*x4-120*x2+10
        F(-2,5)=720*x7-1680*x5+1200*x3-240*x 
        F(-1,5)=5040*x6-8400*x4+3600*x2-240
        F(0,5)=30240*x5-33600*x3+7200*x
        F(1,5)=151200*x4-100800*x2+7200
        F(2,5)=604800*x3-201600*x
        F(3,5)=1814400*x2-201600
        F(4,5)=3628800*x
        F(5,5)=3628800

        F(-6,6)=x12-6*x10+15*x8-20*x6+15*x4-6*x2+1
        F(-5,6)=12*x11-60*x9+120*x7-120*x5+60*x3-12*x
        F(-4,6)=132*x10-540*x8+840*x6-600*x4+180*x2-12
        F(-3,6)=1320*x9-4320*x7+5040*x5-2400*x3+360*x
        F(-2,6)=11880*x8-30240*x6+25200*x4-7200*x2+360
        F(-1,6)=95040*x7-181440*x5+100800*x3-14400*x
        F(0,6)=665280*x6-907200*x4+302400*x2-14400
        F(1,6)=3991680*x5-3628800*x3+604800*x
        F(2,6)=19958400*x4-10886400*x2+604800
        F(3,6)=79833600*x3-21772800*x
        F(4,6)=239500800*x2-21772800
        F(5,6)=479001600*x
        F(6,6)=479001600

        Plm(0,0)=1.d0/dsqrt(4*pi)

        do l=1,6
        call factorial(ifact,l)

        Plm(0,l)=1.d0/2.d0**l/ifact*F(0,l)
        
        do m=1,l
        isign=-1
        if(mod(m,2).eq.0) isign=1
        Plm(m,l)=isign/2.d0**l/ifact*
     &       (1-x**2)**(m/2.d0)*F(m,l)    ! this is the Legendre func. 


        call factorial(ifact1,l-m)
        call factorial(ifact2,l+m)

        Plm(m,l)=dsqrt((2*l+1.d0)*ifact1/(4*pi*ifact2))*
     &      Plm(m,l)            


        Plm(-m,l)=isign*Plm(m,l)
        enddo

        enddo
        
        return
        contains

        subroutine factorial(ifact_tmp,ii)
        implicit double precision (a-h,o-z)
        integer ifact_tmp,ii

        if(ii.eq.0)  then
        ifact_tmp=1
        endif

        if(ii.gt.0) then
        ifact_tmp=1
        do i=1,ii
        ifact_tmp=ifact_tmp*i
        enddo
        endif

        return
        end subroutine factorial

        end

        


     
     
                                                                        
        
