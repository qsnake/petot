      SUBROUTINE PCC_POLY(NR,ICORE,AC,BC,CC,R,CDC)

ccccccccccccccccccccccccccccccccccccccccccc
cccccc   Replace the original full core charge CDC/r^2 inside r(icore), 
cccccc   with a form (ac+bc*r^2+cc*r^3+dc*r^4), match the core charge
cccccc   for its value, slop, second and third derivatives.  
cccccc   This is used for GGA calculation where the second derivative is needed
cccccc   The resulting core charge is smaller than the PCC_EXP.f result, 
cccccc   and it guarentees that the third derivative is matched. 
cccccc   Possible danger: there is no guarantee that the 
cccccc   core charge is positive !! But usually, it does not happen. 
cccccc
cccccc
ccccccccc  Lin-Wang Wang, Feb. 28, 2001

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION R(NR),CDC(NR)

            v0=cdc(icore)/r(icore)**2

            vp2=cdc(icore+2)/r(icore+2)**2
            vp1=cdc(icore+1)/r(icore+1)**2
            vm1=cdc(icore-1)/r(icore-1)**2
            vm2=cdc(icore-2)/r(icore-2)**2

            d0=(vp1-vm1)/(r(icore+1)-r(icore-1))
            v1=d0

            dp1h=(vp1-v0)/(r(icore+1)-r(icore))
            dp2h=(vp2-vp1)/(r(icore+2)-r(icore+1))
            dm1h=(v0-vm1)/(r(icore)-r(icore-1))
            dm2h=(vm1-vm2)/(r(icore-1)-r(icore-2))

            dd0=(dp1h-dm1h)/((r(icore+1)-r(icore-1))/2)
            v2=dd0
 
            ddp1=(dp2h-dp1h)/((r(icore+2)-r(icore))/2)
            ddm1=(dm1h-dm2h)/((r(icore)-r(icore-2))/2)

            ddd0=(ddp1-ddm1)/((r(icore+2)+2*r(icore+1)-
     &         r(icore-2)-2*r(icore-1))/4)
            v3=ddd0

            x=r(icore)


            dc=(x**2*v3-2*x*v2+2*v1)/(8*x**3)
            cc=(v3-24*dc*x)/6
            bc=(v2-12*dc*x**2-6*cc*x)/2
            ac=(v0-bc*x**2-cc*x**3-dc*x**4)


            do j=1,icore
            x=r(j)
            cdc(j)=x**2*(ac+bc*x**2+cc*x**3+dc*x**4)
            enddo

           
           return
           end
  

