        function clebsch_gordan(j1,m1,j2,m2,j,m)
          implicit none
!     calculate a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
!     arguments are integer and twice the true value. 

        real*8    :: cleb,factor,sum,binom
        real*8    :: clebsch_gordan
        integer :: j1,m1,j2,m2,j,m,par,z,zmin,zmax

!     some checks for validity (let's just return zero for bogus arguments)


        if (2*(j1/2)-int(2*(j1/2.0)).ne.2*(abs(m1)/2)-
     &     int(2*(abs(m1)/2.0)) .or. 
     &     2*(j2/2)-int(2*(j2/2.0)).ne.2*(abs(m2)/2)-
     &     int(2*(abs(m2)/2.0)) .or. 
     &     2*(j/2)-int(2*(j/2.0)).ne.2*(abs(m)/2)-
     &     int(2*(abs(m)/2.0)) .or. 
     &     j1.lt.0 .or. j2.lt.0 .or. j.lt.0 .or. abs(m1).gt.j1 .or. 
     &     abs(m2).gt.j2 .or.
     &     abs(m).gt.j .or. j1+j2.lt.j .or. abs(j1-j2).gt.j .or. 
     &     m1+m2.ne.m) then

           cleb= 0.0
        else
        
           factor = 0.0
           factor = binom(j1,(j1+j2-j)/2) /
     &                binom((j1+j2+j+2)/2,(j1+j2-j)/2)
           factor = factor * binom(j2,(j1+j2-j)/2) / 
     &                binom(j1,(j1-m1)/2)
           factor = factor / binom(j2,(j2-m2)/2) / binom(j,(j-m)/2)
           factor = sqrt(factor)
           
           zmin = max(0,j2+(j1-m1)/2-(j1+j2+j)/2,
     &                j1+(j2+m2)/2-(j1+j2+j)/2)
           zmax = min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)
           
           sum=0.0
           do z = zmin,zmax
              par=1
              if(2*(z/2)-int(2*(z/2.0)) /= 0) par=-1
              sum=sum+par*binom((j1+j2-j)/2,z)*
     &          binom((j1-j2+j)/2,(j1-m1)/2-z)*
     &          binom((-j1+j2+j)/2,(j2+m2)/2-z)
           end do
           
           cleb = factor*sum
        end if
           clebsch_gordan=cleb
      end function clebsch_gordan



      recursive function binom(n,r) result(res)
        implicit none
        integer :: n,r
        real*8 :: res
        real*8 :: tmp

        if(n==r .or. r==0) then
           res = 1.0
        else if (r==1) then
           res = real(n)
        else
           res = real(n)/real(n-r)*binom(n-1,r)
        end if
      end function binom
