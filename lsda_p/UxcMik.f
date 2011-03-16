      function  UxcMik(d,Uxc2)
      implicit double precision (a-h,o-z)
      parameter (a0=.4581652932831429d0,a1=2.40875407d0,
     $           a2=.88642404d0,a3=.02600342d0)
      parameter (b1=1.0d0,b2=4.91962865d0,b3=1.34799453d0,
     $           b4=.03120453d0)
      parameter (c1=4.d0*a0*b1/3.0d0,  c2=5.0d0*a0*b2/3.0d0+a1*b1,
     $         c3=2.0d0*a0*b3+4.0d0*a1*b2/3.0d0+2.0d0*a2*b1/3.0d0,
     $      c4=7.0d0*a0*b4/3.0d0+5.0d0*a1*b3/3.0d0+a2*b2+a3*b1/3.0d0,
     $        c5=2.0d0*a1*b4+4.0d0*a2*b3/3.0d0+2.0d0*a3*b2/3.0d0,
     $        c6=5.0d0*a2*b4/3.0d0+a3*b3,c7=4.0d0*a3*b4/3.0d0)
      parameter (pi4=3.1415926d0*4)
      parameter (rsfac=0.62035 04908 994000d0)
      real*8 UxcMik,d

       rh=max(d,1.d-16)
       rs=rsfac*rh**(-.3333333333333333d0)  
       excnum=-(a0+rs*(a1+rs*(a2+rs*a3)))
       excden=rs*(b1+rs*(b2+rs*(b3+rs*b4)))
       exc=excnum/excden 
       vxcnum=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
       vxc=vxcnum/(excden*excden)
       Uxc1=vxc
       Uxc2=exc

       UxcMik=Uxc1

       return
       end


      

