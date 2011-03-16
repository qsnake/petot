      subroutine UxcCA2(rho1,rho2,vxc1,vxc2,uxc1,uxc2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      DATA CEX/-1.969490099D0/
      DATA ONE/1.0D0/,HALF/0.5D0/,THALF/1.5D0/
      DATA ZERO/0.0D0/,THRD/0.333333333333333D0/
      DATA TFT/0.75D0/,TWO/2.0D0/
      DATA RC/0.02258D0/,EPS/1.0D-20/
C-----------------------------------------------------------------------


C.
C. KOHN-SHAM EXCHANGE
C. CEPERLEY-ALDER / RPA CORRELATION
C.
cccc The parameters are for the unit of Ryd. 
cccc beta11=7/6*beta1, beta22=4/3*beta2
cccc B1=B-A/3, C1=2/3*C, D1=(2*D-C)/3

      DATA PI /3.141592654D0/
      DATA gamma_u,beta1_u,beta2_u /-0.2846d0, 1.0529d0, 0.3334d0/
      DATA gamma_p,beta1_p,beta2_p /-0.1686d0, 1.3981d0, 0.2611d0/
      DATA beta11_u,beta22_u /1.228383333D0,0.4445333333D0/
      DATA beta11_p,beta22_p /1.631116667D0,0.3481333333D0/

      DATA A_u,B_u,C_u,D_u /0.0622d0,-0.096d0,0.0040d0,-0.0232d0/
      DATA B1_u,C1_u,D1_u /-0.116733333d0,0.0026666667d0,-0.0168d0/ 

      DATA A_p,B_p,C_p,D_p /0.0311d0,-0.0538d0,0.0014d0,-0.0096d0/
      DATA B1_p,C1_p,D1_p /-0.064166667d0,0.00093333d0,
     &                                            -0.006866667d0/

************************************************
      rho1=max(rho1,1.d-16)
      rho2=max(rho2,1.d-16)

      rho=rho1+rho2



      RS=(PI*RHO/TFT)**(-THRD)
      IF(RS .GE. ONE) THEN

       ROOTRS=SQRT(RS)
       tmp=1.d0/(1.d0+beta1_u*ROOTRS+beta2_u*RS)
       Ec_u=gamma_u*tmp
       Vc_u=gamma_u*(1.d0+beta11_u*ROOTRS+beta22_u*RS)*tmp**2

       tmp=1.d0/(1.d0+beta1_p*ROOTRS+beta2_p*RS)
       Ec_p=gamma_p*tmp
       Vc_p=gamma_p*(1.d0+beta11_p*ROOTRS+beta22_p*RS)*tmp**2

      ELSE

       XLNRS=LOG(RS)
       Ec_u=A_u*XLNRS+B_u+C_u*RS*XLNRS+D_u*RS
       Vc_u=A_u*XLNRS+B1_u+C1_u*RS*XLNRS+D1_u*RS

       Ec_p=A_p*XLNRS+B_p+C_p*RS*XLNRS+D_p*RS
       Vc_p=A_p*XLNRS+B1_p+C1_p*RS*XLNRS+D1_p*RS

      END IF
*************************************
********* interpolation *************
 
        x=(rho1-rho2)/rho
        fx=((1.d0+x)**(4.d0/3)+(1.d0-x)**(4.d0/3)-2)/
     &     ((2.d0)**(4.d0/3)-2)
        dfx=((1.d0+x)**(1.d0/3)-(1.d0-x)**(1.d0/3))*
     &      4.d0/3/((2.d0)**(4.d0/3)-2)

        Ec=Ec_u+fx*(Ec_p-Ec_u)

        Vc1=Vc_u+fx*(Vc_p-Vc_u)+(Ec_p-Ec_u)*(1-x)*dfx
        Vc2=Vc_u+fx*(Vc_p-Vc_u)+(Ec_p-Ec_u)*(-1-x)*dfx

*************************************
****  Add the exchange part
***** 1/TWO: change ryd to Hartree

      Vxc1=(CEX*(2*rho1)**THRD+Vc1)/TWO
      Vxc2=(CEX*(2*rho2)**THRD+Vc2)/TWO
      
      Uxc1=(TFT*CEX*(2*rho1)**THRD+Ec)/TWO
      Uxc2=(TFT*CEX*(2*rho2)**THRD+Ec)/TWO

      RETURN
      END

