      SUBROUTINE PCC_SIN(NR,ICORE,AC,BC,CC,R,CDC)

ccccccccccccccccccccccccccccccccccccccccccc
cccccc   Replace the original full core charge CDC/r inside r(icore), 
cccccc   with a form ac*sin(bc*r), match the original charge at r(icore)
cccccc   for its value and slop.  The purpose is to generate a smoother
cccccc   core charge to be used in plane wave calculation
cccccc
cccccc   Lin-Wang Wang, Feb. 28, 2001
cccccc    Actually, it should be called: PCC_SIN(xxxx)
cccccc
cccccc
ccccccccc  This subroutine  returns a pseudo core  of
ccccccccc  cdc(r)=ac*r*sin(bc*r)   inside r(icore)
ccccccccc  instead of the original
ccccccccc  cdc(r)=r^2 * exp(ac+bc*r^2+cc*r^4)  inside r(icore)
ccccccccc 
ccccccccc  This new pseudo-core is softer than
ccccccccc  the original exp. core, and uses smaller Ecut in 
ccccccccc  PEtot calculation. For PEtot calculation, we suggest
ccccccccc  this pseudo-core for "pe" calculation (for no GGA pseudopotential).  
ccccccccc  Lin-Wang Wang, Feb. 28, 2001

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (NN = 5)

      DIMENSION R(NR),CDC(NR),
     .          DR_K(-NN:NN),DDR_K(-NN:NN),
     .          CDC_SCA(ICORE-NN:ICORE+NN)

            tpfive=2.5d0
            pi=4*datan(1.d0)


            cdcp = (cdc(icore+1)-cdc(icore))/(r(icore+1)-r(icore))
            tanb = cdc(icore)/(r(icore)*cdcp-cdc(icore))
            rbold = tpfive
            do 500 i = 1, 50
               rbnew = pi + datan(tanb*rbold)
               if (abs(rbnew-rbold) .lt. .00001D0) then
                  bc = rbnew/r(icore)
                  ac = cdc(icore)/(r(icore)*sin(rbnew))
                  do 490 j = 1, icore
                     cdc(j) = ac*r(j)*sin(bc*r(j))
  490             continue
           write(*,*) "r(icore),ac,bc", r(icore), ac, bc
           write(*,*) "**************"
           write(*,*) "**** USING sin(bc*r) for core charge ***"
           write(*,*) "**************"
c
                  go to 530
c
               else
                  rbold = rbnew
               end if
  500       continue
            write(*,*) "SHOULD NEVER BE HERE, fail to find core charge"
            Write(*,*) "inside PCC_EXP.f, stop"
            stop

  530       continue
            cc=1.d0
           
           return
           end
  

