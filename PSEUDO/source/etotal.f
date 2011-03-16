      subroutine etotal(itype,zsh,nameat,norb,
     1 no,lo,so,zo,etot,ev,ek,ep)
c 
c  etotal computes the total energy from the 
c  electron charge density.
c  
c  njtj
c  ###  Cray conversions  
c  ###    1)Comment out implicit double precision.
c  ###    2)Switch double precision parameter
c  ###      to single precision parameter statement.
c  ###  Cray conversions
c  njtj
c
c
      implicit double precision (a-h,o-z)
      parameter (zero=0.D0)
Cray      parameter (zero=0.0)
c
      character*1 il(5)
      character*2 nameat
c
      dimension no(norb),lo(norb),so(norb),zo(norb),
     1 etot(10),ev(norb),ek(norb),ep(norb)
c
c      etot(i)    i=1,10 contains various contributions to the total
c                 energy.
c                 (1)   sum of eigenvalues ev
c                 (2)   sum of orbital kinetic energies ek
c                 (3)   el-ion interaction from sum of orbital
c                       potential energies ep
c                 (4)   electrostatic el-el interaction  (from velect)
c                 (5)   vxc (exchange-correlation) correction to sum
c                       of eigenvalues                   (from velect)
c                 (6)   3 * vc - 4 * ec
c                       correction term for virial theorem
c                       when correlation is included     (from velect)
c                 (7)   exchange and correlation energy  (from velect)
c                 (8)   kinetic energy from eigenvalues  (1,3,4,5)
c                 (9)   potential energy
c                 (10)  total energy
c
c
c      sum up eigenvalues ev, kinetic energies ek, and
c      el-ion interaction ep
c
      etot(1) = zero
      etot(2) = zero
      etot(3) = zero
      do 10 i=1,norb
        etot(1) = etot(1) + zo(i)*ev(i)
        etot(2) = etot(2) + zo(i)*ek(i)
        etot(3) = etot(3) + zo(i)*ep(i)
 10   continue
c
c   kinetic energy
c
      etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)
c
c   potential energy
c
      etot(9) = etot(3) + etot(4) + etot(7)
c
c      total energy
c
      etot(10) = etot(1) - etot(4) - etot(5) + etot(7)
c
c   printout
c
      il(1) = 's'
      il(2) = 'p'
      il(3) = 'd'
      il(4) = 'f'
      il(5) = 'g'
      write(6,20) nameat
 20   format(//,1x,a2,' output data for orbitals',/,1x,28('-'),//,
     1 ' nl    s      occ',9x,'eigenvalue',4x,'kinetic energy',
     2 6x,'pot energy',/)
      do 40 i=1,norb
        write(6,30) no(i),il(lo(i)+1),so(i),zo(i),ev(i),ek(i),ep(i)
 30   format(1x,i1,a1,f6.1,f10.4,3f17.8)
 40   continue
      write(6,50) (etot(i),i=1,10)
 50   format(//,' total energies',/,1x,14('-'),/,
     1 /,' sum of eigenvalues        =',f18.8,
     2 /,' kinetic energy from ek    =',f18.8,
     3 /,' el-ion interaction energy =',f18.8,
     4 /,' el-el  interaction energy =',f18.8,
     5 /,' vxc    correction         =',f18.8,
     6 /,' virial correction         =',f18.8,
     7 /,' exchange + corr energy    =',f18.8,
     8 /,' kinetic energy from ev    =',f18.8,
     9 /,' potential energy          =',f18.8,/,1x,45('-'),
     X /,' total energy              =',f18.8)
       if (itype .ge. 4 .or. abs(zsh) .gt. 0.00001) return
c
c   virial theorem
c
       vsum = 2*etot(8) + etot(9) + etot(6)
       write(6,60) 2*etot(8),etot(9),etot(6),vsum
 60    format(//,' virial theorem(nonrelativistic)',/,1x,14('-'),/,
     1 /,' kinetic energy  *  2      =',f18.8,
     2 /,' potential energy          =',f18.8,
     3 /,' virial correction         =',f18.8,/,1x,45('-'),
     4 /,' virial sum                =',f18.8)
       return
       end
