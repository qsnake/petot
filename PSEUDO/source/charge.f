       real*8 function charge(name)
c  
c  ***********************************************************
c  *                                                         *
c  *   Function determines the nuclear charge of an element. *
c  *   All elements from H to Lr are included.               *
c  *                                                         *
c  *  Version dated May 1, 1991                              *
c  *  njtj                                                   *
c  *                                                         *
c  ***********************************************************
c
       character*2 name
c
       if (name .eq. 'H ' .or. name .eq. ' H') then
         charge = 1.0d0
       elseif (name .eq. 'He') then
         charge = 2.0d0
       elseif (name .eq. 'Li') then
         charge = 3.0d0
       elseif (name .eq. 'Be') then
         charge = 4.0d0
       elseif (name .eq. 'B ' .or. name .eq. ' B') then
         charge = 5.0d0
       elseif (name .eq. 'C ' .or. name .eq. ' C') then
         charge = 6.0d0
       elseif (name .eq. 'N ' .or. name .eq. ' N') then
         charge = 7.0d0
       elseif (name .eq. 'O ' .or. name .eq. ' O') then 
         charge = 8.0d0
       elseif (name .eq. 'F ' .or. name .eq. ' F') then 
         charge = 9.0d0
       elseif (name .eq. 'Ne') then 
         charge = 10.0d0
       elseif (name .eq. 'Na') then 
         charge = 11.0d0
       elseif (name .eq. 'Mg') then 
         charge = 12.0d0
       elseif (name .eq. 'Al') then 
         charge = 13.0d0
       elseif (name .eq. 'Si') then 
         charge = 14.0d0
       elseif (name .eq. 'P ' .or. name .eq. ' P') then 
         charge = 15.0d0
       elseif (name .eq. 'S ' .or. name .eq. ' S') then
         charge = 16.0d0
       elseif (name .eq. 'Cl') then 
         charge = 17.0d0
       elseif (name .eq. 'Ar') then 
         charge = 18.0d0
       elseif (name .eq. 'K ' .or. name .eq. ' K') then
         charge = 19.0d0
       elseif (name .eq. 'Ca') then
         charge = 20.0d0
       elseif (name .eq. 'Sc') then 
         charge = 21.0d0
       elseif (name .eq. 'Ti') then 
         charge = 22.0d0
       elseif (name .eq. 'V ' .or. name .eq. ' V') then 
         charge = 23.0d0
       elseif (name .eq. 'Cr') then 
         charge = 24.0d0
       elseif (name .eq. 'Mn') then 
         charge = 25.0d0
       elseif (name .eq. 'Fe') then 
         charge = 26.0d0
       elseif (name .eq. 'Co') then 
         charge = 27.0d0
       elseif (name .eq. 'Ni') then 
         charge = 28.0d0
       elseif (name .eq. 'Cu') then 
         charge = 29.0d0
       elseif (name .eq. 'Zn') then 
         charge = 30.0d0
       elseif (name .eq. 'Ga') then 
         charge = 31.0d0
       elseif (name .eq. 'Ge') then 
         charge = 32.0d0
       elseif (name .eq. 'As') then 
         charge = 33.0d0
       elseif (name .eq. 'Se') then 
         charge = 34.0d0
       elseif (name .eq. 'Br') then 
         charge = 35.0d0
       elseif (name .eq. 'Kr') then 
         charge = 36.0d0
       elseif (name .eq. 'Rb') then 
         charge = 37.0d0
       elseif (name .eq. 'Sr') then 
         charge = 38.0d0
       elseif (name .eq. 'Y ' .or. name .eq. ' Y') then 
         charge = 39.0d0
       elseif (name .eq. 'Zr') then 
         charge = 40.0d0
       elseif (name .eq. 'Nb') then 
         charge = 41.0d0
       elseif (name .eq. 'Mo') then 
         charge = 42.0d0
       elseif (name .eq. 'Tc') then 
         charge = 43.0d0
       elseif (name .eq. 'Ru') then 
         charge = 44.0d0
       elseif (name .eq. 'Rh') then 
         charge = 45.0d0
       elseif (name .eq. 'Pd') then 
         charge = 46.0d0
       elseif (name .eq. 'Ag') then
         charge = 47.0d0
       elseif (name .eq. 'Cd') then 
         charge = 48.0d0
       elseif (name .eq. 'In') then 
         charge = 49.0d0
       elseif (name .eq. 'Sn') then 
         charge = 50.0d0
       elseif (name .eq. 'Sb') then 
         charge = 51.0d0
       elseif (name .eq. 'Te') then 
         charge = 52.0d0
       elseif (name .eq. 'I ' .or. name .eq. ' I') then 
         charge = 53.0d0
       elseif (name .eq. 'Xe') then 
         charge = 54.0d0
       elseif (name .eq. 'Cs') then 
         charge = 55.0d0
       elseif (name .eq. 'Ba') then 
         charge = 56.0d0
       elseif (name .eq. 'La') then 
         charge = 57.0d0
       elseif (name .eq. 'Ce') then 
         charge = 58.0d0
       elseif (name .eq. 'Pr') then 
         charge = 59.0d0
       elseif (name .eq. 'Nd') then 
         charge = 60.0d0
       elseif (name .eq. 'Pm') then 
         charge = 61.0d0
       elseif (name .eq. 'Sm') then 
         charge = 62.0d0
       elseif (name .eq. 'Eu') then 
         charge = 63.0d0
       elseif (name .eq. 'Gd') then
         charge = 64.0d0
       elseif (name .eq. 'Tb') then 
         charge = 65.0d0
       elseif (name .eq. 'Dy') then 
         charge = 66.0d0
       elseif (name .eq. 'Ho') then 
         charge = 67.0d0
       elseif (name .eq. 'Er') then 
         charge = 68.0d0
       elseif (name .eq. 'Tm') then 
         charge = 69.0d0
       elseif (name .eq. 'Yb') then 
         charge = 70.0d0
       elseif (name .eq. 'Lu') then 
         charge = 71.0d0
       elseif (name .eq. 'Hf') then 
         charge = 72.0d0
       elseif (name .eq. 'Ta') then 
         charge = 73.0d0
       elseif (name .eq. 'W ' .or. name .eq. ' W') then 
         charge = 74.0d0
       elseif (name .eq. 'Re') then 
         charge = 75.0d0
       elseif (name .eq. 'Os') then
         charge = 76.0d0
       elseif (name .eq. 'Ir') then 
         charge = 77.0d0
       elseif (name .eq. 'Pt') then 
         charge = 78.0d0
       elseif (name .eq. 'Au') then 
         charge = 79.0d0
       elseif (name .eq. 'Hg') then 
         charge = 80.0d0
       elseif (name .eq. 'Tl') then 
         charge = 81.0d0
       elseif (name .eq. 'Pb') then 
         charge = 82.0d0
       elseif (name .eq. 'Bi') then 
         charge = 83.0d0
       elseif (name .eq. 'Po') then 
         charge = 84.0d0
       elseif (name .eq. 'At') then 
         charge = 85.0d0
       elseif (name .eq. 'Rn') then 
         charge = 86.0d0
       elseif (name .eq. 'Fr') then 
         charge = 87.0d0
       elseif (name .eq. 'Ra') then
         charge = 88.0d0
       elseif (name .eq. 'Ac') then 
         charge = 89.0d0
       elseif (name .eq. 'Th') then 
         charge = 90.0d0
       elseif (name .eq. 'Pa') then 
         charge = 91.0d0
       elseif (name .eq. ' U' .or. name .eq. 'U ') then 
         charge = 92.0d0
       elseif (name .eq. 'Np') then 
         charge = 93.0d0
       elseif (name .eq. 'Pu') then 
         charge = 94.0d0
       elseif (name .eq. 'Am') then 
         charge = 95.0d0
       elseif (name .eq. 'Cm') then 
         charge = 96.0d0
       elseif (name .eq. 'Bk') then 
         charge = 97.0d0
       elseif (name .eq. 'Cf') then 
         charge = 98.0d0
       elseif (name .eq. 'Es') then 
         charge = 99.0d0
       elseif (name .eq. 'Fm') then
         charge = 100.0d0
       elseif (name .eq. 'Md') then 
         charge = 101.0d0
       elseif (name .eq. 'No') then 
         charge = 102.0d0
       elseif (name .eq. 'Lr') then 
         charge = 103.0d0
       else
         write(6,100) name
         stop 'charge one'
       endif
 100   format(//,'Element ',a2,' unknown')
       return
       end
