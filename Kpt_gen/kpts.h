c Parameters for kpts-program
c max number of species and atoms
      integer nsmax,namax
      parameter(nsmax=20,namax=300)      
c max number of kpoints
      integer nkptx
      parameter (nkptx=3000)
      real*8 Pi
      parameter (Pi= 3.1415)
c max number of shells that are analyzed, dimensions of the 
c cube that contains the k-points of the shells
      integer n_shells, i1x,i2x,i3x
      parameter(i1x=30)
      parameter(i2x=30)
      parameter(i3x=30)
      parameter(n_shells = (min(i1x,i2x,i3x)/2)**3)

