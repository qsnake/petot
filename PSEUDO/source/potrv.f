      subroutine potrv(vd,r,nr,k)
c
c ***********************************************************
c *                                                         *
c *    This is a plotting routine; the user should          *
c *  adjust for their own needs.  Prints                    *
c *  out the potential to the current plot.dat              *
c *  file (unit=3) for later ploting.  A marker (marker)    *
c *  is placed at the end of each group of data.            *
c *                                                         *
c ***********************************************************
c  
      implicit real*8 (a-h,o-z)
c
      character*3 marker
c
      dimension vd(nr),r(nr) 
c 
c  Step size of 0.05 is adjustable as seen fit to give 
c  a reasonalble plot.
c
      step=0.0
      do 150,j=2,nr
c        if (r(j) .ge. step) then
          write(3,6000)r(j),vd(j)
c          step=step+0.05
c        endif
 150  continue             
      if (k .eq. 0) then
        marker='vns'
      elseif (k .eq. 1) then
        marker='vnp'
      elseif (k .eq. 2) then
        marker='vnd'
      elseif (k .eq. 3) then
        marker='vnf'
      elseif (k .eq. 4) then
        marker='vng'
      endif
      write(3,6001)marker
      return
c
c  Format statements
c
 6000 format(1x,f12.8,3x,f18.13)
 6001 format(1x,'marker ',a3)
      end
