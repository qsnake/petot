      subroutine potrw(vd,r,nr,k,kj,ist)
c
c ***********************************************************
c *                                                         *
c *    This is a plotting routine; the user should          *
c *  adjust/eliminatebfor their own needs.  Prints          *
c *  out the wave functions to the current plot.dat         *
c *  file (unit=3) for later ploting.  A marker (marker)    *
c *  is placed at the end of each group of data.            *
c *                                                         *
c ***********************************************************
c  
      implicit real*8 (a-h,o-z)
c
      parameter (zero=0.0d0,pzf=0.01)
c
      character*3 marker
c
      dimension vd(nr),r(nr)
c
c  Step size of 0.01 is adjustable as seen fit to give 
c  a reasonable plot.
c
      step=zero
      do 150,j=2,nr
        if (r(j) .ge. step) then
          write(3,6000)r(j),vd(j)*ist
          step=step+pzf
        endif
        if (r(j) .gt. 5.0) goto 151
 150  continue   
 151  continue          
      if (kj .eq. 0) then
        if (k .eq. 0) then
          marker='wsp'
        elseif (k .eq. 1) then
          marker='wpp'
        elseif (k .eq. 2) then
          marker='wdp'
        elseif (k .eq. 3) then
          marker='wfp'
        elseif (k .eq. 4) then
          marker='wgp'
        endif
      else
        if (k .eq. 0) then
          marker='wst'
        elseif (k .eq. 1) then
          marker='wpt'
        elseif (k .eq. 2) then
          marker='wdt'
        elseif (k .eq. 3) then
          marker='wft'
        elseif (k .eq. 4) then
          marker='wgt'
        endif
      endif
      write(3,6001)marker
      return
c
c  Format statements
c
 6000 format(1x,f7.4,3x,f18.14)
 6001 format(1x,'marker ',a3)
      end
