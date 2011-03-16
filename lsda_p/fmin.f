c-----------------------------------------------------------------------
      function fmin(vec,nel)
c-----------------------------------------------------------------------
      implicit none
c
      integer fmin,iel,nel
      integer vec(nel)
c-----------------------------------------------------------------------
      fmin=1
      do iel=1,nel
        if(vec(iel).lt.vec(fmin)) fmin=iel
      enddo
c
      return
      end
