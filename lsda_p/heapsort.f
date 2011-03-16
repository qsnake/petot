c-----------------------------------------------------------------------
      subroutine heapsort(vref,n,v1,v2,v3)
c-----------------------------------------------------------------------
      implicit none
c
c     .. Scalars
      integer i,ir,j,l,n
      integer tmpref,tmp1,tmp2,tmp3
c
c     .. Arrays
      integer vref(n),v1(n),v2(n),v3(n)
c-----------------------------------------------------------------------
      l=n/2+1
      ir=n
   50 continue
      if (l.gt.1) then
         l=l-1
         tmpref=vref(l)
         tmp1=v1(l)
         tmp2=v2(l)
         tmp3=v3(l)
      else
         tmpref=vref(ir)
         tmp1=v1(ir)
         tmp2=v2(ir)
         tmp3=v3(ir)
         vref(ir)=vref(1)
         v1(ir)=v1(1)
         v2(ir)=v2(1)
         v3(ir)=v3(1)
         ir=ir-1
         if (ir.eq.1) then
            vref(ir)=tmpref
            v1(ir)=tmp1
            v2(ir)=tmp2
            v3(ir)=tmp3
            goto 70
         endif
      endif
      i=l
      j=l+l
   60 if (j.le.ir) then
         if (j.lt.ir) then
            if (vref(j).lt.vref(j+1)) j=j+1
         endif
         if (tmpref.lt.vref(j)) then
            vref(i)=vref(j)
            v1(i)=v1(j)
            v2(i)=v2(j)
            v3(i)=v3(j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
         goto 60
      endif
      vref(i)=tmpref
      v1(i)=tmp1
      v2(i)=tmp2
      v3(i)=tmp3
      goto 50
   70 continue
c
      return
      end            
