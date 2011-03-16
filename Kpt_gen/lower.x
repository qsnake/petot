      subroutine lower (n, string)
c     
c     make string lower case
c     $Id: lower.x,v 1.1 89/03/11 19:08:58 sverre Exp $
c
c     $Log:	lower.x,v $
c     Revision 1.1  89/03/11  19:08:58  sverre
c     Initial revision
c     
c     Revision 1.1  89/03/11  19:08:58  sverre
c     Initial revision
c     
c     
      integer n
      character string*(*)
c     
c     input:
c     
c     n        number of characters in string
c     string   mixed case string
c     
c     output:
c     
c     string   lower case string
c     
c     local variables:
c     
      integer i
c
c     rcs id string - allows use of ident to identify binaries
c
      character rcsid*50
      rcsid = '$RCSfile: lower.x,v $$Revision: 1.1 $'
c     
c     convert string to lower case
c     
      do 100 i=1,n
         if (string(i:i) .ge. 'A' .and.
     +        string(i:i) .le. 'Z') then
            string(i:i) = char(ichar(string(i:i)) + 32)
         end if
  100 continue
      return
      end
