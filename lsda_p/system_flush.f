      subroutine system_flush(ifile)
      integer ifile
ccccccc flush the buffer of output file: ifile

ccccc for T3E
c      call flush(ifile)
ccccc for IBM SP
      call flush_(ifile)

      return
      end
