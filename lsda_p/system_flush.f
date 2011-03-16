      subroutine system_flush(ifile)
      integer ifile
ccccccc flush the buffer of output file: ifile

ccccc for T3E
      call flush(ifile)
ccccc for IBM SP
cc      call flush_(ifile)

      return
      end
