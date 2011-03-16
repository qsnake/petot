      function system_time()
      real*8 system_time
cccccccc wall aclock time

ccccccc for T3E
c       system_time=rtc()

cccccc  for IBM SP2
       system_time=rtc()

      return
      end
