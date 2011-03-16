      function system_time()
      real*8 system_time
cccccccc wall aclock time

ccccccc for T3E
       system_time=rtc()

cccccc  for IBM SP2
cc       system_time=rtc()

      return
      end
