subroutine timest(Line)
  implicit none
  character*(*) Line
  character(10) :: dateinfo, dateinfo1
  character(4) :: year*4, month*2, day*2
  character(10) :: timeinfo, timeinfo1
  character(2) :: hour, minute, second*2
  call date_and_time(dateinfo,timeinfo)
  hour = timeinfo(1:2)
  minute = timeinfo(3:4)
  second = timeinfo(5:6)
  year = dateinfo(1:4)
  month = dateinfo(5:6)
  day = dateinfo(7:8)
  timeinfo1 = hour // ':' // minute // ':' // second
  dateinfo1 = day // '.' // month // '.' // year
  write(*,*) Line,' at ', timeinfo1, ' on ', dateinfo1
  call flush(6)
end subroutine timest 
