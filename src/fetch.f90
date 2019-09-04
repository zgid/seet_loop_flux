subroutine fetch_r(x,aa,iplt) ! routine to read double precision variable
  ! Also read character string and
  ! checks that it equals to input string.
  implicit none
  integer iplt
  double precision x
  character*8 aa, bb
  
  write(6,*)' now reading ',aa,' from ',iplt; call flush(6)
  read(iplt, *) bb, x
  if(bb/=aa) then
     write(6,*)' problem; supposed read ',aa,' instead reads ',bb; 
     call flush(6)
     stop
  endif
  write(6,*)' = ',x
end subroutine fetch_r

subroutine fetch_i(ii,aa,iplt) ! same for integer.

  integer ii, iplt
  character*8 aa, bb
  
  write(6,*)' now reading ',aa,' from ',iplt; call flush(6)
  read(iplt, *) bb, ii
  if(bb/=aa) then
     write(6,*)' problem; supposed read ',aa,' instead reads ',bb
     call flush(6)
     stop
  endif
  write(6,*)' = ',ii
end subroutine fetch_i


subroutine fetch_i_bath(ii,jj,aa,iplt) ! same for integer.

  integer ii, jj,iplt
  character*8 aa, bb
  
  write(6,*)' now reading ',aa,' from ',iplt; call flush(6)
  read(iplt, *) bb, ii, jj
  if(bb/=aa) then
     write(6,*)' problem; supposed read ',aa,' instead reads ',bb
     call flush(6)
     stop
  endif
  write(6,*)' = ',ii,jj
end subroutine fetch_i_bath

subroutine fetch_L(L,aa,iplt) 
  implicit none
  logical L
  integer iplt
  character*8 aa, bb
  
  write(6,*)' now reading ',aa,' from ',iplt
  call flush(6)

  read(iplt, *) bb, L
  if(bb/=aa) then
     write(6,*)' problem; supposed read ',aa,' instead reads ',bb
     call flush(6)
     stop
  endif
  write(6,*)' = ',L
end subroutine fetch_L
