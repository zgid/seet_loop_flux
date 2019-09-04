subroutine dffcll
  use gf2, only: ncell, indices_cell, diffcell
  implicit none
  integer :: ic, jc, g0
  integer,dimension(3) :: iindex, jindex, tmpindex 

  allocate(diffcell(ncell,ncell))
  call read_in_indices_cell
  do ic = 1, ncell
    iindex = indices_cell(ic,:)
    do jc = 1, ncell
      jindex = indices_cell(jc,:)
      tmpindex(:) = jindex(:) - iindex(:)
      call ind_search(tmpindex,1,ncell,g0)
      diffcell(ic,jc) = g0
      if(g0.lt.0) then
         diffcell(ic,jc) = 0
      endif
    enddo
  enddo
end subroutine dffcll

subroutine addcll
  use gf2, only: ncellgf2, indices_cell, sumcell, i_cell_p, n_cell_p
  implicit none
  integer :: ic, jc, g0, g1, ncell
  integer,dimension(3) :: iindex, jindex, tmpindex

  ncell = ncellgf2

  allocate(sumcell(ncell,ncell))
  allocate(i_cell_p(ncell,ncell,2))
  allocate(n_cell_p(ncell))
  i_cell_p = 0
  n_cell_p = 0
  call read_in_indices_cell
  do ic = 1, ncell
    iindex = indices_cell(ic,:)
    do jc = 1, ncell
      jindex = indices_cell(jc,:)
      tmpindex(:) = jindex(:) + iindex(:)
      call ind_search(tmpindex,1,ncell,g0)
      sumcell(ic,jc) = g0
      if(g0.lt.0) then
         sumcell(ic,jc) = 0
      endif
      if(g0.gt.0) then
        n_cell_p(g0) = n_cell_p(g0) + 1
        g1 = n_cell_p(g0)
        i_cell_p(g0,g1,1) = ic
        i_cell_p(g0,g1,2) = jc
      endif
    enddo
  enddo
end subroutine addcll

subroutine read_in_indices_cell
  use gf2, only: ncell,indices_cell
  implicit none
  integer:: i,ic,ix,iy,iz

  if (.not. allocated(indices_cell)) then
    allocate(indices_cell(ncell,3))
  end if

  open(45,file='Integer_indices.txt',status='old')
  do i=1,ncell
    read(45,*)ic,ix,iy,iz
    indices_cell(ic,1)=ix
    indices_cell(ic,2)=iy
    indices_cell(ic,3)=iz
  end do
  close(45)
end subroutine read_in_indices_cell


subroutine ind_search(index,non_zero_cell_s,non_zero_cell_e,ind)
  use gf2, only: nao,ncell,indices_cell
  implicit none
  integer::non_zero_cell_s,non_zero_cell_e,ind
  integer,dimension(3)::index
  integer::kdiff,ic
  ind=-1000
  do ic=non_zero_cell_s,non_zero_cell_e
    kdiff=abs(index(1)-indices_cell(ic,1))+abs(index(2)-indices_cell(ic,2))+abs(index(3)-indices_cell(ic,3))
    if (kdiff==0) then
      ind=ic
      exit
    end if
 end do

end subroutine ind_search
