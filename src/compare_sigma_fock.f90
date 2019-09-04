subroutine comp_sig_fock(F_kk,nts)
use gf2, only : nk, nao, ncell, nel_cell, ncellgf2    
use gf2, only : Sigma_tau
implicit none
integer nts
complex*16, dimension(nao,nao,2*nk) :: F_kk
complex*16, dimension(:,:,:),allocatable:: Xr 
integer it,i,j,ic,itt, icc, ii, jj

allocate(Xr(nao,nao,ncell))           

Xr=dcmplx(0.0d0,0.0d0) 
do ic=1, ncell
   call FT_k_to_r_ic(F_kk, Xr(:,:,ic),ic)
end do


     open(1,file="SE_nonsym.txt")
!     open(1,file="SE_nonsym.txt",status="replace")
     do it = 1, nts
        do ic = 1, ncellgf2
           do i = 1, nao
              do j = 1, nao
                read(1,*) itt, icc, ii, jj, Sigma_tau(i,j,ic,it)
!                 write(unit=1,fmt='(4I4XX,F15.10X)') it, ic, i, j, Sigma_tau(i,j,ic,it)
              enddo
           enddo
        enddo
     enddo
     close(1)

     !If(symse.gt.0) then
     !  symmetrize the central cell, then other cells
     write(*,*)'symmetrizing the self-energy'
     call flush(6)
     do it = 1, nts
        do i = 1, nao
           do j = 1, i
              Sigma_tau(i,j,1,it) = 5.0d-01*(Sigma_tau(i,j,1,it)+Sigma_tau(j,i,1,it))
              Sigma_tau(j,i,1,it) = Sigma_tau(i,j,1,it)
           enddo
        enddo
     enddo
     do it = 1, nts
        do ic = 1, ncellgf2/2
           do i = 1, nao
              do j = 1, i
                 Sigma_tau(i,j,ic*2,it) = 5.0d-01*(Sigma_tau(i,j,ic*2,it)+Sigma_tau(j,i,ic*2+1,it))
                 Sigma_tau(j,i,ic*2+1,it) = Sigma_tau(i,j,ic*2,it)
                 
                 Sigma_tau(j,i,ic*2,it) = 5.0d-01*(Sigma_tau(j,i,ic*2,it)+Sigma_tau(i,j,ic*2+1,it))
                 Sigma_tau(i,j,ic*2+1,it) = Sigma_tau(j,i,ic*2,it)
              enddo
           enddo
        enddo
     enddo

!Now compare with Fock matrix

     do it = nts,nts
        do ic = 1, ncell
           do i = 1, nao
              do j = 1, nao
                write(*,*) ic, i, j, dble(Sigma_tau(i,j,ic,it)),dble(Xr(i,j,ic))
!                 write(unit=1,fmt='(4I4XX,F15.10X)') it, ic, i, j, Sigma_tau(i,j,ic,it)
              enddo
           enddo
        enddo
     enddo


deallocate(Xr)


end subroutine comp_sig_fock
