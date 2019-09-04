subroutine read_in_v_integrals_short_n_1
   use gf2, only : nao, ncell, ncellgf2, vi_snc, vi_sne, i_vi_snc, i_vi_sne, i_len_c, i_len_e
   use gf2, only : vi_sncl, vi_snel, i_vi_sncl, i_vi_snel
   use gf2, only : vi_s, v_ind, v_ind_1, v_ind_2
   integer :: i0, j0, k0, l0, nf, nf1, nf2, np
   integer*8 :: i, j, k, l
   integer :: nmax_c, nmax_e, nmax_ct, nmax_et, nc, ne
   character*256 :: file_name
   double precision :: val
   nf = nao*ncellgf2 !ncell
   allocate(i_len_c(nf))
   allocate(i_len_e(nf))
   i_len_c = 0
   i_len_e = 0
   nmax_ct = 0
   nmax_et = 0
     open(950,file='2EI_2',status='old',form='formatted')
     rewind(950)
     do
       read(950,*,end=12)
       nmax_ct = nmax_ct + 1
     enddo
12   continue
     close(950)
!     i_len_c(nf1) = np

     open(950,file='2EI_3',status='old',form='formatted')
     rewind(950)
     do
       read(950,*,end=13)
       nmax_et = nmax_et + 1
     enddo
13   continue
     close(950)
!     i_len_e(nf1) = np
!   nmax_c = maxval(i_len_c)
!   nmax_e = maxval(i_len_e)
!   nmax_ct = sum(i_len_c)
!   nmax_et = sum(i_len_e) 
   write(*,*)'nmax_ct, nmax_et = ', nmax_ct, nmax_et
!   write(*,*) nmax_ct, nmax_et 
!   stop
   allocate(vi_sncl(nmax_ct))
   allocate(vi_snel(nmax_et))
   allocate(i_vi_sncl(nmax_ct,4))
   allocate(i_vi_snel(nmax_et,4)) 
   vi_sncl = 0.0d0
   vi_snel = 0.0d0
   i_vi_sncl = 0
   i_vi_snel = 0
   nc = 0
   ne = 0
      open(950,file='2EI_2',status='old',form='formatted') 
      rewind(950)
      do nc = 1, nmax_ct
        read(950,*) i, j, k, l, val
        
        vi_sncl(nc) = val
        i_vi_sncl(nc,1) = i
        i_vi_sncl(nc,2) = j
        i_vi_sncl(nc,3) = k
        i_vi_sncl(nc,4) = l
        i_len_c(j) = i_len_c(j) + 1
      enddo
      close(950)

      open(950,file='2EI_3',status='old',form='formatted')
      rewind(950)
      do ne = 1, nmax_et 
        read(950,*) i, j, k, l, val

        vi_snel(ne) = val
        i_vi_snel(ne,1) = i
        i_vi_snel(ne,2) = j
        i_vi_snel(ne,3) = k
        i_vi_snel(ne,4) = l
        i_len_e(k) = i_len_e(k) + 1
      enddo
      close(950)
!   stop
   allocate(v_ind(nc,4))
   allocate(vi_s(nc))
   allocate(v_ind_1(nao+1))
   v_ind_1 = 0 

   do i = 1, nao
     write(file_name,fmt='(A5,I0)') '2EI0_',i
     open(950+i,file=file_name,action='write',status='replace',form='unformatted') 
     !open(950+i,file=file_name,action='write',status='replace',form='formatted')
   enddo
   do np = 1, nmax_ct
     i = i_vi_sncl(np,1)
     j = i_vi_sncl(np,2)
     k = i_vi_sncl(np,3)
     l = i_vi_sncl(np,4) 
     val  = vi_sncl(np)
     write(950+i) i, j, k, l, val
     !write(950+i,fmt='(4I3XX,F14.12)') i, j, k, l, val
     v_ind_1(i) = v_ind_1(i) + 1
   enddo
   do i = 1, nao
     close(950+i)
   enddo
!   stop
   write(*,*) v_ind_1
   nf1 = 0
   do nc = 1, nao
     write(file_name,fmt='(A5,I0)') '2EI0_',nc
     open(950+nc,file=file_name,action='read',status='old',form='unformatted')
     rewind(950+nc)
     do np = 1, v_ind_1(nc)
       read(950+nc) i, j, k, l, val
       nf1 = nf1 + 1
       vi_s(nf1) = val
       v_ind(nf1,1) = i
       v_ind(nf1,2) = j
       v_ind(nf1,3) = k
       v_ind(nf1,4) = l 
     enddo
     close(950+nc)
   enddo
   allocate(v_ind_2(nao,2))
   v_ind_2 = 0
   do nc = 1, nao
     v_ind_2(nc,1) = sum(v_ind_1(1:nc-1)) + 1 
     v_ind_2(nc,2) = v_ind_2(nc,1) + v_ind_1(nc) - 1
   enddo
end subroutine read_in_v_integrals_short_n_1
