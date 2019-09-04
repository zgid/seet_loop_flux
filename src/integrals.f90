subroutine read_all_integrals(vint,nf)
   integer :: nf
   integer*8 :: i, j, k, l
   integer :: ii, jj, kk, ll
   real*8 :: val
   double precision, dimension(nf,nf,nf,nf) :: vint

   vint = 0.0d0
   open(950,file='2EIFC',status='old',form='unformatted')
   rewind(950)
   do
     read(950,end=12) i, j, k, l, val
     ii = int(i)
     jj = int(j)
     kk = int(k)
     ll = int(l)
     vint(ii,jj,kk,ll) = dble(val)
   enddo
12 continue
   close(950)
   write(*,*)'integrals successfully read'
end subroutine read_all_integrals
subroutine read_in_v_integrals_short_n
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
   do nf1 = 1, nf
     write(file_name,fmt='(A4,I0,A1,I0)') '2EI_',nf1,'_',12
     open(950,file=file_name,status='old',form='unformatted')
     rewind(950)
     np = 0
     do
       read(950,end=12)
       np = np + 1
     enddo
12   continue
     close(950)
     i_len_c(nf1) = np

     write(file_name,fmt='(A4,I0,A1,I0)') '2EI_',nf1,'_',13
     open(950,file=file_name,status='old',form='unformatted')
     rewind(950)
     np = 0
     do
       read(950,end=13)
       np = np + 1
     enddo
13   continue
     close(950)
     i_len_e(nf1) = np
   enddo
   nmax_c = maxval(i_len_c)
   nmax_e = maxval(i_len_e)
   nmax_ct = sum(i_len_c)
   nmax_et = sum(i_len_e) 
   write(*,*)'nmax_c, nmax_e = ', nmax_c, nmax_e
   write(*,*) nmax_ct, nmax_et 
   
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
   do nf1 = 1, nf 
      write(file_name,fmt='(A4,I0,A1,I0)') '2EI_',nf1,'_',12 
      open(950,file=file_name,status='old',form='unformatted') 
      rewind(950)
      do np = 1, i_len_c(nf1)
        nc = nc + 1
        read(950) i, j, k, l, val
        
        vi_sncl(nc) = val
        i_vi_sncl(nc,1) = i
        i_vi_sncl(nc,2) = j
        i_vi_sncl(nc,3) = k
        i_vi_sncl(nc,4) = l
      enddo
      close(950)

      write(file_name,fmt='(A4,I0,A1,I0)') '2EI_',nf1,'_',13
      open(950,file=file_name,status='old',form='unformatted')
      rewind(950)
      do np = 1, i_len_e(nf1)
        ne = ne + 1
        read(950) i, j, k, l, val

        vi_snel(ne) = val
        i_vi_snel(ne,1) = i
        i_vi_snel(ne,2) = j
        i_vi_snel(ne,3) = k
        i_vi_snel(ne,4) = l
      enddo
      close(950)
   enddo

   allocate(v_ind(nc,4))
   allocate(vi_s(nc))
   allocate(v_ind_1(nao+1))
   v_ind_1 = 0 

   do i = 1, nao
     write(file_name,fmt='(A5,I0)') '2EI0_',i
     !open(950+i,file=file_name,action='write',status='replace',form='unformatted') 
     open(950+i,file=file_name,action='write',status='replace',form='formatted')
   enddo
   do np = 1, nc
     i = i_vi_sncl(np,1)
     j = i_vi_sncl(np,2)
     k = i_vi_sncl(np,3)
     l = i_vi_sncl(np,4) 
     val  = vi_sncl(np)
     !write(950+i) i, j, k, l, val
     write(950+i,fmt='(4I4X,F15.12)') i, j, k, l, val
     v_ind_1(i) = v_ind_1(i) + 1
   enddo
   do i = 1, nao
     close(950+i)
   enddo
   stop
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
end subroutine read_in_v_integrals_short_n

subroutine integrals_ce_n_emb
   use gf2, only : nao, ncell, ncellgf2, vi_snc, vi_sne, i_vi_snc, i_vi_sne, i_len_c, i_len_e, n_cell_p, i_cell_p, sumcell, n_vi_c, n_vi_e
   use gf2, only : vi_sncl, vi_snel, i_vi_sncl, i_vi_snel
   integer :: nf, nf1, n, n0, i0, j0, k0, l0, g_j1, g_k1, g_l1, g_j, g_k, g_l, g_n, ic, jc, np, np1, n1, i1, j1, k1, l1
   integer :: len1, nmaxf_c, nmaxf_e, nmax_ct, nmax_et
   integer*8 :: i, j, k, l
   !double precision :: val
   real*8 :: val
!   real*8 :: val
   character*256 :: file_name
   logical :: exist
   nf = nao*ncellgf2 !ncell
   nf2 = nf !nao*73
   allocate(n_vi_c(nf))
   allocate(n_vi_e(nf))

   n_vi_c = 0
   n_vi_e = 0

   nmax_ct = 0
   nmax_et = 0
     open(950,file='2EIFC',status='old',form='unformatted')
     rewind(950)
     do
       read(950,end=12) i, j, k, l, val
       nmax_ct = nmax_ct + 1
     enddo
12   continue
     close(950)
   open(951,file='2EIFC',status='old',form='unformatted')
   rewind(951)
   do n = 1, nmax_ct
     read(951) i, j, k, l, val
     n_vi_c(j) = n_vi_c(j) + 1
   enddo
   close(951)

   open(951,file='2EIFE',status='old',form='unformatted')
   do n = 1, nmax_ct
     read(951) i, j, k, l, val
     n_vi_e(k) = n_vi_e(k) + 1
   enddo
   close(951)
   nmaxf_c = maxval(n_vi_c)
   nmaxf_e = maxval(n_vi_e)
   write(*,*)'nmaxf_c, nmaxf_e = ', nmaxf_c, nmaxf_e
end subroutine integrals_ce_n_emb

subroutine integrals_ce_n_emb1
   use gf2, only : nao, ncell, ncellgf2, vi_snc, vi_sne, i_vi_snc, i_vi_sne, i_len_c, i_len_e, n_cell_p, i_cell_p, sumcell, n_vi_c, n_vi_e
   use gf2, only : vi_sncl, vi_snel, i_vi_sncl, i_vi_snel 
   integer :: nf, nf1, i, j, k, l, n, n0, i0, j0, k0, l0, g_j1, g_k1, g_l1, g_j, g_k, g_l, g_n, ic, jc, np, np1, n1, i1, j1, k1, l1
   integer :: len1, nmaxf_c, nmaxf_e, nmax_ct, nmax_et
   double precision :: val
!   real*8 :: val
   character*256 :: file_name
   logical :: exist
   nf = nao*ncellgf2 !ncell
   nf2 = nf !nao*73
   allocate(n_vi_c(nf))
   allocate(n_vi_e(nf))

   n_vi_c = 0
   n_vi_e = 0

   nmax_ct = 0
   nmax_et = 0
     open(950,file='2EIFC',status='old',form='unformatted')
     rewind(950)
     do
       read(950,end=12)
       nmax_ct = nmax_ct + 1
     enddo
12   continue
     close(950)

   open(950,file="2EIFC",status='old',action='read',form='unformatted')
   do n = 1, nmax_ct
     read(950) i, j, k, l, val
     !if(i.le.nf2) then
     !write(950) i, j, k, l, val
     n_vi_c(j) = n_vi_c(j) + 1
     !endif
   enddo
   close(950)

   open(950,file="2EIFE",status='old',action='read',form='unformatted')
   do n = 1, nmax_ct
     read(950) i, j, k, l, val
     !if(i.le.nf2) then
     !write(950) i, j, k, l, val
     n_vi_e(k) = n_vi_e(k) + 1
     !endif
   enddo
   close(950)

   nmaxf_c = maxval(n_vi_c)
   nmaxf_e = maxval(n_vi_e)
   write(*,*)'nmaxf_c, nmaxf_e = ', nmaxf_c, nmaxf_e
end subroutine integrals_ce_n_emb1


subroutine integrals_ce_n
   use gf2, only : nao, ncell, ncellgf2, vi_snc, vi_sne, i_vi_snc, i_vi_sne, i_len_c, i_len_e, n_cell_p, i_cell_p, sumcell, n_vi_c, n_vi_e
   use gf2, only : vi_sncl, vi_snel, i_vi_sncl, i_vi_snel 
   integer :: nf, nf1, nf2, i, j, k, l, n, n0, i0, j0, k0, l0, g_j1, g_k1, g_l1, g_j, g_k, g_l, g_n, ic, jc, np, np1, n1, i1, j1, k1, l1
   integer :: len1, nmaxf_c, nmaxf_e
   double precision :: val
!   real*8 :: val
   character*256 :: file_name
   logical :: exist
   nf = nao*ncellgf2 !ncell
   nf2 = nf 
   allocate(n_vi_c(nf))
   allocate(n_vi_e(nf))

   n_vi_c = 0
   n_vi_e = 0

   open(950,file="2EIFC",status='replace',action='write',form='unformatted')
   !open(950,file="2EIFC",status='replace',action='write',form='formatted')
   do n = 1, nf
     g_n = (n-1)/nao + 1
     n0  = n - (g_n-1)*nao
     np  = n_cell_p(g_n)
     do np1 = 1, np
       ic = i_cell_p(g_n,np1,1)
       jc = i_cell_p(g_n,np1,2)
       n1 = n0 + (ic-1)*nao
       len1 = i_len_c(n1)
       ic1 = 0
       ic1 = sum(i_len_c(1:n1-1)) + 1
       ic2 = ic1 + len1 - 1
       do nf1 = ic1, ic2 

         i = i_vi_sncl(nf1,1)
         k = i_vi_sncl(nf1,3)
         l = i_vi_sncl(nf1,4)
         val = vi_sncl(nf1)

         i1 = i + (jc-1)*nao


         g_k = (k-1)/nao + 1
         g_l = (l-1)/nao + 1

         g_k1 = sumcell(jc,g_k)
         g_l1 = sumcell(jc,g_l)
      
         if(g_k1*g_l1.ne.0) then
           i1 = i + (jc-1)*nao


           k0 = k - (g_k-1)*nao
           k1 = k0 + (g_k1-1)*nao

           l0 = l - (g_l-1)*nao
           l1 = l0 + (g_l1-1)*nao
 
           !if(i1.gt.nf.or.n.gt.nf.or.k1.gt.nf.or.l1.gt.nf) then
           !  write(*,*) i1, n, k1, l1
           !  stop
           !endif
           if(i1.le.nf2) then 
           write(950) i1, n, k1, l1, val
           n_vi_c(n) = n_vi_c(n) + 1
           endif
         endif 

       enddo 
     enddo 
   enddo
   close(950)
   open(950,file="2EIFE",status='replace',action='write',form='unformatted')
   do n = 1, nf
     g_n = (n-1)/nao + 1
     n0  = n - (g_n-1)*nao
     np  = n_cell_p(g_n)
     do np1 = 1, np
       ic = i_cell_p(g_n,np1,1)
       jc = i_cell_p(g_n,np1,2)
       n1 = n0 + (ic-1)*nao
       len1 = i_len_e(n1)
       ic1 = 0
       ic1 = sum(i_len_e(1:n1-1)) + 1
       ic2 = ic1 + len1 - 1
       do nf1 = ic1, ic2 

         i = i_vi_snel(nf1,1)
         j = i_vi_snel(nf1,2)
         k = i_vi_snel(nf1,3)
         l = i_vi_snel(nf1,4)
         val = vi_snel(nf1)

         i1 = i + (jc-1)*nao

         g_j = (j-1)/nao + 1

         g_l = (l-1)/nao + 1

         g_l1 = sumcell(jc,g_l)
         g_j1 = sumcell(jc,g_j)

         if(g_j1*g_l1.ne.0) then
           i1 = i + (jc-1)*nao

           j0 = j - (g_j-1)*nao
           j1 = j0 + (g_j1-1)*nao


           l0 = l - (g_l-1)*nao
           l1 = l0 + (g_l1-1)*nao

           if(i1.le.nf2) then
           write(950) i1, j1, n, l1, val 
           n_vi_e(n) = n_vi_e(n) + 1
           endif
         endif

       enddo
     enddo
   enddo 
   close(950)
   nmaxf_c = maxval(n_vi_c)
   nmaxf_e = maxval(n_vi_e)
   write(*,*)'nmaxf_c, nmaxf_e = ', nmaxf_c, nmaxf_e
end subroutine integrals_ce_n

