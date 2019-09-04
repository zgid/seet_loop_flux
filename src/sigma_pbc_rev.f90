subroutine sigma_outer(Gr_full_tau,nts)
  use gf2, only: nao, ncell, ncellgf2, Sigma_tau, n_vi_c, n_vi_e, n_threads, diffcell
  use OMP_LIB
  implicit none
  integer :: nts, t, tt, nf, nfgf2, inmax_c, inmax_e, nmax_c, nmax_e, ifit, ntot_c, ntot_e
  double precision :: t1, t2, td
  double precision, dimension(nao,nao,ncell,nts) :: Gr_full_tau
  integer, parameter :: nmax = 500000000
  double precision, allocatable, dimension(:,:) :: vi_c1, vi_e1
  integer, allocatable, dimension(:,:,:) :: i_vi_c1, i_vi_e1
  integer :: i, j, k, l, i1, j1, k1, l1, np, np1, np2, nbatch, nbatch1, ncycle, nres, nc, n1, n2
  double precision :: val, val1
  double precision, allocatable, dimension(:,:) :: timings
  double precision, dimension(10) :: timings_t

  double precision, allocatable, dimension(:,:,:) :: A_cont, B_cont
  double precision, allocatable, dimension(:,:) :: G_tmp, Gb_tmp
  
  double precision, allocatable, dimension(:,:) :: Aval
  integer, allocatable, dimension(:,:) :: Arow, Acol
  integer :: p, q, p0, q0, g_p, g_q, g_q0

  inmax_c = maxval(n_vi_c)
  inmax_e = maxval(n_vi_e)
  nf = nao*ncell
  nfgf2 = nao*ncellgf2
!  nmax_c = inmax_c*nf
!  nmax_e = inmax_e*nf
  ifit = (inmax_c+inmax_e)*24
  nbatch = min(nf,nmax/ifit)

  nbatch = n_threads

  write(*,*)'inmax_c, inmax_e = ', inmax_c, inmax_e
!  write(*,*)'ifit, ratio, nbatch = ', ifit, nmax/ifit, nbatch

  allocate(vi_c1(inmax_c,nbatch))
  allocate(vi_e1(inmax_e,nbatch))
  allocate(i_vi_c1(4,inmax_c,nbatch))
  allocate(i_vi_e1(4,inmax_e,nbatch))
  allocate(timings(10,n_threads))

  allocate(A_cont(nfgf2,nfgf2,n_threads))
  allocate(B_cont(nfgf2,nfgf2,n_threads))
  allocate(G_tmp(nfgf2,nfgf2))
  allocate(Gb_tmp(nfgf2,nfgf2))

  allocate(Aval(nfgf2*nfgf2,n_threads))
  allocate(Arow(nfgf2*nfgf2,n_threads))
  allocate(Acol(nfgf2*nfgf2,n_threads))

  A_cont = 0.0d0
!  B_cont = 0.0d0
  G_tmp  = 0.0d0
  Gb_tmp = 0.0d0

  ncycle = nfgf2/nbatch
  nres = nfgf2 - ncycle*nbatch
  write(*,*)'ncycle, nres = ', ncycle, nres
  if(nres.gt.0) ncycle = ncycle + 1

  if (.not. allocated(Sigma_tau)) then
     allocate(Sigma_tau(nao,nao,ncell,nts))
  end if
  Sigma_tau(:,:,:,:)=0.0d0

! do one tau-point at a time, parallelization is over i and n in sigma_inner

!  write(*,*)'nts = ', nts
  call timest('Begin GF2 contractions') 
  open(950,file="2EIFC",status='old',form='unformatted')
  open(951,file="2EIFE",status='old',form='unformatted')
  rewind(950)
  rewind(951)
  do nc = 1, ncycle
    n1 = (nc-1)*nbatch + 1
    n2 = nc*nbatch
    nbatch1 = nbatch
    if(nc.eq.ncycle.and.nres.gt.0) then
      n2 = n1 + nres - 1
      nbatch1 = nres
    endif
    write(*,*)'n1, n2, nbatch1 = ', n1, n2, nbatch1
    
    
    t1 = omp_get_wtime()
    nbatch1 = 0
    do np = n1, n2
      nbatch1 = nbatch1 + 1
      do np1 = 1, n_vi_c(np)
        read(950) i, j, k, l, val
        vi_c1(np1,nbatch1) = val
        i_vi_c1(1,np1,nbatch1) = i
        i_vi_c1(2,np1,nbatch1) = j
        i_vi_c1(3,np1,nbatch1) = k
        i_vi_c1(4,np1,nbatch1) = l 
      enddo
    enddo 

    nbatch1 = 0
    do np = n1, n2
      nbatch1 = nbatch1 + 1
      do np1 = 1, n_vi_e(np)
        read(951) i, j, k, l, val
        vi_e1(np1,nbatch1) = val
        i_vi_e1(1,np1,nbatch1) = i
        i_vi_e1(2,np1,nbatch1) = j
        i_vi_e1(3,np1,nbatch1) = k
        i_vi_e1(4,np1,nbatch1) = l
      enddo
    enddo       
    t2 = omp_get_wtime()
    write(*,*) 'time copying:', t2-t1
     
  timings = 0.0d0
  do t=1, nts
     tt=nts-t+1
  do p = 1, nfgf2
   g_p = (p-1)/nao + 1
    p0 = p - (g_p-1)*nao
    do q = 1, nfgf2
      g_q = (q-1)/nao + 1
      q0 = q - (g_q-1)*nao
      g_q0 = diffcell(g_p,g_q)
      if(g_q0.gt.0) then 
        Gb_tmp(p,q) = Gr_full_tau(p0,q0,g_q0,tt)
        G_tmp(p,q) = Gr_full_tau(p0,q0,g_q0,t)
      endif
    enddo
 enddo
     call sigma_inner(Sigma_tau(:,:,:,t),Gr_full_tau(1,1,1,tt),Gr_full_tau(1,1,1,t),vi_c1,vi_e1,i_vi_c1,i_vi_e1,n_vi_c,n_vi_e,inmax_c,inmax_e,nfgf2,nbatch,n1,n2,timings,A_cont,B_cont,G_tmp,Gb_tmp,Aval,Arow,Acol)!, &
!     A_ind,B_ind,A_tmp,B_tmp)
  end do
  enddo
  close(950)
  close(951)
  write(*,*)'Sigma formation timing breakdown:'
  timings_t = 0.0d0
  do n1 = 1, 10
    do n2 = 1, n_threads
      timings_t(n1) = timings_t(n1) + timings(n1,n2)
    enddo
  enddo
  write(*,*)'Overhead: ', timings_t(1)
  write(*,*)'1st contractions: ', timings_t(2)
  write(*,*)'1st DGEMM: ', timings_t(3)
  write(*,*)'2nd DGEMM: ', timings_t(4)
  write(*,*)'Coulomb contraction: ', timings_t(5)
  write(*,*)'Exchange contraction: ', timings_t(6)
  write(*,*)'Total time contracting: ', sum(timings_t(1:6))

  deallocate(A_cont)
  deallocate(B_cont)
  deallocate(G_tmp)
  deallocate(Gb_tmp)

  deallocate(vi_c1)
  deallocate(vi_e1)
  deallocate(i_vi_c1)
  deallocate(i_vi_e1)
  deallocate(timings)

  deallocate(Aval)
  deallocate(Arow)
  deallocate(Acol)
  

end subroutine sigma_outer

subroutine sigma_inner(sig,Gbt,Gt,vi_c,vi_e,i_vi_c,i_vi_e,n_vi_c,n_vi_e,nmax_c,nmax_e,nf,nbatch,n1,n2,timings,A_cont,B_cont,G_tmp,Gb_tmp,Aval,Arow,Acol)!,A_ind,B_ind,A_tmp,B_tmp) 
  use gf2, only: nao, ncell, ncellgf2, n_threads
  double precision, dimension(nao,nao,ncell) :: sig, Gbt, Gt
  double precision, dimension(nf,nf,n_threads) :: A_cont,B_cont
  double precision, dimension(nf,nf) :: G_tmp,Gb_tmp
  double precision, allocatable, dimension(:,:,:,:) :: sig_tmpc, sig_tmpe

  double precision, dimension(nf*nf,n_threads) :: Aval
  integer, dimension(nf*nf,n_threads) :: Arow, Acol

  integer nf, i_n, tn, nmax_c, nmax_e, nbatch, n1, n2
  double precision, dimension(nmax_c,nbatch) :: vi_c
  double precision, dimension(nmax_e,nbatch) :: vi_e  
  integer, dimension(4,nmax_c,nbatch) :: i_vi_c
  integer, dimension(4,nmax_e,nbatch) :: i_vi_e
  integer, dimension(nf) :: n_vi_c
  integer, dimension(nf) :: n_vi_e
  double precision, dimension(10,n_threads) :: timings
     
  INTEGER TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM 

  allocate(sig_tmpc(nao,nao,ncell,n_threads))
  allocate(sig_tmpe(nao,nao,ncell,n_threads)) 
  call omp_set_num_threads(n_threads)
  sig_tmpc(:,:,:,:) = 0.0d0
  sig_tmpe(:,:,:,:) = 0.0d0
!  do i = 1, nao
!$OMP PARALLEL DO DEFAULT(none) & 
!$OMP SHARED(sig_tmpc,sig_tmpe,Gbt,Gt,nao,nf,ncell,ncellgf2,vi_c,vi_e,i_vi_c,i_vi_e,n_vi_c,n_vi_e,nmax_c,nmax_e,n1,n2,nbatch,timings,A_cont,B_cont,G_tmp,Gb_tmp,Aval,Arow,Acol) &
!$OMP PRIVATE(n,i_n,tn,i)   
   do n = n1, n2
   do i = 1, nao
     tn = omp_get_thread_num()
     call contract_in(i,n,nao,nf,ncell,ncellgf2,tn,Gt,Gbt,sig_tmpc(:,:,:,tn+1),sig_tmpe(:,:,:,tn+1),vi_c,vi_e,i_vi_c,i_vi_e,n_vi_c,n_vi_e,nmax_c,nmax_e,n1,n2,nbatch, &
     timings(:,tn+1),A_cont(:,:,tn+1),B_cont(:,:,tn+1),G_tmp,Gb_tmp,Aval(:,tn+1),Arow(:,tn+1),Acol(:,tn+1))!,A_ind(:,:,tn+1),B_ind(:,:,tn+1),A_tmp(:,tn+1),B_tmp(:,tn+1))
     enddo
  enddo
  do tn = 1, n_threads
    sig(:,:,:) = sig(:,:,:) + 2.0d0*sig_tmpc(:,:,:,tn) + sig_tmpe(:,:,:,tn)
  enddo
  deallocate(sig_tmpc)
  deallocate(sig_tmpe)
end subroutine sigma_inner


subroutine contract_in(i,n,nao,nf,ncell,ncellgf2,tn,Gt,Gbt,sigmac,sigmae,vi_c,vi_e,i_vi_c,i_vi_e,n_vi_c,n_vi_e,nmax_c,nmax_e,ni,ne,nbatch,timings,A_cont,B_cont,G_tmp,Gb_tmp,Aval,Arow,Acol)
  use gf2, only : v_ind, v_ind_1, vi_s, diffcell, v_ind_2
  use OMP_LIB
  integer :: i_n, nao, nf, tn, ncell, ncellgf2, nbatch, nc
  integer :: i, n, j, g, g0, j0
  integer :: g_n, n0, g_m, m0, g_n0
  integer :: m1, m2, mi, m, k, l
  integer :: p, q, g_p, g_q, g_q0, p0, q0, g_p0
  integer :: g_k, g_l, g_l0, k0, l0
  integer :: itmp, iitmp, nmax_c, nmax_e
  integer :: l1, n1, p1, nnz
  double precision :: t1, t2, nsparse
  character*256 :: file_name, file_name1
  double precision :: alpha, beta, val
  double precision, dimension(nao,nao,ncell) :: Gbt, Gt, sigmac, sigmae
  double precision, dimension(nf,nf) :: A_cont, B_cont, G_tmp, Gb_tmp
  double precision, dimension(nmax_c,nbatch) :: vi_c
  double precision, dimension(nmax_e,nbatch) :: vi_e
  integer, dimension(4,nmax_c,nbatch) :: i_vi_c
  integer, dimension(4,nmax_e,nbatch) :: i_vi_e
  integer, dimension(nf) :: n_vi_c
  integer, dimension(nf) :: n_vi_e

  double precision, dimension(nf*nf) :: Aval
  integer, dimension(nf*nf) :: Arow, Acol
   

  double precision, dimension(10) :: timings
  double precision, parameter :: thrsh = 1.0d-10

  logical :: sparse

  t1 = omp_get_wtime()
  sparse = .false.
  alpha = 1.0d0
  beta  = 0.0d0
  nc = n - ni + 1
  A_cont(:,:) = 0.0d0

  
  m1 = v_ind_2(i,1)
  m2 = v_ind_2(i,2)

  t2 = omp_get_wtime()
  timings(1) = timings(1) + t2 - t1

  t1 = omp_get_wtime()
  do mi = m1, m2
    m = v_ind(mi,2) 
    q = v_ind(mi,3)
    k = v_ind(mi,4) 
    A_cont(q,k) = A_cont(q,k) + vi_s(mi)*G_tmp(m,n) 
  enddo
  t2 = omp_get_wtime()
  timings(2) = timings(2) + t2 - t1

  t1 = omp_get_wtime()
  nnz = 0
  do q = 1, nf
    do k = 1, nf
      If(dabs(A_cont(q,k)).ge.thrsh) then
        nnz = nnz + 1
        Aval(nnz) = A_cont(q,k)
        Arow(nnz) = q
        Acol(nnz) = k
      Endif
    enddo
  enddo
  nsparse = dble(nnz)/dble(nf*nf)
!  write(*,*)'nsparse = ', nsparse
  If(nsparse.lt.5.0d-02) sparse = .true. 
  t2 = omp_get_wtime()
  timings(1) = timings(1) + t2 - t1 

  t1 = omp_get_wtime()
  If(sparse) then
  call mkl_dcoomm('N',nf,nf,nf,alpha,'GF',Aval,Arow,Acol,nnz,G_tmp,nf,beta,B_cont,nf) 
  Else
!!  call DGEMM('N','N',nf,nf,nf,alpha,Gb_tmp,nf,A_cont,nf,beta,B_cont,nf)
  call DGEMM('N','N',nf,nf,nf,alpha,A_cont,nf,G_tmp,nf,beta,B_cont,nf)
  Endif
  t2 = omp_get_wtime()
  timings(3) = timings(3) + t2 - t1
  
  t1 = omp_get_wtime()  
  t2 = omp_get_wtime()
  timings(1) = timings(1) + t2 - t1
  
  t1 = omp_get_wtime()
!  call DGEMM('N','N',nf,nf,nf,alpha,B_cont,nf,G_tmp,nf,beta,A_cont,nf) 
  call DGEMM('N','N',nf,nf,nf,alpha,Gb_tmp,nf,B_cont,nf,beta,A_cont,nf)
  t2 = omp_get_wtime()
  timings(4) = timings(4) + t2 - t1 

  t1 = omp_get_wtime() 
  do iitmp = 1, n_vi_c(n)
    j = i_vi_c(1,iitmp,nc)
    l = i_vi_c(3,iitmp,nc) 
    p = i_vi_c(4,iitmp,nc) 
    g = (j-1)/nao + 1
    j0 = j - (g-1)*nao
    sigmac(i,j0,g) = sigmac(i,j0,g) + A_cont(p,l)*vi_c(iitmp,nc) 
  enddo
  t2 = omp_get_wtime() 
  timings(5) = timings(5) + t2 - t1

  t1 = omp_get_wtime()
  do iitmp = 1, n_vi_e(n)                                                                           
    j = i_vi_e(1,iitmp,nc)
    l = i_vi_e(2,iitmp,nc)                                                                
    p = i_vi_e(4,iitmp,nc)                                                                
    g = (j-1)/nao + 1                                                                                    
    j0 = j - (g-1)*nao                                                                                   
    sigmae(i,j0,g) = sigmae(i,j0,g) - A_cont(p,l)*vi_e(iitmp,nc)
  enddo
  t2 = omp_get_wtime()
  timings(6) = timings(6) + t2 - t1

  t1 = omp_get_wtime()
100  t2 = omp_get_wtime()
  timings(1) = timings(1) + t2 - t1 
end subroutine contract_in

