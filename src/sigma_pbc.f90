subroutine calc_sigma_outer(Gr_full_tau,nts)
  use gf2, only: nao,ncell,Sigma_tau,max_points_mesh,np,n_threads
  implicit none
  integer:: ic, nts
  double precision, dimension(nao,nao,ncell,nts) :: Gr_full_tau
  integer:: t,tn,i,j,tt
  INTEGER TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  call omp_set_num_threads(n_threads) 
  if (.not. allocated(Sigma_tau)) then
     allocate(Sigma_tau(nao,nao,ncell,nts))
  end if
  Sigma_tau(:,:,:,:)=0.0d0
  
!$OMP PARALLEL DO DEFAULT(none) &                                                                                 
!$OMP SHARED(Sigma_tau,Gr_full_tau,nts) &        
!$OMP PRIVATE(t,tt,tn)
  do t=1,nts
     tn=omp_get_thread_num() 
     tt=nts-t+1
     call calc_sigma_inner(Sigma_tau(:,:,:,t),Gr_full_tau(1,1,1,tt),Gr_full_tau(1,1,1,t))
  end do

!  do t=1,nts
!    do ic = 1, ncell
!     do i=1,nao
!        do j=1,nao
!           write(100,*)t,ic,i,j, Sigma_tau(i,j,ic,t) 
!        end do
!     end do
!    enddo
!  end do
  call flush(6)
  

end subroutine calc_sigma_outer


subroutine calc_sigma_inner(sig,Gbt,Gt)
  use gf2, only: nao,ncell,vi,diffcell,vi_s,v_ind,np,ind,v_ind_last,ind
  implicit none
!  integer:: t
  integer ip, jp
  integer i,i0, m,m0, q,q0, k, j, n, l,l0, p,p0, nf, g, h, f, r, g0, g1, g2, g3, k1, l1, m1, n1, q1, p1, n2, p2, l2, j1, nf1, nf2, indc, nf3, indc1
  double precision val
!  double precision, allocatable, dimension(:,:,:,:)::vi_aux,vi_aux1
  double precision, allocatable, dimension(:)::vi_auxs,vi_aux1s
  integer, allocatable, dimension(:)::v_ind_last1
  integer, allocatable, dimension(:,:)::v_ind1,v_ind2
  double precision, dimension(nao,nao,ncell)::sig
  double precision, dimension(nao,nao,ncell) ::Gbt,Gt
  double precision::thr=1.e-10

  nf=nao*ncell
!  write(*,*)'nao, nf = ', nao, nf
!  stop

  !allocate temporary array
!  allocate(vi_auxs(np*nf))
!  allocate(v_ind1(np*nf,4))
  allocate(vi_auxs(nao*nf*nf*nf))
  allocate(vi_aux1s(nao*nf*nf*nf))
  allocate(v_ind1(nao*nf*nf*nf,4))
  allocate(v_ind2(nao*nf*nf*nf,4))


!  vi_aux=0.0
  vi_auxs = 0.0
  vi_aux1s = 0.0
  v_ind1 = 0
  v_ind2 = 0
  indc = 0 
! 
!  write(*,*) vi_s
!  stop 
  do l = 1, nf
      h = (l-1)/nao+1
      l1 = l - (h-1)*nao
  do nf2 = 1, ind-1!np
     val = 0.0d0
     nf3 = v_ind_last(nf2)
  do nf1 = v_ind_last(nf2), v_ind_last(nf2+1)-1
     k = v_ind(nf1,4)
     g = (k-1)/nao+1
     k1 = k - (g-1)*nao
     g0 = diffcell(g,h)
     if(g0.gt.0) then
       val = val + Gt(k1,l1,g0)*vi_s(nf1)
     endif
  enddo
     if(abs(val).gt.1.0d-20) then 
     indc = indc + 1
     i = v_ind(nf3,1)
     m = v_ind(nf3,2)
     q = v_ind(nf3,3)
     vi_auxs(indc) = val
     v_ind1(indc,1) = i
     v_ind1(indc,2) = m
     v_ind1(indc,3) = q
     v_ind1(indc,4) = l
     endif
  enddo
  enddo
!  write(*,*) vi_auxs
!  stop

  allocate(v_ind_last1(indc+1))
!  write(*,*)'indc = ', indc
  v_ind_last1(1)=1
  i0 = v_ind1(1,1)
  m0 = v_ind1(1,2)
  l0 = v_ind1(1,4)
  indc1 = 1
  do nf1 = 2, indc
    i = v_ind1(nf1,1)
    m = v_ind1(nf1,2)
    l = v_ind1(nf1,4)
    if(i.ne.i0.or.m.ne.m0.or.l.ne.l0) then
       indc1 = indc1 + 1
       v_ind_last1(indc1) = nf1
       i0 = i
       m0 = m
       l0 = l
    endif
  enddo
  indc1 = indc1 + 1
  v_ind_last1(indc1) = indc + 1

!  vi_aux1=0.0

!  allocate(vi_aux1s(indc*nf))
!  allocate(v_ind2(indc*nf,4))


  indc = 0
  do p = 1, nf
      g = (p-1)/nao+1
      p1 = p - (g-1)*nao
  do nf2 = 1, indc1-1!np
     val = 0.0d0
     nf3 = v_ind_last1(nf2)
  do nf1 = v_ind_last1(nf2), v_ind_last1(nf2+1)-1
     q = v_ind1(nf1,3)
     h = (q-1)/nao+1
     q1 = q - (h-1)*nao
     g0 = diffcell(g,h)
     if(g0.gt.0) then
       val = val + Gbt(p1,q1,g0)*vi_auxs(nf1)
     endif
  enddo
     if(abs(val).gt.1.0d-20) then
     indc = indc + 1
     i = v_ind1(nf3,1)
     m = v_ind1(nf3,2)
     l = v_ind1(nf3,4)
     vi_aux1s(indc) = val
     v_ind2(indc,1) = i
     v_ind2(indc,2) = m
     v_ind2(indc,3) = p
     v_ind2(indc,4) = l
     endif
  enddo
  enddo
 
!  write(*,*) vi_aux1s
!  stop

  deallocate(v_ind_last1)
!  deallocate(vi_auxs)
!  deallocate(v_ind1)
  allocate(v_ind_last1(indc+1))

  v_ind_last1(1)=1
  i0 = v_ind2(1,1)
  p0 = v_ind2(1,3)
  l0 = v_ind2(1,4)
  indc1 = 1
  do nf1 = 2, indc
    i = v_ind2(nf1,1)
    p = v_ind2(nf1,3)
    l = v_ind2(nf1,4)
    if(i.ne.i0.or.p.ne.p0.or.l.ne.l0) then
       indc1 = indc1 + 1
       v_ind_last1(indc1) = nf1
       i0 = i
       p0 = p
       l0 = l
    endif
  enddo
  indc1 = indc1 + 1
  v_ind_last1(indc1) = indc + 1

!  vi_aux=0.0

!  allocate(vi_auxs(nf*indc))
  vi_auxs = 0.0
!  allocate(v_ind1(indc*nf,4))
  indc = 0
  do n = 1, nf
     h = (n-1)/nao + 1
     n1 = n - (h-1)*nao
  do nf2 = 1, indc1-1!np
     val = 0.0d0
     nf3 = v_ind_last1(nf2)
  do nf1 = v_ind_last1(nf2), v_ind_last1(nf2+1)-1
     m = v_ind2(nf1,2)
     g = (m-1)/nao + 1
     m1 = m - (g-1)*nao
     g0 = diffcell(g,h)
     if(g0.gt.0) then
       val = val + Gt(m1,n1,g0)*vi_aux1s(nf1) 
     endif
  enddo
     if(abs(val).gt.1.0d-20) then
     indc = indc + 1
       i = v_ind2(nf3,1)
       p = v_ind2(nf3,3)
       l = v_ind2(nf3,4)
!       vi_aux(i,n,p,l)=val
       vi_auxs(indc) = val
     v_ind1(indc,1) = i
     v_ind1(indc,2) = n
     v_ind1(indc,3) = p
     v_ind1(indc,4) = l
     endif
  enddo
  enddo
!  write(*,*)'indc = ', indc, nao*nf*nf*nf
!  do nf1 = 1, indc
!     write(*,*) v_ind1(nf1,1), v_ind1(nf1,2), v_ind1(nf1,3), v_ind1(nf1,4)
!  enddo
!  stop
!  write(*,*) vi_auxs
  sig=0.0d0 !cmplx(0.0d0,0.0d0)
  do nf1 = 1, indc
     i = v_ind1(nf1,1)
     n = v_ind1(nf1,2)
     p = v_ind1(nf1,3)
     l = v_ind1(nf1,4)
        h = (n-1)/nao + 1
        n1 = n - (h-1)*nao
           f = (p-1)/nao + 1
           p1 = p - (f-1)*nao 
                 r = (l-1)/nao + 1
                 l1 = l - (r-1)*nao
                 do j=1,nf
                    g = (j-1)/nao + 1
                    j1 = j - (g-1)*nao
                    g1 = diffcell(g,h) ! n
                    g2 = diffcell(g,f) ! p 
                    g3 = diffcell(g,r) ! l
                    if(g1*g2*g3.gt.0) then
                       n2 = n1 + (g1-1)*nao
                       p2 = p1 + (g2-1)*nao 
                       l2 = l1 + (g3-1)*nao
                       sig(i,j1,g) = sig(i,j1,g) + vi_auxs(nf1)*(2.0d0*vi(j1,n2,l2,p2)-vi(j1,l2,n2,p2)) 
!                       write(*,*) vi(j1,n2,l2,p2), vi(j1,l2,n2,p2)
                    endif
                 end do
  end do
!  write(*,*)'sig'
!  write(*,*) sig
!  stop

!  do i=1,nao
!     do j=1,nao
!        write(*,*)"here5",i,j,sig(i,j)
!     end do
!  end do
!  call flush(6)
  
!  deallocate(vi_aux)
!  deallocate(vi_aux1)

  deallocate(vi_auxs)
  deallocate(vi_aux1s)
  deallocate(v_ind1)
  deallocate(v_ind2)


end subroutine calc_sigma_inner


