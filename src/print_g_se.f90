subroutine print_g_se(S_k,F_k,mu,nts,iterc)

  use gf2, only: nk,nao,ncell,iwmax,beta,omegas,nc,logm
  use gf2, only: nleg,power_mesh_small_int,max_points_mesh,uniform_tau, power_tau,power_mesh_small
  use gf2, only: n_threads, Sigma_tau, nuniform
  implicit none
  
  logical k_write

  integer TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  integer tn
  integer lwork,IO
  integer k ,nts!,uniform
  integer w,i,j,ii,jj,t,kk,l,ip,tt,ic,np,ww,icc,ipt,ipk,itau
  integer nzcg,th
  integer iterc
  character*256 file_name, file_name_se

  integer,dimension(:,:),allocatable::file_io 

  double precision:: mu
  double precision::nel_total,nel_per_cell,val
  double precision::taun
  double precision, save ::pi=4.D0*DATAN(1.D0)


  complex*16::muomega,tmp
  complex*16::nel_t,nel_temp,tmpc
  complex*16::alpha,betaz
  complex*16, dimension(nao,nao,2*nk) :: S_k,F_k
  complex*16, dimension(:,:),allocatable :: aux
  complex*16, dimension(:,:,:),allocatable :: G1,G2,Fk_ortho,X_k,Gw
  complex*16, dimension(:,:,:,:),allocatable ::sigl

  complex*16, dimension(:,:), allocatable :: Gk, sigma_w, unity


!  write(*,*) Sigma_tau
!  stop
    write(file_name,fmt='(A12)') 'GF.txt'
  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)

  call omp_set_num_threads(n_threads) 

  nzcg=ncell

  allocate(G1(nao,nao,ncell))
  allocate(G2(nao,nao,ncell))
  allocate(Fk_ortho(nao,nao,2*nk))
  allocate(X_k(nao,nao,2*nk))
  allocate(sigl(nao,nao,0:nleg,2*nk))
  allocate(aux(nao,nao))
  allocate(Gw(nao,nao,iwmax))
  allocate(unity(nao,nao))

  allocate(Gk(nao,nao))
  allocate(sigma_w(nao,nao))
  unity=dcmplx(0.0d0,0.0d0)
  do i=1,nao
     unity(i,i)=dcmplx(1.0d0,0.0d0)
  end do

    open(950,file=file_name,status='replace',action='write',form='formatted')
!    open(951,file="SEkw.txt",status='replace',action='write',form='formatted')

  do k=1,nk
    kk = 2*k - 1
    call get_quantities_per_k_point(kk,F_k,S_k,nts,nuniform,Fk_ortho,X_k,sigl(:,:,:,kk))
    do w = 1, iwmax
      muomega =omegas(w)+mu
      call calc_sigw_using_tnl_inner(w,sigma_w(:,:),sigl(:,:,:,kk))
!      Gk(:,:) = (muomega*S_k(:,:,k)-Fk_ortho(:,:,k)-sigma_w(:,:)) 
      Gk(:,:) = (muomega*S_k(:,:,kk)-F_k(:,:,kk)-sigma_w(:,:))
      call invert_complex_matrix(nao,Gk(:,:))

      do i = 1, nao
       do j = 1, nao
         write(950,fmt='(2I5,2I3,2F15.10)') w, k, i, j, dble(Gk(i,j)), dimag(Gk(i,j))
       enddo
      enddo


!      do i = 1, nao
!       do j = 1, nao
!         write(951,*) w, k, i, j, dble(sigma_w(i,j)), dimag(sigma_w(i,j))
!       enddo
!      enddo



    enddo
  end do
  close(950)
!  close(951)

  deallocate(G1)
  deallocate(G2)
  deallocate(Fk_ortho)
  deallocate(X_k)
  deallocate(sigl)
  deallocate(aux)
  deallocate(Gw)

  deallocate(Gk)
  deallocate(sigma_w)
 deallocate(unity)
end subroutine print_g_se
