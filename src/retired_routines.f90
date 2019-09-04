subroutine calc_ft_gwk_to_gtk1(S_k,F_k,mu,nts,uniform)
!this routine is no longer compatible with 
!calc_gwr_inner(th,w,mu,Xk,Fk,Sk,nts,uniform,G1,G2,file_io,sigl,nzcg)
  use gf2, only: nk,nao,ncell,iwmax,beta,omegas,nc,logm
  use gf2, only: nleg,power_mesh_small_int,max_points_mesh,uniform_tau, power_tau,power_mesh_small
  use gf2, only: Gr_tau
  use gf2, only: n_threads
  implicit none
  
  logical k_write

  integer TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  integer tn
  integer lwork,IO
  integer k,nts,uniform
  integer w,i,j,ii,jj,t,kk,l,ip,tt,ic,np,ww,icc,ipt,ipk,itau
  integer nzcg,th

  integer,dimension(:,:),allocatable::file_io 

  double precision:: mu
  double precision::nel_total,nel_per_cell,val
  double precision::taun
  double precision, save ::pi=4.D0*DATAN(1.D0)
  double precision, allocatable, dimension(:,:,:) :: Gtau


  complex*16::muomega,tmp
  complex*16::nel_t,nel_temp,tmpc
  complex*16::alpha,betaz
  complex*16, dimension(nao,nao,2*nk) :: S_k,F_k
  complex*16, dimension(:,:),allocatable :: aux
  complex*16, dimension(:,:,:),allocatable :: G1,G2,Fk_ortho,X_k,Gw
  complex*16, dimension(:,:,:,:),allocatable ::sigl,Gkmt



  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)

  Gr_tau = 0.0d0

  call omp_set_num_threads(n_threads) 

  nzcg=ncell

  allocate(G1(nao,nao,ncell))
  allocate(G2(nao,nao,ncell))
  allocate(Fk_ortho(nao,nao,2*nk))
  allocate(X_k(nao,nao,2*nk))
  allocate(sigl(nao,nao,0:nleg,2*nk))
  allocate(aux(nao,nao))
  allocate(Gw(nao,nao,iwmax))
  allocate(Gtau(nao,nao,max_points_mesh))
  allocate(Gkmt(nao,nao,nk,1))

  call calc_G1(S_k,Fk_ortho,G1)

  do k=1,2*nk
     call get_quantities_per_k_point(k,F_k,S_k,nts,uniform,Fk_ortho,X_k,sigl(:,:,:,k))
  end do

  call calc_G2(iwmax,mu,X_k,Fk_ortho,nts,uniform,G2,sigl)


  allocate(file_io(n_threads,nzcg))
  ip=3000
  do ic=1,nzcg
     do th=1,n_threads
        file_io(th,ic)=ip
        ip=ip+1
        rewind(ip)
     end do
  end do
  write(*,*)'ip, nzcg*n_threads = ', ip, nzcg*n_threads
  !stop
  


!$OMP PARALLEL DO DEFAULT(none) &                                                                      
!$OMP SHARED(mu,X_k,Fk_ortho,S_k,nts,uniform,G1,G2,file_io,sigl,n_threads,nzcg,iwmax,F_k) &        
!$OMP PRIVATE(w,tn)

  do w= 1,iwmax
     tn=omp_get_thread_num()
     call calc_gwr_inner(tn+1,w,mu,X_k,Fk_ortho,S_k,nts,uniform,G1,G2,file_io,sigl,nzcg)
!     call calc_gwr_inner(tn+1,w,mu,X_k,F_k,S_k,nts,uniform,G1,G2,file_io,sigl,nzcg)
  end do

  !call system('rm grtbin.*')
  call system('rm gktbin.*')

  do t=1,nts
     ipt=2000+t
     if(ipt>2999) stop ' ipt should be below 3000 -- or change 3000 for k '
     call open_file_ip_write(t,ipt)
  end do

  do ic=1,nzcg
     Gw=dcmplx(0.0d0,0.0d0)
     do th=1,n_threads
        ip=file_io(th,ic)
        open(ip,status='old')
        np=0  
        rewind(ip)
        do
           read(ip,end=11)ww,icc,Gw(:,:,ww)
        end do
        np=np+1
11      continue
        close(ip)
     end do
     
     Gtau(:,:,:)=0.0d0
     !call frequency_to_time_ft_from_fortran(nao, iwmax, beta,power_tau, uniform_tau,Gw,Gtau)
     do t=1,nts
        tt=power_mesh_small_int(t)
        
        do i=1,nao
           do j=1,nao
              aux(i,j)=dcmplx(Gtau(i,j,tt),0.0d0)
           end do
        end do
        aux(:,:)=aux(:,:)-0.5d0*G1(:,:,ic)!-0.25d0*G2(:,:,ic)*(beta-2.0d0*power_mesh_small(t))
        Gr_tau(:,:,ic,t) = dble(aux(:,:))
        ipt=2000+t
     end do
  end do
  do t=1,nts
     ipt=2000+t
     if(ipt>2999) stop ' ipt should be below 3000 -- or change 3000 for k '
     call close_file_ip_write(t,ipt)
  end do
  
  do ic=1,nzcg
     do th=1,n_threads
        ip=file_io(th,ic)
        open(ip,status='old')
        close(ip,status='delete')
     enddo
  end do   
  

!  do t=1,nts 
!     do ic=1,nzcg
!        do i=1,nao
!           do j=1,nao
!              write(112,*)t,ic,i,j,Gr_tau(i,j,ic,nts)
!           end do
!        end do
!     end do
!  end do


  deallocate(G1)
  deallocate(G2)
  deallocate(Fk_ortho)
  deallocate(X_k)
  deallocate(sigl)
  deallocate(aux)
  deallocate(Gw)
  deallocate(Gtau)
  deallocate(Gkmt)
  
end subroutine calc_ft_gwk_to_gtk1

