subroutine calc_G1(Sk,aux,G1)
  use gf2, only:nao,nk,ncell
  implicit none
  complex*16, dimension(nao,nao,2*nk) ::Sk,aux
  complex*16, dimension(nao,nao,ncell) ::G1
  integer k

  G1=dcmplx(0.0d0,0.0d0)
  !calculate S^{-1} for the first row 
  aux=dcmplx(0.0d0,0.0d0)
  do k=1,2*nk
     aux(:,:,k)=Sk(:,:,k)
     call invert_complex_matrix(nao,aux(:,:,k))
  end do

  call FT_k_to_r(aux,G1)
end subroutine calc_G1

subroutine calc_S_eigval(Sk)
  use gf2, only:nao,nk,ncell
  implicit none
  complex*16, dimension(nao,nao,2*nk) ::Sk
  complex*16, dimension(:,:,:), allocatable ::aux
  double precision, dimension(:,:), allocatable ::eigval
  integer k,i
  allocate(aux(nao,nao,2*nk))
  allocate(eigval(nao,2*nk))

  aux=dcmplx(0.0d0,0.0d0)

  do k=1,2*nk
     aux(:,:,k)=Sk(:,:,k)
     call compute_ortho_basis_complex(Sk(:,:,k),aux(:,:,k),nao,eigval(:,k))
     !do i=1,nao
     !   write(*,*)"k,i,Seig",k,i,eigval(i,k)
     !end do
  end do

  deallocate(aux)
  deallocate(eigval)

end subroutine calc_S_eigval
  




subroutine calc_G2(w,mu,Xk,Fk,sig_sol_loc_ndc,sigma_inf_loc,nts,uniform,G2,sigl)
  use gf2, only:nao,nk,ncell,iwmax,nc,logm,nleg,beta,omegas,S_k
  implicit none
  
  integer:: w,nts,uniform,ncc
  double precision:: mu
  complex*16, dimension(nao,nao,2*nk) :: Fk,Xk
  complex*16,dimension(nao,nao,iwmax):: sig_sol_loc_ndc
  complex*16,dimension(nao,nao)::sigma_inf_loc
  complex*16, dimension(nao,nao,ncell) :: G2
  complex*16, dimension(:,:,:),allocatable:: Xr
  complex*16, dimension(:,:,:),allocatable :: Gk,auxx,sigma_w,sigma
  complex*16, dimension(nao,nao,0:nleg,2*nk):: sigl
  double precision, dimension(:,:),allocatable::sigma1
  integer:: i,j,k,ic,kk,t
  complex*16::om,muomega
  double precision::taun
  complex*16::alpha,betaz
  double precision, save::pi=4.D0*DATAN(1.D0)
  complex*16,dimension(:,:),allocatable:: unity
  complex*16, dimension(:,:),allocatable::aux1,aux
  complex*16, dimension(:,:,:),allocatable :: S_kk
  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)
  taun=pi

  allocate(Xr(nao,nao,ncell))
  allocate(Gk(nao,nao,2*nk))
  allocate(S_kk(nao,nao,2*nk))
  allocate(sigma_w(nao,nao,2*nk))
  allocate(sigma(nao,nao,nts))
  allocate(unity(nao,nao))
  allocate(aux1(nao,nao),aux(nao,nao))
  
  unity=dcmplx(0.0d0,0.0d0)
  do i=1,nao
     unity(i,i)=dcmplx(1.0d0,0.0d0)
  end do
!
  do kk = 1, nk
    S_kk(:,:,2*kk-1) = S_k(:,:,kk)
    S_kk(:,:,2*kk) = dconjg(S_k(:,:,kk))
  enddo

  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))

     muomega =omegas(w)+mu
     aux1=dcmplx(0.0d0,0.0d0)
     aux(:,:)=sig_sol_loc_ndc(:,:,w)+sigma_inf_loc(:,:)
     call ZGEMM('t','n',nao,nao,nao,alpha,dconjg(Xk(:,:,k)),nao,aux(:,:),nao,betaz,aux1,nao)
     aux=dcmplx(0.0d0,0.0d0)
     call ZGEMM('n','n',nao,nao,nao,alpha,aux1,nao,Xk(:,:,k),nao,betaz,aux(:,:),nao)

!     sigma_w(:,:,k)=sigma_w(:,:,k)+sig_sol_loc_ndc(:,:,w)
!     Gk(:,:,k) = (muomega*unity(:,:)-Fk(:,:,k)-sigma_w(:,:,k))
     Gk(:,:,k) = (muomega*S_kk(:,:,k)-Fk(:,:,k)-aux(:,:)-sigma_w(:,:,k))

     call invert_complex_matrix(nao,Gk(:,:,k))

  end do

  G2=dcmplx(0.0d0,0.0d0)
  Xr(:,:,:) = dcmplx(0.0d0,0.0d0) 
  do ic=1,ncell
     call FT_k_to_r_ic(Gk(:,:,:), Xr(:,:,ic),ic)
     G2(:,:,ic)=dble(Xr(:,:,ic)*omegas(w)*omegas(w))
  end do
  

  deallocate(Xr)
  deallocate(Gk)
  deallocate(sigma_w)
  deallocate(sigma)
  deallocate(unity)
  deallocate(aux1,aux)
  deallocate(S_kk)
!

end subroutine calc_G2

                               

subroutine get_quantities_per_k_point(k,F_k,S_k,nts,uniform,Fk_ortho,xmatk,sigl)
  use gf2, only : nk,nao,iwmax,beta,omegas,nc,logm,nleg,four_coef,power_mesh_small
  use gf2, only : Sigma_tau
  implicit none
  double precision:: mu
  integer:: k,nts,uniform
  complex*16, dimension(nao,nao,2*nk) :: S_k,F_k,Fk_ortho,xmatk
  complex*16, allocatable, dimension(:,:) :: aux,aux1,unity,dmr,Fk,xmat
  complex*16, allocatable, dimension(:,:,:) :: sigma
  complex*16::alpha,betaz
  complex*16, dimension(:,:),allocatable :: sig
  complex*16, dimension(:,:,:),allocatable :: rsig,isig
  complex*16, dimension(nao,nao,0:nleg)::sigl

  complex*16, allocatable, dimension(:,:) :: lvec, rvec
  complex*16, allocatable, dimension(:) :: eigval, work, rwork 

  integer w,i,j,t,kk,l,ip,ii,n, info
  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)

  allocate(sig(nao,nao))
  allocate(sigma(nao,nao,nts))
  allocate(rsig(nao,nao,0:nleg),isig(nao,nao,0:nleg))
  allocate(aux(nao,nao))
  allocate(aux1(nao,nao))
  allocate(unity(nao,nao))
  allocate(xmat(nao,nao))
  allocate(Fk(nao,nao))

  allocate(lvec(nao,nao))
  allocate(rvec(nao,nao)) 
  allocate(eigval(nao))
  allocate(work(3*nao))
  allocate(rwork(2*nao))

  sigl=dcmplx(0.0d0,0.0d0)
  unity(:,:)=0.0d0
  do i=1,nao
     unity(i,i)=dcmplx(1.0d0,0.0d0)
  end do

     !------------------------- k-points ---------------------------
  xmat(:,:)=unity(:,:) !cmplx(0.0d0,0.0d0)
  xmatk(:,:,k)=xmat(:,:)

  Fk(:,:)=F_k(:,:,k)


!  call compute_diag_basis_complex(S_k(:,:,k),Fk(:,:),xmat,nao)
!  xmatk(:,:,k)=xmat(:,:)
!  call ZGEMM('t','n',nao,nao,nao,alpha,conjg(xmat),nao,Fk(:,:),nao,betaz,aux1,nao)
!  call ZGEMM('n','n',nao,nao,nao,alpha,aux1,nao,xmat,nao,betaz,Fk(:,:),nao)

  Fk_ortho(:,:,k)=Fk(:,:)
  sigma=dcmplx(0.0d0,0.0d0)

  do t=1,nts
     aux1=dcmplx(0.0d0,0.0d0)
     call FT_r_to_single_k(sigma(:,:,t),Sigma_tau(:,:,:,t),k)
     sig(:,:)=sigma(:,:,t)
     sigma(:,:,t)=dcmplx(0.0d0,0.0d0)
     call ZGEMM('t','n',nao,nao,nao,alpha,dconjg(xmat),nao,sig(:,:),nao,betaz,aux1,nao)
     call ZGEMM('n','n',nao,nao,nao,alpha,aux1,nao,xmat,nao,betaz,sigma(:,:,t),nao)
  end do

  call calc_sigl_from_sigtau(nts,uniform,rsig,dble(sigma(:,:,:)))
  call calc_sigl_from_sigtau(nts,uniform,isig,dimag(sigma(:,:,:)))

  sigl(:,:,:)=rsig(:,:,:)+dcmplx(0.0d0,1.0d0)*isig(:,:,:)

  deallocate(sig)
  deallocate(sigma)
  deallocate(rsig,isig)
  deallocate(aux)
  deallocate(aux1)
  deallocate(unity)
  deallocate(xmat)
  deallocate(Fk)

end subroutine get_quantities_per_k_point

subroutine calc_sigl_from_sigtau(nts,uniform,Gl,Gt)
  use gf2, only:nao,max_points_mesh,power_mesh,power_mesh_small,power_mesh_small_int,uniform_tau, power_tau,nleg,beta,n_threads
  implicit none
  integer::nts,uniform
  double precision:: x,h,sqrtl,x1,x2,x3
  double precision,dimension(:,:),allocatable:: pn1,pn2,pn3,pd
  double precision,dimension(nao,nao,nts)::Gt
  complex*16,dimension(nao,nao,0:nleg)::Gl
  complex*16,dimension(:,:,:,:),allocatable::Glt
  INTEGER TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  integer tn
  integer:: ntaul,ntaur,n_int
  integer:: ntau1,ntau2,ntau3
  integer:: l,p,t,i,j


  call omp_set_num_threads(n_threads)

  n_int = 2*power_tau + 2


  allocate(pn1(0:nleg,n_threads))
  allocate(pn2(0:nleg,n_threads))
  allocate(pn3(0:nleg,n_threads))
  allocate(pd(0:nleg,n_threads))
  allocate(Glt(nao,nao,0:nleg,n_threads))

  Gl(:,:,:)=dcmplx(0.0d0,0.0d0)
  Glt(:,:,:,:)=dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO DEFAULT(none) &                                                                                 
!$OMP SHARED(nts,uniform,Glt,Gt,n_int,pn1,pn2,pn3,pd) &        
!$OMP PRIVATE(p,tn)
  do p=1,n_int
     tn=omp_get_thread_num()
     call calc_sigl_from_sigtau_inner(p,nts,uniform,Glt(:,:,:,tn+1),Gt,pn1(:,tn+1),pn2(:,tn+1),pn3(:,tn+1),pd(:,tn+1))
  end do


  do p=1,n_threads
     Gl(:,:,:)=Gl(:,:,:)+Glt(:,:,:,p)
  end do

   deallocate(pn1)
   deallocate(pn2)
   deallocate(pn3)
   deallocate(pd)
   deallocate(Glt)
 end subroutine calc_sigl_from_sigtau

subroutine calc_sigl_from_sigtau_inner(p,nts,uniform,Gl,Gt,pn1,pn2,pn3,pd)
  use gf2,only:nao,max_points_mesh,power_mesh,power_mesh_small,power_mesh_small_int,uniform_tau,power_tau,nleg,beta
  implicit none
  integer::p,nts,uniform
  double precision:: x,h,sqrtl,x1,x2,x3
  double precision,dimension(0:nleg):: pn1,pn2,pn3,pd
  double precision,dimension(nao,nao,nts)::Gt
  complex*16,dimension(nao,nao,0:nleg)::Gl
  integer:: ntaul,ntaur,n_int
  integer:: ntau1,ntau2,ntau3
  integer:: l,t,i,j

     ntaul = uniform*(p-1)
     ntaur = uniform*(p)
     h=power_mesh_small(ntaul+2)-power_mesh_small(ntaul+1)
     do t=1,uniform/2
        ntau1 = ntaul + 2*t-1
        x1=(power_mesh_small(ntau1)/beta-0.5d0)*2.0d0
        call lpn(nleg,x1,pn1(0),pd(0))
        ntau2 = ntaul + 2*t
        x2=(power_mesh_small(ntau2)/beta-0.5d0)*2.0d0
        call lpn(nleg,x2,pn2(0),pd(0))
        ntau3 = ntaul + 2*t+1
        x3=(power_mesh_small(ntau3)/beta-0.5d0)*2.0d0
        call lpn(nleg,x3,pn3(0),pd(0))
        do l=0,nleg
           sqrtl=dsqrt(2.0d0*l+1.0d0)
           Gl(:,:,l)=Gl(:,:,l)+h/3.0d0*sqrtl*(dble(Gt(:,:,ntau1))*pn1(l)+4.0d0*dble(Gt(:,:,ntau2))*pn2(l)+dble(Gt(:,:,ntau3))*pn3(l))
        end do
     end do

end subroutine calc_sigl_from_sigtau_inner

subroutine calc_dm_new_inner(mu,S_k,F_k,Fk_ortho,X_k,xmatk,sig_sol_loc_ndc,sigma_inf_loc,G1,G2,sigl,nts,uniform,dmr,nel_t,nodm)
  use gf2, only: nk,nao,iwmax,beta,omegas,nc,logm,nleg,ncell,logm,n_threads,ncs
  implicit none

  double precision:: mu
  integer:: k,nts,uniform
  complex*16,dimension(nao,nao,2*nk) :: S_k,F_k,Fk_ortho,X_k,xmatk
  complex*16,dimension(nao,nao,iwmax):: sig_sol_loc_ndc
  complex*16,dimension(nao,nao)::sigma_inf_loc
  complex*16,dimension(nao,nao,ncell):: dmr,G1,G2
  complex*16,dimension(nao,nao,0:nleg,2*nk)::sigl
  complex*16, dimension(:,:),allocatable::auxd,tmpd
  complex*16,dimension(:,:,:,:),allocatable:: dm
  logical :: nodm
  integer :: ncell_el
  double precision::nel_t, tt1,tt2
  integer::ic,i,j,ip,t,n,w
  complex*16::alpha,betaz
  INTEGER TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
  integer tn

  if (nodm) then
     ncell_el=ncs
  else
     ncell_el=ncell
  end if

  !CALL CPU_TIME(tt1)

  call omp_set_num_threads(n_threads)


  allocate(auxd(nao,nao),tmpd(nao,nao))
  allocate(dm(nao,nao,ncell,n_threads))

  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)

  X_k(:,:,:)=xmatk(:,:,:)
  do k=1,2*nk
!     write(*,*)"k",X_k(:,:,k)
!     call flush(6)
     call invert_complex_matrix(nao,X_k(:,:,k))
  end do


  G2=dcmplx(0.0d0,0.0d0)
  call calc_G2(iwmax,mu,X_k,F_k,sig_sol_loc_ndc,sigma_inf_loc,nts,uniform,G2,sigl)
  dm=dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO DEFAULT(none) &                                                                     
!$OMP SHARED(mu,X_k,Fk_ortho,nts,uniform,G1,G2,dm,sigl,iwmax,F_k,ncell_el,sig_sol_loc_ndc,sigma_inf_loc) &        
!$OMP PRIVATE(w,tn)
  do w=1,iwmax
     tn=omp_get_thread_num()
!     call calc_dm_inner_omp(w,mu,X_k,Fk_ortho,nts,uniform,G1,G2,dm(:,:,:,tn+1),sigl)
     call calc_dm_inner_omp(w,mu,X_k,F_k,sig_sol_loc_ndc(:,:,w),sigma_inf_loc,nts,uniform,G1,G2,dm(:,:,:,tn+1),sigl,ncell_el)  
  end do


  dmr=dcmplx(0.0d0,0.0d0)
  do n=1,n_threads
     dmr(:,:,:)=dmr(:,:,:)+dm(:,:,:,n)
  end do

   ! do ic=1,10
   !     do i=1,nao
   !        do j=1,nao            
   !         write(*,*)"D00",ic,i,j,dmr(i,j,ic),G1(i,j,ic),G2(i,j,ic)
   !         call flush(6)
   !        end do
   !     end do
   !  end do


  do ic=1, ncell
     dmr(:,:,ic)=dmr(:,:,ic)-0.5*G1(:,:,ic)+0.25*G2(:,:,ic)*beta
     dmr(:,:,ic)=-2.0d0*dmr(:,:,ic)
  end do


!    do ic=1,10
!        do i=1,nao
 !          do j=1,nao            
!            write(*,*)"D0",ic,i,j,0.50d0*dble(dmr(i,j,ic))
!            call flush(6)
!           end do
!        end do
!     end do

 ! CALL CPU_TIME(tt2) 
  write(*,*)"total time"!,tt2-tt1
  call flush(6)
  

   nel_t=0.0d0
   do ic=1, ncs
      auxd(:,:)=dcmplx(0.0d0,0.0d0)
      tmpd=dcmplx(0.0d0,0.0d0)
      call FT_k_to_r_ic(S_k, tmpd,ic)
      call ZGEMM('n','t',nao,nao,nao,alpha,tmpd,nao,dmr(:,:,ic),nao,betaz,auxd(:,:),nao)
      do i=1,nao
         nel_t=nel_t+dble(auxd(i,i))
      end do
!      write(*,*)ic,nel_t
   end do



   write(*,*)"mu nel sigma",mu,nel_t
   call flush(6)

  deallocate(auxd,tmpd)
  deallocate(dm)

end subroutine calc_dm_new_inner

subroutine calc_dm_inner_omp(w,mu,Xk,Fk,sig_sol_loc_ndc,sigma_inf_loc,nts,uniform,G1,G2,dm,sigl,ncell_el)
  use gf2, only:nao,nk,ncell,iwmax,nc,logm,nleg,beta,omegas,S_k,Sigma_kw
  implicit none

  integer:: w,nts,uniform,ncc
  integer:: ncell_el
  double precision:: mu
  complex*16, dimension(nao,nao):: sig_sol_loc_ndc,sigma_inf_loc
  complex*16, dimension(nao,nao,2*nk) :: Fk,Xk
  complex*16, dimension(nao,nao,ncell) :: G1,G2,dm
  complex*16, dimension(:,:,:),allocatable:: Xr
  complex*16, dimension(:,:,:),allocatable :: Gk,auxx,sigma_w,sigma
  complex*16, dimension(nao,nao,0:nleg,2*nk):: sigl
  double precision, dimension(:,:),allocatable::sigma1
  integer:: i,j,k,ic,kk,t
  complex*16::om,muomega
  double precision::taun
  complex*16::alpha,betaz
  double precision, save::pi=4.D0*DATAN(1.D0)
  complex*16,dimension(:,:),allocatable:: unity
  complex*16, dimension(:,:),allocatable::aux2,aux1,aux

  complex*16, dimension(:,:,:),allocatable :: S_kk
  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)
  taun=pi

  allocate(Xr(nao,nao,ncell))
  allocate(Gk(nao,nao,2*nk))
  allocate(sigma_w(nao,nao,2*nk))
  allocate(sigma(nao,nao,nts))
  allocate(unity(nao,nao))
  allocate(aux2(nao,nao),aux1(nao,nao),aux(nao,nao))
  allocate(S_kk(nao,nao,2*nk)) 
  do kk = 1, nk
    S_kk(:,:,2*kk-1) = S_k(:,:,kk)
    S_kk(:,:,2*kk) = dconjg(S_k(:,:,kk))
  enddo
 

  unity=dcmplx(0.0d0,0.0d0)
  do i=1,nao
     unity(i,i)=dcmplx(1.0d0,0.0d0)
  end do
  muomega =omegas(w)+mu
  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))
     aux1=dcmplx(0.0d0,0.0d0)
     aux(:,:)=sig_sol_loc_ndc(:,:)+sigma_inf_loc(:,:)
     !aux2(:,:)=Xk(:,:,k)
     !call invert_complex_matrix(nao,aux2(:,:))
!     call flush(6)
!     if (k<5 .and. w<10)  write(*,*)"a",k,w,aux(:,:)
!     if (k<5 .and. w<10)  write(*,*)"X",k,w,Xk(:,:,k)
!     call flush(6) 
     call ZGEMM('t','n',nao,nao,nao,alpha,dconjg(Xk(:,:,k)),nao,aux(:,:),nao,betaz,aux1,nao)
     aux=dcmplx(0.0d0,0.0d0)
     call ZGEMM('n','n',nao,nao,nao,alpha,aux1,nao,Xk(:,:,k),nao,betaz,aux(:,:),nao)

!     aux1=dcmplx(0.0d0,0.0d0)
!     call ZGEMM('n','n',nao,nao,nao,alpha,Xk(:,:,k),nao,Gk(:,:,k),nao,betaz,aux1,nao)
!     Gk(:,:,k)=dcmplx(0.0d0,0.0d0)
!     call ZGEMM('n','t',nao,nao,nao,alpha,aux1,nao,conjg(Xk(:,:,k)),nao,betaz,Gk(:,:,k),nao)


!     sigma_w(:,:,k)=sigma_w(:,:,k)+sig_sol_loc_ndc(:,:)
     !muomega =omegas(w)+mu
!     Gk(:,:,k) = (muomega*unity(:,:)-Fk(:,:,k)-sigma_w(:,:,k))
!     if (k<5 .and. w<10) write(*,*)k,w,aux(:,:)
     Gk(:,:,k) = (muomega*S_kk(:,:,k)-Fk(:,:,k)-aux(:,:)-sigma_w(:,:,k)) 

     call invert_complex_matrix(nao,Gk(:,:,k))

  end do


  Xr(:,:,:) = dcmplx(0.0d0,0.0d0)
  do ic=1, ncell_el
     call FT_k_to_r_ic(Gk(:,:,:), Xr(:,:,ic),ic)
     Xr(:,:,ic)=Xr(:,:,ic)-G1(:,:,ic)/omegas(w)-G2(:,:,ic)/(omegas(w)*omegas(w))
     call ft_w_tau_inner(Xr(:,:,ic),dm(:,:,ic),nao,w,taun,beta)
  end do
  

  deallocate(Xr)
  deallocate(Gk)
  deallocate(sigma_w)
  deallocate(sigma)
  deallocate(unity)
  deallocate(aux2,aux1,aux)
  deallocate(S_kk)
!

end subroutine calc_dm_inner_omp

    
subroutine calc_sigw_using_tnl_inner(w,sigr,sigl)
  use gf2, only:nao,nleg,iwmax,tnl,omegas
  implicit none
  integer::w
  complex*16, dimension(nao,nao) :: sigr
  complex*16, dimension(nao,nao,0:nleg) :: sigl
  integer:: l,i,j,n

  sigr(:,:)=dcmplx(0.0d0,0.0d0)
  do l=1,nleg
           sigr(:,:)=sigr(:,:)+tnl(w,l)*sigl(:,:,l-1)
        end do
end subroutine calc_sigw_using_tnl_inner
       

subroutine calc_gwr_inner(w,mu,Xr,Fk,Sk,sig_sol_loc_ndc,sigma_inf_loc,G1,G2,sigl,nzcg,unum)
  use gf2, only:nao,nk,ncell,iwmax,nc,nleg,omegas
  implicit none
  
  integer:: w,nzcg,unum
  double precision:: mu
  complex*16, dimension(nao,nao) :: sig_sol_loc_ndc,sigma_inf_loc
  complex*16, dimension(nao,nao,2*nk) :: Fk,Sk
  complex*16, dimension(nao,nao,ncell) :: G1,G2
  complex*16, dimension(nao,nao,nzcg) :: Xr
  complex*16, dimension(nao,nao,0:nleg,2*nk):: sigl
  complex*16, dimension(:,:,:),allocatable :: Gk,sigma_w
  integer:: i,j,k,ic,kk,t
  complex*16::om,muomega

  
  allocate(Gk(nao,nao,2*nk))
  allocate(sigma_w(nao,nao,2*nk))  


  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))

     muomega =omegas(w)+mu

     sigma_w(:,:,k)=sigma_w(:,:,k)+sig_sol_loc_ndc(:,:)

     Gk(:,:,k) = (muomega*Sk(:,:,k)-Fk(:,:,k)-sigma_inf_loc(:,:)-sigma_w(:,:,k))

     call invert_complex_matrix(nao,Gk(:,:,k))
     do i=1,nao
        write(unum,fmt='(2I5,2I3,2F15.10)') w, k, i, i, dble(Gk(i,i,k)), dimag(Gk(i,i,k))
     end do
  end do

  Xr=dcmplx(0.0d0,0.0d0)
  do ic=1,nzcg
     call FT_k_to_r_ic(Gk(:,:,:), Xr(:,:,ic),ic)
     Xr(:,:,ic)=Xr(:,:,ic)-G1(:,:,ic)/omegas(w)!-G2(:,:,ic)/(omegas(w)*omegas(w))
  end do
  
  deallocate(Gk)
  deallocate(sigma_w)
 
end subroutine calc_gwr_inner


subroutine calc_gwr_inv_inner(w,mu,Xr,Fk,Sk,sig_sol_loc_ndc,sigma_inf_loc,sigl,nzcg)
  use gf2, only:nao,nk,ncell,iwmax,nc,nleg,omegas
  implicit none
  
  integer:: w,nzcg
  double precision:: mu
  complex*16, dimension(nao,nao) :: sig_sol_loc_ndc,sigma_inf_loc
  complex*16, dimension(nao,nao,2*nk) :: Fk,Sk
  complex*16, dimension(nao,nao,nzcg) :: Xr
  complex*16, dimension(nao,nao,0:nleg,2*nk):: sigl
  complex*16, dimension(:,:,:),allocatable :: Gk,sigma_w
  integer:: i,j,k,ic,kk,t
  complex*16::om,muomega

  
  allocate(Gk(nao,nao,2*nk))
  allocate(sigma_w(nao,nao,2*nk))  


  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))

     muomega =omegas(w)+mu

     sigma_w(:,:,k)=sigma_w(:,:,k)+sig_sol_loc_ndc(:,:)

     Gk(:,:,k) = (muomega*Sk(:,:,k)-Fk(:,:,k)-sigma_inf_loc(:,:)-sigma_w(:,:,k))

!     call invert_complex_matrix(nao,Gk(:,:,k))
     
  end do

  Xr=dcmplx(0.0d0,0.0d0)
  do ic=1,nzcg
     call FT_k_to_r_ic(Gk(:,:,:), Xr(:,:,ic),ic)
!     Xr(:,:,ic)=Xr(:,:,ic)-G1(:,:,ic)/omegas(w)!-G2(:,:,ic)/(omegas(w)*omegas(w))
  end do
  
  deallocate(Gk)
  deallocate(sigma_w)
 
end subroutine calc_gwr_inv_inner


subroutine calc_gwr_inner2(th,w,mu,Xk,Fk,Sk,Gtloc,sig_sol_loc_ndc,sigma_inf_loc,nts,uniform,G1,G2,sigl,nzcg)
  use gf2, only:nao,nk,ncell,iwmax,nc,logm,nleg,beta,omegas,n_threads,S_k,Sigma_kw,power_mesh_small
  implicit none
  
  integer:: th,w,nts,uniform,nzcg
  double precision:: mu
  complex*16, dimension(nao,nao,2*nk) :: Fk,Xk,Sk
  complex*16, dimension(nao,nao,ncell) :: G1,G2
  complex*16, dimension(nao,nao)::sig_sol_loc_ndc,sigma_inf_loc
  integer,dimension(n_threads,ncell)::file_io
  complex*16, dimension(:,:,:),allocatable:: Xr
  complex*16, dimension(:,:,:),allocatable :: Gk,auxx,sigma_w,sigma
  complex*16, dimension(nao,nao,0:nleg,2*nk):: sigl
  complex*16, dimension(nao,nao,nzcg,nts):: Gtloc
  double precision, dimension(:,:),allocatable::sigma1
  integer:: i,j,k,ic,kk,t
  integer:: ip
  complex*16::om,muomega
  double precision::taun
  double precision, save::pi=4.D0*DATAN(1.D0)
  complex*16::alpha,betaz
  complex*16,dimension(:,:),allocatable:: unity
  complex*16, dimension(:,:),allocatable::aux1
  complex*16, dimension(:,:,:),allocatable :: S_kk
  alpha=dcmplx(1.0d0,0.0d0)
  betaz=dcmplx(0.0d0,0.0d0)
  

  allocate(Xr(nao,nao,ncell))
  allocate(Gk(nao,nao,2*nk))
  allocate(sigma_w(nao,nao,2*nk))
  allocate(sigma(nao,nao,nts))
  allocate(unity(nao,nao))
  allocate(aux1(nao,nao))
  allocate(S_kk(nao,nao,2*nk))
  
  unity=dcmplx(0.0d0,0.0d0)
  do i=1,nao
     unity(i,i)=dcmplx(1.0d0,0.0d0)
  end do
!  do kk = 1, nk
!    S_kk(:,:,2*kk-1) = S_k(:,:,kk)
!    S_kk(:,:,2*kk) = dconjg(S_k(:,:,kk))
!  enddo

  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))
     muomega =omegas(w)+mu
     sigma_w(:,:,k)=sigma_w(:,:,k)+sig_sol_loc_ndc(:,:)
     Gk(:,:,k) = (muomega*Sk(:,:,k)-Fk(:,:,k)-sigma_inf_loc(:,:)-sigma_w(:,:,k))
     call invert_complex_matrix(nao,Gk(:,:,k))
  end do

  Xr(:,:,:) = dcmplx(0.0d0,0.0d0)
  do ic=1,nzcg
     call FT_k_to_r_ic(Gk(:,:,:), Xr(:,:,ic),ic)
     Xr(:,:,ic)=Xr(:,:,ic)-G1(:,:,ic)/omegas(w)-G2(:,:,ic)/(omegas(w)*omegas(w))
     do t=1,nts
        taun=power_mesh_small(t)*pi/beta
        call ft_w_tau_inner(Xr(:,:,ic),Gtloc(:,:,ic,t),nao,w,taun,beta)
     end do
  end do
  

  deallocate(Xr)
  deallocate(Gk)
  deallocate(sigma_w)
  deallocate(sigma)
  deallocate(unity)
  deallocate(aux1)
  deallocate(S_kk)
!

end subroutine calc_gwr_inner2

subroutine read_sigma_tau(nts)
  use gf2, only: nao,ncell,Sigma_tau
  implicit none
  integer::nts, ic, icc
  integer::ii,jj,n,m,t,tt
  double precision::val, val1

  open(100,file="fort.100",status='old')
  rewind(100)
  do t=1,nts
    do ic = 1, ncell
     do n=1,nao
        do m=1,nao
        read(100,*)tt,icc,ii,jj,val,val1
        call flush(6)
        Sigma_tau(ii,jj,icc,tt)=val
     end do
    enddo
  end do
end do
close(100)
end subroutine read_sigma_tau

subroutine calc_gw_using_tnl
  use gf2, only:nao,nleg,iwmax,tnl,Sigma_rw,Gl,omegas
  implicit none
  integer:: w,l,i,j

  Sigma_rw=dcmplx(0.0d0,0.0d0)
  do w=1,iwmax
     do l=1,nleg
        Sigma_rw(:,:,:,w)=Sigma_rw(:,:,:,w)+tnl(w,l)*Gl(:,:,:,l-1)
     end do
  end do

  call flush(6)

end subroutine calc_gw_using_tnl

subroutine check_symmetry(mat,symmetric)
  use gf2, only: nao,ncell
  implicit none
  logical symmetric
  double precision,dimension(nao,nao,ncell)::mat
  integer i,j,k,l,ic

  symmetric=.True.

 ! !----checking symmetrizing first cell -----!                                                                                                                                                         
  ic=1
   do i=1,nao
      do j=1,i
         if (abs(mat(j,i,ic)-mat(i,j,ic))>1.e-10) then
            write(*,*)"first cell non-symmetric", i,j,mat(j,i,ic),mat(i,j,ic), abs(mat(j,i,ic)-mat(i,j,ic))
            symmetric=.False.
         end if
      enddo
   enddo




   do ic=2,ncell,2
        !----symmetrizing remaining cells -----!                                                                                                                                                
      do i=1,nao
         do j=1,i
            if (abs(mat(j,i,ic+1)-mat(i,j,ic))>1.e-10) then
               symmetric=.False.
               write(*,*)ic, "cell non-symmetric", i,j,mat(j,i,ic+1),mat(i,j,ic),abs(mat(j,i,ic+1)-mat(i,j,ic))
            end if
            if (abs(mat(i,j,ic+1)-mat(j,i,ic))>1.e-10) then
               symmetric=.False.
               write(*,*)ic, "cell non-symmetric", i,j,mat(i,j,ic+1),mat(j,i,ic),abs(mat(i,j,ic+1)-mat(j,i,ic))
            end if
         enddo
      enddo
   end do



 end subroutine check_symmetry



subroutine restore_symmetry(mat,symmetric)
  use gf2, only: nao,ncell
  implicit none
  logical symmetric
  double precision,dimension(nao,nao,ncell)::mat
  integer i,j,k,l,ic

  

 ! !----symmetrizing first cell -----!                                                                                                                
  if (.not. symmetric) then
     ic=1
     do i=1,nao
        do j=1,i
           mat(i,j,ic)=(mat(i,j,ic)+mat(j,i,ic))/2
           mat(j,i,ic)=mat(i,j,ic)
        enddo
     enddo
  end if



  if(.not. symmetric) then
     do ic=2,ncell,2
        !----symmetrizing remaining cells -----!                                                                                                                                                
        do i=1,nao
           do j=1,i
              mat(i,j,ic)=(mat(i,j,ic)+mat(j,i,ic+1))/2
              mat(j,i,ic+1)=mat(i,j,ic)
              
              mat(j,i,ic)=(mat(j,i,ic)+mat(i,j,ic+1))/2
              mat(i,j,ic+1)=mat(j,i,ic)
           enddo
        enddo
     end do
  end if



end subroutine restore_symmetry


subroutine get_sigma_kl(sigl,nts,uniform)
  use gf2, only:nao,nk,nleg,Sigma_tau
  implicit none
  integer::nts,uniform
  complex*16, dimension(1:nao,1:nao,0:nleg,1:2*nk):: sigl
  complex*16, dimension(:,:,:),allocatable :: rsig,isig,sigma_tau_k
  integer k,t

  write(*,*)"I am here 2",nao,nts,nleg
  call flush(6)

  allocate(sigma_tau_k(nao,nao,nts))
  write(*,*)"I am here 3"
  call flush(6)
  allocate(rsig(nao,nao,0:nleg),isig(nao,nao,0:nleg))

  write(*,*)"I am here 4"
  call flush(6)
  sigma_tau_k=dcmplx(0.0d0,0.0d0)
  rsig=dcmplx(0.0d0,0.0d0)
  isig=dcmplx(0.0d0,0.0d0)

  write(*,*)"I am here 5"
  call flush(6)

  do k=1,2*nk

     do t = 1, nts
       call FT_r_to_single_k(sigma_tau_k(:,:,t),Sigma_tau(:,:,:,t),k)
     enddo
     call calc_sigl_from_sigtau(nts,uniform,rsig,dble(sigma_tau_k(:,:,:)))
     call calc_sigl_from_sigtau(nts,uniform,isig,dimag(sigma_tau_k(:,:,:)))


     sigl(:,:,:,k)=rsig(:,:,:)+dcmplx(0.0d0,1.0d0)*isig(:,:,:)
  end do

  deallocate(rsig,isig)
  deallocate(sigma_tau_k)

end subroutine get_sigma_kl
