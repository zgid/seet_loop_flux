subroutine FT_k_to_r(Xk, Xr)
  use gf2, only: nk,nao,ncell,four_coef,k_weights
  implicit none
  complex*16,dimension(nao,nao,2*nk) ::Xk
  complex*16,dimension(nao,nao,ncell) ::Xr
  integer:: k,kk,ic,i,j
  double precision:: nrmlz,nkpw
  double precision,dimension(:),allocatable::w_func
  double precision,dimension(:),allocatable::k_w


  do ic=1,ncell
     call FT_k_to_r_ic(Xk, Xr(:,:,ic),ic)
  end do


end subroutine FT_k_to_r

subroutine FT_k_to_r_ic(Xk, Xr,ic1)
  use gf2, only: nk,nao,ncell,four_coef,k_weights
  use gf2, only: nkpw, w_func
  implicit none
  complex*16,dimension(nao,nao,2*nk) ::Xk
  complex*16,dimension(nao,nao) ::Xr
  integer:: k,kk,ic,i,j,ic1
  double precision:: nrmlz!,nkpw

     do k=1,nk

        Xr(:,:)=Xr(:,:)+Xk(:,:,2*k-1)*dconjg(four_coef(k,ic1))
        Xr(:,:)=Xr(:,:)+Xk(:,:,2*k)*four_coef(k,ic1)*w_func(k)

     end do
  Xr(:,:) = Xr(:,:)*nkpw

end subroutine FT_k_to_r_ic

subroutine FT_r_to_single_k(Xk,Xr,k)
  use gf2, only: nk, nao, ncell, four_coef
  implicit none
  complex*16, dimension(nao,nao) :: Xk
  double precision, dimension(nao,nao,ncell) :: Xr
  integer:: k,kk,k1,ic

  Xk(:,:)=dcmplx(0.0d0,0.0d0)
  kk = modulo(k,2)
  k1 = k/2+kk
  if(kk.ne.0) then
     do ic=1,ncell
       Xk(:,:)=Xk(:,:)+Xr(:,:,ic)*four_coef(k1,ic)
     enddo
  else
     do ic=1,ncell             
        Xk(:,:)=Xk(:,:)+Xr(:,:,ic)*dconjg(four_coef(k1,ic))
     end do
  endif

end subroutine FT_r_to_single_k



subroutine calc_ft_gwk_to_gtk2(S_k,F_k,sig_sol_loc_ndc,sigma_inf_loc,mu,nts,uniform,dm_r)
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
  integer w,i,j,n,ii,jj,t,kk,l,ip,tt,ic,np,ww,icc,ipt,ipk,itau
  integer nzcg,th

  integer,dimension(:,:),allocatable::file_io 

  double precision:: mu
  double precision::nel_total,nel_per_cell,val
  double precision::taun
  double precision, save ::pi=4.D0*DATAN(1.D0)



  complex*16::muomega,tmp
  complex*16::nel_t,nel_temp,tmpc
  complex*16::alpha,betaz
  complex*16, dimension(nao,nao,2*nk) :: S_k,F_k
  complex*16, dimension(nao,nao,iwmax)::sig_sol_loc_ndc
  complex*16, dimension(nao,nao)::sigma_inf_loc
  complex*16,dimension(nao,nao,ncell):: dm_r
  complex*16, dimension(:,:),allocatable :: aux
  complex*16, dimension(:,:,:),allocatable :: G1,G2,Fk_ortho,X_k,Gw
  complex*16, dimension(:,:,:,:),allocatable ::sigl,Gkmt
  complex*16, allocatable, dimension(:,:,:,:,:) :: Gtau


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
  allocate(Gtau(nao,nao,nzcg,nts,n_threads))
  allocate(Gkmt(nao,nao,nk,1))
  write(*,*)"allocated"
  call flush(6)
  call calc_G1(S_k,Fk_ortho,G1)

  do k=1,2*nk
     call get_quantities_per_k_point(k,F_k,S_k,nts,uniform,Fk_ortho,X_k,sigl(:,:,:,k))
  end do

  
  call calc_G2(iwmax,mu,X_k,Fk_ortho,sig_sol_loc_ndc,sigma_inf_loc,nts,uniform,G2,sigl)

  
Gtau=dcmplx(0.0d0,0.0d0)

!$OMP PARALLEL DO DEFAULT(none) &                                                                      
!$OMP SHARED(mu,X_k,S_k,Fk_ortho,nts,uniform,G1,G2,sigl,nzcg,iwmax,Gtau,sig_sol_loc_ndc,sigma_inf_loc) &        
!$OMP PRIVATE(w,tn)

  do w= 1,iwmax
     tn=omp_get_thread_num()
     call calc_gwr_inner2(tn+1,w,mu,X_k,Fk_ortho,S_k,Gtau(:,:,:,:,tn+1),sig_sol_loc_ndc(:,:,w),sigma_inf_loc,nts,uniform,G1,G2,sigl,nzcg)
  end do


  Gr_tau(:,:,:,:)=0.0d0
  

     do t=1,nts
        do ic=1, ncell
           do n=1,n_threads
              Gr_tau(:,:,ic,t)=Gr_tau(:,:,ic,t)+dble(Gtau(:,:,ic,t,n))
           end do
        end do
     end do


   ! do ic=1,10
   !     do i=1,nao
   !        do j=1,nao            
   !         write(*,*)"D11",ic,i,j,Gr_tau(i,j,ic,nts),G1(i,j,ic),G2(i,j,ic)
   !         call flush(6)
   !        end do
   !     end do
   !  end do

     do ic=1, ncell
        do t=1,nts
           Gr_tau(:,:,ic,t)=Gr_tau(:,:,ic,t)-0.5*G1(:,:,ic)+0.25*G2(:,:,ic)*(-1.0d0*beta+2.0d0*power_mesh_small(t))
        end do
     end do


   ! do ic=1,10
   !     do i=1,nao
   !        do j=1,nao            
   !         write(*,*)"D111",ic,i,j,Gr_tau(i,j,ic,nts)
   !         call flush(6)
   !        end do
   !     end do
   !  end do
     dm_r(:,:,:)=dcmplx(0.0d0,0.0d0)
     do ic=1, ncell
        t=nts
        dm_r(:,:,ic)=Gr_tau(:,:,ic,t)!+0.25*G2(:,:,ic)*beta
        dm_r(:,:,ic)=-2.0d0*dm_r(:,:,ic)
     end do


!    do ic=1,10
!        do i=1,nao
!           do j=1,nao            
!            write(*,*)"D1",ic,i,j,0.50d0*dble(dm_r(i,j,ic)),dble(Gr_tau(i,j,ic,1)),G1(i,j,ic)!,G2(i,j,ic)
!            call flush(6)
!           end do
!        end do
!     end do

  
  deallocate(G1)
  deallocate(G2)
  deallocate(Fk_ortho)
  deallocate(X_k)
  deallocate(sigl)
  deallocate(aux)
  deallocate(Gw)
  deallocate(Gtau)
  deallocate(Gkmt)
  
end subroutine calc_ft_gwk_to_gtk2

subroutine FT_r_to_k(Xk, Xr)                                                                             
  use gf2, only: nk,nao,ncell,four_coef                                                                  
  implicit none                                                                                          
  complex*16,dimension(nao,nao,2*nk) ::Xk                                                                
  double precision,dimension(nao,nao,ncell) ::Xr                                                                 
  integer:: k,kk,ic                                                                                      
                                                                                                         
  Xk(:,:,:)=dcmplx(0.0d0,0.0d0)                                                                           
  kk=1                                                                                            
  do k=1,nk                                                                                          
     do ic=1,ncell                                                                                       
        Xk(:,:,kk)=Xk(:,:,kk)+Xr(:,:,ic)*four_coef(k,ic)                                                 
        Xk(:,:,kk+1)=Xk(:,:,kk+1)+Xr(:,:,ic)*dconjg(four_coef(k,ic))                                      
     end do                                                                                              
     kk=kk+2                                                                                             
  end do                                                                                                 
                                                                                                         
end subroutine FT_r_to_k

subroutine FT_r_to_k_cmplx(Xk, Xr)
  use gf2, only: nk,nao,ncell,four_coef
  implicit none
  complex*16,dimension(nao,nao,2*nk) ::Xk
  complex*16,dimension(nao,nao,ncell) ::Xr
  integer:: k,kk,ic

  Xk(:,:,:)=dcmplx(0.0d0,0.0d0)
  kk=1
  do k=1,nk
     do ic=1,ncell
        Xk(:,:,kk)=Xk(:,:,kk)+Xr(:,:,ic)*four_coef(k,ic)
        Xk(:,:,kk+1)=Xk(:,:,kk+1)+Xr(:,:,ic)*dconjg(four_coef(k,ic))
     end do
     kk=kk+2
  end do

end subroutine FT_r_to_k_cmplx
