subroutine mu_dens(F_kk,S_kk,xmatk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,nts,print_sigl,nodm)
  use seet, only : cells_seet
  use gf2, only : nao, nk, ncell,iwmax
  use gf2, only : mu
  use gf2, only : nuniform
  complex*16, dimension(nao,nao,2*nk) :: F_kk, S_kk,xmatk
  complex*16, dimension(nao,nao,ncell) :: dm_r
  complex*16, dimension(nao,nao,iwmax) :: sig_sol_loc_ndc
  complex*16, dimension(nao,nao) ::sigma_inf_loc 
  double precision mu_max, mu_min
  integer nts
  logical :: print_sigl
  logical :: nodm


  mu_max = 5.0d0 !mu + 1.0d0
  mu_min = -5.0d0 !mu - 1.0d0

  call FalsiMethod(mu,mu_max,mu_min,F_kk,S_kk,xmatk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,nts,nuniform,print_sigl,nodm)
end subroutine mu_dens
subroutine FalsiMethod(mu,mu_max,mu_min,F_kk,S_kk,xmatk,dm_r,sig_sol_loc_ndc,sigma_inf_loc,nts,nuniform,print_sigl,nodm)
  use gf2, only: nk,nao,F_k,S_k,Sigma,iwmax,beta,omegas,nel_cell,ncell,nleg,ncellgf2
  implicit none
  double precision:: mu
  double precision:: mu_max,mu_min
  double precision:: nel_cell_now
  complex*16, dimension(nao,nao,iwmax) :: sig_sol_loc_ndc
  complex*16, dimension(nao,nao) ::sigma_inf_loc 
  complex*16, dimension(nao,nao,2*nk) :: S_kk,F_kk,xmatk
  complex*16, dimension(nao,nao,ncell) :: dm_r
  integer:: nts,nuniform
  logical :: print_sigl
  logical :: nodm
  complex*16,dimension(:,:,:),allocatable:: G1,G2,F_kk_ortho,X_kk,Xk1,Xk2,Xk3,Xk4
  complex*16,dimension(:,:,:,:),allocatable:: sigl,YY
  double precision,save:: tol=1.d-6
  double precision:: mu_upper,mu_lower,el_upper,el_lower
  integer::l,k,j,i,ic,w
  integer::step_max=300

  integer:: n, side=0
! starting values at endpoints of interval 
   double precision:: fs,ft,s,t,r,fr
   complex*16::diff

  allocate(G1(nao,nao,ncell))
  allocate(G2(nao,nao,ncell))
  allocate(F_kk_ortho(nao,nao,2*nk))
  allocate(X_kk(nao,nao,2*nk))
  allocate(sigl(nao,nao,0:nleg,2*nk))
  G1=dcmplx(0.0d0,0.0d0)
  call calc_G1(S_kk,F_kk_ortho,G1)

  

  do k=1,2*nk
     call get_quantities_per_k_point(k,F_kk,S_kk,nts,nuniform,F_kk_ortho,X_kk,sigl(:,:,:,k))
  end do

  if (print_sigl) then
     open(950,file='GF.bin',status='replace',action='write',form='unformatted')
     do k=1,2*nk 
           do i=1,nao
              do j=1,nao

                 do l=0,2
!                 write(950)k,l,i,j,dble(sigl(i,j,l,k)),dimag(sigl(i,j,l,k))
!                    write(*,*)"sigl",k,l,i,j,dble(sigl(i,j,l,k)),dimag(sigl(i,j,l,k))
                 end do

                 

           end do
        end do
     end do
     close(950)

     do k=1,2*nk 
           do i=1,nao
              do j=1,nao

!                 write(950)k,l,i,j,dble(sigl(i,j,l,k)),dimag(sigl(i,j,l,k))
!                    write(*,*)"S/F",k,i,j,F_kk(i,j,k),S_kk(i,j,k)


                 

           end do
        end do
     end do


end if
  

   
   call el_difference(mu,fr,F_kk,S_kk,F_kk_ortho,X_kk,xmatk,sig_sol_loc_ndc,sigma_inf_loc,G1,G2,sigl,dm_r,nts,nuniform,nodm)
   call flush(6)

   r=mu
   if (abs(fr)>tol) then

     s=mu_min

     call el_difference(s,fs,F_kk,S_kk,F_kk_ortho,X_kk,xmatk,sig_sol_loc_ndc,sigma_inf_loc,G1,G2,sigl,dm_r,nts,nuniform,nodm)

     t=mu_max

     call el_difference(t,ft,F_kk,S_kk,F_kk_ortho,X_kk,xmatk,sig_sol_loc_ndc,sigma_inf_loc,G1,G2,sigl,dm_r,nts,nuniform,nodm)

      do n=1,step_max
         r = (fs*t - ft*s)/(fs - ft)
         call el_difference(r,fr,F_kk,S_kk,F_kk_ortho,X_kk,xmatk,sig_sol_loc_ndc,sigma_inf_loc,G1,G2,sigl,dm_r,nts,nuniform,nodm)

         if (abs(fr)>tol) then
            if (fr*ft > 0) then
               t = r
               ft = fr
               if (side==-1) fs=fs/2
               side = -1
            else if (fs * fr > 0) then
               s = r
               fs = fr
               if (side==+1) ft=ft/2
               side = +1
            else
               exit
            end if
         else
            exit
         end if
      end do
      call flush(6)
   end if
   mu=r
   call flush(6)

  if (print_sigl) then
     open(950,file='mu.txt',status='replace',action='write',form='formatted')
     write(950,*)mu
     close(950)
  end if

  deallocate(G1)
  deallocate(G2)
  deallocate(F_kk_ortho)
  deallocate(X_kk)
  deallocate(sigl)
 end subroutine FalsiMethod

subroutine el_difference(mu,el_diff,F_kk,S_kk,F_kk_ortho,X_kk,xmatk,sig_sol_loc_ndc,sigma_inf_loc,G1,G2,sigl,dm_r,nts,nuniform,nodm)
  use seet, only: cells_seet
  use gf2, only: nk,nao,F_k, S_k,Sigma,iwmax,beta,omegas,nel_cell,ncell,nleg
  implicit none
  double precision:: mu
  double precision:: mu_max,mu_min
  double precision::el_diff
  double precision:: nel_cell_now
  integer:: nts,nuniform
  complex*16, dimension(nao,nao,2*nk) :: S_kk,F_kk,F_kk_ortho,X_kk,xmatk
  complex*16,dimension(nao,nao,iwmax):: sig_sol_loc_ndc
  complex*16,dimension(nao,nao):: sigma_inf_loc   
  complex*16, dimension(nao,nao,ncell) :: dm_r,G1,G2
  complex*16, dimension(nao,nao,0:nleg,2*nk) :: sigl
  logical :: nodm
  integer k


  call calc_dm_new_inner(mu,S_kk,F_kk,F_kk_ortho,X_kk,xmatk,sig_sol_loc_ndc,sigma_inf_loc,G1,G2,sigl,nts,nuniform,dm_r,nel_cell_now,nodm)
  el_diff=nel_cell-nel_cell_now
end subroutine el_difference


