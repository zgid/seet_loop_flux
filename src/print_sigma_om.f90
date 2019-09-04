subroutine plot_sigma_om_real(w,sigl,sigmar)
  use gf2, only:nao,nk,ncell,iwmax,nc,logm,nleg,beta,omegas,ncellgf2
  implicit none
  integer:: w
  complex*16, dimension(nao,nao,0:nleg,2*nk):: sigl
  complex*16, dimension(nao,nao,ncellgf2):: sigmar
  complex*16, dimension(:,:,:),allocatable:: sigma_w
  integer i,k,ic
  allocate(sigma_w(nao,nao,2*nk))

  do k=1,2*nk
     call calc_sigw_using_tnl_inner(w,sigma_w(:,:,k),sigl(:,:,:,k))
  end do
  sigmar(:,:,:) = dcmplx(0.0d0,0.0d0)
  do ic=1, ncellgf2
     call FT_k_to_r_ic(sigma_w(:,:,:), sigmar(:,:,ic),ic)
  end do

  deallocate(sigma_w)

end subroutine plot_sigma_om_real
