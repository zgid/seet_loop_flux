subroutine pbc_sao(S_kk,wannier,U_SAO)
  use gf2, only: nao, ncell, diffcell, nk
  use gf2, only: four_coef, k_weights, nkpw, w_func
  implicit none
  complex*16, dimension(nao,nao,2*nk) :: S_kk, U_SAO
  double precision, dimension(ncell*nao,nao) :: wannier

  double precision, allocatable, dimension(:,:,:) :: U_SAOR


  complex*16, allocatable, dimension(:,:) :: BB

  integer :: i, j, k, l, p, q, r, g, h, gg, hh, g1, nf, nfk, k1, kk, ic

  double precision :: alpha, beta


  alpha = 1.0d0
  beta =  0.0d0

  allocate(BB(nao,nao))
  allocate(U_SAOR(nao,nao,ncell))

  U_SAO(:,:,:) = dcmplx(0.0d0,0.0d0)


  write(*,*)'Sk eigenvectors'
  do k = 1, 2*nk
     write(*,*)"k",k
    call symmetric_ortho_basis(S_kk(:,:,k),U_SAO(:,:,k),nao)
  enddo




  do gg = 1, ncell
     BB=dcmplx(0.0d0,0.0d0)
     call FT_k_to_r_ic(U_SAO,BB,gg)
     U_SAOR(:,:,gg) = dble(BB(:,:))
  enddo

  do i = 1, nao
    do g = 1, ncell
      do j = 1, nao
        wannier((g-1)*nao+j,i) = U_SAOR(j,i,g)
        write(200,*)i,j,g,wannier((g-1)*nao+j,i)
      enddo
    enddo
  enddo
  call flush(6)
  stop
  deallocate(U_SAOR,BB)
end subroutine pbc_sao

